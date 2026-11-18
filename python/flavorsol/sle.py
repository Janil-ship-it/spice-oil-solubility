"""
Solid-liquid equilibrium of ONE flavour compound in an oil
(Member 2 - Single-solute solver).

Equilibrium condition (pure solid solute in contact with the liquid solution):

    ln(x_s * gamma_s) = ln a_sat(T)
    ln a_sat(T) = -(dHfus/R)(1/T - 1/Tm) + (dCp/R)[Tm/T - 1 - ln(Tm/T)]

a_sat is the "ideal solubility" (the solubility if gamma = 1). Because gamma
depends on x, the equation is implicit; we solve it with Brent's method in the
variable u = ln x, which keeps x positive and handles solubilities spanning
many orders of magnitude.

Edge cases handled
------------------
* Several roots / liquid-liquid split: UNIFAC can predict a miscibility gap
  between a molten (or sub-cooled) solute and oil. We scan ln a(x) on a grid;
  if it is not monotonic, the binary liquid-liquid equilibrium (LLE) is solved.
  If the LLE activity is below a_sat, the solid "oils out" and the oil-rich
  liquid composition is the solubility limit.
* T >= Tm: no stable solid exists. The solute is either fully miscible with
  the oil or limited by LLE.
"""
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq, root

from .unifac import UNIFAC

R_GAS = 8.314462618  # J/(mol K)


def ln_a_sat(solute, T):
    """Ideal-solubility activity of the pure solid solute at T (K)."""
    T = np.asarray(T, dtype=float)
    lna = -(solute.dHfus / R_GAS) * (1.0 / T - 1.0 / solute.Tm)
    if solute.dCp:
        lna = lna + (solute.dCp / R_GAS) * (solute.Tm / T - 1.0 - np.log(solute.Tm / T))
    return lna


def ideal_solubility(solute, T):
    """Mole-fraction solubility for an ideal solution (gamma = 1), capped at 1."""
    return np.minimum(np.exp(ln_a_sat(solute, T)), 1.0)


@dataclass
class SLEResult:
    solute: str
    solvent: str
    T: float
    x: float                    # solubility, mole fraction of solute
    gamma: float                # activity coefficient at the solution
    ln_gamma_C: float
    ln_gamma_R: float
    x_ideal: float
    status: str                 # 'SLE', 'LLE-limited', 'miscible'
    iterations: int = 0
    function_calls: int = 0
    n_roots: int = 1
    lle: tuple = None           # (x_oil_rich, x_solute_rich, a_LLE) if a gap exists
    message: str = ""
    history: list = field(default_factory=list)


def _x_grid(n=500):
    lo = np.logspace(-12, np.log10(0.5), n)
    hi = 1.0 - np.logspace(-9, np.log10(0.5), n // 2)
    return np.unique(np.concatenate([lo, hi]))


def _ln_activity(model, x, T):
    lg = model.ln_gamma([x, 1.0 - x], T)
    return np.log(x) + lg[0], np.log(1.0 - x) + lg[1]


def _lower_hull(xs, gs):
    """Indices of the lower convex hull of points (xs sorted ascending)."""
    hull = []
    for i in range(len(xs)):
        while len(hull) >= 2:
            i1, i2 = hull[-2], hull[-1]
            cross = (xs[i2] - xs[i1]) * (gs[i] - gs[i1]) - (gs[i2] - gs[i1]) * (xs[i] - xs[i1])
            if cross <= 0:
                hull.pop()
            else:
                break
        hull.append(i)
    return hull


def binary_lle(model, T, grid=None):
    """Solve the binary liquid-liquid equilibrium solute(0)+solvent(1).

    Returns (x_a, x_b, ln_a_solute) with x_a < x_b, or None if no gap.
    Initial guesses come from the convex hull of dG_mix/RT; Newton-type
    refinement (scipy 'hybr') enforces equal activities of both components.
    """
    xs = _x_grid(400) if grid is None else grid
    g = np.array([model.g_mix([x, 1 - x], T) for x in xs])
    hull = _lower_hull(xs, g)
    best = None
    for a, b in zip(hull[:-1], hull[1:]):
        if b - a > 1:
            seg = slice(a + 1, b)
            chord = g[a] + (g[b] - g[a]) * (xs[seg] - xs[a]) / (xs[b] - xs[a])
            excess = np.max(g[seg] - chord)
            if excess > 1e-7 and (best is None or excess > best[2]):
                best = (a, b, excess)
    if best is None:
        return None
    xa0, xb0 = xs[best[0]], xs[best[1]]

    # unknowns: v0 = ln x_a (solute-lean phase), v1 = ln(1 - x_b) (solute-rich phase);
    # this keeps both compositions well conditioned even when x_b -> 1
    def unpack(v):
        return np.exp(v[0]), 1.0 - np.exp(v[1])

    def eqs(v):
        xa, xb = unpack(v)
        a1 = _ln_activity(model, xa, T)
        a2 = _ln_activity(model, xb, T)
        return [a1[0] - a2[0], a1[1] - a2[1]]

    v0 = [np.log(xa0), np.log(1.0 - xb0)]
    best_sol = None
    for method in ("hybr", "lm"):
        sol = root(eqs, v0, method=method, options={"xtol": 1e-14})
        if np.max(np.abs(eqs(sol.x))) < 1e-9:
            best_sol = sol
            break
    if best_sol is None:
        xa, xb = xa0, xb0                          # fall back to the hull estimate
    else:
        xa, xb = sorted(unpack(best_sol.x))
        if abs(xb - xa) < 1e-6:                    # collapsed onto the trivial solution
            xa, xb = xa0, xb0
    return float(xa), float(xb), float(_ln_activity(model, xa, T)[0])


def solubility(solute, solvent, T, model=None, keep_history=False, ln_gamma_offset=0.0):
    """Solubility of a solid solute in a solvent (oil) at T (K).

    Returns an SLEResult. `model` can be a pre-built UNIFAC([solute, solvent]).
    `ln_gamma_offset` (delta) is added to UNIFAC's ln gamma of the solute; it is
    how measured data "anchor" the model (see engine.py, routes 2-3).
    """
    model = model or UNIFAC([solute, solvent])
    T = float(T)
    lnas = float(ln_a_sat(solute, T))
    x_id = float(min(np.exp(lnas), 1.0))
    history = []

    def f(u):
        x = np.exp(u)
        val = u + model.ln_gamma([x, 1 - x], T)[0] + ln_gamma_offset - lnas
        if keep_history:
            history.append((x, val))
        return val

    # 1) scan ln a(x) to detect non-monotonic behaviour (=> liquid-liquid gap)
    xs = _x_grid(300)
    lna = np.array([_ln_activity(model, x, T)[0] for x in xs]) + ln_gamma_offset
    nonmono = np.any(np.diff(lna) < -1e-10)
    lle = binary_lle(model, T) if nonmono else None

    # 2) roots of ln a(x) = ln a_sat on the grid
    resid = lna - lnas
    idx = np.where(np.sign(resid[:-1]) != np.sign(resid[1:]))[0]
    n_roots = len(idx)

    def finish(x, status, it=0, fc=0, msg=""):
        lnC, lnR = model.ln_gamma_parts([x, 1 - x], T) if x < 1 else (np.zeros(2), np.zeros(2))
        return SLEResult(solute.name, solvent.name, T, float(x),
                         float(np.exp(lnC[0] + lnR[0] + ln_gamma_offset)), float(lnC[0]), float(lnR[0]),
                         x_id, status, it, fc, n_roots, lle, msg, history)

    below_melting = T < solute.Tm
    if lle is not None and lnas >= lle[2] + ln_gamma_offset:
        # the liquid-liquid split is reached before solid saturation (or no solid exists)
        msg = ("solute forms a separate liquid phase; oil-rich phase composition "
               "is the solubility limit")
        return finish(lle[0], "LLE-limited", msg=msg)

    if n_roots == 0:
        if not below_melting:
            return finish(1.0, "miscible", msg="above Tm and fully miscible with the oil")
        if ln_gamma_offset < 0:
            # a negative anchoring offset can push the activity below a_sat everywhere
            return finish(1.0, "miscible", msg="with this gamma correction the solute is fully miscible")
        raise RuntimeError(f"No SLE root found for {solute.name} in {solvent.name} at {T} K")

    i = idx[0]                                   # oil-rich (smallest) root
    u_lo, u_hi = np.log(xs[i]), np.log(xs[i + 1])
    u, r = brentq(f, u_lo, u_hi, xtol=1e-14, rtol=1e-14, full_output=True)
    msg = "" if n_roots == 1 else f"{n_roots} roots; oil-rich root selected (stable branch)"
    if not below_melting:
        msg = "above Tm: liquid solute with a miscibility gap"
    return finish(np.exp(u), "SLE", r.iterations, r.function_calls, msg)


def temperature_sweep(solute, solvent, temperatures):
    """Solubility at each temperature (K). Returns a list of SLEResult."""
    model = UNIFAC([solute, solvent])
    return [solubility(solute, solvent, T, model) for T in temperatures]


def ln_a_sat_numeric(solute, T, dCp_of_T=None):
    """ln a_sat by numerical integration (for a temperature-dependent dCp).

    Gibbs-Helmholtz along the subcooled-liquid path:
        d ln a_sat / dT = dH_fus(T) / (R T^2),  dH_fus(T) = dHfus + int_Tm^T dCp(T') dT'
    integrated from Tm (where ln a_sat = 0) to T with scipy quad.
    With a constant dCp this reproduces the closed form in ln_a_sat() exactly,
    so numerical integration is only worth it when dCp really varies with T.
    """
    from scipy.integrate import quad
    f = dCp_of_T if dCp_of_T is not None else (lambda t: solute.dCp)

    def dH(t):
        return solute.dHfus + quad(f, solute.Tm, t)[0]

    val, _ = quad(lambda t: dH(t) / (R_GAS * t ** 2), solute.Tm, float(T))
    return val                      # signed: quad from Tm to T
