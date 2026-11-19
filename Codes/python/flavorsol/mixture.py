"""
Several flavour compounds dissolving in the same oil (Member 3 - Mixture solver).

All solutes share one liquid phase, so each one's activity coefficient depends
on how much of the OTHERS has dissolved. For every solute i that is in excess
(solid present) we need

    F_i = ln x_i + ln gamma_i(x_1, ..., x_n, x_oil; T) - ln a_sat,i(T) = 0

which is a coupled nonlinear system. We solve it with a damped Newton-Raphson
method in the variables u_i = ln n_i (n_i = moles dissolved per mole of oil):

    J du = -F,   J_ij = dF_i/du_j  (central finite differences)
    u <- u + lambda * du,  lambda halved until ||F|| decreases (line search)

Mass-limited dissolution (real cooking): a recipe adds a fixed amount of each
compound. A compound either saturates (solid left over, F_i = 0) or dissolves
completely (n_i = n_i,total, x_i gamma_i <= a_sat,i). Which case applies is
found with an active-set loop around the Newton solver.
"""
from dataclasses import dataclass, field

import numpy as np

from .sle import ln_a_sat, solubility
from .unifac import UNIFAC
from .units import mg_per_100g_oil, mg_per_100g_to_moles_per_mol_oil


@dataclass
class MixtureResult:
    T: float
    oil: str
    solutes: list
    x: np.ndarray                 # liquid mole fractions (solutes..., oil)
    n: np.ndarray                 # mol dissolved per mol oil
    gamma: np.ndarray
    saturated: list               # True if solid remains (solubility-limited)
    mg_per_100g: np.ndarray
    extracted_fraction: np.ndarray  # dissolved / added (nan if no loading)
    newton_iterations: int
    residual_history: list = field(default_factory=list)
    notes: list = field(default_factory=list)


def _x_from_n(n):
    tot = 1.0 + n.sum()
    return np.append(n, 1.0) / tot


def _residuals(model, n, T, lnas, active):
    x = _x_from_n(n)
    lg = model.ln_gamma(x, T)
    return (np.log(x[:-1]) + lg[:-1] - lnas)[active], x, lg


def _newton(model, n0, T, lnas, active, tol=1e-10, max_iter=50, history=None):
    """Damped Newton-Raphson on u = ln n for the active (saturated) solutes."""
    n = n0.copy()
    idx = np.where(active)[0]
    if len(idx) == 0:
        return n, 0
    u = np.log(n[idx])

    def F_of(u_vec):
        nn = n.copy()
        nn[idx] = np.exp(u_vec)
        return _residuals(model, nn, T, lnas, active)[0]

    F = F_of(u)
    for it in range(1, max_iter + 1):
        if history is not None:
            history.append(float(np.max(np.abs(F))))
        if np.max(np.abs(F)) < tol:
            n[idx] = np.exp(u)
            return n, it - 1
        h = 1e-6
        J = np.empty((len(idx), len(idx)))
        for j in range(len(idx)):
            e = np.zeros(len(idx))
            e[j] = h
            J[:, j] = (F_of(u + e) - F_of(u - e)) / (2 * h)
        du = np.linalg.solve(J, -F)
        step = np.max(np.abs(du))
        if step > 2.0:                       # trust region in log space
            du *= 2.0 / step
        lam, normF = 1.0, np.linalg.norm(F)
        while lam > 1e-6:
            F_new = F_of(u + lam * du)
            if np.linalg.norm(F_new) < normF:
                break
            lam *= 0.5
        u = u + lam * du
        F = F_new
    raise RuntimeError(f"Newton-Raphson did not converge at T = {T:.2f} K "
                       f"(max|F| = {np.max(np.abs(F)):.2e})")


def solve_mixture(solutes, oil, T, loadings_mg_per_100g=None, model=None):
    """Solve multi-solute SLE in one oil at temperature T (K).

    Parameters
    ----------
    solutes : list of Compound (solid solutes)
    oil     : Oil pseudo-component
    loadings_mg_per_100g : None (every solute in excess) or a list with the
              amount of each compound added, in mg per 100 g oil.
    """
    k = len(solutes)
    model = model or UNIFAC([*solutes, oil])
    lnas = np.array([float(ln_a_sat(s, T)) for s in solutes])
    notes, history = [], []

    if loadings_mg_per_100g is None:
        n_tot = np.full(k, np.inf)
    else:
        n_tot = np.array([mg_per_100g_to_moles_per_mol_oil(L, s, oil)
                          for L, s in zip(loadings_mg_per_100g, solutes)], dtype=float)

    # initial guess and starting active set from single-solute solubilities
    n0 = np.empty(k)
    active = np.zeros(k, dtype=bool)
    for i, s in enumerate(solutes):
        if T >= s.Tm:
            if not np.isfinite(n_tot[i]):
                raise ValueError(f"{s.name} is molten at {T:.1f} K (Tm = {s.Tm} K); "
                                 "give a loading to treat it as a dissolved liquid")
            r = solubility(s, oil, T)
            n_cap = np.inf if r.status == "miscible" else r.x / (1.0 - r.x)
            n0[i] = min(n_tot[i], n_cap)
            if n_tot[i] > n_cap:
                notes.append(f"{s.name} is molten and exceeds its liquid-liquid limit; "
                             "the excess forms a separate liquid phase (binary estimate)")
            else:
                notes.append(f"{s.name} above its melting point: dissolved as a liquid")
            continue
        r = solubility(s, oil, T)
        n_single = r.x / (1.0 - r.x)
        if n_single < n_tot[i]:
            n0[i], active[i] = n_single, True
        else:
            n0[i] = n_tot[i]
    n0 = np.maximum(n0, 1e-30)

    total_iter = 0
    for _ in range(4 * k + 2):          # active-set iterations
        n, it = _newton(model, n0, T, lnas, active, history=history)
        total_iter += it
        x = _x_from_n(n)
        lg = model.ln_gamma(x, T)
        changed = False
        for i in range(k):
            if active[i] and n[i] > n_tot[i]:             # would dissolve more than added
                active[i], n[i], changed = False, n_tot[i], True
            elif (not active[i]) and T < solutes[i].Tm and \
                    np.log(x[i]) + lg[i] - lnas[i] > 1e-9:  # supersaturated
                active[i], changed = True, True
        n0 = n
        if not changed:
            break

    x = _x_from_n(n)
    gam = np.exp(model.ln_gamma(x, T))
    mg = np.array([mg_per_100g_oil(x[i], s.M, x[-1], oil.M) for i, s in enumerate(solutes)])
    frac = np.where(np.isfinite(n_tot), n / np.where(np.isfinite(n_tot), n_tot, 1.0), np.nan)
    return MixtureResult(T, oil.name, [s.name for s in solutes], x, n, gam,
                         [bool(a) for a in active], mg, frac, total_iter, history, notes)
