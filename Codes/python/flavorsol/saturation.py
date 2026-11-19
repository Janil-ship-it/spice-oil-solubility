"""
Saturation point: at which temperature / amount the oil can take no more
of a compound (Member 3 / Member 5).

How the saturation point is found from the data we have
-------------------------------------------------------
Inputs used (all already in the package):
    * fusion data of the compound: Tm and dHfus (NIST / DSC literature)
    * UNIFAC groups of the compound and of the oil (oil built from its
      fatty-acid composition)
    * the recipe: grams of oil, grams of compound added, temperature

Step 1  Saturation solubility S(T).
        Solid and dissolved compound are in equilibrium when
            ln(x * gamma(x, T)) = ln a_sat(T)
                                = -(dHfus/R)(1/T - 1/Tm)
        This is solved for x with Brent's method (sle.solubility), then
        converted to mg of compound per 100 g of oil.  S(T) is the most the
        oil can hold at T.
Step 2  Capacity of the pan.  capacity (mg) = S(T) * oil mass (g) / 100.
Step 3  Compare with the amount added.
            added >= capacity  ->  SATURATION IS ACHIEVED at this T.
                                   dissolved = capacity, the rest cannot mix
                                   in and stays as solid (or, above Tm, as a
                                   separate liquid layer).
            added <  capacity  ->  not saturated; everything dissolves and
                                   (capacity - added) more could still mix in.
Step 4  Saturation temperature T_sat for this loading.
        Because dissolving is endothermic (dHfus > 0, and the apparent
        enthalpy of solution from the model is 25-42 kJ/mol), S(T) increases
        monotonically with T.  T_sat is the single root of
            S(T) - loading = 0
        found by scanning 25-200 C and refining with Brent's method.
        Below T_sat the oil is saturated; at and above T_sat everything
        dissolves.
Mixtures (several compounds in one oil)
        All compounds share one liquid, so each gamma depends on the others.
        mixture.solve_mixture solves the coupled equations with Newton-Raphson
        and an active-set loop; a compound is saturated when its equilibrium
        equation is active (solid left) and its dissolved amount is below what
        was added.  The extra capacity of a non-saturated compound is found by
        solving the mixture again with that compound in excess.  The
        temperature at which the whole recipe dissolves is found by bisection
        on T (saturation can only disappear as T rises).
"""
import math

import numpy as np
from scipy.optimize import brentq

from .data.compounds import FLAVOURS
from .data.oils import get_oil
from .mixture import solve_mixture
from .sle import solubility
from .unifac import UNIFAC
from .units import binary_mg_per_100g, to_kelvin

T_MIN_C, T_MAX_C = 25.0, 200.0
PRETTY = {"curcumin": "curcumin", "vanillin": "vanillin", "piperine": "piperine"}


def _fmt_mg(mg):
    if mg is None or not math.isfinite(mg):
        return "no limit"
    if mg >= 1000:
        return f"{mg / 1000:.2f} g"
    if mg >= 100:
        return f"{mg:.0f} mg"
    if mg >= 10:
        return f"{mg:.1f} mg"
    return f"{mg:.2f} mg"


def _fmt_T(t):
    return f"{t:.1f} °C"


# ---------------------------------------------------------------- one compound
def solubility_mg(solute, oil, T_C, model=None):
    """(mg per 100 g oil, status) at T_C. inf for full miscibility."""
    r = solubility(solute, oil, float(to_kelvin(T_C)), model)
    mg = binary_mg_per_100g(r.x, solute, oil) if r.status != "miscible" else math.inf
    return float(mg), r.status


def saturation_temperature(solute, oil, loading_mg_per_100g, model=None,
                           t_lo=T_MIN_C, t_hi=T_MAX_C, step=1.0):
    """Temperature (deg C) at which `loading` mg per 100 g oil is exactly saturated.

    Returns (T_sat, kind):
        kind = 'root'        T_sat found inside [t_lo, t_hi]
        kind = 'below_range' the oil is not saturated even at t_lo (T_sat < t_lo)
        kind = 'above_range' the oil is still saturated at t_hi (T_sat > t_hi)
        kind = 'none'        loading <= 0
    """
    if loading_mg_per_100g <= 0:
        return None, "none"
    model = model or UNIFAC([solute, oil])
    L = float(loading_mg_per_100g)

    def g(tc):
        mg, _ = solubility_mg(solute, oil, tc, model)
        # work in log space: S spans five decades; inf (miscible) -> large positive
        return (math.log(mg) if math.isfinite(mg) else 50.0) - math.log(L)

    grid = np.arange(t_lo, t_hi + 1e-9, step)
    vals = [g(t) for t in grid]
    if vals[0] >= 0:
        return None, "below_range"
    for t0, t1, v0, v1 in zip(grid[:-1], grid[1:], vals[:-1], vals[1:]):
        if v0 < 0 <= v1:
            if v1 >= 49.0:          # jump to full miscibility at Tm
                return float(t1 if solute.Tm - 273.15 > t1 else solute.Tm - 273.15), "root"
            return float(brentq(g, t0, t1, xtol=1e-6)), "root"
    return None, "above_range"


def single_saturation(solute_name, oil_name, T_C, oil_g, added_g):
    """Saturation check for one compound in one oil (the website 'Single compound' tab).

    Returns a dict with capacity, dissolved, excess, saturated flag,
    saturation temperature and a printable statement.
    """
    s = FLAVOURS[solute_name]
    oil = get_oil(oil_name)
    model = UNIFAC([s, oil])
    S, status = solubility_mg(s, oil, T_C, model)
    added_mg = 1000.0 * added_g
    loading = 1e5 * added_g / oil_g
    cap_mg = S * oil_g / 100.0 if math.isfinite(S) else math.inf
    saturated = math.isfinite(cap_mg) and added_mg >= cap_mg * (1 - 1e-9) and added_mg > 0
    dissolved = min(added_mg, cap_mg)
    excess = max(added_mg - dissolved, 0.0)
    T_sat, kind = saturation_temperature(s, oil, loading, model)

    name, o = PRETTY[s.name], oil_name
    if status == "miscible":
        text = (f"At {_fmt_T(T_C)} {name} is above its melting point ({s.Tm - 273.15:.0f} °C) and mixes "
                f"with {o} oil in all proportions, so there is no saturation point: "
                f"all {_fmt_mg(added_mg)} mixes in.")
    elif saturated:
        where = ("stays as undissolved solid" if status == "SLE"
                 else "separates as a second liquid layer (molten compound)")
        text = (f"At {_fmt_T(T_C)} saturation is achieved: {oil_g:g} g of {o} oil can hold only "
                f"{_fmt_mg(cap_mg)} of {name} ({_fmt_mg(S)} per 100 g oil). More mixing cannot be done; "
                f"the extra {_fmt_mg(excess)} {where}.")
        if kind == "root":
            text += f" To dissolve all {_fmt_mg(added_mg)}, heat to at least {_fmt_T(T_sat)}."
        elif kind == "above_range":
            text += f" Even at {T_MAX_C:.0f} °C this amount does not fully dissolve."
    elif added_mg <= 0:
        text = (f"At {_fmt_T(T_C)} the oil would be saturated at {_fmt_mg(cap_mg)} of {name} "
                f"({_fmt_mg(S)} per 100 g oil). Nothing has been added yet.")
    else:
        text = (f"At {_fmt_T(T_C)} saturation is not reached: all {_fmt_mg(added_mg)} of {name} dissolves. "
                f"The oil can take {_fmt_mg(cap_mg - added_mg)} more before it saturates "
                f"(saturation point {_fmt_mg(cap_mg)} in {oil_g:g} g oil).")
        if kind == "root":
            text += f" On cooling, this mixture becomes saturated at {_fmt_T(T_sat)}."
        elif kind == "below_range":
            text += f" It stays unsaturated down to {T_MIN_C:.0f} °C."

    return {
        "compound": s.name, "oil": oil_name, "T_C": float(T_C), "oil_g": float(oil_g),
        "added_mg": added_mg, "loading_mg_per_100g": loading,
        "solubility_mg_per_100g": S, "status": status,
        "capacity_mg": cap_mg, "dissolved_mg": dissolved, "excess_mg": excess,
        "saturated": bool(saturated),
        "saturation_T_C": T_sat, "saturation_T_kind": kind,
        "statement": text,
    }


# ---------------------------------------------------------------- mixtures
def _solve(names, oil, T_C, loads):
    sols = [FLAVOURS[n] for n in names]
    return solve_mixture(sols, oil, float(to_kelvin(T_C)), loadings_mg_per_100g=loads)


def _mixture_flags(names, oil, T_C, loads):
    mix = _solve(names, oil, T_C, loads)
    flags = []
    for i, n in enumerate(names):
        short = mix.mg_per_100g[i] < loads[i] * (1 - 1e-7)
        flags.append(bool(mix.saturated[i] or short))
    return mix, flags


def mixture_all_dissolved_temperature(names, oil, loads, t_lo=T_MIN_C, t_hi=T_MAX_C):
    """Lowest T (deg C) at which every compound of the recipe is fully dissolved."""
    def any_sat(tc):
        return any(_mixture_flags(names, oil, tc, loads)[1])

    if not any_sat(t_lo):
        return None, "below_range"
    if any_sat(t_hi):
        return None, "above_range"
    grid = np.arange(t_lo, t_hi + 1e-9, 5.0)
    prev = grid[0]
    for t in grid[1:]:
        if not any_sat(t):
            lo, hi = prev, t
            for _ in range(25):
                mid = 0.5 * (lo + hi)
                lo, hi = (mid, hi) if any_sat(mid) else (lo, mid)
            return float(hi), "root"
        prev = t
    return None, "above_range"


def mixture_saturation(oil_name, T_C, oil_g, added_g: dict, with_temperature=True):
    """Saturation check for several compounds added together to one oil.

    added_g : {'curcumin': g, 'piperine': g, 'vanillin': g} (zeros are dropped)
    """
    oil = get_oil(oil_name)
    names = [n for n in ("vanillin", "curcumin", "piperine") if added_g.get(n, 0) > 0]
    if not names:
        raise ValueError("add at least one compound")
    loads = [1e5 * added_g[n] / oil_g for n in names]
    mix, flags = _mixture_flags(names, oil, T_C, loads)

    rows, lines = [], []
    for i, n in enumerate(names):
        s = FLAVOURS[n]
        added_mg = 1000.0 * added_g[n]
        dissolved = float(mix.mg_per_100g[i]) * oil_g / 100.0
        dissolved = min(dissolved, added_mg)
        excess = max(added_mg - dissolved, 0.0)
        if flags[i]:
            cap_mg = dissolved
            molten = T_C + 273.15 >= s.Tm
            where = ("separates as a second liquid layer" if molten else "stays as undissolved solid")
            line = (f"{n.capitalize()}: at {_fmt_T(T_C)} saturation is achieved in the mixture at "
                    f"{_fmt_mg(cap_mg)} in {oil_g:g} g oil. More {n} cannot be mixed in; "
                    f"the extra {_fmt_mg(excess)} {where}.")
        else:
            # capacity with this compound in excess and the others as added
            l2 = list(loads)
            l2[i] = math.inf
            binary_note = ""
            if T_C + 273.15 >= s.Tm:
                # molten: the solver needs a finite amount, so use the binary
                # liquid-liquid limit of this compound in the oil (estimate)
                S_bin, st = solubility_mg(s, oil, T_C)
                cap_mg = S_bin * oil_g / 100.0 if math.isfinite(S_bin) else math.inf
                binary_note = " (molten; binary liquid-liquid estimate)"
            else:
                m2 = _solve(names, oil, T_C, l2)
                cap_mg = float(m2.mg_per_100g[i]) * oil_g / 100.0
            if math.isfinite(cap_mg):
                line = (f"{n.capitalize()}: at {_fmt_T(T_C)} not saturated, all {_fmt_mg(added_mg)} dissolves; "
                        f"about {_fmt_mg(max(cap_mg - added_mg, 0))} more could mix in before saturation{binary_note}.")
            else:
                line = (f"{n.capitalize()}: at {_fmt_T(T_C)} it is molten and mixes in all proportions "
                        "(no saturation limit).")
        rows.append({
            "compound": n, "added_mg": added_mg, "loading_mg_per_100g": loads[i],
            "dissolved_mg": dissolved, "excess_mg": excess,
            "capacity_in_mixture_mg": cap_mg, "gamma": float(mix.gamma[i]),
            "saturated": bool(flags[i]), "statement": line,
        })
        lines.append(line)

    T_all, kind = (mixture_all_dissolved_temperature(names, oil, loads)
                   if with_temperature else (None, "skipped"))
    sat_names = [r["compound"] for r in rows if r["saturated"]]
    if sat_names:
        listed = sat_names[0] if len(sat_names) == 1 else ", ".join(sat_names[:-1]) + " and " + sat_names[-1]
        head = (f"At {_fmt_T(T_C)} saturation is achieved for {listed} in {oil_name} oil: "
                "more mixing cannot be done for "
                f"{'this compound' if len(sat_names) == 1 else 'these compounds'}.")
        if kind == "root":
            head += f" The whole recipe dissolves only at {_fmt_T(T_all)} or above."
        elif kind == "above_range":
            head += f" The recipe does not fully dissolve even at {T_MAX_C:.0f} °C."
    else:
        head = f"At {_fmt_T(T_C)} no compound is saturated: the whole recipe dissolves in {oil_g:g} g of {oil_name} oil."
        if kind == "below_range":
            head += f" It stays fully dissolved down to {T_MIN_C:.0f} °C."
    return {"oil": oil_name, "T_C": float(T_C), "oil_g": float(oil_g), "rows": rows,
            "all_dissolved_T_C": T_all, "all_dissolved_kind": kind,
            "headline": head, "statements": [head] + lines,
            "newton_iterations": mix.newton_iterations, "notes": mix.notes}
