"""
Analysis of the predictions (Member 5 - Analysis).

* temperature sweeps in kitchen units (mg per 100 g oil)
* apparent enthalpy of solution (van't Hoff) and what it says about gamma
* the temperature at which vanillin "oils out" instead of staying solid
* competition between compounds when they dissolve together
* how much of a recipe's spice compound actually ends up in the oil
* flavour index = dissolved concentration / sensory threshold
* sensitivity of the prediction to each input
"""
from dataclasses import replace

import numpy as np

from .data.compounds import CURCUMIN, FLAVOURS, PIPERINE, VANILLIN
from .data.oils import OILS, make_oil
from .data.sensory import THRESHOLDS_MG_PER_KG, recipe_loading_mg_per_100g
from .mixture import solve_mixture
from .sle import R_GAS, solubility
from .unifac import UNIFAC
from .units import binary_mg_per_100g, to_celsius, to_kelvin

MAIN_OILS = ("sunflower", "coconut")


def sweep(solutes=None, oils=MAIN_OILS, t_celsius=np.arange(25, 201, 5)):
    """Rows of solubility results for every solute/oil/temperature."""
    rows = []
    for s in (solutes or FLAVOURS.values()):
        for oname in oils:
            oil = OILS[oname]
            model = UNIFAC([s, oil])
            for tc in t_celsius:
                T = float(to_kelvin(tc))
                r = solubility(s, oil, T, model)
                rows.append(dict(
                    solute=s.name, oil=oname, T_C=float(tc), T_K=T, x=r.x,
                    x_ideal=r.x_ideal, gamma=r.gamma, ln_gamma_C=r.ln_gamma_C,
                    ln_gamma_R=r.ln_gamma_R, status=r.status,
                    mg_per_100g=binary_mg_per_100g(r.x, s, oil),
                    mg_per_100g_ideal=binary_mg_per_100g(min(r.x_ideal, 0.999999), s, oil)
                    if r.x_ideal < 1 else np.inf,
                    brent_iterations=r.iterations))
    return rows


def vant_hoff_enthalpy(solute, oil, t1_c=25.0, t2_c=60.0, n=8):
    """Apparent enthalpy of solution from the slope of ln x vs 1/T (J/mol).

    dH_sol = dHfus + partial molar excess enthalpy of the solute in solution;
    dH_sol > dHfus means mixing is endothermic (gamma falls as T rises).
    """
    T = to_kelvin(np.linspace(t1_c, t2_c, n))
    model = UNIFAC([solute, oil])
    lnx = np.array([np.log(solubility(solute, oil, t, model).x) for t in T])
    slope = np.polyfit(1.0 / T, lnx, 1)[0]
    return -R_GAS * slope


def oiling_out_temperature(solute, oil, t_lo_c=25.0):
    """Lowest temperature where the solute forms a separate liquid phase
    in the oil instead of a solid (None if it never happens below Tm + 50 K)."""
    model = UNIFAC([solute, oil])

    def lle_flag(T):
        return solubility(solute, oil, T, model).status == "LLE-limited"

    Ts = to_kelvin(np.arange(t_lo_c, to_celsius(solute.Tm) + 50, 2.0))
    flags = [lle_flag(T) for T in Ts]
    if not any(flags):
        return None
    k = flags.index(True)
    if k == 0:
        return float(Ts[0])
    lo, hi = Ts[k - 1], Ts[k]
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if lle_flag(mid) else (mid, hi)
    return float(hi)


def cross_effects(oils=MAIN_OILS, t_celsius=(25, 40, 60, 75)):
    """Solubility of each compound with the other two also saturated,
    relative to its solubility alone (ratio > 1: co-solvency)."""
    out = []
    sol = [VANILLIN, CURCUMIN, PIPERINE]
    for oname in oils:
        oil = OILS[oname]
        model = UNIFAC([*sol, oil])
        for tc in t_celsius:
            T = float(to_kelvin(tc))
            mix = solve_mixture(sol, oil, T, model=model)
            for i, s in enumerate(sol):
                alone = binary_mg_per_100g(solubility(s, oil, T).x, s, oil)
                out.append(dict(oil=oname, T_C=tc, solute=s.name,
                                mg_alone=alone, mg_mixture=mix.mg_per_100g[i],
                                ratio=mix.mg_per_100g[i] / alone,
                                newton_iterations=mix.newton_iterations))
    return out


def mixture_matrix(oils=MAIN_OILS, t_celsius=(25, 40, 60, 75)):
    """Every pair and the triple, all solids in excess: solubility of each
    compound in the mixture divided by its solubility alone (1 = no interaction)."""
    from itertools import combinations
    base = [VANILLIN, CURCUMIN, PIPERINE]
    combos = [c for k in (2, 3) for c in combinations(base, k)]
    out = []
    for oname in oils:
        oil = OILS[oname]
        for combo in combos:
            model = UNIFAC([*combo, oil])
            label = " + ".join(s.name for s in combo)
            for tc in t_celsius:
                T = float(to_kelvin(tc))
                mix = solve_mixture(list(combo), oil, T, model=model)
                for i, s in enumerate(combo):
                    alone = binary_mg_per_100g(solubility(s, oil, T).x, s, oil)
                    out.append(dict(oil=oname, mixture=label, T_C=tc, solute=s.name,
                                    mg_alone=alone, mg_in_mixture=float(mix.mg_per_100g[i]),
                                    ratio=float(mix.mg_per_100g[i] / alone),
                                    gamma=float(mix.gamma[i]),
                                    newton_iterations=mix.newton_iterations))
    return out


def recipe_extraction(oils=MAIN_OILS, t_celsius=np.arange(25, 181, 5)):
    """Fraction of each compound added in the recipe that dissolves in the oil."""
    sol = [VANILLIN, CURCUMIN, PIPERINE]
    loads = [recipe_loading_mg_per_100g(s.name) for s in sol]
    out = []
    for oname in oils:
        oil = OILS[oname]
        model = UNIFAC([*sol, oil])
        for tc in t_celsius:
            T = float(to_kelvin(tc))
            mix = solve_mixture(sol, oil, T, loadings_mg_per_100g=loads, model=model)
            for i, s in enumerate(sol):
                thr = THRESHOLDS_MG_PER_KG[s.name]["value"]
                mgkg = mix.mg_per_100g[i] * 10.0
                out.append(dict(oil=oname, T_C=float(tc), solute=s.name,
                                loading_mg_per_100g=loads[i],
                                dissolved_mg_per_100g=mix.mg_per_100g[i],
                                extracted_pct=100 * mix.extracted_fraction[i],
                                saturated=mix.saturated[i],
                                flavour_index=(mgkg / thr) if thr else np.nan,
                                newton_iterations=mix.newton_iterations))
    return out


def sensitivity(solute, oil_name="sunflower", t_c=60.0):
    """Change in predicted solubility (%) when each input is perturbed."""
    oil = OILS[oil_name]
    T = float(to_kelvin(t_c))
    base = binary_mg_per_100g(solubility(solute, oil, T).x, solute, oil)

    def pct(s=solute, o=oil, model=None):
        v = binary_mg_per_100g(solubility(s, o, T, model).x, s, o)
        return 100.0 * (v - base) / base

    cases = {
        "dHfus -10%": pct(replace(solute, dHfus=0.9 * solute.dHfus)),
        "dHfus +10%": pct(replace(solute, dHfus=1.1 * solute.dHfus)),
        "Tm -5 K": pct(replace(solute, Tm=solute.Tm - 5)),
        "Tm +5 K": pct(replace(solute, Tm=solute.Tm + 5)),
        "dCp = dHfus/Tm": pct(replace(solute, dCp=solute.dHfus / solute.Tm)),
    }
    # all UNIFAC interaction parameters scaled by +/-10 %
    for lab, sc in (("a_mn -10%", 0.9), ("a_mn +10%", 1.1)):
        model = UNIFAC([solute, oil])
        model.A = model.A * sc
        model._pure_cache.clear()
        cases[lab] = pct(model=model)
    if oil_name == "sunflower":
        hi_oleic = make_oil("sunflower", {"C16:0": 4.0, "C18:0": 4.0, "C18:1": 82.0, "C18:2": 10.0})
        b2 = binary_mg_per_100g(solubility(solute, hi_oleic, T).x, solute, hi_oleic)
        cases["high-oleic sunflower"] = 100.0 * (b2 - base) / base
    if solute.name == "curcumin":
        enol = replace(solute, groups={"ACH": 6, "AC": 4, "ACOH": 2, "CH3O": 2, "CH=CH": 2,
                                       "CH=C": 1, "OH": 1, "C=O*": 1})
        cases["keto-enol groups"] = pct(enol)
    return base, cases
