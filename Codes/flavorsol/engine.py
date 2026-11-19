"""
Prediction engine: mandatory vs optional inputs, decision pathways,
fallbacks, boundary checks and uncertainty propagation.

MANDATORY inputs                      OPTIONAL inputs (skip or replace a step)
------------------------------------  ----------------------------------------------
solute  (library name, or M + Tm +    dHfus, dCp       (else Walden's rule / dCp = 0)
         groups for a new compound)   groups           (needed only for UNIFAC routes)
oil     (name, or fatty-acid comp.)   gamma (+ its T)  -> skips UNIFAC          (route R1)
T_C     (one value or a list)         measurements     -> anchor the model:
                                         same oil, any T             (route R2)
                                         another solvent, any T      (route R3)
                                         headspace pressure          (route R4)
                                      uncertainties    (else defaults below)
                                      oil_grade, loading, monte_carlo

Route R5 (nothing optional given) = pure UNIFAC prediction.

Decision order (see Figure 1 of the manuscript):
  0 check inputs -> stop on errors
  1 fill missing properties with documented fallbacks
  2 choose route:  gamma given? -> R1 ; measurements? -> R2/R3/R4 ; else R5
  3 solve SLE at each target temperature
  4 propagate input uncertainties (first-order + optional Monte Carlo)
  5 post-checks (degradation, smoke point, melting, extrapolation distance)
"""
import json
from dataclasses import dataclass, field, replace

import numpy as np

from . import checks
from .data.compounds import ALL_SOLUTES, SOLVENTS, Compound
from .data.oils import Oil, get_oil, make_oil
from .data.spices import EXTRA_COMPOUNDS, FLAVORDB, flavor_profile
from .derivatives import (dissolution_properties, dlnx_dT, extrapolate_ln_gamma,
                          partial_molar_excess_enthalpy, partial_molar_excess_entropy)
from .sle import R_GAS, ln_a_sat, solubility
from .unifac import UNIFAC
from .units import binary_mg_per_100g, mg_per_100g_to_moles_per_mol_oil

WALDEN_DS = 56.5            # J/(mol K), Walden's rule for the entropy of fusion
# Default 1-sigma uncertainties when the user gives none
DEFAULT_SIGMA = {"T": 1.0, "Tm": 1.0, "dHfus_rel": 0.05, "dHfus_walden_rel": 0.25,
                 # 1-sigma of ln gamma, from validation: mean |ln error| = 0.798 sigma (normal)
                 "model_lngamma": 0.96,   # pure UNIFAC in oils: mean |ln error| 0.77
                 "transfer_lngamma": 0.51,  # R3 solvent transfer: mean |ln error| 0.41
                 "meas_rel": 0.05}


# ----------------------------------------------------------------------------
# input containers
# ----------------------------------------------------------------------------
@dataclass
class Measurement:
    """A measured data point used to anchor the model (optional input).

    kind 'solubility': value + unit at T_C in `solvent` (saturated solution)
         units: 'mole fraction', 'mg/100 g', 'mg/g', 'g/kg', 'mg/mL', 'g/L'
    kind 'headspace' : value = solute partial pressure (Pa) above a solution of
         known mole fraction `x_solution`; `p_sub_Pa` = sublimation pressure of
         the pure solid at the same T. Henry's law gives x_sat = x_solution * p_sub / p.
    """
    value: float
    unit: str
    T_C: float
    solvent: str = None            # default: the target oil
    kind: str = "solubility"
    sigma_rel: float = None        # relative 1-sigma uncertainty of `value`
    sigma_T: float = 0.5           # K
    x_solution: float = None
    p_sub_Pa: float = None
    source: str = ""


@dataclass
class Request:
    solute: object                 # name or Compound           (mandatory)
    oil: object                    # name or Oil                (mandatory)
    T_C: object                    # float or list of floats    (mandatory)
    M: float = None
    Tm: float = None
    dHfus: float = None
    dCp: float = None
    groups: dict = None
    gamma: float = None
    gamma_T_C: float = None
    measurements: list = field(default_factory=list)
    oil_composition: dict = None
    oil_grade: str = "refined"
    uncertainties: dict = field(default_factory=dict)
    loading_mg_per_100g: float = None
    monte_carlo: int = 0
    T_extrapolation: str = "unifac"   # 'unifac' | 'regular' | 'constant'  (routes R1/R2/R4)


@dataclass
class Prediction:
    T_C: float
    route: str
    status: str
    x: float
    mg_per_100g: float
    gamma: float
    sigma_lnx: float = np.nan
    ci95_mg_per_100g: tuple = (np.nan, np.nan)
    variance_share: dict = field(default_factory=dict)
    mc_ci95_mg_per_100g: tuple = None
    derivatives: dict = field(default_factory=dict)
    assumptions: list = field(default_factory=list)
    issues: list = field(default_factory=list)
    anchor: dict = None

    def summary(self):
        ci = self.ci95_mg_per_100g
        s = (f"T = {self.T_C:.1f} C | route {self.route} | {self.status} | "
             f"{self.mg_per_100g:.4g} mg/100 g")
        if np.isfinite(self.sigma_lnx):
            s += f" (95 % range {ci[0]:.3g}-{ci[1]:.3g})"
        return s


# ----------------------------------------------------------------------------
def _resolve_solute(req: Request, issues, assumptions):
    base = None
    if isinstance(req.solute, Compound):
        base = req.solute
    else:
        key = str(req.solute).strip().lower()
        base = ALL_SOLUTES.get(key) or EXTRA_COMPOUNDS.get(key)
        if base is None:
            base = Compound(name=key, M=req.M, groups=req.groups or {})
    over = {k: getattr(req, k) for k in ("M", "Tm", "dHfus", "dCp", "groups")
            if getattr(req, k) is not None}
    s = replace(base, **over) if over else base
    if not s.M:
        issues.append(checks.Issue("error", "M_MISSING", "Molar mass is mandatory for a new compound."))
    # Identity cross-check against the FlavorDB2 record (IIIT-Delhi), when we have one.
    rec = FLAVORDB.get(s.name.lower(), {})
    if rec:
        words, src = flavor_profile(s.name)
        assumptions.append(f"identity: PubChem CID {rec.get('pubchem_cid')}; FlavorDB status "
                           f"{rec.get('status')}; flavour {', '.join(words)} [{src}]")
        try:
            mw_db = float(rec.get("molecular_weight") or "nan")
        except ValueError:
            mw_db = float("nan")
        if s.M and mw_db == mw_db and abs(mw_db - s.M) / s.M > 0.005:
            issues.append(checks.Issue("warning", "MW_MISMATCH",
                                       f"Molar mass {s.M} g/mol differs from FlavorDB ({mw_db} g/mol): "
                                       "check that this is the same compound."))
    issues.extend(checks.check_fusion_data(s.Tm, s.dHfus, s.dCp))
    if checks.has_errors(issues):
        return s
    if s.dHfus is None:                                   # fallback 1: Walden's rule
        s = replace(s, dHfus=WALDEN_DS * s.Tm)
        assumptions.append(f"dHfus not given: Walden's rule dSfus = {WALDEN_DS} J/(mol K) "
                           f"-> dHfus = {s.dHfus / 1e3:.1f} kJ/mol (+/-25 %)")
        req.uncertainties.setdefault("dHfus", DEFAULT_SIGMA["dHfus_walden_rel"] * s.dHfus)
    if s.dCp is None:
        s = replace(s, dCp=0.0)
    if s.dCp == 0.0:
        assumptions.append("dCp = 0 (usual first approximation; its effect is in the uncertainty)")
    return s


def _resolve_oil(req: Request, issues, assumptions):
    if isinstance(req.oil, Oil):
        return req.oil
    name = str(req.oil).replace(" oil", "").strip()
    if req.oil_composition:
        issues.extend(checks.check_composition(req.oil_composition))
        assumptions.append(f"oil built from the fatty-acid composition you supplied")
        return make_oil(name, req.oil_composition)
    try:
        return get_oil(name)
    except KeyError as e:
        issues.append(checks.Issue("error", "OIL_UNKNOWN", str(e) + " Or give oil_composition."))
        return None


def _x_from_measurement(m: Measurement, solute, solvent_obj):
    """Convert a measured value to a saturation mole fraction."""
    if m.kind == "headspace":
        if not (m.x_solution and m.p_sub_Pa and m.value):
            raise ValueError("headspace measurement needs value (Pa), x_solution and p_sub_Pa")
        return m.x_solution * m.p_sub_Pa / m.value
    M_s, M_v = solute.M, solvent_obj.M
    u = m.unit.replace(" ", "").lower()
    if u in ("molefraction", "x"):
        return m.value
    per_g = {"mg/100g": m.value / 100.0, "mg/g": m.value, "g/kg": m.value}.get(u)
    if per_g is None:
        rho = getattr(solvent_obj, "rho", None) or getattr(solvent_obj, "rho20", None)
        if rho is None:
            raise ValueError(f"unit '{m.unit}' needs the solvent density")
        per_mL = {"mg/ml": m.value, "g/l": m.value}.get(u)
        if per_mL is None:
            raise ValueError(f"unknown unit '{m.unit}'")
        per_g = per_mL / rho                                    # mg solute per g solvent
    n_s, n_v = per_g / 1000.0 / M_s, 1.0 / M_v
    return n_s / (n_s + n_v)


def _solvent_object(name, oil):
    if name is None:
        return oil
    key = name.replace(" oil", "").strip()
    if key.lower() in SOLVENTS:
        return SOLVENTS[key.lower()]
    return get_oil(key)


def _build_anchors(req, solute, oil, issues, assumptions):
    """Turn measurements into ln-gamma offsets (delta) relative to UNIFAC,
    or into plain ln-gamma values when UNIFAC groups are missing."""
    anchors = []
    for k, m in enumerate(req.measurements):
        label = f"measurement {k + 1}"
        issues.extend(checks.check_measurement(m.value, m.unit, label=label))
        if checks.has_errors(issues):
            return anchors
        solv = _solvent_object(m.solvent, oil)
        T_ref = m.T_C + 273.15
        x_ref = _x_from_measurement(m, solute, solv)
        x_id = float(np.exp(ln_a_sat(solute, T_ref)))
        issues.extend(checks.check_measurement(x_ref, "mole fraction", x_meas=x_ref,
                                               x_ideal=x_id, label=label))
        ln_g_exp = float(ln_a_sat(solute, T_ref) - np.log(x_ref))
        same_oil = solv is oil or getattr(solv, "name", "") == oil.name
        a = dict(label=label, solvent=solv.name, T_ref=T_ref, x_ref=x_ref, ln_gamma_exp=ln_g_exp,
                 same_oil=same_oil, kind=m.kind,
                 sigma=m.sigma_rel if m.sigma_rel is not None else DEFAULT_SIGMA["meas_rel"])
        if solute.groups:
            mdl = UNIFAC([solute, solv])
            a["delta"] = ln_g_exp - mdl.ln_gamma([x_ref, 1 - x_ref], T_ref)[0]
        elif not same_oil:
            issues.append(checks.Issue("error", "TRANSFER_NEEDS_GROUPS",
                                       f"{label} is in {solv.name}; moving it to {oil.name} needs "
                                       "UNIFAC groups for the solute."))
        anchors.append(a)
    return anchors


def _route(req, solute, anchors):
    if req.gamma is not None:
        return "R1"
    if anchors:
        if any(a["kind"] == "headspace" for a in anchors):
            return "R4"
        return "R2" if all(a["same_oil"] for a in anchors) else "R3"
    return "R5"


def _core_lnx(route, solute, oil, T, req, anchors, model, delta_extra=0.0):
    """ln x and status at temperature T for one route (no uncertainty)."""
    use_unifac = solute.groups and (route == "R3" or req.T_extrapolation == "unifac")
    if route == "R5" or (route in ("R2", "R3", "R4") and use_unifac):
        delta = np.mean([a["delta"] for a in anchors]) if anchors else 0.0
        r = solubility(solute, oil, T, model, ln_gamma_offset=delta + delta_extra)
        return np.log(r.x), r.status, r.gamma, r
    # explicit routes: gamma known -> x = a_sat / gamma
    if route == "R1":
        T_g = (req.gamma_T_C + 273.15) if req.gamma_T_C is not None else T
        lng_ref = np.log(req.gamma)
    else:                                             # R2/R4 without groups
        T_g = np.mean([a["T_ref"] for a in anchors])
        lng_ref = np.mean([a["ln_gamma_exp"] for a in anchors])
    if abs(T - T_g) < 1e-9 or req.T_extrapolation == "constant":
        lng = lng_ref
    elif solute.groups and req.T_extrapolation == "unifac":
        hE = partial_molar_excess_enthalpy(model, [1e-6, 1 - 1e-6], 0.5 * (T + T_g))
        lng = extrapolate_ln_gamma(lng_ref, T_g, T, hE=hE)
    else:
        lng = extrapolate_ln_gamma(lng_ref, T_g, T, hE=None)
    lng += delta_extra
    lnx = float(ln_a_sat(solute, T)) - lng
    status = "SLE" if lnx < 0 else "miscible"
    return min(lnx, 0.0), status, float(np.exp(lng)), None


def predict(req: Request):
    """Run the full workflow. Returns (list_of_Prediction, global_issues)."""
    issues, assumptions = [], []
    issues.extend(checks.check_uncertainties(req.uncertainties))
    solute = _resolve_solute(req, issues, assumptions)
    oil = _resolve_oil(req, issues, assumptions)
    if checks.has_errors(issues):
        return [], issues
    if req.gamma is not None and req.gamma <= 0:
        issues.append(checks.Issue("error", "GAMMA_SIGN", "gamma must be positive."))
        return [], issues
    if not solute.groups and req.gamma is None and not req.measurements:
        issues.append(checks.Issue("error", "NO_ROUTE",
                                   "Without UNIFAC groups you must give gamma or a measurement."))
        return [], issues

    model = UNIFAC([solute, oil]) if solute.groups else None
    if model is not None and model.missing_pairs:
        issues.append(checks.Issue("warning", "UNIFAC_MISSING_PAIRS",
                                   f"no published parameters for main groups {sorted(model.missing_pairs)}; set to 0"))
    anchors = _build_anchors(req, solute, oil, issues, assumptions)
    if checks.has_errors(issues):
        return [], issues
    route = _route(req, solute, anchors)
    route_text = {"R1": "activity coefficient supplied: UNIFAC skipped",
                  "R2": "anchored to measured solubility in the same oil",
                  "R3": "measured in another solvent, transferred to the oil with UNIFAC",
                  "R4": "headspace (partial-pressure) data converted with Henry's law",
                  "R5": "pure UNIFAC prediction"}[route]
    assumptions.append(f"route {route}: {route_text}")
    if route in ("R1", "R2", "R4") and (not solute.groups or req.T_extrapolation == "regular"):
        assumptions.append("temperature extrapolation assumes a regular solution "
                           "(excess entropy = 0, ln gamma proportional to 1/T)")
    elif route in ("R1", "R2", "R4") and req.T_extrapolation == "constant":
        assumptions.append("temperature extrapolation assumes gamma is independent of T")
    elif route in ("R1", "R2", "R4"):
        assumptions.append("temperature dependence of gamma from UNIFAC (Gibbs-Helmholtz)")
    if route == "R3":
        assumptions.append("solvent transfer assumes the measured-vs-UNIFAC error (delta) of the "
                           "solute is the same in the oil")

    sig = dict(req.uncertainties or {})
    sig.setdefault("T", DEFAULT_SIGMA["T"])
    sig.setdefault("Tm", DEFAULT_SIGMA["Tm"])
    sig.setdefault("dHfus", DEFAULT_SIGMA["dHfus_rel"] * solute.dHfus)
    sig.setdefault("dCp", 0.5 * solute.dHfus / solute.Tm)
    if route == "R5":
        sig.setdefault("ln_gamma_model", DEFAULT_SIGMA["model_lngamma"])
    elif route == "R1":
        sig.setdefault("ln_gamma_model", 0.1)
    else:
        spread = np.std([a.get("delta", a["ln_gamma_exp"]) for a in anchors]) if len(anchors) > 1 else 0.0
        meas = np.sqrt(np.mean([a["sigma"] ** 2 for a in anchors]))
        transfer = DEFAULT_SIGMA["transfer_lngamma"] if route == "R3" else 0.0
        sig.setdefault("ln_gamma_model", float(np.sqrt(spread ** 2 + meas ** 2 + transfer ** 2)))

    temps = req.T_C if isinstance(req.T_C, (list, tuple, np.ndarray)) else [req.T_C]
    preds = []
    for tc in temps:
        p_issues = checks.check_temperature(tc, solute.name, oil.name, req.oil_grade, solute.Tm)
        if checks.has_errors(p_issues):
            preds.append(Prediction(tc, route, "rejected", np.nan, np.nan, np.nan, issues=p_issues))
            continue
        T = tc + 273.15
        lnx, status, gam, raw = _core_lnx(route, solute, oil, T, req, anchors, model)
        x = float(np.exp(lnx))
        mg = binary_mg_per_100g(x, solute, oil) if status != "miscible" else np.inf
        pr = Prediction(tc, route, status, x, mg, gam, assumptions=list(assumptions), issues=p_issues)
        if anchors:
            far = max(abs(T - a["T_ref"]) for a in anchors)
            if far > 50:
                p_issues.append(checks.Issue("warning", "LONG_EXTRAPOLATION",
                                             f"{far:.0f} K away from the nearest measurement."))
            pr.anchor = {"n": len(anchors), "delta": [a.get("delta") for a in anchors],
                         "T_ref_C": [a["T_ref"] - 273.15 for a in anchors]}
        if req.loading_mg_per_100g is not None and np.isfinite(mg):
            if req.loading_mg_per_100g > mg:
                p_issues.append(checks.Issue("info", "EXCESS_SOLID",
                                             f"loading {req.loading_mg_per_100g:g} mg/100 g exceeds solubility: "
                                             f"{100 * mg / req.loading_mg_per_100g:.0f} % dissolves, the rest stays solid."))
        if status == "SLE" and model is not None:
            xs = [x, 1 - x]
            pr.derivatives = {
                "dlnx_dT_per_K": dlnx_dT(solute, model, x, T),
                "hE_J_mol": partial_molar_excess_enthalpy(model, xs, T),
                "sE_J_molK": partial_molar_excess_entropy(model, xs, T),
                **dissolution_properties(solute, model, x, T)}
        if status == "SLE":
            _propagate(pr, route, solute, oil, T, req, anchors, model, sig)
        preds.append(pr)
    return preds, issues


def _propagate(pr, route, solute, oil, T, req, anchors, model, sig):
    """First-order (delta-method) propagation with numerical derivatives,
    plus optional Monte Carlo check.  sigma^2(ln x) = sum_i (d ln x/d p_i)^2 sigma_i^2"""
    base = {"T": T, "Tm": solute.Tm, "dHfus": solute.dHfus, "dCp": solute.dCp, "ln_gamma_model": 0.0}
    steps = {"T": 0.05, "Tm": 0.05, "dHfus": 1e-3 * solute.dHfus, "dCp": 0.5, "ln_gamma_model": 1e-3}

    def lnx_of(p):
        s2 = replace(solute, Tm=p["Tm"], dHfus=p["dHfus"], dCp=p["dCp"])
        # anchors were back-calculated with the base fusion data: re-derive them so that
        # errors in Tm/dHfus partly cancel, exactly as they do in reality
        a2 = []
        for a in anchors:
            shift = float(ln_a_sat(s2, a["T_ref"]) - ln_a_sat(solute, a["T_ref"]))
            b = dict(a, ln_gamma_exp=a["ln_gamma_exp"] + shift)
            if "delta" in a:
                b["delta"] = a["delta"] + shift
            a2.append(b)
        return _core_lnx(route, s2, oil, p["T"], req, a2, model, p["ln_gamma_model"])[0]

    grads, var = {}, {}
    for k, h in steps.items():
        if sig.get(k, 0) == 0:
            continue
        up, dn = dict(base), dict(base)
        up[k] += h
        dn[k] -= h
        grads[k] = (lnx_of(up) - lnx_of(dn)) / (2 * h)
        var[k] = (grads[k] * sig[k]) ** 2
    total = sum(var.values())
    s = float(np.sqrt(total))
    pr.sigma_lnx = s
    pr.ci95_mg_per_100g = (pr.mg_per_100g * np.exp(-1.96 * s), pr.mg_per_100g * np.exp(1.96 * s))
    pr.variance_share = {k: v / total for k, v in var.items()} if total > 0 else {}
    if req.monte_carlo:
        rng = np.random.default_rng(1)
        vals = []
        for _ in range(int(req.monte_carlo)):
            p = {k: base[k] + (rng.normal(0, sig[k]) if k in grads else 0.0) for k in base}
            vals.append(lnx_of(p))
        vals = np.array(vals)
        to_mg = lambda lnx: binary_mg_per_100g(float(np.exp(lnx)), solute, oil)
        pr.mc_ci95_mg_per_100g = (to_mg(np.percentile(vals, 2.5)), to_mg(np.percentile(vals, 97.5)))


# ----------------------------------------------------------------------------
def load_request(path):
    """Read a JSON input file (see the inputs/ folder for examples)."""
    with open(path) as f:
        d = json.load(f)
    d = {k: v for k, v in d.items() if not k.startswith("_")}       # "_comment" keys allowed
    d["measurements"] = [Measurement(**{k: v for k, v in m.items() if not k.startswith("_")})
                         for m in d.get("measurements", [])]
    return Request(**d)


def format_report(preds, issues):
    lines = []
    for i in issues:
        lines.append(str(i))
    if preds:
        lines.append("Assumptions:")
        lines += [f"  - {a}" for a in preds[0].assumptions]
    for p in preds:
        lines.append(p.summary())
        if p.variance_share:
            share = ", ".join(f"{k} {100 * v:.0f}%" for k, v in
                              sorted(p.variance_share.items(), key=lambda kv: -kv[1]))
            lines.append(f"    uncertainty comes from: {share}")
        if p.mc_ci95_mg_per_100g:
            lines.append(f"    Monte Carlo 95 % range: {p.mc_ci95_mg_per_100g[0]:.3g}-{p.mc_ci95_mg_per_100g[1]:.3g}")
        if p.derivatives:
            d = p.derivatives
            lines.append(f"    d ln x/dT = {d['dlnx_dT_per_K']:.4f} 1/K; dH_sol = {d['dH_sol'] / 1e3:.1f} kJ/mol; "
                         f"hE = {d['hE_J_mol'] / 1e3:.2f} kJ/mol (Gibbs-Helmholtz)")
        for i in p.issues:
            lines.append(f"    {i}")
    return "\n".join(lines)
