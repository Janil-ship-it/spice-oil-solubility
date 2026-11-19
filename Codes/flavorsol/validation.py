"""
Validation of the model against literature and lab data (Member 4 - Validation).

Level 1  UNIFAC implementation vs the DDBST reference solution
Level 2  UNIFAC + SLE vs measured solubility of vanillin analogues in pure
         solvents that carry the same functional groups as triglycerides
         (esters, ketones, alcohols) - ethyl vanillin and vanillyl alcohol
Level 3  Partitioning check: vanillin octanol-water partition coefficient
Level 4  Predictions vs published solubilities in edible oils
Level 5  The group's own lab measurements (data/lab_measurements.csv)
"""
import csv
import os

import numpy as np

from .data.compounds import SOLVENTS, get_compound
from .data.oils import get_oil
from .sle import ideal_solubility, solubility
from .unifac import UNIFAC
from .units import binary_mg_per_100g, to_kelvin

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(os.path.dirname(HERE), "data")


def aard(calc, exp):
    """Average absolute relative deviation in percent."""
    calc, exp = np.asarray(calc, float), np.asarray(exp, float)
    return 100.0 * np.mean(np.abs(calc - exp) / exp)


def unifac_benchmark():
    """n-hexane(1)/2-butanone(2), 333.15 K, x1 = 0.5 (DDBST reference solution)."""
    m = UNIFAC([SOLVENTS["n-hexane"], SOLVENTS["2-butanone"]])
    calc = m.gamma([0.5, 0.5], 333.15)
    ref = np.array([1.4276025835, 1.3646545010])
    return {"system": "n-hexane / 2-butanone, 333.15 K, x1 = 0.5",
            "gamma_calc": calc.tolist(), "gamma_ref": ref.tolist(),
            "max_abs_error": float(np.max(np.abs(calc - ref)))}


def load_literature(path=None):
    path = path or os.path.join(DATA_DIR, "literature_solubility.csv")
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def validate_pure_solvents(rows=None):
    """Return per-point and per-system results for analogues in pure solvents."""
    rows = [r for r in (rows or load_literature()) if r["group"] == "pure solvent"]
    points, systems = [], {}
    models = {}
    for r in rows:
        s, v = get_compound(r["solute"]), get_compound(r["solvent"])
        key = (s.name, v.name)
        models.setdefault(key, UNIFAC([s, v]))
        T, x_exp = float(r["T_K"]), float(r["x_exp"])
        res = solubility(s, v, T, models[key])
        p = dict(solute=s.name, solvent=v.name, T=T, x_exp=x_exp, x_calc=res.x,
                 x_ideal=float(ideal_solubility(s, T)), gamma=res.gamma, status=res.status)
        points.append(p)
        systems.setdefault(key, []).append(p)
    summary = []
    for (s, v), pts in systems.items():
        xe = [p["x_exp"] for p in pts]
        summary.append(dict(solute=s, solvent=v, n=len(pts),
                            T_range=(min(p["T"] for p in pts), max(p["T"] for p in pts)),
                            AARD_unifac=aard([p["x_calc"] for p in pts], xe),
                            AARD_ideal=aard([p["x_ideal"] for p in pts], xe)))
    overall = dict(n=len(points),
                   AARD_unifac=aard([p["x_calc"] for p in points], [p["x_exp"] for p in points]),
                   AARD_ideal=aard([p["x_ideal"] for p in points], [p["x_exp"] for p in points]))
    return points, summary, overall


def vanillin_water():
    """Hard case: vanillin in water (UNIFAC is known to struggle with water)."""
    s, w = get_compound("vanillin"), SOLVENTS["water"]
    out = []
    for r in load_literature():
        if r["group"] != "water":
            continue
        T, c = float(r["T_K"]), float(r["value"])           # g per L ~ g per kg water
        x_exp = (c / s.M) / (c / s.M + 1000.0 / w.M)
        res = solubility(s, w, T)
        out.append(dict(T=T, x_exp=x_exp, x_calc=res.x, status=res.status,
                        x_ideal=float(ideal_solubility(s, T))))
    return out


def log_kow(solute, T=295.15):
    """Octanol-water partition coefficient from UNIFAC infinite-dilution gammas.

    Equal activity of the solute in both phases: x_o gamma_o = x_w gamma_w, so
    K_ow = c_o / c_w = (gamma_w / gamma_o) (V_w / V_o).
    Water-saturated octanol: x_octanol = 0.725, x_water = 0.275,
    V_o = 0.725*157.6 + 0.275*18.07 = 119.2 mL/mol; V_w = 18.07 mL/mol.
    """
    w, o = SOLVENTS["water"], SOLVENTS["1-octanol"]
    m = UNIFAC([solute, w, o])
    eps = 1e-9
    ln_gw = m.ln_gamma([eps, 1 - eps, 0.0], T)[0]
    ln_go = m.ln_gamma([eps, 0.275 * (1 - eps), 0.725 * (1 - eps)], T)[0]
    V_w, V_o = 18.07, 0.725 * 157.6 + 0.275 * 18.07
    return float(np.log10(np.exp(ln_gw - ln_go) * V_w / V_o))


def database_check():
    """Cross-check our compound library against FlavorDB2 (IIIT-Delhi) records.

    Molar mass must agree (identity check). XLogP is FlavorDB2's atom-based
    estimate, independent of UNIFAC, so the difference from our UNIFAC log Kow
    shows where the group assignment is least reliable.
    """
    from .data.spices import FLAVORDB, EXTRA_COMPOUNDS
    out = []
    for name, r in FLAVORDB.items():
        if r.get("status") != "retrieved":
            continue
        c = EXTRA_COMPOUNDS.get(name) or get_compound(name)
        lk = log_kow(c)
        xl = float(r["xlogp"])
        out.append({"compound": name, "M_library": round(c.M, 3),
                    "M_flavordb": float(r["molecular_weight"]),
                    "logKow_unifac": round(lk, 2), "xlogp_flavordb": xl,
                    "difference": round(lk - xl, 2), "source": r["flavordb_source"]})
    return out


def validate_oils(rows=None):
    """Compare predicted solubility in oils (as mg per mL oil) with literature."""
    rows = [r for r in (rows or load_literature()) if r["group"].startswith("oil")]
    out = []
    for r in rows:
        s, oil = get_compound(r["solute"]), get_oil(r["solvent"])
        T = float(r["T_K"])
        res = solubility(s, oil, T)
        mg100 = binary_mg_per_100g(res.x, s, oil)
        pred_mg_mL = mg100 / 100.0 * oil.rho20        # literature volumes at room T
        exp = float(r["value"])
        out.append(dict(solute=s.name, oil=r["solvent"], T=T, exp_mg_mL=exp,
                        pred_mg_mL=pred_mg_mL, ratio=pred_mg_mL / exp,
                        kind=r["group"], source=r["source"], note=r["note"],
                        ideal_mg_mL=binary_mg_per_100g(float(ideal_solubility(s, T)), s, oil)
                        / 100.0 * oil.rho20, status=res.status))
    return out


def validate_lab(path=None):
    """Compare with the group's own measurements, if the file has data.

    Expected columns (see data/lab_measurements_TEMPLATE.csv):
    solute, oil, T_C, solubility_mg_per_100g, uncertainty_mg_per_100g, ...
    """
    path = path or os.path.join(DATA_DIR, "lab_measurements.csv")
    if not os.path.exists(path):
        return None
    with open(path, newline="") as f:
        rows = [r for r in csv.DictReader(f) if r.get("solute")]
    if not rows:
        return None
    out = []
    for r in rows:
        s, oil = get_compound(r["solute"]), get_oil(r["oil"])
        T = float(to_kelvin(float(r["T_C"])))
        res = solubility(s, oil, T)
        pred = binary_mg_per_100g(res.x, s, oil)
        exp = float(r["solubility_mg_per_100g"])
        unc = float(r["uncertainty_mg_per_100g"] or "nan")
        out.append(dict(solute=s.name, oil=r["oil"], T_C=float(r["T_C"]), exp=exp, unc=unc,
                        pred=pred, rel_dev=100 * (pred - exp) / exp, status=res.status))
    return out


def extrapolation_validation(rows=None):
    """Test the data-anchored routes on the analogue data (leave-out style).

    (a) temperature extrapolation, route R2: anchor each solute-solvent system at
        its lowest temperature and predict the remaining temperatures with
        anchored UNIFAC, with the regular-solution fallback (no groups) and with
        a constant gamma; compare with pure UNIFAC.
    (b) solvent transfer, route R3: anchor at 298.15 K (or nearest) in one solvent
        and predict every point in every other solvent.
    Returns dicts of AARD (%) and mean |ln(x_calc/x_exp)|.
    """
    from .derivatives import extrapolate_ln_gamma
    from .sle import ln_a_sat
    rows = [r for r in (rows or load_literature()) if r["group"] == "pure solvent"]
    systems = {}
    for r in rows:
        systems.setdefault((r["solute"], r["solvent"]), []).append((float(r["T_K"]), float(r["x_exp"])))
    res_a = {"anchored UNIFAC (R2)": [], "regular solution, no groups": [],
             "constant gamma": [], "pure UNIFAC (R5)": []}
    for (sn, vn), pts in systems.items():
        s, v = get_compound(sn), get_compound(vn)
        m = UNIFAC([s, v])
        pts = sorted(pts)
        T0, x0 = pts[0]
        lng0 = float(ln_a_sat(s, T0) - np.log(x0))
        delta = lng0 - m.ln_gamma([x0, 1 - x0], T0)[0]
        for T, xe in pts[1:]:
            res_a["anchored UNIFAC (R2)"].append((solubility(s, v, T, m, ln_gamma_offset=delta).x, xe))
            res_a["pure UNIFAC (R5)"].append((solubility(s, v, T, m).x, xe))
            lnas = float(ln_a_sat(s, T))
            res_a["regular solution, no groups"].append(
                (min(np.exp(lnas - extrapolate_ln_gamma(lng0, T0, T)), 1.0), xe))
            res_a["constant gamma"].append((min(np.exp(lnas - lng0), 1.0), xe))
    res_b = {"anchored transfer (R3)": [], "pure UNIFAC (R5)": []}
    by_solute = {}
    for (sn, vn), pts in systems.items():
        by_solute.setdefault(sn, {})[vn] = sorted(pts)
    for sn, solv in by_solute.items():
        s = get_compound(sn)
        for ref, pts in solv.items():
            v_ref = get_compound(ref)
            T0, x0 = min(pts, key=lambda p: abs(p[0] - 298.15))
            m_ref = UNIFAC([s, v_ref])
            delta = float(ln_a_sat(s, T0) - np.log(x0)) - m_ref.ln_gamma([x0, 1 - x0], T0)[0]
            for tgt, tpts in solv.items():
                if tgt == ref:
                    continue
                v = get_compound(tgt)
                m = UNIFAC([s, v])
                for T, xe in tpts:
                    res_b["anchored transfer (R3)"].append((solubility(s, v, T, m, ln_gamma_offset=delta).x, xe))
                    res_b["pure UNIFAC (R5)"].append((solubility(s, v, T, m).x, xe))

    def stats(pairs):
        c, e = np.array(pairs).T
        return {"n": len(pairs), "AARD": aard(c, e), "mean_abs_ln": float(np.mean(np.abs(np.log(c / e))))}
    return ({k: stats(v) for k, v in res_a.items()}, {k: stats(v) for k, v in res_b.items()})


def uncertainty_check(n_mc=2000):
    """Compare first-order (delta-method) and Monte Carlo 95 % ranges."""
    from .engine import Request, predict
    out = []
    for sname, oil, tc in (("curcumin", "sunflower", 60.0), ("piperine", "coconut", 60.0),
                           ("vanillin", "sunflower", 40.0)):
        p = predict(Request(sname, oil, [tc], monte_carlo=n_mc))[0][0]
        out.append(dict(solute=sname, oil=oil, T_C=tc, mg=p.mg_per_100g,
                        delta_lo=p.ci95_mg_per_100g[0], delta_hi=p.ci95_mg_per_100g[1],
                        mc_lo=p.mc_ci95_mg_per_100g[0], mc_hi=p.mc_ci95_mg_per_100g[1],
                        share=p.variance_share))
    return out


def metrics(pred, exp):
    """Error metrics for predicted vs measured solubility (any one unit).

    Solubilities span orders of magnitude, so the log-based metrics are the
    meaningful headline numbers; absolute ones are dominated by the largest values.
    """
    p, e = np.asarray(pred, float), np.asarray(exp, float)
    lr = np.log10(p / e)
    return {
        "n": int(len(e)),
        "MAE": float(np.mean(np.abs(p - e))),
        "RMSE": float(np.sqrt(np.mean((p - e) ** 2))),
        "AARD_percent": float(100 * np.mean(np.abs(p - e) / e)),
        "mean_abs_log10": float(np.mean(np.abs(lr))),           # 0.30 = typical factor 2
        "rmse_log10": float(np.sqrt(np.mean(lr ** 2))),
        "bias_log10": float(np.mean(lr)),                        # > 0 = over-prediction
        "typical_factor": float(10 ** np.mean(np.abs(lr))),
        "within_factor_2_percent": float(100 * np.mean(np.abs(lr) <= np.log10(2))),
    }
