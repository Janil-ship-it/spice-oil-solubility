"""
Build the submission Data folder (../../Data) from the model.

    python export_data.py          (run after run_all.py; about 3 minutes)

Every file is CSV or TXT. Data/README.txt lists what each file contains.
Inputs  : everything the code reads (compound data, oil compositions, UNIFAC
          parameters, literature data, database records, scenario definitions).
Outputs : everything the code produces, for every scenario handled
          (single compound, saturation point, mixtures, recipes, routes R1-R5,
          validation, sensitivity, uncertainty).
"""
import contextlib
import csv
import glob
import io
import json
import math
import os
import shutil
import warnings

import numpy as np

from flavorsol.data.compounds import FLAVOURS
from flavorsol.data.oils import FATTY_ACIDS, OIL_COMPOSITIONS, OILS
from flavorsol.data.unifac_params import SUBGROUPS, _PAIRS, MAIN_GROUP_NAMES
from flavorsol.saturation import mixture_saturation, saturation_temperature, single_saturation
from flavorsol.sle import solubility
from flavorsol.unifac import UNIFAC
from flavorsol.units import binary_mg_per_100g, to_kelvin

warnings.filterwarnings("ignore", message="UNIFAC: no published parameters")
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "..", "Data"))
IN, OUT = os.path.join(DATA, "inputs"), os.path.join(DATA, "outputs")
MAIN_OILS = ("sunflower", "coconut")
ORDER = ("vanillin", "curcumin", "piperine")
TADKA_MG = {"vanillin": 15.0, "curcumin": 47.0, "piperine": 50.0}   # in 15 g oil
OIL_G = 15.0


def num(v, sig=6):
    if v is None:
        return ""
    if isinstance(v, (bool, np.bool_)):
        return str(bool(v))
    if isinstance(v, (int, np.integer)):
        return int(v)
    if isinstance(v, (float, np.floating)):
        if not math.isfinite(v):
            return "inf" if v > 0 else ""
        return float(f"{v:.{sig}g}")
    return v


def write_csv(path, rows, header_note=None):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        if not rows:
            f.write("no rows\n")
            return
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow({k: num(v) for k, v in r.items()})


def write_txt(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(text.rstrip() + "\n")


# ------------------------------------------------------------------ inputs
def export_inputs():
    rows = [{"compound": k, "formula": c.formula, "M_g_per_mol": c.M, "Tm_K": c.Tm,
             "Tm_C": c.Tm - 273.15, "dHfus_J_per_mol": c.dHfus, "dCp_J_per_mol_K": c.dCp,
             "unifac_groups": "; ".join(f"{g} x{n:g}" for g, n in c.groups.items()),
             "source": c.source} for k, c in FLAVOURS.items()]
    write_csv(os.path.join(IN, "compound_properties.csv"), rows)

    rows = []
    for o, comp in OIL_COMPOSITIONS.items():
        oil = OILS[o]
        r = {"oil": o, "used_for": "main study" if o in MAIN_OILS else "literature validation only",
             "pseudo_triglyceride_M_g_per_mol": oil.M, "density_20C_g_per_mL": oil.rho20}
        for fa in FATTY_ACIDS:
            r[f"{fa}_wt_pct"] = comp.get(fa, 0.0)
        r["unifac_groups"] = "; ".join(f"{g} x{n:.4f}" for g, n in oil.groups.items())
        rows.append(r)
    write_csv(os.path.join(IN, "oil_fatty_acid_compositions.csv"), rows)

    write_csv(os.path.join(IN, "unifac_subgroups.csv"),
              [{"subgroup": g, "subgroup_id": s.sid, "main_group": s.main,
                "main_group_name": MAIN_GROUP_NAMES.get(s.main, ""), "R": s.R, "Q": s.Q,
                "note": "derived by group additivity" if s.sid > 900 else "DDBST / Hansen 1991"}
               for g, s in sorted(SUBGROUPS.items(), key=lambda kv: kv[1].sid)])
    write_csv(os.path.join(IN, "unifac_interaction_parameters.csv"),
              [{"main_n": n, "main_m": m, "name_n": MAIN_GROUP_NAMES[n], "name_m": MAIN_GROUP_NAMES[m],
                "a_nm_K": a, "a_mn_K": b} for (n, m), (a, b) in sorted(_PAIRS.items())])

    for f in ("literature_solubility.csv", "flavordb_records.csv", "flavordb_spice_entities.csv",
              "spicerx_links.csv", "lab_measurements_TEMPLATE.csv"):
        shutil.copy(os.path.join(HERE, "data", f), os.path.join(IN, f))

    rdir = os.path.join(IN, "route_inputs_R1_to_R5")
    os.makedirs(rdir, exist_ok=True)
    for p in sorted(glob.glob(os.path.join(HERE, "inputs", "*.json"))):
        shutil.copy(p, os.path.join(rdir, os.path.basename(p).replace(".json", ".txt")))

    sc = []
    for c in ORDER:
        for o in MAIN_OILS:
            for T in (25, 60, 100, 140, 180):
                for label, mg in (("tadka", TADKA_MG[c]), ("high", 1000.0)):
                    sc.append({"scenario": f"S-{c}-{o}-{T}C-{label}", "compound": c, "oil": o,
                               "T_C": T, "oil_g": OIL_G, "added_mg": mg})
    write_csv(os.path.join(IN, "scenarios_saturation_single.csv"), sc)
    ms = []
    recipes = {"tadka_V+C+P": dict(TADKA_MG),
               "turmeric_pepper_C+P": {"vanillin": 0.0, "curcumin": 47.0, "piperine": 50.0},
               "high_loading_V+C+P": {"vanillin": 1000.0, "curcumin": 500.0, "piperine": 2000.0}}
    for o in MAIN_OILS:
        for T in (25, 60, 100, 140):
            for name, mg in recipes.items():
                ms.append({"scenario": f"M-{name}-{o}-{T}C", "oil": o, "T_C": T, "oil_g": OIL_G,
                           **{f"{c}_mg": mg[c] for c in ORDER}})
    write_csv(os.path.join(IN, "scenarios_saturation_mixture.csv"), ms)
    return sc, ms


# ------------------------------------------------------------------ outputs
def export_single_curves():
    d = os.path.join(OUT, "01_single_compound_solubility")
    for c in ORDER:
        s = FLAVOURS[c]
        for o in MAIN_OILS:
            oil = OILS[o]
            m = UNIFAC([s, oil])
            rows = []
            for tc in range(25, 201):
                r = solubility(s, oil, float(to_kelvin(tc)), m)
                mg = binary_mg_per_100g(r.x, s, oil) if r.status != "miscible" else math.inf
                mg_id = binary_mg_per_100g(min(r.x_ideal, 0.999999), s, oil) if r.x_ideal < 1 else math.inf
                rows.append({"T_C": tc, "T_K": tc + 273.15, "status": r.status, "x_solubility": r.x,
                             "x_ideal": r.x_ideal, "gamma": r.gamma, "ln_gamma_combinatorial": r.ln_gamma_C,
                             "ln_gamma_residual": r.ln_gamma_R, "solubility_mg_per_100g_oil": mg,
                             "ideal_mg_per_100g_oil": mg_id,
                             "saturation_point_mg_in_15g_oil": mg * OIL_G / 100 if math.isfinite(mg) else math.inf})
            write_csv(os.path.join(d, f"solubility_{c}_in_{o}_oil_25-200C.csv"), rows)


def export_saturation_single(scenarios):
    d = os.path.join(OUT, "02_saturation_point_single_compound")
    rows, lines = [], []
    for sc in scenarios:
        r = single_saturation(sc["compound"], sc["oil"], sc["T_C"], sc["oil_g"], sc["added_mg"] / 1000)
        rows.append({"scenario": sc["scenario"], "compound": r["compound"], "oil": r["oil"], "T_C": r["T_C"],
                     "oil_g": r["oil_g"], "added_mg": r["added_mg"], "loading_mg_per_100g": r["loading_mg_per_100g"],
                     "status": r["status"], "solubility_mg_per_100g": r["solubility_mg_per_100g"],
                     "saturation_point_mg_in_pan": r["capacity_mg"], "dissolved_mg": r["dissolved_mg"],
                     "cannot_mix_in_mg": r["excess_mg"], "saturated": r["saturated"],
                     "saturation_temperature_C": r["saturation_T_C"],
                     "saturation_temperature_note": r["saturation_T_kind"], "statement": r["statement"]})
        lines.append(f"[{sc['scenario']}]\n{r['statement']}\n")
    write_csv(os.path.join(d, "saturation_single_scenarios.csv"), rows)
    write_txt(os.path.join(d, "saturation_single_statements.txt"),
              "Saturation statements, one compound in one oil (15 g oil).\n"
              "Generated by flavorsol.saturation.single_saturation.\n\n" + "\n".join(lines))
    # saturation temperature for a range of loadings
    tab = []
    for c in ORDER:
        for o in MAIN_OILS:
            s, oil = FLAVOURS[c], OILS[o]
            m = UNIFAC([s, oil])
            for L in (10, 30, 100, 300, 1000, 3000, 10000):
                T, kind = saturation_temperature(s, oil, L, m)
                tab.append({"compound": c, "oil": o, "loading_mg_per_100g_oil": L,
                            "loading_mg_in_15g_oil": L * OIL_G / 100, "saturation_temperature_C": T,
                            "note": {"root": "saturated below this T, fully dissolved above",
                                     "below_range": "fully dissolved at all T >= 25 C",
                                     "above_range": "still saturated at 200 C"}[kind]})
    write_csv(os.path.join(d, "saturation_temperature_vs_loading.csv"), tab)


def export_saturation_mixture(scenarios):
    d = os.path.join(OUT, "03_mixtures")
    rows, lines = [], []
    for sc in scenarios:
        added = {c: sc[f"{c}_mg"] / 1000 for c in ORDER if sc[f"{c}_mg"] > 0}
        r = mixture_saturation(sc["oil"], sc["T_C"], sc["oil_g"], added)
        for x in r["rows"]:
            rows.append({"scenario": sc["scenario"], "oil": sc["oil"], "T_C": sc["T_C"], "oil_g": sc["oil_g"],
                         "compound": x["compound"], "added_mg": x["added_mg"], "dissolved_mg": x["dissolved_mg"],
                         "cannot_mix_in_mg": x["excess_mg"],
                         "saturation_point_in_mixture_mg": x["capacity_in_mixture_mg"],
                         "gamma": x["gamma"], "saturated": x["saturated"],
                         "recipe_fully_dissolved_above_C": r["all_dissolved_T_C"],
                         "recipe_note": r["all_dissolved_kind"], "statement": x["statement"]})
        lines.append(f"[{sc['scenario']}]\n" + "\n".join(r["statements"]) + "\n")
    write_csv(os.path.join(d, "saturation_mixture_scenarios.csv"), rows)
    write_txt(os.path.join(d, "saturation_mixture_statements.txt"),
              "Saturation statements, several compounds in one oil (15 g oil).\n"
              "Generated by flavorsol.saturation.mixture_saturation.\n\n" + "\n".join(lines))


def export_routes():
    from flavorsol.engine import format_report, load_request, predict
    d = os.path.join(OUT, "05_routes_R1_to_R5")
    for p in sorted(glob.glob(os.path.join(HERE, "inputs", "*.json"))):
        name = os.path.basename(p).replace(".json", "")
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf):
                preds, issues = predict(load_request(p))
            text = format_report(preds, issues)
        except Exception as e:                          # bad inputs are meant to stop
            text = f"stopped with an input error (expected for {name}):\n{type(e).__name__}: {e}"
        write_txt(os.path.join(d, f"output_{name}.txt"),
                  f"python -m flavorsol predict inputs/{name}.json\n\n{buf.getvalue()}{text}")


def copy_results():
    res = os.path.join(HERE, "results")
    plan = {
        "03_mixtures": ["mixture_pairs_and_triple.csv", "mixture_cross_effects.csv", "recipe_extraction.csv",
                        "batch_example_sunflower_15g.csv"],
        "04_validation": ["validation_pure_solvents.csv", "validation_pure_solvents_summary.csv",
                          "validation_oils.csv", "validation_flavordb.csv"],
        "06_analysis": ["thermo_summary.csv", "uncertainty_check.csv", "solubility_vs_temperature.csv"],
    }
    for sub, files in plan.items():
        os.makedirs(os.path.join(OUT, sub), exist_ok=True)
        for f in files:
            src = os.path.join(res, f)
            if os.path.exists(src):
                shutil.copy(src, os.path.join(OUT, sub, f))
    sp = os.path.join(res, "summary.json")
    if os.path.exists(sp):
        with open(sp) as f:
            summ = json.load(f)
        write_txt(os.path.join(OUT, "06_analysis", "summary_headline_numbers.txt"), json.dumps(summ, indent=2))
        sens = summ.get("sensitivity_60C_sunflower", {})
        rows = []
        for comp, d in sens.items():
            base = d.get("base_mg_per_100g")
            for k, v in d.items():
                if k != "base_mg_per_100g":
                    rows.append({"compound": comp, "oil": "sunflower", "T_C": 60, "base_mg_per_100g": base,
                                 "perturbation": k, "change_percent": v})
        write_csv(os.path.join(OUT, "06_analysis", "sensitivity_60C_sunflower.csv"), rows)
    # validation report (text)
    from flavorsol import validation
    lines = []
    b = validation.unifac_benchmark()
    lines.append(f"UNIFAC benchmark (n-hexane/2-butanone, 333.15 K) max abs error: {b['max_abs_error']:.1e}")
    _, summ, _ = validation.validate_pure_solvents()
    for s in summ:
        lines.append(f"{s['solute']:17s} in {s['solvent']:15s} n={s['n']:2d}  AARD UNIFAC {s['AARD_unifac']:6.1f} %"
                     f"   ideal {s['AARD_ideal']:6.1f} %")
    for r in validation.validate_oils():
        lines.append(f"{r['solute']:9s} {r['oil']:9s} {r['T'] - 273.15:5.0f} C  published {r['exp_mg_mL']:6.2f}"
                     f"  predicted {r['pred_mg_mL']:7.3f}  ideal {r['ideal_mg_mL']:7.2f} mg/mL  ({r['kind']})")
    write_txt(os.path.join(OUT, "04_validation", "validation_report.txt"), "\n".join(lines))


def thermo_checks():
    """Second-law and van't Hoff checks on every predicted point (written to 06_analysis)."""
    rows = []
    for c in ORDER:
        for o in MAIN_OILS:
            s, oil = FLAVOURS[c], OILS[o]
            m = UNIFAC([s, oil])
            xs, ok_ideal, ok_mono, gmin = [], True, True, 1e9
            for tc in range(25, 201):
                r = solubility(s, oil, float(to_kelvin(tc)), m)
                if r.status == "SLE":
                    ok_ideal &= r.x <= r.x_ideal * (1 + 1e-9)
                    gmin = min(gmin, r.gamma)
                xs.append(r.x)
            ok_mono = all(b >= a - 1e-15 for a, b in zip(xs, xs[1:]))
            pure = m.ln_gamma([1.0, 1e-12], 330.0)[0]
            rows.append({"compound": c, "oil": o, "min_gamma_at_saturation_SLE": gmin,
                         "x_below_ideal_everywhere": ok_ideal, "solubility_increases_with_T": ok_mono,
                         "ln_gamma_pure_solute_limit": pure})
    write_csv(os.path.join(OUT, "06_analysis", "thermodynamic_consistency_checks.csv"), rows)
    return rows


def main():
    os.makedirs(DATA, exist_ok=True)
    sc, ms = export_inputs()
    export_single_curves()
    export_saturation_single(sc)
    export_saturation_mixture(ms)
    export_routes()
    copy_results()
    rows = thermo_checks()
    print("Data folder written to", DATA)
    for r in rows:
        print(r)


if __name__ == "__main__":
    main()
