"""
Reproduce every number, table and figure in the report.

    python run_all.py            (about 5 minutes; also fills ../../Data via export_data.py)

Outputs
    results/*.csv            raw tables (open in Excel)
    results/summary.json     headline numbers
    results/tables/*.tex     LaTeX tables used by the manuscript
    figures/*.png, *.pdf     figures
"""
import csv
import json
import os
import time
import warnings

import numpy as np

from flavorsol import analysis, plotting, validation
from flavorsol.data.compounds import CURCUMIN, FLAVOURS, PIPERINE, VANILLIN
from flavorsol.data.oils import OILS
from flavorsol.data.sensory import THRESHOLDS_MG_PER_KG, recipe_loading_mg_per_100g
from flavorsol.mixture import solve_mixture
from flavorsol.sle import solubility
from flavorsol.unifac import UNIFAC
from flavorsol.units import binary_mg_per_100g, to_kelvin

warnings.filterwarnings("ignore", message="UNIFAC: no published parameters")
ROOT = os.path.dirname(os.path.abspath(__file__))
RES, FIG, TAB = (os.path.join(ROOT, d) for d in ("results", "figures", "results/tables"))
for d in (RES, FIG, TAB):
    os.makedirs(d, exist_ok=True)


def write_csv(name, rows):
    if not rows:
        return
    with open(os.path.join(RES, name), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow({k: (float(v) if isinstance(v, np.floating) else v) for k, v in r.items()})


def fmt(v, sig=3):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return "--"
    if abs(v) >= 100:                      # 3 significant figures, written as an integer
        return f"{float(f'{v:.{sig}g}'):.0f}"
    return f"{v:.{sig}g}"


def main():
    t0 = time.time()
    summary = {}

    # ---------------- Level 1: UNIFAC benchmark
    summary["unifac_benchmark"] = validation.unifac_benchmark()

    # ---------------- Level 2-3: analogues, water, log Kow
    points, per_system, overall = validation.validate_pure_solvents()
    write_csv("validation_pure_solvents.csv", points)
    write_csv("validation_pure_solvents_summary.csv",
              [{**s, "T_range": f"{s['T_range'][0]:.2f}-{s['T_range'][1]:.2f}"} for s in per_system])
    aprotic = [p for p in points if p["solvent"] not in ("ethanol", "1-propanol", "2-propanol")]
    alcohol = [p for p in points if p["solvent"] in ("ethanol", "1-propanol", "2-propanol")]
    summary["pure_solvent_validation"] = {
        "overall": overall,
        "aprotic": dict(n=len(aprotic),
                        AARD_unifac=validation.aard([p["x_calc"] for p in aprotic], [p["x_exp"] for p in aprotic]),
                        AARD_ideal=validation.aard([p["x_ideal"] for p in aprotic], [p["x_exp"] for p in aprotic])),
        "alcohols": dict(n=len(alcohol),
                         AARD_unifac=validation.aard([p["x_calc"] for p in alcohol], [p["x_exp"] for p in alcohol]),
                         AARD_ideal=validation.aard([p["x_ideal"] for p in alcohol], [p["x_exp"] for p in alcohol])),
    }
    summary["vanillin_water"] = validation.vanillin_water()
    dbc = validation.database_check()
    write_csv("validation_flavordb.csv", dbc)
    summary["flavordb_check"] = dbc
    summary["log_kow_vanillin"] = {"calc": validation.log_kow(VANILLIN), "exp": 1.23,
                                   "source": "OECD SIDS vanillin (measured, 22 C)"}

    # ---------------- Level 4: oils
    oil_rows = validation.validate_oils()
    write_csv("validation_oils.csv", oil_rows)
    eq = [r for r in oil_rows if r["kind"] == "oil"]
    logdev = [abs(np.log10(r["ratio"])) for r in eq]
    summary["oil_validation_metrics"] = validation.metrics([r["pred_mg_mL"] for r in eq],
                                                           [r["exp_mg_mL"] for r in eq])
    summary["oil_validation"] = {
        "n_equilibrium_points": len(eq),
        "mean_abs_log10_dev": float(np.mean(logdev)),
        "max_factor": float(10 ** np.max(logdev)),
        "AARD_unifac": validation.aard([r["pred_mg_mL"] for r in eq], [r["exp_mg_mL"] for r in eq]),
        "AARD_ideal": validation.aard([r["ideal_mg_mL"] for r in eq], [r["exp_mg_mL"] for r in eq]),
    }

    # ---------------- Level 5: lab data (if the group has filled it in)
    lab = validation.validate_lab()
    summary["lab_validation"] = "no lab file found (data/lab_measurements.csv)" if lab is None else lab
    if lab:
        write_csv("validation_lab.csv", lab)
        with open(os.path.join(TAB, "lab.tex"), "w") as f:
            for r in lab:
                f.write(f"{r['solute']} & {r['oil']} & {r['T_C']:.0f} & {fmt(r['exp'])} & "
                        f"{fmt(r['pred'])} & {r['rel_dev']:+.0f} \\\\\n")

    # ---------------- predictions
    rows = analysis.sweep()
    write_csv("solubility_vs_temperature.csv", rows)

    table_T = [25, 60, 100, 140, 180]
    pred_table = {}
    for s in FLAVOURS:
        for o in analysis.MAIN_OILS:
            pred_table[(s, o)] = {r["T_C"]: r for r in rows
                                  if r["solute"] == s and r["oil"] == o and r["T_C"] in table_T}

    thermo = []
    for s in FLAVOURS.values():
        for o in analysis.MAIN_OILS:
            oil = OILS[o]
            m = UNIFAC([s, oil])
            lnC, lnR = m.ln_gamma_parts([1e-9, 1 - 1e-9], 298.15)
            out_T = analysis.oiling_out_temperature(s, oil)
            thermo.append(dict(solute=s.name, oil=o, lnC=float(lnC[0]), lnR=float(lnR[0]),
                               gamma_inf_25C=float(np.exp(lnC[0] + lnR[0])),
                               dH_sol_kJ=analysis.vant_hoff_enthalpy(s, oil) / 1e3,
                               dHfus_kJ=s.dHfus / 1e3,
                               oiling_out_C=None if out_T is None else out_T - 273.15))
    write_csv("thermo_summary.csv", thermo)

    mm = analysis.mixture_matrix()
    write_csv("mixture_pairs_and_triple.csv", mm)
    with open(os.path.join(TAB, "mixtures.tex"), "w") as f:
        keyed = {}
        for r in mm:
            keyed.setdefault((r["mixture"], r["solute"]), {})[(r["oil"], r["T_C"])] = r["ratio"]
        for (mix_name, sol), v in keyed.items():
            short = mix_name.replace("vanillin", "V").replace("curcumin", "C").replace("piperine", "P")
            cells = [f"{v[(o, t)]:.2f}" for o in analysis.MAIN_OILS for t in (25, 60, 75)]
            f.write(f"{short} & {sol} & " + " & ".join(cells) + r" \\" + "\n")
    summary["mixture_matrix_max_ratio"] = {
        f"{r['oil']} | {r['mixture']} | {r['solute']}": r["ratio"] for r in mm if r["T_C"] == 60}
    from flavorsol.operating import batch
    bt = [r for tc in (25, 60, 100, 140) for r in
          batch("sunflower", tc, 15.0, {"curcumin": 0.047, "piperine": 0.05, "vanillin": 0.015})]
    write_csv("batch_example_sunflower_15g.csv", bt)

    cross = analysis.cross_effects()
    write_csv("mixture_cross_effects.csv", cross)
    recipe = analysis.recipe_extraction()
    write_csv("recipe_extraction.csv", recipe)

    def full_T(oil):
        for r in sorted([r for r in recipe if r["solute"] == "curcumin" and r["oil"] == oil],
                        key=lambda r: r["T_C"]):
            if r["extracted_pct"] >= 99.999:
                return r["T_C"]
    summary["recipe"] = {
        "loadings_mg_per_100g": {s: recipe_loading_mg_per_100g(s) for s in FLAVOURS},
        "curcumin_full_extraction_T_C": {o: full_T(o) for o in analysis.MAIN_OILS},
        "curcumin_extracted_pct": {o: {int(r["T_C"]): r["extracted_pct"] for r in recipe
                                       if r["solute"] == "curcumin" and r["oil"] == o
                                       and r["T_C"] in (25, 60, 100)} for o in analysis.MAIN_OILS},
    }

    sens = {s.name: analysis.sensitivity(s, "sunflower", 60.0) for s in FLAVOURS.values()}
    summary["sensitivity_60C_sunflower"] = {k: {"base_mg_per_100g": v[0], **v[1]} for k, v in sens.items()}

    # convergence histories for the demo figure
    sun = OILS["sunflower"]
    b1 = solubility(CURCUMIN, sun, to_kelvin(100), keep_history=True)
    b2 = solubility(VANILLIN, sun, to_kelvin(25), keep_history=True)
    n1 = solve_mixture([VANILLIN, CURCUMIN, PIPERINE], sun, float(to_kelvin(60)))
    n2 = solve_mixture([VANILLIN, CURCUMIN, PIPERINE], OILS["coconut"], float(to_kelvin(60)))
    summary["solver_stats"] = {
        "brent_curcumin_sunflower_100C": dict(iterations=b1.iterations, calls=b1.function_calls),
        "brent_vanillin_sunflower_25C": dict(iterations=b2.iterations, calls=b2.function_calls),
        "newton_mix_sunflower_60C": n1.newton_iterations,
        "newton_mix_coconut_60C": n2.newton_iterations,
    }

    # ---------------- data-anchored routes, extrapolation and uncertainty (meeting feedback)
    from flavorsol.engine import Measurement, Request, predict
    ext_T, ext_S = validation.extrapolation_validation()
    summary["route_validation"] = {"temperature_extrapolation": ext_T, "solvent_transfer": ext_S}
    unc = validation.uncertainty_check(1000)
    summary["uncertainty_check"] = unc
    write_csv("uncertainty_check.csv", [{k: v for k, v in r.items() if k != "share"} for r in unc])
    t_grid = list(np.arange(25, 141, 5.0))
    r5 = predict(Request("curcumin", "corn", t_grid))[0]
    r2 = predict(Request("curcumin", "corn", t_grid,
                         measurements=[Measurement(0.07, "mg/mL", 25, sigma_rel=0.10)]))[0]
    corn = OILS["corn"]
    lit_cc = {"anchor": (25, 100 * 0.07 / corn.rho20),
              "hot": (100, 100 * 3.0 / corn.rho20, 100 * 4.0 / corn.rho20)}
    p100 = {lab: next(p for p in pr if p.T_C == 100) for lab, pr in (("R5", r5), ("R2", r2))}
    summary["curcumin_corn_100C"] = {k: {"mg_per_100g": v.mg_per_100g, "ci95": v.ci95_mg_per_100g}
                                     for k, v in p100.items()}
    summary["curcumin_corn_100C"]["published_mg_per_100g"] = lit_cc["hot"][1:]
    for meth in ("regular", "constant"):
        pr = predict(Request("curcumin", "corn", [100.0], T_extrapolation=meth,
                             measurements=[Measurement(0.07, "mg/mL", 25, sigma_rel=0.10)]))[0][0]
        summary["curcumin_corn_100C"]["R2_" + meth] = pr.mg_per_100g

    # ---------------- figures
    tm_c = {s.name: s.Tm - 273.15 for s in FLAVOURS.values()}
    lit = {"curcumin": [dict(T_C=25, mg_per_100g=100 * 0.10 / OILS["coconut"].rho20,
                             label="coconut, lit.")]}
    plotting.fig_validation(points, oil_rows, FIG)
    plotting.fig_solubility_T(rows, tm_c, lit, FIG)
    plotting.fig_gamma(thermo, rows, FIG)
    plotting.fig_recipe(recipe, cross, FIG)
    plotting.fig_sensitivity(sens, 60, FIG)
    plotting.fig_convergence(
        [("curcumin, sunflower, 100 °C", [h[1] for h in b1.history], plotting.COLORS["curcumin"]),
         ("vanillin, sunflower, 25 °C", [h[1] for h in b2.history], plotting.COLORS["vanillin"])],
        [("sunflower oil, 60 °C", n1.residual_history, "-"),
         ("coconut oil, 60 °C", n2.residual_history, plotting.OIL_STYLE["coconut"])], FIG)
    plotting.fig_extrapolation(r5, r2, lit_cc, ext_T, ext_S, FIG)

    # ---------------- saturation point (flavorsol.saturation)
    from flavorsol.saturation import saturation_temperature
    sat_curve = [dict(oil=r["oil"], T_C=r["T_C"], mg=r["mg_per_100g"]) for r in rows
                 if r["solute"] == "curcumin" and np.isfinite(r["mg_per_100g"])]
    L_tadka = 47.0 * 100 / 15.0
    tsat_example = {o: saturation_temperature(CURCUMIN, OILS[o], L_tadka)[0] for o in analysis.MAIN_OILS}
    tsat_rows = []
    for s in FLAVOURS.values():
        for o in analysis.MAIN_OILS:
            m = UNIFAC([s, OILS[o]])
            for L in np.logspace(1, 4.5, 15):
                T_sat, kind = saturation_temperature(s, OILS[o], float(L), m)
                tsat_rows.append(dict(solute=s.name, oil=o, loading=float(L), T_sat=T_sat, kind=kind))
    write_csv("saturation_temperature_vs_loading.csv", tsat_rows)
    summary["saturation_tadka_curcumin_T_C"] = tsat_example
    plotting.fig_saturation(sat_curve, tsat_rows, {"loading": L_tadka, "T_sat": tsat_example}, FIG)

    # ---------------- LaTeX tables
    with open(os.path.join(TAB, "predictions.tex"), "w") as f:
        for (s, o), d in pred_table.items():
            cells = []
            for t in table_T:
                r = d[t]
                if r["status"] == "miscible":
                    cells.append("miscible")
                elif r["status"] == "LLE-limited":
                    cells.append(fmt(r["mg_per_100g"]) + "$^{b}$")
                else:
                    cells.append(fmt(r["mg_per_100g"]))
            f.write(f"{s} & {o} & " + " & ".join(cells) + r" \\" + "\n")
    with open(os.path.join(TAB, "validation_solvents.tex"), "w") as f:
        for srow in per_system:
            ab = {"ethyl vanillin": "EV", "vanillyl alcohol": "VA"}.get(srow["solute"], srow["solute"])
            f.write(f"{ab} & {srow['solvent']} & {srow['n']} & "
                    f"{srow['T_range'][0]:.0f}--{srow['T_range'][1]:.0f} & "
                    f"{srow['AARD_unifac']:.0f} & {srow['AARD_ideal']:.0f} \\\\\n")
    with open(os.path.join(TAB, "validation_oils.tex"), "w") as f:
        for r in oil_rows:
            kind = {"oil": "", "oil-lower-bound": "$^{c}$", "oil-high-cluster": "$^{d}$"}[r["kind"]]
            f.write(f"{r['solute']} & {r['oil']} & {r['T'] - 273.15:.0f} & "
                    f"{r['exp_mg_mL']:g}{kind} & {fmt(r['pred_mg_mL'], 2)} & "
                    f"{fmt(r['ideal_mg_mL'], 2)} \\\\\n")
    with open(os.path.join(TAB, "thermo.tex"), "w") as f:
        for r in thermo:
            oo = "--" if r["oiling_out_C"] is None else f"{r['oiling_out_C']:.0f}"
            f.write(f"{r['solute']} & {r['oil']} & {r['lnC']:.2f} & {r['lnR']:.2f} & "
                    f"{r['gamma_inf_25C']:.1f} & {r['dH_sol_kJ']:.1f} & {oo} \\\\\n")

    summary["runtime_s"] = time.time() - t0
    with open(os.path.join(RES, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2, default=lambda o: float(o) if isinstance(o, np.floating) else str(o))
    print(json.dumps({k: summary[k] for k in ("unifac_benchmark", "pure_solvent_validation",
                                              "log_kow_vanillin", "oil_validation", "recipe",
                                              "solver_stats", "route_validation",
                                              "curcumin_corn_100C", "runtime_s")},
                     indent=1, default=float))


if __name__ == "__main__":
    main()
    # fill the submission Data folder (inputs and outputs of every scenario)
    import export_data
    export_data.main()
