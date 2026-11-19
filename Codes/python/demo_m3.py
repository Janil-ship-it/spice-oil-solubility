"""
Milestone 3 - proof-of-concept demo (run this live in the meeting).

    python demo_m3.py
    python demo_m3.py --solute curcumin --oil coconut --T 100

It walks through the calculation in the same order as our flowchart,
printing what each team member's module does.
"""
import argparse
import warnings

import numpy as np

from flavorsol.data.compounds import FLAVOURS
from flavorsol.data.oils import OILS
from flavorsol.mixture import solve_mixture
from flavorsol.sle import ln_a_sat, solubility
from flavorsol.unifac import UNIFAC
from flavorsol.units import binary_mg_per_100g, to_kelvin
from flavorsol.validation import log_kow, unifac_benchmark, validate_lab, validate_oils

warnings.filterwarnings("ignore", message="UNIFAC: no published parameters")


def banner(step, who, title):
    print(f"\n{'=' * 72}\nSTEP {step}  [{who}]  {title}\n{'=' * 72}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--solute", default="curcumin", choices=sorted(FLAVOURS))
    ap.add_argument("--oil", default="sunflower", choices=["sunflower", "coconut"])
    ap.add_argument("--T", type=float, default=100.0, help="temperature in deg C")
    a = ap.parse_args()
    s, oil, T = FLAVOURS[a.solute], OILS[a.oil], float(to_kelvin(a.T))

    banner(1, "Member 1 - data", "Inputs")
    print(f"Solute : {s.name} ({s.formula}), M = {s.M} g/mol, Tm = {s.Tm:.2f} K, "
          f"dHfus = {s.dHfus / 1e3:.2f} kJ/mol")
    print(f"         UNIFAC groups: {s.groups}")
    print(f"Oil    : {oil.name}, M = {oil.M:.1f} g/mol (average triglyceride)")
    print("         UNIFAC groups: " + ", ".join(f"{g} {n:.2f}" for g, n in oil.groups.items()))

    banner(2, "Member 2 - UNIFAC", "Check the activity-coefficient code against DDBST")
    b = unifac_benchmark()
    print(f"{b['system']}: gamma = {np.round(b['gamma_calc'], 6)}  reference = "
          f"{np.round(b['gamma_ref'], 6)}  max error = {b['max_abs_error']:.1e}")

    banner(3, "Member 2 - SLE", f"Ideal solubility at {a.T:.0f} C")
    lnas = float(ln_a_sat(s, T))
    print(f"ln a_sat = -(dHfus/R)(1/T - 1/Tm) = {lnas:.4f}  ->  x_ideal = {np.exp(lnas):.4e}")

    banner(4, "Member 2 - single-solute solver", "Brent's method on ln x + ln gamma(x) = ln a_sat")
    r = solubility(s, oil, T, keep_history=True)
    print(f"{'call':>4} {'x':>14} {'residual':>14}")
    for k, (x, f) in enumerate(r.history, 1):
        print(f"{k:4d} {x:14.6e} {f:14.3e}")
    print(f"status = {r.status}; x = {r.x:.5e}; gamma = {r.gamma:.3f} "
          f"(ln gC = {r.ln_gamma_C:.3f}, ln gR = {r.ln_gamma_R:.3f}); "
          f"{r.iterations} iterations")
    print(f"=> {binary_mg_per_100g(r.x, s, oil):.1f} mg per 100 g {oil.name}  "
          f"(ideal solution would give {binary_mg_per_100g(min(r.x_ideal, .999999), s, oil):.0f})")

    banner(5, "Member 3 - mixture solver", "Vanillin + curcumin + piperine together, Newton-Raphson")
    T_mix = float(to_kelvin(60.0))
    sols = [FLAVOURS["vanillin"], FLAVOURS["curcumin"], FLAVOURS["piperine"]]
    mix = solve_mixture(sols, oil, T_mix)
    print("All three in excess at 60 C. max|F| per Newton iteration:",
          ", ".join(f"{v:.1e}" for v in mix.residual_history))
    for i, sn in enumerate(mix.solutes):
        alone = binary_mg_per_100g(solubility(sols[i], oil, T_mix).x, sols[i], oil)
        print(f"  {sn:9s} {mix.mg_per_100g[i]:9.1f} mg/100 g in mixture vs {alone:9.1f} alone")

    banner(6, "Member 4 - validation", "Predictions vs published data")
    print(f"vanillin log Kow: UNIFAC {log_kow(FLAVOURS['vanillin']):.2f}, measured 1.23")
    for v in validate_oils():
        if v["kind"] == "oil":
            print(f"  {v['solute']:9s} {v['oil']:8s} {v['T'] - 273.15:4.0f} C: published "
                  f"{v['exp_mg_mL']:5.2f} mg/mL, predicted {v['pred_mg_mL']:6.3f} mg/mL "
                  f"(x{v['ratio']:.2f})")
    print("Lab data:", "none yet - fill data/lab_measurements.csv" if validate_lab() is None
          else f"{len(validate_lab())} points compared")

    banner(7, "Member 5 - analysis", "Solubility vs temperature (mg per 100 g oil)")
    print(f"{'T (C)':>6} " + " ".join(f"{n:>12s}" for n in FLAVOURS))
    for tc in (25, 60, 100, 140, 180):
        cells = []
        for sv in FLAVOURS.values():
            rr = solubility(sv, oil, float(to_kelvin(tc)))
            v = binary_mg_per_100g(rr.x, sv, oil)
            cells.append("miscible" if rr.status == "miscible" else
                         f"{v:.3g}{'*' if rr.status == 'LLE-limited' else ''}")
        print(f"{tc:6d} " + " ".join(f"{c:>12s}" for c in cells))
    print("* molten solute forms a second liquid phase; value = oil-rich phase limit")

    banner(8, "All - engine", "Decision pathways: optional inputs change the route")
    from flavorsol.engine import format_report, load_request, predict
    for f in ("inputs/03_anchor_same_oil.json", "inputs/07_bad_inputs.json"):
        print(f"\n--- {f}")
        req = load_request(f)
        req.monte_carlo = 0
        print(format_report(*predict(req)))


if __name__ == "__main__":
    main()
