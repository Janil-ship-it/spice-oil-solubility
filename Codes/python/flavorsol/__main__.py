"""
Command-line interface.

    python -m flavorsol solubility curcumin sunflower --T 25 60 100 140
    python -m flavorsol mixture --oil coconut --T 60 --loadings 100 314 333
    python -m flavorsol validate
    python -m flavorsol database
    python -m flavorsol predict inputs/03_anchor_same_oil.json
    python -m flavorsol batch --oil sunflower --T 100 --oil-g 15 --curcumin 0.047 --piperine 0.05
    python -m flavorsol saturation curcumin coconut --T 60 --oil-g 15 --added-g 0.047
    python -m flavorsol mixsat --oil coconut --T 60 --oil-g 15 --curcumin 0.047 --piperine 0.05
"""
import argparse
import warnings

from .data.compounds import FLAVOURS
from .data.oils import OILS
from .mixture import solve_mixture
from .sle import solubility
from .units import binary_mg_per_100g, to_kelvin
from . import validation

warnings.filterwarnings("ignore", message="UNIFAC: no published parameters")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="flavorsol",
                                 description="Solubility of spice flavour compounds in cooking oils")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p1 = sub.add_parser("solubility", help="one compound in one oil")
    p1.add_argument("solute", choices=sorted(FLAVOURS))
    p1.add_argument("oil", choices=sorted(OILS))
    p1.add_argument("--T", type=float, nargs="+", default=[25, 60, 100], help="deg C")
    p2 = sub.add_parser("mixture", help="vanillin + curcumin + piperine together")
    p2.add_argument("--oil", default="sunflower", choices=sorted(OILS))
    p2.add_argument("--T", type=float, default=60.0, help="deg C")
    p2.add_argument("--loadings", type=float, nargs=3, metavar=("VAN", "CUR", "PIP"),
                    help="mg added per 100 g oil (omit = all in excess)")
    sub.add_parser("validate", help="run the validation checks")
    sub.add_parser("database", help="what FlavorDB/SpiceRx (IIIT-Delhi) records are loaded")
    p4 = sub.add_parser("predict", help="full workflow from a JSON input file (see inputs/)")
    p4.add_argument("json_file")
    p5 = sub.add_parser("batch", help="grams of oil + grams of each compound -> dissolved and undissolved grams")
    p5.add_argument("--oil", default="sunflower", choices=sorted(OILS))
    p5.add_argument("--T", type=float, default=100.0, help="deg C")
    p5.add_argument("--oil-g", type=float, required=True)
    for n in ("curcumin", "vanillin", "piperine"):
        p5.add_argument(f"--{n}", type=float, default=None, help=f"grams of {n} added")
    p6 = sub.add_parser("saturation", help="is the oil saturated at T? prints the saturation statement")
    p6.add_argument("solute", choices=sorted(FLAVOURS))
    p6.add_argument("oil", choices=sorted(OILS))
    p6.add_argument("--T", type=float, required=True, help="deg C")
    p6.add_argument("--oil-g", type=float, default=15.0, help="grams of oil")
    p6.add_argument("--added-g", type=float, required=True, help="grams of compound added")
    p7 = sub.add_parser("mixsat", help="saturation statements for several compounds in one oil")
    p7.add_argument("--oil", default="coconut", choices=sorted(OILS))
    p7.add_argument("--T", type=float, required=True, help="deg C")
    p7.add_argument("--oil-g", type=float, default=15.0)
    for n in ("curcumin", "vanillin", "piperine"):
        p7.add_argument(f"--{n}", type=float, default=0.0, help=f"grams of {n} added")
    a = ap.parse_args(argv)

    if a.cmd == "saturation":
        from .saturation import single_saturation
        r = single_saturation(a.solute, a.oil, a.T, a.oil_g, a.added_g)
        print(r["statement"])
        return
    if a.cmd == "mixsat":
        from .saturation import mixture_saturation
        added = {n: getattr(a, n) for n in ("vanillin", "curcumin", "piperine") if getattr(a, n) > 0}
        if not added:
            ap.error("give at least one of --curcumin, --vanillin, --piperine (grams)")
        r = mixture_saturation(a.oil, a.T, a.oil_g, added)
        print("\n".join(r["statements"]))
        return

    if a.cmd == "solubility":
        s, oil = FLAVOURS[a.solute], OILS[a.oil]
        print(f"{'T (C)':>6} {'status':>12} {'x':>11} {'gamma':>8} {'mg/100 g oil':>14}")
        for tc in a.T:
            if tc < -50 or tc > 300:
                print(f"{tc:6.1f}  skipped: outside -50..300 C")
                continue
            r = solubility(s, oil, float(to_kelvin(tc)))
            mg = binary_mg_per_100g(r.x, s, oil)
            print(f"{tc:6.1f} {r.status:>12} {r.x:11.4e} {r.gamma:8.3f} {mg:14.4g}")
    elif a.cmd == "mixture":
        sols = [FLAVOURS[n] for n in ("vanillin", "curcumin", "piperine")]
        mix = solve_mixture(sols, OILS[a.oil], float(to_kelvin(a.T)), a.loadings)
        for i, n in enumerate(mix.solutes):
            state = "saturated (solid left)" if mix.saturated[i] else "fully dissolved"
            print(f"{n:9s} {mix.mg_per_100g[i]:10.2f} mg/100 g  gamma = {mix.gamma[i]:7.3f}  {state}")
        for note in mix.notes:
            print("note:", note)
        print(f"Newton iterations: {mix.newton_iterations}")
    elif a.cmd == "batch":
        from .operating import batch
        added = {n: getattr(a, n) for n in ("curcumin", "vanillin", "piperine") if getattr(a, n)}
        if not added:
            ap.error("give at least one of --curcumin, --vanillin, --piperine (grams)")
        print(f"{'compound':9s} {'added g':>9s} {'dissolved g':>12s} {'left solid g':>13s} "
              f"{'% dissolved':>12s} {'ppm':>9s} {'mg/100 g':>9s}  limited by")
        for r in batch(a.oil, a.T, a.oil_g, added):
            print(f"{r['compound']:9s} {r['available_g']:9.4f} {r['dissolved_g']:12.4f} "
                  f"{r['undissolved_g']:13.4f} {r['dissolved_percent']:12.1f} {r['ppm_w']:9.0f} "
                  f"{r['mg_per_100g_oil']:9.1f}  {r['limited_by']}")
    elif a.cmd == "predict":
        from .engine import format_report, load_request, predict
        preds, issues = predict(load_request(a.json_file))
        print(format_report(preds, issues))
    elif a.cmd == "database":
        from .data.spices import database_status
        print(database_status())
        print("Refresh from a machine with internet access: python tools/fetch_flavordb.py")
    else:
        b = validation.unifac_benchmark()
        print("UNIFAC benchmark max error:", f"{b['max_abs_error']:.1e}")
        _, summ, ov = validation.validate_pure_solvents()
        for srow in summ:
            print(f"{srow['solute']:17s} in {srow['solvent']:15s} AARD UNIFAC "
                  f"{srow['AARD_unifac']:6.1f} %   ideal {srow['AARD_ideal']:6.1f} %")
        for r in validation.validate_oils():
            print(f"{r['solute']:9s} {r['oil']:9s} {r['T'] - 273.15:5.0f} C  published "
                  f"{r['exp_mg_mL']:6.2f}  predicted {r['pred_mg_mL']:7.3f} mg/mL  ({r['kind']})")


if __name__ == "__main__":
    main()
