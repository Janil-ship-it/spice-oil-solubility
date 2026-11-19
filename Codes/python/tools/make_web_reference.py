"""
Reference results from the Python model, used by the website test
(Codes/website/tests/engine-vs-python.test.ts) to prove that the browser
solver gives the same numbers.   Run from Codes/python:
    python tools/make_web_reference.py
"""
import json
import math
import os
import sys
import warnings

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from flavorsol.data.compounds import FLAVOURS                      # noqa: E402
from flavorsol.data.oils import OILS                               # noqa: E402
from flavorsol.mixture import solve_mixture                        # noqa: E402
from flavorsol.saturation import mixture_saturation, single_saturation   # noqa: E402
from flavorsol.units import to_kelvin                              # noqa: E402

warnings.filterwarnings("ignore", message="UNIFAC: no published parameters")
OUT = os.path.join(os.path.dirname(__file__), "..", "..", "website", "tests", "python_reference.json")


def clean(v):
    if isinstance(v, float) and not math.isfinite(v):
        return None
    return v


def main():
    mix = []
    order = ["vanillin", "curcumin", "piperine"]
    cases = []
    for oil in ("sunflower", "coconut"):
        for T in (25, 40, 60, 75):
            cases.append((oil, T, order, None))                                  # all in excess
            cases.append((oil, T, ["vanillin", "curcumin"], None))
            cases.append((oil, T, ["curcumin", "piperine"], None))
        for T in (25, 60, 100, 120, 140, 160):
            cases.append((oil, T, order, [100.0, 314.0, 333.0]))                 # tadka recipe
            cases.append((oil, T, ["curcumin", "piperine"], [313.3, 333.3]))
    for oil, T, names, loads in cases:
        r = solve_mixture([FLAVOURS[n] for n in names], OILS[oil], float(to_kelvin(T)), loads)
        mix.append({"oil": oil, "T_C": T, "names": names, "loads": loads,
                    "mg": [float(v) for v in r.mg_per_100g], "saturated": r.saturated})
    single = []
    for oil in ("sunflower", "coconut"):
        for c in order:
            for T, g in ((25, 0.05), (60, 0.047), (100, 0.5), (120, 0.047), (150, 1.0)):
                s = single_saturation(c, oil, T, 15.0, g)
                single.append({k: clean(s[k]) for k in ("compound", "oil", "T_C", "oil_g", "added_mg",
                                                        "capacity_mg", "saturated", "saturation_T_C",
                                                        "saturation_T_kind")})
    msat = []
    for oil, T, added in (("coconut", 60, {"curcumin": 0.047, "piperine": 0.05}),
                          ("sunflower", 60, {"curcumin": 0.0471, "piperine": 0.05, "vanillin": 0.015}),
                          ("sunflower", 100, {"curcumin": 0.0471, "piperine": 0.05, "vanillin": 0.015}),
                          ("coconut", 140, {"curcumin": 0.0471, "piperine": 0.05, "vanillin": 0.015})):
        r = mixture_saturation(oil, T, 15.0, added)
        msat.append({"oil": oil, "T_C": T, "added_mg": {k: 1000 * v for k, v in added.items()},
                     "rows": [{"compound": x["compound"], "dissolved_mg": x["dissolved_mg"],
                               "saturated": x["saturated"]} for x in r["rows"]],
                     "all_dissolved_T_C": r["all_dissolved_T_C"], "kind": r["all_dissolved_kind"]})
    with open(OUT, "w") as f:
        json.dump({"mixture": mix, "single": single, "mixture_saturation": msat}, f, indent=1)
    print("wrote", os.path.normpath(OUT))


if __name__ == "__main__":
    main()
