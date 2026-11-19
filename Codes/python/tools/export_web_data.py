"""
Export the model data used by the website (Codes/website/src/data/model.json).

The website runs without a server, so it needs:
  * the single-compound solubility curves of the Python model
    (sunflower and coconut oil, 25-200 C in 1 C steps): mg/100 g, x, gamma, status
  * the UNIFAC data (subgroups R, Q, main group; interaction parameters) and the
    compound / oil group counts, so the browser can solve the MIXTURE equations
    with the same Newton-Raphson + active-set method as flavorsol/mixture.py

Run from Codes/python:
    python tools/export_web_data.py
"""
import json
import math
import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from flavorsol.data.compounds import FLAVOURS                      # noqa: E402
from flavorsol.data.oils import OILS                               # noqa: E402
from flavorsol.data.unifac_params import SUBGROUPS, interaction     # noqa: E402
from flavorsol.sle import solubility                               # noqa: E402
from flavorsol.unifac import UNIFAC                                # noqa: E402
from flavorsol.units import binary_mg_per_100g, to_kelvin          # noqa: E402

warnings.filterwarnings("ignore", message="UNIFAC: no published parameters")
OUT = os.path.join(os.path.dirname(__file__), "..", "..", "website", "src", "data", "model.json")
WEB_OILS = ("sunflower", "coconut")
T_GRID = list(range(25, 201))


def r6(v):
    return None if v is None or not math.isfinite(v) else float(f"{v:.7g}")


def main():
    comps = {}
    for k, c in FLAVOURS.items():
        comps[k] = {"M": c.M, "Tm": c.Tm, "dHfus": c.dHfus, "dCp": c.dCp,
                    "formula": c.formula, "groups": c.groups}
    oils, curves = {}, {}
    for o in WEB_OILS:
        oil = OILS[o]
        oils[o] = {"M": oil.M, "groups": oil.groups, "composition": oil.composition}
        curves[o] = {}
        for k, s in FLAVOURS.items():
            m = UNIFAC([s, oil])
            mg, x, g, st = [], [], [], []
            for tc in T_GRID:
                r = solubility(s, oil, float(to_kelvin(tc)), m)
                mg.append(r6(binary_mg_per_100g(r.x, s, oil)) if r.status != "miscible" else None)
                x.append(r6(r.x))
                g.append(r6(r.gamma))
                st.append({"SLE": "S", "LLE-limited": "L", "miscible": "M"}[r.status])
            curves[o][k] = {"mg": mg, "x": x, "gamma": g, "status": "".join(st)}
    used = sorted({g for c in [*FLAVOURS.values(), *(OILS[o] for o in WEB_OILS)] for g in c.groups},
                  key=lambda g: SUBGROUPS[g].sid)
    sub = {g: {"sid": SUBGROUPS[g].sid, "main": SUBGROUPS[g].main, "R": SUBGROUPS[g].R,
               "Q": SUBGROUPS[g].Q} for g in used}
    mains = sorted({SUBGROUPS[g].main for g in used})
    inter = {}
    for a in mains:
        for b in mains:
            v = interaction(a, b)
            inter[f"{a},{b}"] = 0.0 if v is None else v
    data = {"generated_by": "Codes/python/tools/export_web_data.py",
            "T_grid_C": [T_GRID[0], T_GRID[-1]], "compounds": comps, "oils": oils,
            "subgroups": sub, "interactions": inter, "curves": curves}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(data, f, separators=(",", ":"))
    print("wrote", os.path.normpath(OUT), os.path.getsize(OUT), "bytes")


if __name__ == "__main__":
    main()
