"""Unit tests. Run:  python -m pytest -q   (or  python tests/test_flavorsol.py)"""
import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
warnings.filterwarnings("ignore", message="UNIFAC: no published parameters")

from flavorsol.data.compounds import (ALL_SOLUTES, CURCUMIN, PIPERINE, SOLVENTS, VANILLIN,
                                      Compound, atom_check)
from flavorsol.data.oils import OILS
from flavorsol.mixture import solve_mixture
from flavorsol.sle import _ln_activity, ideal_solubility, ln_a_sat, solubility
from flavorsol.unifac import UNIFAC
from flavorsol.units import binary_mg_per_100g, x_from_mg_per_100g

SUN, COCO = OILS["sunflower"], OILS["coconut"]


def test_unifac_matches_ddbst_reference():
    m = UNIFAC([SOLVENTS["n-hexane"], SOLVENTS["2-butanone"]])
    assert np.allclose(m.gamma([0.5, 0.5], 333.15), [1.4276025835, 1.3646545010], atol=1e-9)


def test_pure_component_gamma_is_one():
    m = UNIFAC([CURCUMIN, SUN])
    lg = m.ln_gamma([1 - 1e-12, 1e-12], 350.0)
    assert abs(lg[0]) < 1e-9


def test_gibbs_duhem():
    m = UNIFAC([VANILLIN, COCO])
    T, x, h = 330.0, 0.2, 1e-6
    d0 = (m.ln_gamma([x + h, 1 - x - h], T) - m.ln_gamma([x - h, 1 - x + h], T)) / (2 * h)
    assert abs(x * d0[0] + (1 - x) * d0[1]) < 1e-7


def test_group_formulas_have_right_carbon_and_oxygen():
    expect = {"vanillin": (8, 3), "curcumin": (21, 6), "piperine": (17, 3),
              "ethyl vanillin": (9, 3), "vanillyl alcohol": (8, 3)}
    for name, c in ALL_SOLUTES.items():
        assert atom_check(c) == expect[name], name


def test_oil_pseudocomponent():
    assert 860 < SUN.M < 890 and 640 < COCO.M < 680
    assert SUN.groups["CH2COO"] == 3 and SUN.groups["CH"] == 1


def test_ideal_limit_recovered_when_solvent_is_chemically_identical():
    twin = Compound("liquid twin", VANILLIN.M, dict(VANILLIN.groups))
    r = solubility(VANILLIN, twin, 310.0)
    assert abs(r.x - float(ideal_solubility(VANILLIN, 310.0))) < 1e-10


def test_solution_satisfies_equilibrium_equation():
    for s in (VANILLIN, CURCUMIN, PIPERINE):
        r = solubility(s, SUN, 333.15)
        m = UNIFAC([s, SUN])
        resid = np.log(r.x) + m.ln_gamma([r.x, 1 - r.x], 333.15)[0] - ln_a_sat(s, 333.15)
        assert abs(resid) < 1e-10


def test_solubility_increases_with_temperature_below_melting():
    for s in (VANILLIN, CURCUMIN, PIPERINE):
        xs = [solubility(s, SUN, T).x for T in (298.15, 318.15, 338.15)]
        assert xs[0] < xs[1] < xs[2]


def test_above_melting_point_edge_cases():
    assert solubility(PIPERINE, SUN, 423.15).status == "miscible"
    r = solubility(VANILLIN, SUN, 393.15)
    assert r.status == "LLE-limited"
    xa, xb, _ = r.lle
    m = UNIFAC([VANILLIN, SUN])
    a1, a2 = _ln_activity(m, xa, 393.15), _ln_activity(m, xb, 393.15)
    assert np.allclose(a1, a2, atol=1e-8)        # equal activities in both liquids


def test_mixture_with_one_solute_equals_single_solver():
    single = solubility(CURCUMIN, COCO, 343.15).x
    mix = solve_mixture([CURCUMIN], COCO, 343.15)
    assert abs(mix.x[0] - single) / single < 1e-8


def test_mixture_mass_limited_and_saturated():
    mix = solve_mixture([VANILLIN, CURCUMIN, PIPERINE], SUN, 333.15,
                        loadings_mg_per_100g=[100, 314, 333])
    assert mix.saturated == [False, True, False]
    assert abs(mix.mg_per_100g[0] - 100) < 1e-6 and abs(mix.mg_per_100g[2] - 333) < 1e-6
    m = UNIFAC([VANILLIN, CURCUMIN, PIPERINE, SUN])
    lg = m.ln_gamma(mix.x, 333.15)
    assert abs(np.log(mix.x[1]) + lg[1] - ln_a_sat(CURCUMIN, 333.15)) < 1e-9


def test_unit_round_trip():
    x = 3.2e-3
    mg = binary_mg_per_100g(x, CURCUMIN, SUN)
    assert abs(x_from_mg_per_100g(mg, CURCUMIN, SUN) - x) < 1e-15


if __name__ == "__main__":
    tests = [v for k, v in dict(globals()).items() if k.startswith("test_")]
    for t in tests:
        t()
        print("PASS", t.__name__)
    print(f"{len(tests)} tests passed")
