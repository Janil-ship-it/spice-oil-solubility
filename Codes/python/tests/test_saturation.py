"""Saturation point: statements, capacity, saturation temperature (flavorsol.saturation)."""
import math
import warnings

import pytest

from flavorsol.data.compounds import FLAVOURS
from flavorsol.data.oils import get_oil
from flavorsol.saturation import (mixture_saturation, saturation_temperature,
                                  single_saturation, solubility_mg)

warnings.filterwarnings("ignore", message="UNIFAC: no published parameters")


def test_saturated_when_added_exceeds_capacity():
    r = single_saturation("curcumin", "coconut", 60, 15.0, 0.047)
    assert r["saturated"]
    assert "saturation is achieved" in r["statement"]
    assert "More mixing cannot be done" in r["statement"]
    # mass balance
    assert r["dissolved_mg"] + r["excess_mg"] == pytest.approx(r["added_mg"])
    assert r["dissolved_mg"] == pytest.approx(r["capacity_mg"])


def test_not_saturated_when_capacity_is_larger():
    r = single_saturation("curcumin", "coconut", 120, 15.0, 0.047)
    assert not r["saturated"]
    assert r["excess_mg"] == 0
    assert "not reached" in r["statement"]


def test_saturation_temperature_is_root_of_solubility_curve():
    s, oil = FLAVOURS["curcumin"], get_oil("sunflower")
    L = 313.3
    T_sat, kind = saturation_temperature(s, oil, L)
    assert kind == "root"
    mg, _ = solubility_mg(s, oil, T_sat)
    assert mg == pytest.approx(L, rel=1e-5)
    # below T_sat the oil is saturated, above it is not
    assert solubility_mg(s, oil, T_sat - 2)[0] < L < solubility_mg(s, oil, T_sat + 2)[0]


def test_solubility_increases_with_temperature():
    """Dissolution is endothermic, so S(T) must rise monotonically (needed for a unique T_sat)."""
    for name in FLAVOURS:
        for o in ("sunflower", "coconut"):
            s, oil = FLAVOURS[name], get_oil(o)
            vals = [solubility_mg(s, oil, t)[0] for t in range(25, 201, 5)]
            assert all(b >= a for a, b in zip(vals, vals[1:])), (name, o)


def test_molten_piperine_has_no_saturation():
    r = single_saturation("piperine", "sunflower", 150, 15.0, 5.0)
    assert not r["saturated"] and math.isinf(r["capacity_mg"])
    assert "no saturation point" in r["statement"]


def test_mixture_saturation_flags_curcumin_only():
    r = mixture_saturation("coconut", 60, 15.0, {"curcumin": 0.047, "piperine": 0.05},
                           with_temperature=False)
    sat = {row["compound"]: row["saturated"] for row in r["rows"]}
    assert sat == {"curcumin": True, "piperine": False}
    assert "saturation is achieved for curcumin" in r["headline"]
    for row in r["rows"]:
        assert row["dissolved_mg"] + row["excess_mg"] == pytest.approx(row["added_mg"])


def test_mixture_all_dissolved_temperature_close_to_single():
    r = mixture_saturation("coconut", 60, 15.0, {"curcumin": 0.047, "piperine": 0.05})
    single = single_saturation("curcumin", "coconut", 60, 15.0, 0.047)["saturation_T_C"]
    assert r["all_dissolved_kind"] == "root"
    # co-solvency of dissolved piperine lowers it slightly, never raises it much
    assert single - 5 < r["all_dissolved_T_C"] <= single + 0.5
