"""Tests for the input/decision layer, checks, derivatives and uncertainty."""
import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
warnings.filterwarnings("ignore", message="UNIFAC: no published parameters")

from flavorsol import checks
from flavorsol.data.compounds import CURCUMIN, PIPERINE, VANILLIN
from flavorsol.data.oils import OILS
from flavorsol.derivatives import (dlnx_dT, extrapolate_ln_gamma,
                                   partial_molar_excess_enthalpy)
from flavorsol.engine import Measurement, Request, load_request, predict
from flavorsol.sle import solubility
from flavorsol.unifac import UNIFAC
from flavorsol.units import binary_mg_per_100g

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SUN = OILS["sunflower"]


def test_route_selection():
    assert predict(Request("curcumin", "sunflower", 60))[0][0].route == "R5"
    assert predict(Request("curcumin", "sunflower", 60, gamma=50))[0][0].route == "R1"
    m_same = Measurement(5.0, "mg/100 g", 60)
    assert predict(Request("curcumin", "sunflower", 80, measurements=[m_same]))[0][0].route == "R2"
    m_other = Measurement(0.181, "mole fraction", 25, solvent="ethanol")
    assert predict(Request("ethyl vanillin", "sunflower", 40, measurements=[m_other]))[0][0].route == "R3"


def test_anchor_reproduces_measurement_at_its_temperature():
    p = predict(Request("curcumin", "corn", [25.0], measurements=[Measurement(0.07, "mg/mL", 25)]))[0][0]
    assert abs(p.mg_per_100g - 100 * 0.07 / OILS["corn"].rho20) < 1e-6


def test_gamma_route_skips_unifac_and_matches_formula():
    p = predict(Request("piperine", "coconut", [60.0], gamma=2.8, gamma_T_C=60.0))[0][0]
    from flavorsol.sle import ln_a_sat
    x = np.exp(ln_a_sat(PIPERINE, 333.15)) / 2.8
    assert abs(p.x - x) < 1e-12


def test_headspace_route_recovers_saturation():
    T = 313.15
    x_true = solubility(VANILLIN, SUN, T).x
    x_sol, p_sub = 1e-3, 0.05                      # synthetic consistency test
    p_meas = p_sub * x_sol / x_true                # Henry's law
    m = Measurement(p_meas, "Pa", 40.0, kind="headspace", x_solution=x_sol, p_sub_Pa=p_sub)
    p = predict(Request("vanillin", "sunflower", [40.0], measurements=[m]))[0][0]
    assert p.route == "R4" and abs(p.x - x_true) / x_true < 1e-8


def test_walden_fallback_used_when_dHfus_missing():
    preds, _ = predict(Request("thymol", "coconut", [30.0]))
    assert any("Walden" in a for a in preds[0].assumptions)


def test_boundary_checks():
    _, issues = predict(Request("vanillin", "sunflower", [60], dHfus=22.4))
    assert any(i.code == "DHFUS_UNITS" for i in issues)
    preds, _ = predict(Request("vanillin", "sunflower", [350]))
    assert preds[0].status == "rejected"
    codes = [i.code for i in checks.check_temperature(180, "curcumin", "coconut oil", "unrefined", 456.0)]
    assert "SMOKE_POINT" in codes and "THERMAL_STABILITY" in codes
    _, issues = predict(Request("curcumin", "sunflower", [60],
                                measurements=[Measurement(5.4, "mg/mL", 25)]))
    assert any(i.code == "ABOVE_IDEAL" for i in issues)


def test_gibbs_helmholtz_extrapolation_matches_unifac():
    m = UNIFAC([CURCUMIN, SUN])
    x = [1e-6, 1 - 1e-6]
    hE = lambda T: partial_molar_excess_enthalpy(m, x, T)
    lng = extrapolate_ln_gamma(m.ln_gamma(x, 300.0)[0], 300.0, 380.0, hE=hE)
    assert abs(lng - m.ln_gamma(x, 380.0)[0]) < 1e-5


def test_implicit_derivative_matches_finite_difference():
    m = UNIFAC([CURCUMIN, SUN])
    T, h = 350.0, 0.01
    x = solubility(CURCUMIN, SUN, T, m).x
    fd = (np.log(solubility(CURCUMIN, SUN, T + h, m).x) - np.log(solubility(CURCUMIN, SUN, T - h, m).x)) / (2 * h)
    assert abs(dlnx_dT(CURCUMIN, m, x, T) - fd) < 1e-5


def test_delta_method_agrees_with_monte_carlo_for_small_errors():
    req = Request("vanillin", "sunflower", [40.0], gamma=9.0, gamma_T_C=40.0, monte_carlo=3000,
                  uncertainties={"T": 0.5, "Tm": 0.5, "dHfus": 300.0, "dCp": 1.0, "ln_gamma_model": 0.05})
    p = predict(req)[0][0]
    lo, hi = p.ci95_mg_per_100g
    mlo, mhi = p.mc_ci95_mg_per_100g
    assert abs(lo - mlo) / lo < 0.03 and abs(hi - mhi) / hi < 0.03


def test_example_input_files_run():
    for f in sorted(os.listdir(os.path.join(HERE, "inputs"))):
        preds, issues = predict(load_request(os.path.join(HERE, "inputs", f)))
        assert preds or checks.has_errors(issues), f
