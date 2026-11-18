"""Tests for the audit additions: HSP layer, numerical dCp, batch, metrics, mixtures."""
import os
import sys
import warnings
from dataclasses import replace

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
warnings.filterwarnings("ignore", message="UNIFAC: no published parameters")

from flavorsol.data.compounds import CURCUMIN, PIPERINE, VANILLIN
from flavorsol.data.oils import OILS
from flavorsol.hsp import HSP, hansen_distance, ln_gamma_hansen, red
from flavorsol.mixture import solve_mixture
from flavorsol.operating import batch, convert
from flavorsol.sle import ln_a_sat, ln_a_sat_numeric, solubility
from flavorsol.validation import metrics


def test_numeric_dcp_integration_equals_closed_form():
    c = replace(CURCUMIN, dCp=80.0)
    for T in (298.15, 373.15, 440.0):
        assert abs(ln_a_sat_numeric(c, T) - float(ln_a_sat(c, T))) < 1e-9


def test_hsp_distance_red_and_regular_solution():
    a = HSP(18.0, 8.0, 10.0, V=250.0, R0=8.0)
    b = HSP(16.0, 2.0, 4.0, V=900.0)
    assert hansen_distance(a, a) == 0.0
    assert abs(hansen_distance(a, b) - np.sqrt(16 + 36 + 36)) < 1e-12
    assert abs(red(a, b) - np.sqrt(88) / 8.0) < 1e-12
    assert ln_gamma_hansen(a, HSP(18.0, 8.0, 10.0, V=250.0), 300.0) == 0.0


def test_batch_mass_balance_and_limits():
    rows = batch("sunflower", 25, 15.0, {"curcumin": 0.047, "vanillin": 0.015})
    for r in rows:
        assert abs(r["dissolved_g"] + r["undissolved_g"] - r["available_g"]) < 1e-9
    cur = [r for r in rows if r["compound"] == "curcumin"][0]
    van = [r for r in rows if r["compound"] == "vanillin"][0]
    assert cur["limited_by"] == "solubility" and van["limited_by"] == "amount added"


def test_unit_conversions_consistent():
    x = solubility(PIPERINE, OILS["coconut"], 313.15).x
    u = convert(x, PIPERINE, OILS["coconut"])
    assert abs(u["ppm_w"] / 1e4 - u["wt_percent"]) < 1e-9
    assert abs(u["mg_per_g_oil"] * 100 - u["mg_per_100g_oil"]) < 1e-9


def test_metrics():
    m = metrics([2.0, 0.5], [1.0, 1.0])
    assert abs(m["mean_abs_log10"] - np.log10(2)) < 1e-12 and abs(m["bias_log10"]) < 1e-12


def test_binary_mixture_with_absent_partner_limit():
    # a vanishing loading of the partner leaves the other solute's solubility unchanged
    oil, T = OILS["sunflower"], 333.15
    alone = solubility(CURCUMIN, oil, T).x
    mix = solve_mixture([CURCUMIN, VANILLIN], oil, T, loadings_mg_per_100g=[1e6, 1e-9])
    assert abs(mix.x[0] - alone) / alone < 1e-6


def test_flavordb_records_loaded_and_consistent():
    """FlavorDB2 (IIIT-Delhi) record for vanillin is used, and agrees with our library."""
    from flavorsol.data import spices
    from flavorsol.data.compounds import get_compound
    rec = spices.FLAVORDB["vanillin"]
    assert rec["status"] == "retrieved" and int(rec["pubchem_cid"]) == 1183
    words, src = spices.flavor_profile("vanillin")
    assert "vanilla" in words and src.startswith("FlavorDB2")
    assert abs(float(rec["molecular_weight"]) - get_compound("vanillin").M) < 0.1
    # compounds without a retrieved record say so instead of posing as FlavorDB data
    assert "not in FlavorDB2" in spices.flavor_profile("curcumin")[1]
    # SpiceRx supplies curcumin's spice link that FlavorDB2 lacks
    assert "Turmeric" in spices.spices_containing("curcumin")["SpiceRx"]
    assert "Pepper" in spices.spices_containing("piperine")["FlavorDB2"]
    # every retrieved FlavorDB2 record has the same molar mass as our library
    from flavorsol.validation import database_check
    rows = database_check()
    assert {r["compound"] for r in rows} == {"vanillin", "thymol", "piperine"}
    assert all(abs(r["M_library"] - r["M_flavordb"]) < 0.01 for r in rows)


def test_fetcher_parses_entity_record():
    import importlib.util, os
    path = os.path.join(os.path.dirname(__file__), "..", "tools", "fetch_flavordb.py")
    spec = importlib.util.spec_from_file_location("fetch_flavordb", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    m = {"pubchem_id": 1183, "cas_id": "121-33-5", "fema_flavor_profile": "Vanilla",
         "functional_groups": "aldehyde@phenol", "molecular_weight": 152.149}
    row = mod.molecule_row("vanillin", m, 284, "Beans")
    assert row["cas"] == "121-33-5" and row["functional_groups"] == "aldehyde; phenol"
    assert row["status"] == "retrieved" and "284" in row["flavordb_source"]
