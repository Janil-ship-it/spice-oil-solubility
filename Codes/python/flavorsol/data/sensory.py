"""
Sensory thresholds and typical recipe loadings (Member 5 - Analysis).

Thresholds are reported in different matrices (mostly water); oil changes
volatility and release, so they are used as order-of-magnitude reference
lines, not as exact perception limits in oil.

Recipe loadings are illustrative assumptions for one tadka (tempering):
15 g oil (1 tbsp) with 1.5 g turmeric and 1 g black pepper; vanillin at a
typical 0.1 % flavouring level in a fat-based dessert. Change them to match
your own experiments.
"""

THRESHOLDS_MG_PER_KG = {
    # vanillin orthonasal odour threshold in water ~20 ug/L (retronasal ~30 ug/L),
    # as compiled by Plotto et al., Flavour Fragr. J. 2004, 19, 491.
    "vanillin": {"value": 0.020, "kind": "odour (orthonasal, water)",
                 "source": "Plotto et al., Flavour Fragr. J. 2004, 19, 491"},
    # piperine detection threshold 0.55 +/- 0.15 ppm (Moss et al. 2023),
    # cited by Spence, Int. J. Gastron. Food Sci. 2024, 35, 100900.
    "piperine": {"value": 0.55, "kind": "pungency detection",
                 "source": "Moss et al. 2023, cited in Spence, Int. J. Gastron. Food Sci. 2024, 35, 100900"},
    # curcumin is mainly a colourant; no flavour threshold was found.
    "curcumin": {"value": None, "kind": "colourant (no flavour threshold found)",
                 "source": "-"},
}

# Spice content of the active compound (mass fraction) and amount used
SPICE_CONTENT = {
    "curcumin": {"spice": "turmeric", "fraction": 0.0314,
                 "source": "Tayyem et al., Nutr. Cancer 2006, 55, 126 (mean 3.14 % in turmeric powder)"},
    "piperine": {"spice": "black pepper", "fraction": 0.05,
                 "source": "Gorgani et al., Compr. Rev. Food Sci. Food Saf. 2017, 16, 124 (about 2-7 %)"},
    "vanillin": {"spice": "vanillin flavouring", "fraction": 1.0, "source": "pure flavouring"},
}

RECIPE = {"oil_g": 15.0, "spice_g": {"curcumin": 1.5, "piperine": 1.0, "vanillin": 0.015}}


def recipe_loading_mg_per_100g(name: str) -> float:
    """mg of compound added per 100 g oil in the illustrative recipe."""
    mg = RECIPE["spice_g"][name] * SPICE_CONTENT[name]["fraction"] * 1000.0
    return mg * 100.0 / RECIPE["oil_g"]
