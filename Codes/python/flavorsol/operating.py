"""
Operating conditions: from equilibrium solubility to what actually ends up in
the oil (Member 3 / Member 5).

Three different quantities - keep them apart:

1. Saturation solubility x_sat(T, oil)           -> needs only T and the oil
   (property of the system; independent of how much you add)
2. Dissolved amount at a given loading            -> needs masses added
   dissolved = min(added, capacity), solved jointly for mixtures
   (capacity of m_oil grams = x_sat converted to grams for that oil mass)
3. Extraction yield from a real spice             -> needs spice mass, compound
   content, and a release/accessibility factor eta (0-1) that thermodynamics
   cannot give: diffusion out of plant tissue, binding to fibre, time.
   eta = 1 gives the thermodynamic UPPER BOUND; measure eta experimentally.
"""
import numpy as np

from .data.compounds import FLAVOURS
from .data.oils import get_oil
from .mixture import solve_mixture
from .sle import solubility
from .units import binary_mg_per_100g, to_kelvin


def convert(x_solute, solute, oil, x_oil=None):
    """One saturation mole fraction in all the units people ask for."""
    x_oil = (1.0 - x_solute) if x_oil is None else x_oil
    if x_oil <= 0 or x_solute >= 1:
        return {"mole_fraction": x_solute, "wt_percent": 100.0, "ppm_w": 1e6,
                "mg_per_100g_oil": np.inf, "mg_per_g_oil": np.inf}
    m_s = x_solute * solute.M
    m_o = x_oil * oil.M
    return {
        "mole_fraction": float(x_solute),
        "wt_percent": float(100 * m_s / (m_s + m_o)),          # of the solution
        "ppm_w": float(1e6 * m_s / (m_s + m_o)),               # mg per kg solution
        "mg_per_100g_oil": float(1e5 * m_s / m_o),             # per 100 g of oil (solvent basis)
        "mg_per_g_oil": float(1e3 * m_s / m_o),
    }


def batch(oil, T_C, oil_mass_g, added_mass_g: dict, spice=None):
    """Dissolved and undissolved mass for a real batch.

    oil          : oil name or Oil
    T_C          : temperature, deg C
    oil_mass_g   : grams of oil
    added_mass_g : {'curcumin': g, 'piperine': g, ...} pure compound added
    spice        : optional {'curcumin': (spice_g, content_fraction, eta)} to
                   compute the amount AVAILABLE from a spice instead of pure compound
    Returns a list of dicts, one per compound.
    """
    oil = get_oil(oil) if isinstance(oil, str) else oil
    names = list(added_mass_g) if added_mass_g else list(spice)
    solutes = [FLAVOURS[n] for n in names]
    avail = {}
    for n in names:
        if spice and n in spice:
            g, frac, eta = spice[n]
            avail[n] = g * frac * eta
        else:
            avail[n] = added_mass_g[n]
    loads = [1e5 * avail[n] / oil_mass_g for n in names]      # g/g -> mg per 100 g oil
    T = float(to_kelvin(T_C))
    mix = solve_mixture(solutes, oil, T, loadings_mg_per_100g=loads)
    out = []
    for i, (n, s) in enumerate(zip(names, solutes)):
        single = solubility(s, oil, T)
        cap_alone_g = binary_mg_per_100g(single.x, s, oil) / 1e5 * oil_mass_g
        dissolved_g = mix.mg_per_100g[i] / 1e5 * oil_mass_g
        out.append({
            "compound": n, "T_C": T_C, "oil": oil.name, "oil_mass_g": oil_mass_g,
            "available_g": avail[n], "oil_to_solute_ratio": oil_mass_g / avail[n] if avail[n] else np.inf,
            "capacity_alone_g": cap_alone_g,
            "dissolved_g": dissolved_g, "undissolved_g": max(avail[n] - dissolved_g, 0.0),
            "dissolved_percent": 100 * dissolved_g / avail[n] if avail[n] else np.nan,
            "limited_by": "solubility" if mix.saturated[i] else "amount added",
            **{k: v for k, v in convert(mix.x[i], s, oil, mix.x[-1]).items()},
        })
    return out
