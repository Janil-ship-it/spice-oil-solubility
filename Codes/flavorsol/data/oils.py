"""
Cooking oils as triglyceride pseudo-components (Member 1 - Data).

An oil is a mixture of many triglycerides. We represent it as ONE average
triglyceride whose UNIFAC groups are the composition-weighted average of its
fatty-acid chains. This is the usual approach for vegetable oils in
group-contribution models and keeps the solver binary-like.

Fatty-acid compositions (wt% of total fatty acids) are typical values chosen
inside the ranges of Codex Alimentarius CXS 210-1999 (Standard for Named
Vegetable Oils). MCT oil (Miglyol 812 type) follows the Ph. Eur. definition
(caprylic 50-80 %, capric 20-50 %). Edit the dictionaries below to model a
specific brand or your lab's oil.

Groups for one acyl chain Cn:d  (R-CH2-C(=O)-O-):
    CH2COO 1, CH3 1, CH=CH d, CH2 (n - 3 - 2d)
Glycerol backbone: CH2 2, CH 1
"""
from dataclasses import dataclass

from .compounds import Compound

# fatty acid: (carbons, double bonds, molar mass g/mol)
FATTY_ACIDS = {
    "C6:0":  (6, 0, 116.16),  "C8:0":  (8, 0, 144.21),  "C10:0": (10, 0, 172.26),
    "C12:0": (12, 0, 200.32), "C14:0": (14, 0, 228.37), "C16:0": (16, 0, 256.42),
    "C16:1": (16, 1, 254.41), "C18:0": (18, 0, 284.48), "C18:1": (18, 1, 282.46), "C18:2": (18, 2, 280.45),
    "C18:3": (18, 3, 278.43),
}

OIL_COMPOSITIONS = {
    # high-linoleic sunflower oil
    "sunflower": {"C16:0": 6.5, "C18:0": 4.0, "C18:1": 25.0, "C18:2": 64.5},
    # coconut oil (lauric)
    "coconut": {"C6:0": 0.5, "C8:0": 7.0, "C10:0": 6.0, "C12:0": 48.0,
                "C14:0": 18.0, "C16:0": 9.0, "C18:0": 3.0, "C18:1": 6.5,
                "C18:2": 2.0},
    # used only to compare with literature data
    "corn": {"C16:0": 11.0, "C18:0": 2.0, "C18:1": 28.0, "C18:2": 58.0, "C18:3": 1.0},
    "olive": {"C16:0": 12.5, "C16:1": 1.0, "C18:0": 2.5, "C18:1": 73.0, "C18:2": 10.0,
              "C18:3": 1.0},
    "MCT": {"C8:0": 57.0, "C10:0": 42.0, "C12:0": 1.0},
}

# density at 20 C (g/mL); vegetable oils lose ~0.00064 g/mL per K
# (Noureddini, Teoh & Clements, JAOCS 1992, 69, 1184).
DENSITY_20C = {"sunflower": 0.919, "coconut": 0.924, "corn": 0.919, "olive": 0.912,
               "MCT": 0.945}
DRHO_DT = -0.00064

GLYCEROL_M, WATER_M = 92.094, 18.015


@dataclass(frozen=True)
class Oil(Compound):
    composition: dict = None
    rho20: float = 0.92

    def density(self, T_kelvin: float) -> float:
        """Liquid density in g/mL at temperature T (linear in T)."""
        return self.rho20 + DRHO_DT * (T_kelvin - 293.15)


def chain_groups(fa: str) -> dict:
    n, d, _ = FATTY_ACIDS[fa]
    return {"CH2COO": 1.0, "CH3": 1.0, "CH=CH": float(d), "CH2": float(n - 3 - 2 * d)}


def make_oil(name: str, composition: dict = None) -> Oil:
    """Build the triglyceride pseudo-component for a named oil."""
    comp = composition or OIL_COMPOSITIONS[name]
    total = sum(comp.values())
    # weight % of fatty acids -> mole fraction of acyl chains
    moles = {fa: w / total / FATTY_ACIDS[fa][2] for fa, w in comp.items()}
    ntot = sum(moles.values())
    y = {fa: n / ntot for fa, n in moles.items()}

    groups = {"CH2": 2.0, "CH": 1.0}               # glycerol backbone
    for fa, yi in y.items():
        for g, nu in chain_groups(fa).items():
            groups[g] = groups.get(g, 0.0) + 3.0 * yi * nu
    groups = {g: round(v, 6) for g, v in groups.items() if v > 1e-12}

    M_fa = sum(yi * FATTY_ACIDS[fa][2] for fa, yi in y.items())
    M = GLYCEROL_M + 3.0 * M_fa - 3.0 * WATER_M
    return Oil(name=f"{name} oil", M=M, groups=groups, composition=dict(comp),
               rho20=DENSITY_20C.get(name, 0.92),
               source="Codex Alimentarius CXS 210-1999 typical composition")


OILS = {k: make_oil(k) for k in OIL_COMPOSITIONS}


def get_oil(name: str) -> Oil:
    key = name.replace(" oil", "").strip()
    for k in OILS:
        if k.lower() == key.lower():
            return OILS[k]
    raise KeyError(f"Unknown oil '{name}'. Known: {sorted(OILS)}")
