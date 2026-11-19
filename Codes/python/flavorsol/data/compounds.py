"""
Pure-component data for the flavour solutes, validation solutes and small
solvents (Member 1 - Data).

Every solute needs four things for solid-liquid equilibrium (SLE):
    M      molar mass, g/mol
    Tm     normal melting temperature, K
    dHfus  enthalpy of fusion, J/mol
    groups UNIFAC subgroup counts {subgroup name: count}
dCp (liquid minus solid heat capacity) is set to 0 by default, the usual
first approximation; the sensitivity study tests dCp = dHfus/Tm.

Replace a value here if your group measured or chose a different one in M1;
everything downstream recomputes automatically.
"""
from dataclasses import dataclass, field


@dataclass(frozen=True)
class Compound:
    name: str
    M: float                    # g/mol
    groups: dict                # UNIFAC subgroups
    Tm: float = None            # K   (solids only)
    dHfus: float = None         # J/mol (solids only)
    dCp: float = 0.0            # J/(mol K)
    formula: str = ""
    source: str = ""
    note: str = ""
    rho: float = None           # liquid density at 25 C, g/mL (solvents only)

    @property
    def is_solid_solute(self):
        return self.Tm is not None and self.dHfus is not None


# --------------------------------------------------------------------------
# Flavour solutes (the three compounds of the project)
# --------------------------------------------------------------------------
VANILLIN = Compound(
    name="vanillin", formula="C8H8O3", M=152.15,
    Tm=355.1, dHfus=22.4e3,
    # 4-hydroxy-3-methoxybenzaldehyde
    groups={"ACH": 3, "AC": 2, "ACOH": 1, "CHO": 1, "CH3O": 1},
    source=("Tm: NIST WebBook average of 7 values (355.1 +/- 0.9 K); "
            "dHfus: Temprado, Roux & Chickos, J. Therm. Anal. Calorim. 2008, 94, 257."),
)

CURCUMIN = Compound(
    name="curcumin", formula="C21H20O6", M=368.38,
    Tm=456.0, dHfus=37.19e3,
    # diketo form: Ar-CH=CH-C(=O)-CH2-C(=O)-CH=CH-Ar, Ar = 4-hydroxy-3-methoxyphenyl
    groups={"ACH": 6, "AC": 4, "ACOH": 2, "CH3O": 2, "CH=CH": 2,
            "CH2CO": 1, "C=O*": 1},
    source=("Tm: 183 C (form I), PubChem; dHfus: Ranjbar et al., "
            "CrystEngComm 2024, 26, 1077 (value used for pure curcumin)."),
    note=("beta-diketone fragmented as CH2CO + derived C=O* so atom counts are "
          "exact; the keto-enol form's intramolecular H-bond is not resolved "
          "by UNIFAC. Literature dHfus spans roughly 37-49 kJ/mol."),
)

PIPERINE = Compound(
    name="piperine", formula="C17H19NO3", M=285.34,
    Tm=406.05, dHfus=135.32 * 285.34,      # 135.32 J/g x M = 38.61 kJ/mol
    # 1,3-benzodioxole ring + penta-2,4-dienoyl chain + piperidine amide
    groups={"ACH": 3, "AC": 3, "CH2O": 1, "O*": 1, "CH=CH": 2,
            "HCON(CH2)2": 1, "CH2": 3},
    source=("Tm 132.90 C and dHfus 135.32 J/g by DSC: Fitriani, Simbolon & "
            "Zaini, J. Sains Farm. Klin. 2023, 10, 353."),
    note=("Tertiary amide represented by the DMF main group (HCON(CH2)2, one "
          "surplus H) because the amide main group has no published "
          "parameters with ACH, CH2O or CCOO."),
)

FLAVOURS = {c.name: c for c in (VANILLIN, CURCUMIN, PIPERINE)}

# --------------------------------------------------------------------------
# Structural analogues with measured SLE data in pure solvents (validation)
# --------------------------------------------------------------------------
ETHYL_VANILLIN = Compound(
    name="ethyl vanillin", formula="C9H10O3", M=166.17,
    Tm=352.98, dHfus=23.75e3,
    groups={"ACH": 3, "AC": 2, "ACOH": 1, "CHO": 1, "CH2O": 1, "CH3": 1},
    source="Wu et al., J. Chem. Thermodyn. 2017, 105, 345 (via NIST ThermoML).",
)

VANILLYL_ALCOHOL = Compound(
    name="vanillyl alcohol", formula="C8H10O3", M=154.16,
    Tm=389.0, dHfus=30.74e3,
    groups={"ACH": 3, "ACCH2": 1, "AC": 1, "ACOH": 1, "CH3O": 1, "OH": 1},
    source="Guo et al., J. Chem. Thermodyn. 2017, 106, 276 (via NIST ThermoML).",
)

ANALOGUES = {c.name: c for c in (ETHYL_VANILLIN, VANILLYL_ALCOHOL)}

# --------------------------------------------------------------------------
# Small solvents (validation)
# --------------------------------------------------------------------------
SOLVENTS = {
    "ethanol":        Compound("ethanol", 46.07, {"CH3": 1, "CH2": 1, "OH": 1}, rho=0.785),
    "1-propanol":     Compound("1-propanol", 60.10, {"CH3": 1, "CH2": 2, "OH": 1}, rho=0.8),
    "2-propanol":     Compound("2-propanol", 60.10, {"CH3": 2, "CH": 1, "OH": 1}, rho=0.781),
    "1-octanol":      Compound("1-octanol", 130.23, {"CH3": 1, "CH2": 7, "OH": 1}, rho=0.824),
    "acetone":        Compound("acetone", 58.08, {"CH3": 1, "CH3CO": 1}, rho=0.784),
    "methyl acetate": Compound("methyl acetate", 74.08, {"CH3": 1, "CH3COO": 1}, rho=0.932),
    "ethyl acetate":  Compound("ethyl acetate", 88.11, {"CH3": 1, "CH2": 1, "CH3COO": 1}, rho=0.894),
    "butyl acetate":  Compound("butyl acetate", 116.16, {"CH3": 1, "CH2": 3, "CH3COO": 1}, rho=0.876),
    "water":          Compound("water", 18.015, {"H2O": 1}, rho=0.997),
    "n-hexane":       Compound("n-hexane", 86.18, {"CH3": 2, "CH2": 4}, rho=0.655),
    "2-butanone":     Compound("2-butanone", 72.11, {"CH3": 1, "CH2": 1, "CH3CO": 1}, rho=0.8),
}

ALL_SOLUTES = {**FLAVOURS, **ANALOGUES}


def get_compound(name: str) -> Compound:
    key = name.strip().lower()
    for table in (FLAVOURS, ANALOGUES, SOLVENTS):
        if key in table:
            return table[key]
    raise KeyError(f"Unknown compound '{name}'. Known: "
                   f"{sorted([*FLAVOURS, *ANALOGUES, *SOLVENTS])}")


def atom_check(c: Compound):
    """Count C and O atoms implied by the UNIFAC groups (used in tests)."""
    C = {"CH3": 1, "CH2": 1, "CH": 1, "C": 1, "CH=CH": 2, "CH=C": 2, "ACH": 1,
         "AC": 1, "ACCH3": 2, "ACCH2": 2, "OH": 0, "H2O": 0, "ACOH": 1,
         "CH3CO": 2, "CH2CO": 2, "CHO": 1, "CH3COO": 2, "CH2COO": 2,
         "CH3O": 1, "CH2O": 1, "HCON(CH2)2": 3, "C=O*": 1, "O*": 0}
    O = {"OH": 1, "H2O": 1, "ACOH": 1, "CH3CO": 1, "CH2CO": 1, "CHO": 1,
         "CH3COO": 2, "CH2COO": 2, "CH3O": 1, "CH2O": 1, "HCON(CH2)2": 1,
         "C=O*": 1, "O*": 1}
    nC = sum(C[g] * n for g, n in c.groups.items())
    nO = sum(O.get(g, 0) * n for g, n in c.groups.items())
    return nC, nO
