"""
Boundary checks: flag unrealistic inputs and conditions where the compound or
the oil breaks down instead of simply dissolving (Member 4 - Validation).

Every check returns Issue objects with a level:
    'error'   - the calculation cannot be trusted; the engine stops
    'warning' - result is computed but needs a caveat in the report
    'info'    - context the user should know
"""
from dataclasses import dataclass

# Temperatures above which simple dissolution is no longer the whole story.
THERMAL_LIMITS_C = {
    "curcumin": [
        (150.0, "warning", "curcumin degrades in hot oil: 38.9 % loss was measured in corn oil "
                           "heated up to 200 C (Siddiqui, Pak. J. Pharm. Sci. 2015, 28, 299). "
                           "Dissolved amount will fall with time; prediction is an upper bound."),
    ],
    "piperine": [
        (100.0, "info", "piperine was stable as a solid at 100 C for 4 days (Kotte et al., "
                        "Anal. Methods 2014); stability above 100 C in oil is not verified."),
        (180.0, "warning", "no published stability data for piperine in oil at this temperature."),
    ],
    "vanillin": [
        (147.0, "warning", "above vanillin's flash point (147 C; b.p. 285 C): evaporation to air "
                           "removes vanillin from the pan; equilibrium solubility ignores this loss. "
                           "Vanillin itself is chemically stable to 250 C in water "
                           "(Doctor et al., Molecules 2020, 25, 1061)."),
    ],
}

# Smoke points (deg C) by oil and grade; above these the oil itself breaks down.
# Values from the smoke-point table in Wikipedia "Smoke point" (accessed Nov 2026):
# high-oleic refined sunflower 232 C, unrefined sunflower 107 C,
# refined coconut 204 C, unrefined (virgin) coconut 177 C. Add rows for other oils.
SMOKE_POINT_C = {
    ("sunflower", "refined"): 232.0, ("sunflower", "unrefined"): 107.0,
    ("coconut", "refined"): 204.0, ("coconut", "unrefined"): 177.0,
}
MELTING_RANGE_OIL_C = {"coconut": 24.0}     # coconut oil is partly solid below ~24 C

UNIFAC_FIT_RANGE_K = (275.0, 425.0)          # most original-UNIFAC data lie in this window


@dataclass
class Issue:
    level: str
    code: str
    message: str

    def __str__(self):
        return f"[{self.level.upper()}] {self.code}: {self.message}"


def check_temperature(T_C, solute_name=None, oil_name=None, oil_grade="refined", Tm=None):
    out = []
    if not (-40.0 <= T_C <= 300.0):
        out.append(Issue("error", "T_RANGE", f"T = {T_C} C is outside -40..300 C (kitchen and lab range)."))
        return out
    T_K = T_C + 273.15
    if not (UNIFAC_FIT_RANGE_K[0] <= T_K <= UNIFAC_FIT_RANGE_K[1]):
        out.append(Issue("info", "UNIFAC_EXTRAPOLATION",
                         f"T = {T_C:.0f} C is outside the 2-152 C window where most UNIFAC "
                         "parameters were fitted; treat predictions as extrapolations."))
    if Tm is not None and T_K >= Tm:
        out.append(Issue("info", "ABOVE_MELTING",
                         f"T is above the solute melting point ({Tm - 273.15:.0f} C): no solid can "
                         "exist; the result is a miscibility (liquid-liquid) limit, not a solubility."))
    for lim, level, msg in THERMAL_LIMITS_C.get(solute_name, []):
        if T_C >= lim:
            out.append(Issue(level, "THERMAL_STABILITY", msg))
    if oil_name:
        key = oil_name.replace(" oil", "")
        sp = SMOKE_POINT_C.get((key, oil_grade))
        if sp is not None and T_C >= sp:
            out.append(Issue("warning", "SMOKE_POINT",
                             f"T is at/above the smoke point of {oil_grade} {key} oil (~{sp:.0f} C): the oil "
                             "oxidises and decomposes, so its composition (the solvent) changes."))
        mr = MELTING_RANGE_OIL_C.get(key)
        if mr is not None and T_C < mr:
            out.append(Issue("warning", "OIL_PARTLY_SOLID",
                             f"{key} oil is partly solid below ~{mr:.0f} C; prediction assumes liquid oil."))
    return out


def check_fusion_data(Tm, dHfus, dCp=0.0):
    out = []
    if Tm is None:
        out.append(Issue("error", "TM_MISSING",
                         "Melting temperature Tm is mandatory (DSC, PubChem or supplier datasheet)."))
        return out
    if not (150.0 <= Tm <= 700.0):
        out.append(Issue("error", "TM_RANGE", f"Tm = {Tm} K is not a realistic melting point for a flavour compound."))
    if dHfus is not None:
        if dHfus <= 0:
            out.append(Issue("error", "DHFUS_SIGN", "Enthalpy of fusion must be positive (J/mol)."))
        else:
            dS = dHfus / Tm
            if dHfus < 1000:
                out.append(Issue("error", "DHFUS_UNITS",
                                 f"dHfus = {dHfus} looks like kJ/mol; the program expects J/mol."))
            elif not (20.0 <= dS <= 200.0):
                out.append(Issue("warning", "WALDEN",
                                 f"Entropy of fusion dHfus/Tm = {dS:.0f} J/(mol K) is far from the usual "
                                 "40-120 range (Walden's rule ~56.5); check the value."))
    if dCp is not None and abs(dCp) > 500:
        out.append(Issue("warning", "DCP_RANGE", f"|dCp| = {abs(dCp)} J/(mol K) is unusually large."))
    return out


def check_measurement(value, unit, x_meas=None, x_ideal=None, label="measurement"):
    out = []
    if value is None or value <= 0:
        out.append(Issue("error", "MEAS_VALUE", f"{label}: value must be positive."))
        return out
    if unit == "mole fraction" and value >= 1:
        out.append(Issue("error", "MEAS_FRACTION", f"{label}: a mole fraction must be below 1."))
    if x_meas is not None and x_ideal is not None and x_meas > 1.05 * x_ideal:
        out.append(Issue("warning", "ABOVE_IDEAL",
                         f"{label}: measured solubility exceeds the ideal-solution limit "
                         f"({x_meas:.3g} > {x_ideal:.3g}), which requires gamma < 1. Possible causes: "
                         "supersaturation after heating, a less stable solid form, impurities, or "
                         "colloidal particles counted as dissolved."))
    return out


def check_uncertainties(unc: dict):
    out = []
    for k, v in (unc or {}).items():
        if v is not None and v < 0:
            out.append(Issue("error", "UNCERTAINTY_SIGN", f"uncertainty '{k}' must be >= 0."))
    return out


def check_composition(comp: dict):
    out = []
    if comp is None:
        return out
    total = sum(comp.values())
    if any(v < 0 for v in comp.values()):
        out.append(Issue("error", "COMPOSITION_NEGATIVE", "fatty-acid percentages must be >= 0."))
    elif abs(total - 100.0) > 2.0:
        out.append(Issue("warning", "COMPOSITION_SUM",
                         f"fatty-acid composition sums to {total:.1f} %, renormalised to 100 %."))
    return out


def has_errors(issues):
    return any(i.level == "error" for i in issues)
