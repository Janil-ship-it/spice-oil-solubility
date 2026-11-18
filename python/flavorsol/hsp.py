"""
Hansen solubility parameters (HSP): a SCREENING layer, kept separate from UNIFAC
(Member 1 / Member 4).

Why separate. UNIFAC already computes the full activity coefficient
(size + group interactions). HSP describes the same non-ideality with three
cohesive-energy terms. Multiplying a UNIFAC gamma by an HSP "correction" would
count the solute-oil interaction twice, so the two are never combined into one
gamma. HSP is used for two jobs only:

1. Screening / consistency check: the Hansen distance Ra and the relative
   energy difference RED = Ra / R0 say whether solute and oil are "alike"
   (RED < 1) - a qualitative cross-check of the UNIFAC ranking.
2. An independent ALTERNATIVE gamma (regular-solution form with Hansen terms),
   reported beside UNIFAC so the two models can be compared with data:

       ln gamma_s = (V_s / R T) * phi_oil^2 * [(dD_s-dD_o)^2 + 0.25 (dP_s-dP_o)^2 + 0.25 (dH_s-dH_o)^2]

   (Hansen, Hansen Solubility Parameters: A User's Handbook, 2nd ed., 2007),
   optionally plus a Flory-Huggins size term.

No HSP values are shipped for curcumin, vanillin, piperine or the oils: we found
no measured values we could verify. Enter your own (measured by the solvent-
sphere method, or estimated with a group-contribution method such as
Stefanis & Panayiotou, Int. J. Thermophys. 2008) in a JSON/dict and pass them in.
"""
from dataclasses import dataclass

import numpy as np

R_GAS = 8.314462618


@dataclass(frozen=True)
class HSP:
    dD: float          # dispersion, MPa^0.5
    dP: float          # polar, MPa^0.5
    dH: float          # hydrogen bonding, MPa^0.5
    V: float = None    # molar volume, cm3/mol
    R0: float = None   # interaction radius (solute sphere), MPa^0.5
    source: str = ""


def hansen_distance(a: HSP, b: HSP) -> float:
    """Ra = sqrt(4 (dD1-dD2)^2 + (dP1-dP2)^2 + (dH1-dH2)^2)   [MPa^0.5]"""
    return float(np.sqrt(4 * (a.dD - b.dD) ** 2 + (a.dP - b.dP) ** 2 + (a.dH - b.dH) ** 2))


def red(solute: HSP, solvent: HSP) -> float:
    """Relative energy difference Ra/R0; < 1 = inside the solubility sphere."""
    if not solute.R0:
        raise ValueError("RED needs the solute's interaction radius R0")
    return hansen_distance(solute, solvent) / solute.R0


def ln_gamma_hansen(solute: HSP, solvent: HSP, T: float, phi_solvent: float = 1.0,
                    flory_huggins: bool = True) -> float:
    """Alternative (NOT combined with UNIFAC) activity coefficient of the solute.

    Regular-solution enthalpic term with Hansen's three-part distance; with
    flory_huggins=True the athermal size term ln(phi_s/x_s) + 1 - phi_s/x_s is
    added in its infinite-dilution form ln(V_s/V_o) + 1 - V_s/V_o.
    """
    if not (solute.V and solvent.V):
        raise ValueError("molar volumes V are needed for the Hansen activity coefficient")
    A = (solute.dD - solvent.dD) ** 2 + 0.25 * (solute.dP - solvent.dP) ** 2 \
        + 0.25 * (solute.dH - solvent.dH) ** 2               # MPa = J/cm3
    lng = solute.V * A * phi_solvent ** 2 / (R_GAS * T)
    if flory_huggins:
        r = solute.V / solvent.V
        lng += np.log(r) + 1.0 - r
    return float(lng)


def mix_hsp(components, volume_fractions):
    """Volume-fraction average HSP of a mixture (e.g. an oil from its triglycerides)."""
    phi = np.asarray(volume_fractions, float)
    phi = phi / phi.sum()
    dD = sum(p * c.dD for p, c in zip(phi, components))
    dP = sum(p * c.dP for p, c in zip(phi, components))
    dH = sum(p * c.dH for p, c in zip(phi, components))
    V = sum(p * c.V for p, c in zip(phi, components)) if all(c.V for c in components) else None
    return HSP(dD, dP, dH, V, source="volume-fraction average")
