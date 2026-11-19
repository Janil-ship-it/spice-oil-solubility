"""
Thermodynamic derivatives evaluated numerically (Member 2 / Member 5).

The Gibbs-Helmholtz equation, which follows from G = H - TS and the
Maxwell-type relation (d(G/T)/dT)_P = -H/T^2, links an activity coefficient
to the partial molar excess enthalpy of the solute:

    (d ln gamma_i / dT)_{P,x} = -hE_i / (R T^2)

With it we can (1) split non-ideality into enthalpy and entropy parts,
(2) carry a measured activity coefficient from one temperature to another,
and (3) get the slope of the solubility curve without re-solving it.

Solubility slope (implicit-function theorem on F = ln x + ln gamma - ln a_sat = 0):

    d ln x / dT = [dHfus(T)/(R T^2) - (d ln gamma/dT)_x] / [1 + (d ln gamma / d ln x)_T]

dHfus(T) = dHfus + dCp (T - Tm) is the Kirchhoff-corrected fusion enthalpy and
the denominator is the "thermodynamic factor" (1 for an ideal solution).
"""
import numpy as np
from scipy.integrate import quad

from .sle import R_GAS


def central_diff(f, x0, h):
    """df/dx at x0 by a second-order central difference."""
    return (f(x0 + h) - f(x0 - h)) / (2.0 * h)


def dlngamma_dT(model, x, T, i=0, h=0.01):
    return central_diff(lambda t: model.ln_gamma(x, t)[i], T, h)


def partial_molar_excess_enthalpy(model, x, T, i=0):
    """hE_i = -R T^2 (d ln gamma_i / dT)  [J/mol]  (Gibbs-Helmholtz)."""
    return -R_GAS * T ** 2 * dlngamma_dT(model, x, T, i)


def partial_molar_excess_entropy(model, x, T, i=0):
    """sE_i = (hE_i - gE_i)/T with gE_i = R T ln gamma_i  [J/(mol K)]."""
    gE = R_GAS * T * model.ln_gamma(x, T)[i]
    return (partial_molar_excess_enthalpy(model, x, T, i) - gE) / T


def thermodynamic_factor(model, x_solute, T, offset=0.0):
    """1 + (d ln gamma_s / d ln x_s)_T for a binary solute(0) + solvent(1)."""
    def lng(lnx):
        xs = np.exp(lnx)
        return model.ln_gamma([xs, 1.0 - xs], T)[0]
    lnx0 = np.log(x_solute)
    h = 1e-4 if x_solute < 0.5 else 1e-6
    return 1.0 + central_diff(lng, lnx0, h)


def dlnx_dT(solute, model, x_solute, T, offset_dT=0.0):
    """Slope of the solubility curve d ln x / dT (1/K) from the implicit-function theorem.

    offset_dT is d(delta)/dT if a temperature-dependent correction is applied.
    """
    dH_T = solute.dHfus + solute.dCp * (T - solute.Tm)
    dlnasat_dT = dH_T / (R_GAS * T ** 2)
    x = [x_solute, 1.0 - x_solute]
    return (dlnasat_dT - dlngamma_dT(model, x, T) - offset_dT) / thermodynamic_factor(model, x_solute, T)


def dissolution_properties(solute, model, x_solute, T):
    """Apparent standard dissolution functions at the solubility limit.

    dH_sol = R T^2 (d ln x/dT);  dG_sol = -R T ln x;  dS_sol = (dH - dG)/T
    """
    dH = R_GAS * T ** 2 * dlnx_dT(solute, model, x_solute, T)
    dG = -R_GAS * T * np.log(x_solute)
    return {"dH_sol": dH, "dG_sol": dG, "dS_sol": (dH - dG) / T}


def extrapolate_ln_gamma(ln_gamma_ref, T_ref, T, hE=None):
    """Carry ln gamma from T_ref to T with Gibbs-Helmholtz.

    hE : None    -> regular-solution assumption (sE = 0): hE = R T_ref ln gamma_ref,
                    so ln gamma(T) = ln gamma_ref * T_ref / T   (fallback, no groups needed)
         float   -> constant hE (J/mol)
         callable-> hE(T) in J/mol, integrated numerically
    """
    if hE is None:
        return ln_gamma_ref * T_ref / T
    if callable(hE):
        integral, _ = quad(lambda t: hE(t) / (R_GAS * t ** 2), T_ref, T)
        return ln_gamma_ref - integral
    return ln_gamma_ref + (hE / R_GAS) * (1.0 / T - 1.0 / T_ref)
