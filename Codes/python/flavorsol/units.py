"""Unit conversions between mole fraction and kitchen/lab units."""
import numpy as np


def mg_per_100g_oil(x_solute, M_solute, x_oil, M_oil):
    """mg of solute per 100 g of oil, from mole fractions in the liquid."""
    if x_oil <= 0:
        return np.inf
    return 1e5 * x_solute * M_solute / (x_oil * M_oil)


def binary_mg_per_100g(x, solute, oil):
    """Binary solute + oil: mg solute per 100 g oil."""
    if x >= 1.0:
        return np.inf
    return mg_per_100g_oil(x, solute.M, 1.0 - x, oil.M)


def mg_per_100g_to_moles_per_mol_oil(mg100, solute, oil):
    """Loading (mg per 100 g oil) -> mol solute per mol oil."""
    return (mg100 / 1e3 / solute.M) / (100.0 / oil.M)


def mg_per_mL_to_mg_per_100g(c_mg_mL, rho_g_mL):
    """Literature values in mg/mL of oil -> mg per 100 g oil."""
    return 100.0 * c_mg_mL / rho_g_mL


def x_from_mg_per_100g(mg100, solute, oil):
    """Binary mole fraction from mg solute per 100 g oil."""
    n_s = mg100 / 1e3 / solute.M
    n_o = 100.0 / oil.M
    return n_s / (n_s + n_o)


def to_kelvin(t_celsius):
    return np.asarray(t_celsius, dtype=float) + 273.15


def to_celsius(T_kelvin):
    return np.asarray(T_kelvin, dtype=float) - 273.15
