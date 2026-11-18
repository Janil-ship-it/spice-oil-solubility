"""
Original UNIFAC parameters (Member 1 - Data).

Source: DDBST "Parameters of the Original UNIFAC Model" (published parameter set,
https://www.ddbst.com/published-parameters-unifac.html), which reproduces the
Hansen et al. (1991) revision:
    Hansen, H. K.; Rasmussen, P.; Fredenslund, Aa.; Schiller, M.; Gmehling, J.
    Ind. Eng. Chem. Res. 1991, 30, 2352-2355.

Only the subgroups and main-group pairs needed for this project are stored.

Two "derived" subgroups are added by group additivity (documented in the
manuscript, Section 2.3). UNIFAC subgroups inside one main group differ by CH2
increments, so subtracting the CH2 contribution gives a consistent value:
    C=O  (main group 9)  = CH2CO - CH2  ->  R = 1.4457 - 0.6744, Q = 1.180 - 0.540
    O    (main group 13) = CH2O  - CH2  ->  R = 0.9183 - 0.6744, Q = 0.780 - 0.540
They are needed for curcumin's central beta-diketone and piperine's
methylenedioxy bridge, which have no exact original-UNIFAC subgroup.
"""
from dataclasses import dataclass


@dataclass(frozen=True)
class Subgroup:
    sid: int      # subgroup number (DDBST numbering; 900+ = derived here)
    main: int     # main group number
    R: float      # van der Waals volume parameter
    Q: float      # van der Waals surface parameter


MAIN_GROUP_NAMES = {
    1: "CH2", 2: "C=C", 3: "ACH", 4: "ACCH2", 5: "OH", 7: "H2O", 8: "ACOH",
    9: "CH2CO", 10: "CHO", 11: "CCOO", 13: "CH2O", 39: "DMF",
}

SUBGROUPS = {
    # name        sid main   R       Q
    "CH3":      Subgroup(1, 1, 0.9011, 0.8480),
    "CH2":      Subgroup(2, 1, 0.6744, 0.5400),
    "CH":       Subgroup(3, 1, 0.4469, 0.2280),
    "C":        Subgroup(4, 1, 0.2195, 0.0000),
    "CH=CH":    Subgroup(6, 2, 1.1167, 0.8670),
    "CH=C":     Subgroup(8, 2, 0.8886, 0.6760),
    "ACH":      Subgroup(9, 3, 0.5313, 0.4000),
    "AC":       Subgroup(10, 3, 0.3652, 0.1200),
    "ACCH3":    Subgroup(11, 4, 1.2663, 0.9680),
    "ACCH2":    Subgroup(12, 4, 1.0396, 0.6600),
    "ACCH":     Subgroup(13, 4, 0.8121, 0.3480),
    "OH":       Subgroup(14, 5, 1.0000, 1.2000),
    "H2O":      Subgroup(16, 7, 0.9200, 1.4000),
    "ACOH":     Subgroup(17, 8, 0.8952, 0.6800),
    "CH3CO":    Subgroup(18, 9, 1.6724, 1.4880),
    "CH2CO":    Subgroup(19, 9, 1.4457, 1.1800),
    "CHO":      Subgroup(20, 10, 0.9980, 0.9480),   # aldehyde
    "CH3COO":   Subgroup(21, 11, 1.9031, 1.7280),
    "CH2COO":   Subgroup(22, 11, 1.6764, 1.4200),
    "CH3O":     Subgroup(24, 13, 1.1450, 1.0880),
    "CH2O":     Subgroup(25, 13, 0.9183, 0.7800),
    "HCON(CH2)2": Subgroup(73, 39, 2.6322, 2.1200),
    # derived by group additivity (see module docstring)
    "C=O*":     Subgroup(901, 9, round(1.4457 - 0.6744, 4), round(1.180 - 0.540, 4)),
    "O*":       Subgroup(902, 13, round(0.9183 - 0.6744, 4), round(0.780 - 0.540, 4)),
}

# Interaction parameters a_nm in kelvin, stored as (n, m): (a_nm, a_mn).
# Psi_nm = exp(-a_nm / T).  Copied from the DDBST published table.
_PAIRS = {
    (1, 2): (86.02, -35.36),     (1, 3): (61.13, -11.12),
    (1, 4): (76.50, -69.70),     (1, 5): (986.5, 156.4),
    (1, 7): (1318.0, 300.0),     (1, 8): (1333.0, 275.8),
    (1, 9): (476.4, 26.76),      (1, 10): (677.0, 505.7),
    (1, 11): (232.1, 114.8),     (1, 13): (251.5, 83.36),
    (1, 39): (485.3, -31.95),
    (2, 3): (38.81, 3.446),      (2, 4): (74.15, -113.6),
    (2, 5): (524.1, 457.0),      (2, 7): (270.6, 496.1),
    (2, 8): (526.1, 217.5),      (2, 9): (182.6, 42.92),
    (2, 10): (448.75, 56.30),    (2, 11): (37.85, 132.1),
    (2, 13): (214.5, 26.51),     (2, 39): (-70.45, 249.0),
    (3, 4): (167.0, -146.8),     (3, 5): (636.1, 89.60),
    (3, 7): (903.8, 362.3),      (3, 8): (1329.0, 25.34),
    (3, 9): (25.77, 140.1),      (3, 10): (347.3, 23.39),
    (3, 11): (5.994, 85.84),     (3, 13): (32.14, 52.13),
    (3, 39): (245.6, -133.9),
    (4, 5): (803.2, 25.82),      (4, 7): (5695.0, 377.6),
    (4, 8): (884.9, 244.2),      (4, 9): (-52.10, 365.8),
    (4, 10): (586.8, 106.0),     (4, 11): (5688.0, -170.0),
    (4, 13): (213.1, 65.69),     (4, 39): (5629.0, -240.2),
    (5, 7): (353.5, -229.1),     (5, 8): (-259.7, -451.6),
    (5, 9): (84.00, 164.5),      (5, 10): (-203.6, 529.0),
    (5, 11): (101.1, 245.4),     (5, 13): (28.06, 237.7),
    (5, 39): (-143.9, 64.16),
    (7, 8): (324.5, -601.8),     (7, 9): (-195.4, 472.5),
    (7, 10): (-116.0, 480.8),    (7, 11): (72.87, 200.8),
    (7, 13): (540.5, -314.7),    (7, 39): (319.0, -287.1),
    (8, 9): (-356.1, -133.1),    (8, 10): (-271.1, -155.6),
    (8, 11): (-449.4, -36.72),   (8, 13): (-162.8742, -178.5461),
    (9, 10): (-37.36, 128.0),    (9, 11): (-213.7, 372.2),
    (9, 13): (-103.6, 191.1),    (9, 39): (-61.70, 97.04),
    (10, 11): (-110.3, 185.1),   (10, 13): (304.1, -7.838),
    (10, 39): (-268.8, 13.89),
    (11, 13): (-235.7, 461.3),   (11, 39): (85.33, -82.12),
    (13, 39): (254.8, -158.2),
}


def interaction(n: int, m: int):
    """Return a_nm (K) for main groups n -> m, or None if not published."""
    if n == m:
        return 0.0
    if (n, m) in _PAIRS:
        return _PAIRS[(n, m)][0]
    if (m, n) in _PAIRS:
        return _PAIRS[(m, n)][1]
    return None
