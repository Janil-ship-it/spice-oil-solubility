"""
Original UNIFAC activity-coefficient model (Member 2 - Single-solute solver).

    ln gamma_i = ln gamma_i^C (combinatorial: size and shape)
               + ln gamma_i^R (residual: group-group energy interactions)

Combinatorial part (Staverman-Guggenheim, z = 10):
    J_i = r_i / sum_j r_j x_j          L_i = q_i / sum_j q_j x_j
    ln gamma_i^C = 1 - J_i + ln J_i - 5 q_i [1 - J_i/L_i + ln(J_i/L_i)]

Residual part (solution of groups):
    ln gamma_i^R = sum_k nu_ki [ln Gamma_k - ln Gamma_k^(i)]
    ln Gamma_k  = Q_k [1 - ln(sum_m theta_m Psi_mk) - sum_m theta_m Psi_km / sum_n theta_n Psi_nm]
    Psi_mn      = exp(-a_mn / T)

Reference: Fredenslund, Jones & Prausnitz, AIChE J. 1975, 21, 1086;
parameters from Hansen et al., Ind. Eng. Chem. Res. 1991, 30, 2352.
"""
import warnings

import numpy as np

from .data.unifac_params import SUBGROUPS, interaction


class UNIFAC:
    """UNIFAC model for a fixed set of components.

    Parameters
    ----------
    components : list of Compound (anything with a `.groups` dict) or dicts
    missing    : what to do if a main-group pair has no published parameter:
                 'zero' (set a_mn = 0 and warn) or 'raise'.
    """

    def __init__(self, components, missing="zero"):
        groups_list = [c if isinstance(c, dict) else c.groups for c in components]
        self.names = [getattr(c, "name", f"comp{i}") for i, c in enumerate(components)]
        for gl in groups_list:
            for g in gl:
                if g not in SUBGROUPS:
                    raise KeyError(f"Subgroup '{g}' is not in the parameter table")
        self.subgroups = sorted({g for gl in groups_list for g in gl},
                                key=lambda g: SUBGROUPS[g].sid)
        ng, nc = len(self.subgroups), len(groups_list)
        self.nu = np.array([[gl.get(g, 0.0) for g in self.subgroups] for gl in groups_list],
                           dtype=float)                       # (nc, ng)
        self.R = np.array([SUBGROUPS[g].R for g in self.subgroups])
        self.Q = np.array([SUBGROUPS[g].Q for g in self.subgroups])
        self.main = [SUBGROUPS[g].main for g in self.subgroups]

        self.A = np.zeros((ng, ng))
        self.missing_pairs = set()
        for a, ma in enumerate(self.main):
            for b, mb in enumerate(self.main):
                val = interaction(ma, mb)
                if val is None:
                    self.missing_pairs.add(tuple(sorted((ma, mb))))
                    if missing == "raise":
                        raise ValueError(f"No UNIFAC parameter for main groups {ma}-{mb}")
                    val = 0.0
                self.A[a, b] = val
        if self.missing_pairs:
            warnings.warn(f"UNIFAC: no published parameters for main-group pairs "
                          f"{sorted(self.missing_pairs)}; set to 0.", stacklevel=2)

        self.r = self.nu @ self.R                              # (nc,)
        self.q = self.nu @ self.Q
        self._pure_cache = {}
        self.nc, self.ng = nc, ng

    # ------------------------------------------------------------------
    def _ln_Gamma(self, X, Psi):
        """Group residual activity coefficients for group mole fractions X."""
        theta = X * self.Q
        theta = theta / theta.sum()
        S = theta @ Psi                         # S_k = sum_m theta_m Psi_mk
        return self.Q * (1.0 - np.log(S) - (theta / S) @ Psi.T)

    def _pure_ln_Gamma(self, T):
        key = round(float(T), 10)
        if key not in self._pure_cache:
            Psi = np.exp(-self.A / T)
            out = np.zeros((self.nc, self.ng))
            for i in range(self.nc):
                Xi = self.nu[i] / self.nu[i].sum()
                out[i] = self._ln_Gamma(Xi, Psi)
            if len(self._pure_cache) > 2000:
                self._pure_cache.clear()
            self._pure_cache[key] = out
        return self._pure_cache[key]

    def ln_gamma_parts(self, x, T):
        """Return (ln gamma^C, ln gamma^R) for mole fractions x at T (K)."""
        x = np.asarray(x, dtype=float)
        x = x / x.sum()
        # combinatorial
        J = self.r / (x @ self.r)
        L = self.q / (x @ self.q)
        lnC = 1.0 - J + np.log(J) - 5.0 * self.q * (1.0 - J / L + np.log(J / L))
        # residual
        Psi = np.exp(-self.A / T)
        Xg = x @ self.nu
        Xg = Xg / Xg.sum()
        lnG = self._ln_Gamma(Xg, Psi)
        lnG_pure = self._pure_ln_Gamma(T)
        lnR = (self.nu * (lnG[None, :] - lnG_pure)).sum(axis=1)
        return lnC, lnR

    def ln_gamma(self, x, T):
        lnC, lnR = self.ln_gamma_parts(x, T)
        return lnC + lnR

    def gamma(self, x, T):
        return np.exp(self.ln_gamma(x, T))

    def ln_gamma_inf(self, i, T, eps=1e-10):
        """ln gamma of component i at infinite dilution in the others (equal parts)."""
        x = np.full(self.nc, (1.0 - eps) / (self.nc - 1))
        x[i] = eps
        return self.ln_gamma(x, T)[i]

    def g_mix(self, x, T):
        """Dimensionless Gibbs energy of mixing  dG_mix/(RT) = sum x_i ln(x_i gamma_i)."""
        x = np.asarray(x, dtype=float)
        x = x / x.sum()
        lg = self.ln_gamma(x, T)
        with np.errstate(divide="ignore", invalid="ignore"):
            terms = np.where(x > 0, x * (np.log(x) + lg), 0.0)
        return terms.sum()
