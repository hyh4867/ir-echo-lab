"""Order-of-magnitude relations used in docs/THEORY.md (sources cited there)."""

from __future__ import annotations

import numpy as np

PC_IN_LIGHT_DAYS = 3.0856775814913673e16 / (2.99792458e8 * 86400.0)  # 1191.3
L_SUN_ERG_S = 3.828e33  # IAU 2015 nominal solar luminosity


def sublimation_radius_pc(l_agn_erg_s: float, grain: str = "graphite",
                          t_sub_k: float | None = None, f_theta: float = 1.0) -> float:
    """Mean dust sublimation radius, Netzer (2015, ARA&A 53, 365) Eq. 1:

        R_sub,C  = 0.5 L46^(1/2) (1800 K / T_sub)^2.6 f(theta) pc   (graphite)
        R_sub,Si = 1.3 L46^(1/2) (1500 K / T_sub)^2.6 f(theta) pc   (silicate)

    with L46 = L_AGN / 1e46 erg/s, for a mean grain size of 0.05 micron.
    """
    coeff, t_ref = {"graphite": (0.5, 1800.0), "silicate": (1.3, 1500.0)}[grain]
    t_sub = t_ref if t_sub_k is None else t_sub_k
    return coeff * np.sqrt(l_agn_erg_s / 1e46) * (t_ref / t_sub) ** 2.6 * f_theta


def light_travel_days(radius_pc: float) -> float:
    """Light-crossing time R/c in days."""
    return radius_pc * PC_IN_LIGHT_DAYS


def lyu2019_lag_days(l_agn_erg_s: float, band: str = "w1") -> float:
    """Empirical mid-IR lag of PG quasars, Lyu, Rieke & Smith (2019, ApJ 886, 33)
    Eqs. 15-16: lag/day = 10^2.10 (L/1e11 Lsun)^0.47 (W1), 10^2.20 (...)^0.45 (W2).
    Calibrated on quasars with L_AGN ~ 1e11-1e13 Lsun; lower luminosities are an
    extrapolation."""
    log_norm, slope = {"w1": (2.10, 0.47), "w2": (2.20, 0.45)}[band]
    return 10.0**log_norm * (l_agn_erg_s / L_SUN_ERG_S / 1e11) ** slope
