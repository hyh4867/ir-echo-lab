"""Physical quantities of the variable component (docs/THEORY.md section 7).

The variable flux dF = F(visit) - F(baseline) in W1 and W2 is modelled as a
single-temperature blackbody. With N = dF_nu / B_nu(T_obs) (an effective solid
angle), the bolometric flux is N sigma T_obs^4 / pi and

    L_dust = 4 pi D_L^2 * N sigma T_obs^4 / pi
    R_bb   = sqrt(L_dust / (4 pi sigma T_rest^4)),   T_rest = (1 + z) T_obs.

R_bb is the radius of a sphere that would emit L_dust as a perfect blackbody.
Real dust has emissivity and covering factor below one, so R_bb is a lower
limit on the size of the emitting region.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad

from .photometry import C_LIGHT, blackbody_temperature, planck_nu
from .config import LAMBDA_ISO_UM

SIGMA_SB = 5.670374419e-8  # W m^-2 K^-4 (CODATA 2018)
MPC_M = 3.0856775814913673e22
LIGHT_DAY_M = C_LIGHT * 86400.0
JY = 1e-26  # W m^-2 Hz^-1
ERG_S_PER_W = 1e7

# Cosmology used by Son et al. (2022), so that luminosities are comparable.
H0_KM_S_MPC = 67.4
OMEGA_M = 0.315

# Flux corrections f_c for blackbodies, Wright et al. (2010, AJ 140, 1868) Table 1:
# the monochromatic flux density at the isophotal wavelength is F*_nu / f_c.
WRIGHT2010_FC = {
    # T (K): (f_c W1, f_c W2)
    400.0: (1.1316, 1.0229),
    566.0: (1.0263, 0.9919),
    800.0: (0.9884, 0.9853),
    1131.0: (0.9801, 0.9877),
}


def luminosity_distance_mpc(z: float, h0: float = H0_KM_S_MPC, omega_m: float = OMEGA_M) -> float:
    """Flat LambdaCDM luminosity distance (no peculiar-velocity correction)."""
    integrand = lambda x: 1.0 / np.sqrt(omega_m * (1.0 + x) ** 3 + 1.0 - omega_m)  # noqa: E731
    d_c = (C_LIGHT / 1e3 / h0) * quad(integrand, 0.0, z)[0]
    return (1.0 + z) * d_c


def flux_correction(temperature_k: float, band: str) -> float:
    """f_c interpolated linearly in log T between the tabulated blackbodies
    (clamped to 400-1131 K)."""
    temps = np.array(sorted(WRIGHT2010_FC))
    col = 0 if band == "w1" else 1
    vals = np.array([WRIGHT2010_FC[t][col] for t in temps])
    return float(np.interp(np.log(np.clip(temperature_k, temps[0], temps[-1])),
                           np.log(temps), vals))


def blackbody_quantities(dflux_w1_jy: float, dflux_w2_jy: float, z: float,
                         d_l_mpc: float, color_correct: bool = False) -> dict:
    """T, L and R_bb of a blackbody that reproduces the two excess fluxes.
    Returns NaNs when the excess is not positive in both bands or its colour is
    outside the blackbody range."""
    nan = {"t_obs_k": np.nan, "t_rest_k": np.nan, "l_dust_erg_s": np.nan, "r_bb_light_days": np.nan}
    if not (np.isfinite(dflux_w1_jy) and np.isfinite(dflux_w2_jy)
            and dflux_w1_jy > 0 and dflux_w2_jy > 0):
        return nan
    f1, f2 = dflux_w1_jy, dflux_w2_jy
    t_obs = blackbody_temperature(f1 / f2)
    if color_correct and np.isfinite(t_obs):
        for _ in range(4):  # fixed-point iteration; converges in 2-3 steps
            c1, c2 = flux_correction(t_obs, "w1"), flux_correction(t_obs, "w2")
            t_obs = blackbody_temperature((dflux_w1_jy / c1) / (dflux_w2_jy / c2))
            if not np.isfinite(t_obs):
                return nan
        f1, f2 = dflux_w1_jy / c1, dflux_w2_jy / c2
    if not np.isfinite(t_obs):
        return nan
    solid = f2 * JY / planck_nu(LAMBDA_ISO_UM["w2"], t_obs)  # N, steradian
    f_bol = solid * SIGMA_SB * t_obs**4 / np.pi  # W m^-2
    l_w = 4.0 * np.pi * (d_l_mpc * MPC_M) ** 2 * f_bol
    t_rest = (1.0 + z) * t_obs
    r_m = np.sqrt(l_w / (4.0 * np.pi * SIGMA_SB * t_rest**4))
    return {"t_obs_k": float(t_obs), "t_rest_k": float(t_rest),
            "l_dust_erg_s": float(l_w * ERG_S_PER_W), "r_bb_light_days": float(r_m / LIGHT_DAY_M)}


def monte_carlo(flux: dict[str, float], flux_err: dict[str, float],
                base: dict[str, float], base_err: dict[str, float], z: float,
                d_l_mpc: float, n_draws: int, rng: np.random.Generator,
                base_draws: dict[str, np.ndarray] | None = None) -> dict:
    """Propagate visit and baseline flux errors. ``base_draws`` lets several
    visits share the same baseline realisations (their errors are correlated
    through the baseline)."""
    f = {b: rng.normal(flux[b], flux_err[b], n_draws) for b in ("w1", "w2")}
    if base_draws is None:
        base_draws = {b: rng.normal(base[b], base_err[b], n_draws) for b in ("w1", "w2")}
    keys = ("t_rest_k", "l_dust_erg_s", "r_bb_light_days")
    draws = {k: np.full(n_draws, np.nan) for k in keys}
    for i in range(n_draws):
        q = blackbody_quantities(f["w1"][i] - base_draws["w1"][i],
                                 f["w2"][i] - base_draws["w2"][i], z, d_l_mpc)
        for k in keys:
            draws[k][i] = q[k]
    out = {"valid_draw_fraction": float(np.isfinite(draws["t_rest_k"]).mean())}
    for k in keys:
        good = draws[k][np.isfinite(draws[k])]
        lo, med, hi = np.percentile(good, [16, 50, 84]) if good.size else (np.nan,) * 3
        out[f"{k}_p16"], out[f"{k}_p50"], out[f"{k}_p84"] = float(lo), float(med), float(hi)
    return out
