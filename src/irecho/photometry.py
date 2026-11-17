"""Magnitude/flux conversion, colour, and two-band blackbody temperature."""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from .config import F0_JY, LAMBDA_ISO_UM

LN10_OVER_2P5 = np.log(10.0) / 2.5  # d(flux)/flux = this * d(mag)

# CODATA 2018 exact values (SI).
H_PLANCK = 6.62607015e-34
C_LIGHT = 2.99792458e8
K_BOLTZ = 1.380649e-23


def mag_to_flux_jy(mag, band: str):
    """Vega magnitude -> flux density in Jy (F_nu = const zero point, no colour
    correction). All-Sky Explanatory Supplement IV.4.h Eq. 1."""
    return F0_JY[band] * 10.0 ** (-0.4 * np.asarray(mag, dtype=float))


def flux_jy_to_mag(flux, band: str):
    flux = np.asarray(flux, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        return -2.5 * np.log10(flux / F0_JY[band])


def mag_err_to_flux_err(flux, mag_err):
    return np.asarray(flux, dtype=float) * LN10_OVER_2P5 * np.asarray(mag_err, dtype=float)


def flux_err_to_mag_err(flux, flux_err):
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.asarray(flux_err, dtype=float) / np.asarray(flux, dtype=float) / LN10_OVER_2P5


def color_w1_w2(w1mag, w1err, w2mag, w2err):
    """W1-W2 and its uncertainty (independent errors added in quadrature)."""
    color = np.asarray(w1mag, dtype=float) - np.asarray(w2mag, dtype=float)
    err = np.hypot(np.asarray(w1err, dtype=float), np.asarray(w2err, dtype=float))
    return color, err


def planck_nu(wavelength_um: float, temperature_k: float) -> float:
    """B_nu(T) in SI units at the frequency c/lambda."""
    nu = C_LIGHT / (wavelength_um * 1e-6)
    x = H_PLANCK * nu / (K_BOLTZ * temperature_k)
    return 2.0 * H_PLANCK * nu**3 / C_LIGHT**2 / np.expm1(x)


def blackbody_flux_ratio(temperature_k: float) -> float:
    """F_nu(W1)/F_nu(W2) of a blackbody, evaluated at the isophotal wavelengths."""
    return planck_nu(LAMBDA_ISO_UM["w1"], temperature_k) / planck_nu(
        LAMBDA_ISO_UM["w2"], temperature_k
    )


def blackbody_temperature(flux_ratio_w1_w2: float, t_min: float = 150.0,
                          t_max: float = 20000.0) -> float:
    """Colour temperature from F_nu(W1)/F_nu(W2).

    The ratio rises monotonically with T, so the root is unique. Returns NaN
    when the ratio is outside the range reachable between t_min and t_max
    (for example bluer than the Rayleigh-Jeans limit).
    This is a monochromatic estimate: no bandpass integration, no dust
    emissivity law, observed frame.
    """
    if not np.isfinite(flux_ratio_w1_w2) or flux_ratio_w1_w2 <= 0:
        return float("nan")
    lo, hi = blackbody_flux_ratio(t_min), blackbody_flux_ratio(t_max)
    if not lo < flux_ratio_w1_w2 < hi:
        return float("nan")
    return float(brentq(lambda t: blackbody_flux_ratio(t) - flux_ratio_w1_w2, t_min, t_max))


def color_to_flux_ratio(color_w1_w2_mag: float) -> float:
    """F_nu(W1)/F_nu(W2) for a Vega colour W1-W2."""
    return F0_JY["w1"] / F0_JY["w2"] * 10.0 ** (-0.4 * color_w1_w2_mag)
