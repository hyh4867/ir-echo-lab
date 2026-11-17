"""Small statistics helpers: robust scatter, reduced chi-square, error floor."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq
from scipy.stats import chi2 as chi2_dist

MAD_TO_SIGMA = 1.4826  # Gaussian: sigma = 1.4826 * MAD


def robust_std(values) -> float:
    """1.4826 * median absolute deviation (NaN for fewer than 2 values)."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size < 2:
        return float("nan")
    return float(MAD_TO_SIGMA * np.median(np.abs(v - np.median(v))))


@dataclass(frozen=True)
class Chi2Result:
    n: int
    dof: int
    weighted_mean: float
    chi2: float
    chi2_red: float
    p_value: float  # P(chi2 >= observed | constant source, errors correct)
    rms: float  # unweighted standard deviation (ddof=1)
    median_err: float


def constant_chi2(values, errors) -> Chi2Result:
    """Chi-square of the hypothesis 'the source is constant' given the errors.

    The constant is the inverse-variance weighted mean, so dof = n - 1.
    """
    v = np.asarray(values, dtype=float)
    e = np.asarray(errors, dtype=float)
    ok = np.isfinite(v) & np.isfinite(e) & (e > 0)
    v, e = v[ok], e[ok]
    n = int(v.size)
    if n < 2:
        nan = float("nan")
        return Chi2Result(n, max(n - 1, 0), nan, nan, nan, nan, nan, nan)
    w = 1.0 / e**2
    mean = float(np.sum(w * v) / np.sum(w))
    chi2 = float(np.sum(w * (v - mean) ** 2))
    dof = n - 1
    return Chi2Result(
        n=n, dof=dof, weighted_mean=mean, chi2=chi2, chi2_red=chi2 / dof,
        p_value=float(chi2_dist.sf(chi2, dof)),
        rms=float(np.std(v, ddof=1)), median_err=float(np.median(e)),
    )


def error_floor(values, errors, max_floor: float = 1.0) -> float:
    """Extra error sigma_sys (added in quadrature) that makes chi2_red = 1.

    Returns 0 when the quoted errors already give chi2_red <= 1.
    """
    v = np.asarray(values, dtype=float)
    e = np.asarray(errors, dtype=float)
    ok = np.isfinite(v) & np.isfinite(e) & (e > 0)
    v, e = v[ok], e[ok]
    if v.size < 2:
        return float("nan")

    def excess(floor: float) -> float:
        return constant_chi2(v, np.hypot(e, floor)).chi2_red - 1.0

    if excess(0.0) <= 0:
        return 0.0
    if excess(max_floor) > 0:
        return float("nan")
    return float(brentq(excess, 0.0, max_floor, xtol=1e-6))


def pooled_error_floor(value_sets, error_sets, max_floor: float = 1.0) -> float:
    """One sigma_sys for many sources: total chi2 equals total dof."""
    pairs = []
    for v, e in zip(value_sets, error_sets):
        v = np.asarray(v, dtype=float)
        e = np.asarray(e, dtype=float)
        ok = np.isfinite(v) & np.isfinite(e) & (e > 0)
        if ok.sum() >= 2:
            pairs.append((v[ok], e[ok]))
    if not pairs:
        return float("nan")

    def excess(floor: float) -> float:
        chi2 = dof = 0.0
        for v, e in pairs:
            r = constant_chi2(v, np.hypot(e, floor))
            chi2 += r.chi2
            dof += r.dof
        return chi2 / dof - 1.0

    if excess(0.0) <= 0:
        return 0.0
    if excess(max_floor) > 0:
        return float("nan")
    return float(brentq(excess, 0.0, max_floor, xtol=1e-6))
