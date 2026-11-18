"""Stage 8: validation of the hysteresis statistic (docs/VALIDATION_8.md).

1. Null light curves without any colour loop, built on the real visit times,
   errors and W2 light curve of every object (section 3.1).
2. A direct comparison of the two branches that uses no reference visits, no
   shared peak and no blackbody: the offset in W1 flux between the rising and
   the declining branch at the same W2 flux (section 3.2).
3. Driver prescriptions for the lag-model re-test (section 3.3).
"""

from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq, least_squares

from . import breath

DUST_FRACTION = 0.5          # share of dust in the reference W2 flux (null b)
MIN_EACH_DIRECT = 3          # visits per branch inside the common W2 range
MIN_WIDTH_SIGMA = 3.0        # common range >= this many median W2 errors
MAX_SIGMA_COLOUR = 0.03      # mag; 'directly comparable'
LN10_OVER_2P5 = math.log(10.0) / 2.5


# --- branches ------------------------------------------------------------------------------
def branch_rising(f1, f2) -> np.ndarray:
    """Rising (True) or declining visit: sign of the central difference of the
    normalised sum of the two bands. Same definition as breath._prepare."""
    f1, f2 = np.asarray(f1, dtype=float), np.asarray(f2, dtype=float)
    z = f1 / np.median(f1) + f2 / np.median(f2)
    dz = np.empty(z.size)
    dz[1:-1] = z[2:] - z[:-2]
    dz[0], dz[-1] = z[1] - z[0], z[-1] - z[-2]
    return dz > 0


# --- straight line with errors in both variables -------------------------------------------
def eiv_fit(x, y, sx, sy, r: float = 0.0, extra=None) -> dict | None:
    """y = c + m x + sum_j g_j extra_j with Gaussian errors in x and y.

    Effective-variance chi2: residual variance sy^2 + m^2 sx^2 - 2 r m sx sy
    (r = correlation of the x and y errors of one point). ``extra`` is a list
    of indicator columns measured without error. Returns the parameters in the
    order [c, m, g_1, ...], their covariance, chi2 and dof.
    """
    x, y, sx, sy = (np.asarray(v, dtype=float) for v in (x, y, sx, sy))
    cols = [np.ones(x.size), x] + [np.asarray(e, dtype=float) for e in (extra or [])]
    a = np.vstack(cols).T
    n_par = a.shape[1]
    if x.size <= n_par:
        return None
    xs, ys = float(np.median(np.abs(x))) or 1.0, float(np.median(np.abs(y))) or 1.0
    w0 = 1.0 / sy
    try:
        theta0 = np.linalg.lstsq(a * w0[:, None], y * w0, rcond=None)[0]
    except np.linalg.LinAlgError:
        return None

    def resid(theta):
        m = theta[1]
        var = sy ** 2 + (m * sx) ** 2 - 2.0 * r * m * sx * sy
        return (y - a @ theta) / np.sqrt(np.maximum(var, 1e-300))

    scale = np.array([ys, ys / xs] + [ys] * (n_par - 2))
    try:
        sol = least_squares(resid, theta0, x_scale=np.maximum(np.abs(scale), 1e-300))
        cov = np.linalg.inv(sol.jac.T @ sol.jac)
    except (ValueError, np.linalg.LinAlgError, FloatingPointError):
        return None
    if not np.all(np.isfinite(cov)) or np.any(np.diag(cov) <= 0):
        return None
    return {"theta": sol.x, "cov": cov, "chi2": float(2.0 * sol.cost), "dof": int(x.size - n_par)}


# --- null light curves ---------------------------------------------------------------------
def truth_linear(f1, s1, f2, s2) -> np.ndarray | None:
    """Null (a): W1 on a straight line of W2, fitted to the observed visits."""
    fit = eiv_fit(f2, f1, s2, s1)
    if fit is None:
        return None
    c, m = fit["theta"][:2]
    return c + m * np.asarray(f2, dtype=float)


def truth_single_temperature(f1, s1, f2, s2, dust_fraction: float = DUST_FRACTION) -> np.ndarray | None:
    """Null (b): the varying component is a blackbody of fixed area whose
    temperature alone changes (heating in place), so W1 is a single-valued
    function of W2. The reference temperature is set so that the summed W1
    excess of the brightened visits equals the observed one (700 K if no root)."""
    prep = breath._prepare(f1, s1, f2, s2, breath.SETTINGS)
    if prep is None:
        return None
    a1, a2 = prep["a"]
    f1, f2 = np.asarray(f1, dtype=float), np.asarray(f2, dtype=float)
    dust = dust_fraction * a2
    e2 = np.maximum(f2 - a2, -0.9 * dust)
    pos = e2 > 0
    if pos.sum() < 3:
        return None
    target = float(np.sum(f1[pos] - a1))

    def miss(t_ref: float) -> float:
        return float(np.sum(breath.heating_in_place_excess(e2[pos], dust, t_ref))) - target

    t_ref = 700.0
    try:
        if miss(250.0) < 0 < miss(3000.0):
            t_ref = brentq(miss, 250.0, 3000.0, xtol=0.5)
    except (ValueError, RuntimeError):
        pass
    return a1 + breath.heating_in_place_excess(e2, dust, t_ref)


def smooth3(x) -> np.ndarray:
    """Three-visit running mean (end points: two-visit mean)."""
    x = np.asarray(x, dtype=float)
    out = x.copy()
    out[1:-1] = (x[:-2] + x[1:-1] + x[2:]) / 3.0
    out[0], out[-1] = 0.5 * (x[0] + x[1]), 0.5 * (x[-1] + x[-2])
    return out


def add_noise(rng, t1, s1, t2, s2, r: float = 0.0, bands: str = "both"):
    """Bivariate Gaussian visit noise with correlation r between the bands."""
    n = np.asarray(t1).size
    z1, z2 = rng.standard_normal(n), rng.standard_normal(n)
    n1, n2 = z1, r * z1 + math.sqrt(max(1.0 - r * r, 0.0)) * z2
    k1 = 0.0 if bands == "w2" else 1.0
    k2 = 0.0 if bands == "w1" else 1.0
    return np.asarray(t1) + k1 * np.asarray(s1) * n1, np.asarray(t2) + k2 * np.asarray(s2) * n2


# --- the 5b statistic with switches (cause isolation) ----------------------------------------
def d_statistic(f1, s1, f2, s2, pivot: str = "peak", reference=None, rising=None,
                settings: breath.SlopeSettings = breath.SETTINGS) -> dict:
    """D = beta_decline - beta_rise as in breath.measure, with three switches.

    pivot 'peak' (5b: temperature shared at the largest W2 excess) or 'mean'.
    reference: (a1, a2) to use instead of the three faintest visits.
    rising: branch labels to use instead of the ones from the noisy data.
    With the defaults the result equals breath.measure.
    """
    out = {"d": np.nan, "sigma_d": np.nan, "measurable": False, "n_rising": 0, "n_declining": 0}
    prep = breath._prepare(f1, s1, f2, s2, settings)
    if prep is None:
        return out
    if reference is not None:
        a1, a2 = reference
        sa1, sa2 = 1e-3 * abs(a1), 1e-3 * abs(a2)
        e1, e2 = prep["f1"] - a1, prep["f2"] - a2
        snr = (e1 / np.hypot(prep["s1"], sa1) + e2 / np.hypot(prep["s2"], sa2)) / np.sqrt(2.0)
        prep = dict(prep, a=(a1, a2), sa=(sa1, sa2), used=(snr >= settings.min_snr) & (e1 > 0) & (e2 > 0))
    if rising is not None:
        prep = dict(prep, rising=np.asarray(rising, dtype=bool))
    u = prep["used"]
    if u.sum() < settings.min_used:
        return out
    lab = prep["rising"][u]
    out["n_rising"], out["n_declining"] = int(lab.sum()), int((~lab).sum())
    if min(out["n_rising"], out["n_declining"]) < settings.hyst_min_each:
        return out
    e2 = prep["f2"][u] - prep["a"][1]
    x_pivot = float(np.log(e2.max())) if pivot == "peak" else float(np.mean(np.log(e2)))
    fit = breath._fit(prep, np.where(lab, 0, 1), 2, x_pivot, settings)
    if fit is None:
        return out
    c = fit["cov"]
    out["d"] = float(fit["theta"][2] - fit["theta"][1])
    out["sigma_d"] = float(settings.sigma_scale * np.sqrt(max(c[1, 1] + c[2, 2] - 2.0 * c[1, 2], 0.0)))
    out["measurable"] = bool(0 < out["sigma_d"] <= settings.hyst_max_sigma)
    e2r, e2d = e2[lab], e2[~lab]
    out["range_ratio"] = float((e2r.max() - e2r.min()) / max(e2d.max() - e2d.min(), 1e-300))
    return out


# --- direct comparison (section 3.2) ----------------------------------------------------------
def direct_offset(f1, s1, f2, s2, rising=None, r: float = 0.0, parity=None, trend=None) -> dict:
    """Offset of W1 flux between the rising and the declining branch at the same W2 flux.

    Only visits inside the common W2 range of the two branches are used (no
    extrapolation). Model: F1 = c + m F2 + delta [rising] (+ gamma [parity]).
    delta > 0: the rising branch is bluer; this is the sign of D > 0.
    """
    f1, s1, f2, s2 = (np.asarray(v, dtype=float) for v in (f1, s1, f2, s2))
    rising = branch_rising(f1, f2) if rising is None else np.asarray(rising, dtype=bool)
    out = {"usable": False, "comparable": False, "n_rise_in": 0, "n_fade_in": 0, "lo": np.nan, "hi": np.nan,
           "delta": np.nan, "sigma_delta": np.nan, "d_colour": np.nan, "sigma_colour": np.nan,
           "slope": np.nan, "gamma": np.nan, "sigma_gamma": np.nan, "delta_two_lines": np.nan,
           "sigma_two_lines": np.nan, "n_rise": int(rising.sum()), "n_fade": int((~rising).sum())}
    if rising.sum() < 1 or (~rising).sum() < 1:
        return out
    lo = max(f2[rising].min(), f2[~rising].min())
    hi = min(f2[rising].max(), f2[~rising].max())
    inside = (f2 >= lo) & (f2 <= hi)
    out.update(lo=float(lo), hi=float(hi), n_rise_in=int((inside & rising).sum()), n_fade_in=int((inside & ~rising).sum()))
    if min(out["n_rise_in"], out["n_fade_in"]) < MIN_EACH_DIRECT or hi <= lo:
        return out
    if hi - lo < MIN_WIDTH_SIGMA * float(np.median(s2[inside])):
        return out
    x, y, sx, sy, lab = f2[inside], f1[inside], s2[inside], s1[inside], rising[inside].astype(float)
    extra = [lab]
    if parity is not None:
        par = np.asarray(parity, dtype=float)[inside]
        if 0 < par.sum() < par.size:
            extra.append(par)
    if trend is not None:      # a continuous covariate, e.g. time (added after the plan; see the stage 8 report)
        tr = np.asarray(trend, dtype=float)[inside]
        extra.append((tr - tr.mean()) * float(np.median(y)))
        out["mean_trend_rise_minus_fade"] = float(tr[lab > 0.5].mean() - tr[lab < 0.5].mean())
    fit = eiv_fit(x, y, sx, sy, r, extra)
    if fit is None:
        return out
    c, m, delta = fit["theta"][:3]
    sig = float(np.sqrt(fit["cov"][2, 2]))
    x0 = 0.5 * (lo + hi)
    f_fade = c + m * x0
    f_rise = f_fade + delta
    if f_fade <= 0 or f_rise <= 0:
        return out
    out.update(usable=True, delta=float(delta), sigma_delta=sig, slope=float(m),
               d_colour=float(-2.5 * np.log10(f_rise / f_fade)),
               sigma_colour=float(sig / (LN10_OVER_2P5 * 0.5 * (f_rise + f_fade))),
               chi2=fit["chi2"], dof=fit["dof"], frac_delta=float(delta / (0.5 * (f_rise + f_fade))))
    if len(extra) == 2:
        out["gamma"], out["sigma_gamma"] = float(fit["theta"][3]), float(np.sqrt(fit["cov"][3, 3]))
        # with ``trend``: gamma is the fractional change of W1 at fixed W2 per unit of the covariate
    out["comparable"] = bool(out["sigma_colour"] <= MAX_SIGMA_COLOUR)
    # auxiliary: one line per branch, compared at the centre of the common range
    vals = []
    for sel in (lab > 0.5, lab < 0.5):
        fb = eiv_fit(x[sel], y[sel], sx[sel], sy[sel], r) if sel.sum() >= 3 else None
        if fb is None:
            vals = []
            break
        j = np.array([1.0, x0])
        vals.append((float(j @ fb["theta"]), float(j @ fb["cov"] @ j)))
    if len(vals) == 2:
        out["delta_two_lines"] = vals[0][0] - vals[1][0]
        out["sigma_two_lines"] = math.sqrt(vals[0][1] + vals[1][1])
    return out


# --- drivers for the lag-model re-test (section 3.3) -----------------------------------------
MACLEOD = {"sf_inf": (-0.51, -0.479, 0.131, 0.18), "tau": (2.4, 0.17, 0.03, 0.21)}  # MacLeod et al. 2010, Table 1


def absolute_i_mag(log_lbol: float) -> float:
    """K-corrected (z = 0) i-band absolute magnitude from L_bol (erg/s):
    M_i(z=2) = 90 - 2.5 log L_bol (Shen et al. 2009), M_i(z=0) = M_i(z=2) + 0.596."""
    return 90.596 - 2.5 * log_lbol


def macleod_drw(log_lbol: float, log_mbh: float, lambda_rf_angstrom: float = 4000.0) -> tuple[float, float]:
    """(sigma in mag, rest-frame tau in days) of the damped random walk, MacLeod et al. (2010) Eq. 7."""
    m_i = absolute_i_mag(log_lbol)
    out = []
    for key in ("sf_inf", "tau"):
        a, b, c, d = MACLEOD[key]
        out.append(10.0 ** (a + b * math.log10(lambda_rf_angstrom / 4000.0) + c * (m_i + 23.0) + d * (log_mbh - 9.0)))
    return out[0] / math.sqrt(2.0), out[1]


def burke_tau(log_mbh: float) -> float:
    """Rest-frame damping time in days, Burke et al. (2021): 107 d (M_BH / 1e8 Msun)^0.38."""
    return 107.0 * 10.0 ** (0.38 * (log_mbh - 8.0))
