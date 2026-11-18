"""Stage 9: band-to-band transfer W1 -> W2 without a driver (docs/PREDICTIONS_9.md).

Model for one object, at the W2 visits that can be predicted without extrapolation:

    F_W2(t_i) = a + b * (K * F_W1)(t_i) [+ eta * rising_i]

K is a top-hat of width ``w`` centred on the lag ``tau`` (w = 0: a pure shift).
F_W1 between visits is the straight line through the two neighbouring visits
(linear interpolation), so (K * F_W1)(t_i) is a weighted sum of the W1 visit
fluxes, sum_j c_ij F1_j, with weights that are computed exactly. The W1 errors
therefore propagate linearly into the prediction and the covariance of the
residuals is known (``residual_covariance``).

The ladder of fitted models (section 2 of the plan):
    T0  tau = 0, w = 0           (a, b)
    T1  tau free, w = 0          (a, b, tau)
    T2  tau and w free           (a, b, tau, w)
    T3  T1 or T2 plus eta        (a branch offset that a fixed kernel does not give)

The second half of the module makes the synthetic light curves of section 3
(truth models T0, T1, T2, T3 seen with the real visit times and errors).
"""

from __future__ import annotations

import math

import numpy as np
from scipy.signal import lfilter

from . import breath, hyst8

# --- numbers fixed in docs/PREDICTIONS_9.md before any real light curve was fitted -------------
SUPPORT_PAST_DAYS = 350.0      # the kernel may reach this far into the past of a W2 visit
SUPPORT_FUTURE_DAYS = 150.0    # and this far into its future (negative lags are allowed)
MAX_GAP_DAYS = 250.0           # no interpolation across a longer gap between W1 visits
TAU_STEP_DAYS = 5.0
WIDTHS_DAYS = (0.0, 100.0, 200.0, 300.0)
MIN_VISITS_FIT = 10            # W2 visits left after the no-extrapolation rule
N_ITERATIONS = 3               # updates of b in the covariance (b starts at 0)
WELL_DS = 9.0                  # 'well determined lag': T1 better than T0 by this much ...
WELL_HALF_WIDTH_DAYS = 45.0    # ... and the interval of tau no wider than this (half width)
DUST_FRACTION = 0.5            # variable share of the flux in the linear truth models
GRID_STEP_DAYS = 5.0           # time step of the synthetic driver
BURN_IN_DAYS = 3000.0
TAIL_DAYS = 400.0               # the driver grid runs this far past the last visit (negative delays)


# --- linear interpolation and its exact average over a window ------------------------------------
def interp_weights(t, u) -> np.ndarray:
    """Weights W (len(u) x len(t)) with F_lin(u) = W @ F for the linear interpolant.
    ``u`` must lie inside [t[0], t[-1]]; outside, an exception is raised (no extrapolation)."""
    t, u = np.asarray(t, dtype=float), np.atleast_1d(np.asarray(u, dtype=float))
    if np.any(u < t[0] - 1e-9) or np.any(u > t[-1] + 1e-9):
        raise ValueError("extrapolation is not allowed")
    u = np.clip(u, t[0], t[-1])
    k = np.clip(np.searchsorted(t, u, side="right") - 1, 0, t.size - 2)
    x = (u - t[k]) / (t[k + 1] - t[k])
    w = np.zeros((u.size, t.size))
    rows = np.arange(u.size)
    w[rows, k] = 1.0 - x
    w[rows, k + 1] += x
    return w


def antiderivative_weights(t, u) -> np.ndarray:
    """Weights A (len(u) x len(t)) with integral_{t[0]}^{u} F_lin = A @ F (exact)."""
    t, u = np.asarray(t, dtype=float), np.atleast_1d(np.asarray(u, dtype=float))
    if np.any(u < t[0] - 1e-9) or np.any(u > t[-1] + 1e-9):
        raise ValueError("extrapolation is not allowed")
    u = np.clip(u, t[0], t[-1])
    h = np.diff(t)
    full = np.zeros((t.size, t.size))          # full[k] = weights of the integral up to node k
    for k in range(1, t.size):
        full[k] = full[k - 1]
        full[k, k - 1] += 0.5 * h[k - 1]
        full[k, k] += 0.5 * h[k - 1]
    k = np.clip(np.searchsorted(t, u, side="right") - 1, 0, t.size - 2)
    x = (u - t[k]) / h[k]
    a = full[k].copy()
    rows = np.arange(u.size)
    a[rows, k] += h[k] * (x - 0.5 * x * x)
    a[rows, k + 1] += h[k] * 0.5 * x * x
    return a


def kernel_weights(t, t_eval, tau: float, width: float = 0.0) -> np.ndarray:
    """c (len(t_eval) x len(t)): (K * F_lin)(t_eval) = c @ F for a top-hat of ``width`` centred on ``tau``."""
    t_eval = np.atleast_1d(np.asarray(t_eval, dtype=float))
    if width <= 0:
        return interp_weights(t, t_eval - tau)
    hi = antiderivative_weights(t, t_eval - tau + 0.5 * width)
    lo = antiderivative_weights(t, t_eval - tau - 0.5 * width)
    return (hi - lo) / width


def usable_visits(t, past: float = SUPPORT_PAST_DAYS, future: float = SUPPORT_FUTURE_DAYS,
                  max_gap: float = MAX_GAP_DAYS) -> np.ndarray:
    """W2 visits whose whole window [t - past, t + future] is covered by W1 visits
    without a gap longer than ``max_gap``. The same set is used by every model."""
    t = np.asarray(t, dtype=float)
    ok = (t - past >= t[0]) & (t + future <= t[-1])
    gaps = np.diff(t)
    for k in np.nonzero(gaps > max_gap)[0]:
        ok &= ~((t - past < t[k + 1]) & (t + future > t[k]))
    return ok


def grid(widths=WIDTHS_DAYS, step: float = TAU_STEP_DAYS) -> list[tuple[float, float]]:
    """(tau, width) pairs whose support lies inside [-SUPPORT_FUTURE_DAYS, SUPPORT_PAST_DAYS]."""
    out = []
    for w in widths:
        lo, hi = -SUPPORT_FUTURE_DAYS + 0.5 * w, SUPPORT_PAST_DAYS - 0.5 * w
        n = int(math.floor((hi - lo) / step + 1e-9)) + 1
        out.extend((float(lo + i * step), float(w)) for i in range(n))
    return out


# --- generalised least squares with the propagated W1 errors --------------------------------------
def residual_covariance(c, idx, s1, s2, r: float, b):
    """Covariance of R = F2[idx] - a - b c @ F1 for weights ``c`` (.. x m x n).

    Cov = diag(s2^2) + b^2 c S1 c^T - b (X + X^T),  X_ik = c_k,idx_i r s1_idx_i s2_idx_i,
    where r is the correlation of the W1 and W2 errors of one visit. ``b`` may be
    an array with one value per leading index of ``c``.
    """
    c = np.asarray(c, dtype=float)
    b = np.asarray(b, dtype=float)[..., None, None]
    m = c.shape[-2]
    cov = np.einsum("...ij,j,...kj->...ik", c, s1 ** 2, c) * b ** 2
    cross = np.swapaxes(c[..., :, idx], -1, -2) * (r * s1[idx] * s2[idx])[:, None]   # X_ik
    cov = cov - b * (cross + np.swapaxes(cross, -1, -2))
    cov[..., np.arange(m), np.arange(m)] += s2[idx] ** 2
    return cov


def covariance_parts(c, idx, s1, s2, r: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The three parts of ``residual_covariance`` that do not depend on b: diag, Q, X + X^T."""
    c = np.asarray(c, dtype=float)
    q = np.einsum("...ij,j,...kj->...ik", c, s1 ** 2, c)
    cross = np.swapaxes(c[..., :, idx], -1, -2) * (r * s1[idx] * s2[idx])[:, None]
    return s2[idx] ** 2, q, cross + np.swapaxes(cross, -1, -2)


def _solve(cov, a_mat, y):
    cinv = np.linalg.inv(cov)
    ca = cinv @ a_mat                                              # G x m x p
    normal = np.swapaxes(a_mat, -1, -2) @ ca
    rhs = np.einsum("gmp,m->gp", ca, y)
    theta = np.linalg.solve(normal, rhs[..., None])[..., 0]
    chi2 = np.einsum("m,gmk,k->g", y, cinv, y) - np.einsum("gp,gp->g", theta, rhs)
    return theta, chi2, normal


def gls(c, idx, f1, s1, f2, s2, r: float, extra=None, n_iter: int = N_ITERATIONS) -> dict:
    """Fit F2[idx] = a + b c @ F1 for a stack of weight matrices ``c`` (G x m x n).

    The covariance depends on b; b starts at 0 and is updated ``n_iter`` times.
    Returns arrays over the G grid points: the statistic S = chi2 + ln det C, chi2,
    the parameters, the prediction. With ``extra`` (one column, m values) the fit
    is repeated with that column added, using the covariance of the fit without
    it, so the two chi2 are nested: ``chi2_extra`` <= ``chi2`` at every grid point.
    """
    c = np.asarray(c, dtype=float)
    if c.ndim == 2:
        c = c[None]
    y = f2[idx]
    g, m = c.shape[0], c.shape[1]
    diag, q, xs = covariance_parts(c, idx, s1, s2, r)
    eye = np.eye(m) * diag
    pred = c @ f1                                                  # G x m
    a_mat = np.stack([np.ones((g, m)), pred], axis=-1)             # G x m x 2
    b = np.zeros(g)
    for _ in range(n_iter + 1):
        cov = eye + b[:, None, None] ** 2 * q - b[:, None, None] * xs
        theta, chi2, _ = _solve(cov, a_mat, y)
        b = theta[:, 1]
    cov = eye + b[:, None, None] ** 2 * q - b[:, None, None] * xs
    theta, chi2, normal = _solve(cov, a_mat, y)
    logdet = np.linalg.slogdet(cov)[1]
    out = {"s": chi2 + logdet, "chi2": chi2, "theta": theta, "pred": theta[:, :1] + theta[:, 1:2] * pred}
    if extra is not None:
        a3 = np.concatenate([a_mat, np.broadcast_to(np.asarray(extra, dtype=float), (g, m))[..., None]], axis=-1)
        theta3, chi3, normal3 = _solve(cov, a3, y)
        out.update(theta_extra=theta3, chi2_extra=chi3, theta_extra_cov=np.linalg.inv(normal3))
    return out


def prepare(t, grid_points=None) -> dict:
    """What depends on the visit times only: the visits used and the weights of every grid point."""
    t = np.asarray(t, dtype=float)
    use = usable_visits(t)
    idx = np.nonzero(use)[0]
    pts = grid() if grid_points is None else grid_points
    out = {"idx": idx, "taus": np.array([p[0] for p in pts]), "widths": np.array([p[1] for p in pts]), "c": None}
    if idx.size >= MIN_VISITS_FIT:
        out["c"] = np.stack([kernel_weights(t, t[idx], tau, w) for tau, w in pts])
    return out


def fit_object(t, f1, s1, f2, s2, r: float = 0.0, rising=None, prep: dict | None = None) -> dict:
    """The ladder T0, T1, T2, T3 for one object. Returns a flat dict (NaN where not usable)."""
    t, f1, s1, f2, s2 = (np.asarray(v, dtype=float) for v in (t, f1, s1, f2, s2))
    prep = prepare(t) if prep is None else prep
    idx = prep["idx"]
    out = {"usable": False, "n_used": int(idx.size)}
    if prep["c"] is None:
        return out
    rising = hyst8.branch_rising(f1, f2) if rising is None else np.asarray(rising, dtype=bool)
    lab = rising[idx].astype(float)
    taus, widths, c = prep["taus"], prep["widths"], prep["c"]
    scale = float(np.median(f2))
    has_eta = bool(0 < lab.sum() < lab.size)
    plain = gls(c, idx, f1, s1, f2, s2, r, extra=lab * scale if has_eta else None)
    shift = widths == 0
    i0 = int(np.nonzero(shift & (np.abs(taus) < 1e-9))[0][0])
    i1 = int(np.nonzero(shift)[0][np.argmin(plain["s"][shift])])
    i2 = int(np.argmin(plain["s"]))
    m = idx.size
    out.update(usable=True, idx=idx, rising_used=lab > 0.5,
               s_t0=float(plain["s"][i0]), s_t1=float(plain["s"][i1]), s_t2=float(plain["s"][i2]),
               chi2_t0=float(plain["chi2"][i0]), chi2_t1=float(plain["chi2"][i1]), chi2_t2=float(plain["chi2"][i2]),
               tau_t1=float(taus[i1]), tau_t2=float(taus[i2]), width_t2=float(widths[i2]),
               a_t0=float(plain["theta"][i0, 0]), b_t0=float(plain["theta"][i0, 1]),
               a_t1=float(plain["theta"][i1, 0]), b_t1=float(plain["theta"][i1, 1]),
               a_t2=float(plain["theta"][i2, 0]), b_t2=float(plain["theta"][i2, 1]),
               pred_t0=plain["pred"][i0], pred_t1=plain["pred"][i1], pred_t2=plain["pred"][i2],
               ds_01=float(plain["s"][i0] - plain["s"][i1]), ds_02=float(plain["s"][i0] - plain["s"][i2]),
               ds_12=float(plain["s"][i1] - plain["s"][i2]))
    # interval of the pure lag: where S rises by max(1, reduced chi2 at the best lag)
    tau_s, s_s = taus[shift], plain["s"][shift]
    level = max(1.0, out["chi2_t1"] / max(m - 3, 1))
    inside = s_s - s_s.min() <= level
    out.update(tau_lo=float(tau_s[inside].min()), tau_hi=float(tau_s[inside].max()),
               tau_at_edge=bool(out["tau_t1"] <= tau_s.min() + 1e-9 or out["tau_t1"] >= tau_s.max() - 1e-9))
    out["tau_half_width"] = 0.5 * (out["tau_hi"] - out["tau_lo"])
    out["well_determined"] = bool(out["ds_01"] >= WELL_DS and out["tau_half_width"] <= WELL_HALF_WIDTH_DAYS
                                  and not out["tau_at_edge"])
    for key in ("eta_t0", "eta_t1", "eta_t2", "sigma_eta_t0", "sigma_eta_t1", "sigma_eta_t2", "ds_eta_t0",
                "ds_eta_t1", "ds_eta_t2", "eta_colour_t0", "eta_colour_t1", "eta_colour_t2"):
        out[key] = float("nan")
    if has_eta:
        # T3: the branch offset eta added to T0, T1, T2; the lag (and width) is fitted again with eta.
        # S of the model with eta uses the covariance of the fit without it (same ln det).
        s_eta = plain["s"] - plain["chi2"] + plain["chi2_extra"]
        j0 = i0
        j1 = int(np.nonzero(shift)[0][np.argmin(s_eta[shift])])
        j2 = int(np.argmin(s_eta))
        for name, j, base in (("t0", j0, out["s_t0"]), ("t1", j1, out["s_t1"]), ("t2", j2, out["s_t2"])):
            eta = float(plain["theta_extra"][j, 2]) * scale
            out[f"eta_{name}"] = eta
            out[f"sigma_eta_{name}"] = float(np.sqrt(plain["theta_extra_cov"][j, 2, 2])) * scale
            out[f"ds_eta_{name}"] = float(base - s_eta[j])
            ratio = 1.0 + eta / scale
            out[f"eta_colour_{name}"] = float(2.5 * np.log10(ratio)) if ratio > 0 else float("nan")
    return out


def offsets(f1, s1, f2, s2, r: float, fit: dict, rising=None) -> dict:
    """Direct offset of stage 8 on the visits used by the fit: measured (observed W2) and
    predicted (W2 replaced by the fitted model), with the same branch labels and the same
    visits inside the common W2 range. d_colour < 0: the rising branch is bluer."""
    out = {}
    keys = ("delta", "sigma_delta", "d_colour", "sigma_colour", "usable", "comparable")
    for name in ("obs", "t0", "t1", "t2"):
        for k in keys:
            out[f"{k}_{name}"] = float("nan") if k not in ("usable", "comparable") else False
    if not fit.get("usable"):
        return out
    idx = fit["idx"]
    rising = hyst8.branch_rising(f1, f2) if rising is None else np.asarray(rising, dtype=bool)
    lab = rising[idx]
    obs = hyst8.direct_offset(f1[idx], s1[idx], f2[idx], s2[idx], rising=lab, r=r)
    for k in keys:
        out[f"{k}_obs"] = obs[k]
    if not obs["usable"]:
        return out
    inside = (f2[idx] >= obs["lo"]) & (f2[idx] <= obs["hi"])
    for name in ("t0", "t1", "t2"):
        d = hyst8.direct_offset(f1[idx], s1[idx], fit[f"pred_{name}"], s2[idx], rising=lab, r=r, inside=inside)
        for k in keys:
            out[f"{k}_{name}"] = d[k]
    return out


# --- synthetic light curves (section 3 of the plan) -----------------------------------------------
def drw_magnitudes(rng, n: int, dt: float, tau_days: float) -> np.ndarray:
    """Unit-variance damped random walk on an even grid."""
    a = math.exp(-dt / tau_days)
    eps = rng.normal(size=n)
    eps[0] /= math.sqrt(1.0 - a * a)
    return lfilter([math.sqrt(1.0 - a * a)], [1.0, -a], eps)


def boxcar_response(ell, dt: float, lag_days: float) -> np.ndarray:
    """Boxcar of width ``lag_days`` with mean delay ``lag_days`` (as breath._respond)."""
    return breath._respond(np.arange(ell.size) * dt, ell, lag_days)


def shifted(x, dt: float, shift_days: float) -> np.ndarray:
    """x(t - shift) on the same grid (linear interpolation; the first value fills the start)."""
    t = np.arange(x.size) * dt
    return np.interp(t - shift_days, t, x, left=x[0], right=x[-1])


def exponential_tail(x, dt: float, tail_days: float) -> np.ndarray:
    """x convolved with a one-sided exponential of e-folding time ``tail_days``."""
    a = math.exp(-dt / tail_days)
    y = lfilter([1.0 - a], [1.0, -a], x - x[0]) + x[0]
    return y


TRUTHS = {
    # name: (kind, extra W2 delay in days or 'lyu', tau_re in years)
    "T0": ("linear", 0.0, None),
    "T1_20": ("linear", 20.0, None),
    "T1_40": ("linear", 40.0, None),
    "T1_80": ("linear", 80.0, None),
    "T1_160": ("linear", 160.0, None),
    "T1_lyu": ("linear", "lyu", None),
    "T2_40": ("tail", 40.0, None),
    "T2_80": ("tail", 80.0, None),
    "T0b": ("boundary", 0.0, 0.0),
    "T3_1yr": ("boundary", 0.0, 1.0),
    "T3_3yr": ("boundary", 0.0, 3.0),
}
TAIL_SHARE = 0.5     # T2: half of the W2 response has an exponential tail of 2 * delay (mean extra delay = delay)
BOUNDARY_HOST_W1 = 0.4


def truth_fluxes(kind: str, ell, dt: float, lag1_days: float, delay_days: float, tau_re_yr) -> tuple[np.ndarray, np.ndarray]:
    """Noise-free relative W1 and W2 light curves on the driver grid (arbitrary scale)."""
    if kind == "boundary":
        toy = breath.ToyModel("boundary", float(tau_re_yr), host_w1=BOUNDARY_HOST_W1, lag_days=lag1_days)
        flux = breath.toy_fluxes(toy, np.arange(ell.size) * dt, ell)
        return flux["w1"], flux["w2"]
    d1 = boxcar_response(ell, dt, lag1_days)
    if kind == "linear":
        d2 = shifted(d1, dt, delay_days) if delay_days else d1
    elif kind == "tail":
        d2 = (1.0 - TAIL_SHARE) * d1 + TAIL_SHARE * exponential_tail(d1, dt, delay_days / TAIL_SHARE)
    else:
        raise ValueError(kind)
    mean1, mean2 = d1.mean(), d2.mean()
    return 1.0 - DUST_FRACTION + DUST_FRACTION * d1 / mean1, 1.0 - DUST_FRACTION + DUST_FRACTION * d2 / mean2


def mock_truth(mags, visit_mjd, f1_obs, s1, f2_obs, s2, kind: str, lag1_days: float, delay_days: float,
               tau_re_yr=None, dt: float = GRID_STEP_DAYS) -> tuple[np.ndarray, np.ndarray] | None:
    """Noise-free W1 and W2 at the visits for one driver realisation ``mags`` (unit DRW).

    The amplitude of the driver is scaled so that the fractional rms of the noise-free W2
    at the visits equals the observed one after subtracting the measurement variance
    (not below half the median error). Mean levels are the observed medians.
    """
    visit_mjd = np.asarray(visit_mjd, dtype=float)
    grid_t = visit_mjd[-1] + TAIL_DAYS - dt * np.arange(mags.size)[::-1]
    target = math.sqrt(max(np.var(f2_obs) - np.mean(s2 ** 2), (0.5 * np.median(s2)) ** 2)) / float(np.mean(f2_obs))

    def at_visits(amp: float):
        ell = 10.0 ** (-0.4 * amp * mags)
        w1, w2 = truth_fluxes(kind, ell, dt, lag1_days, delay_days, tau_re_yr)
        return np.interp(visit_mjd, grid_t, w1), np.interp(visit_mjd, grid_t, w2)

    amp, best = 0.3, None
    for _ in range(6):
        v1, v2 = at_visits(amp)
        got = float(np.std(v2) / np.mean(v2))
        best = (v1, v2)
        if got <= 0 or abs(got / target - 1.0) < 0.02:
            break
        amp = float(np.clip(amp * target / got, 1e-3, 3.0))
    v1, v2 = best
    if not (np.all(np.isfinite(v1)) and np.all(np.isfinite(v2))):
        return None
    return v1 * float(np.median(f1_obs)) / float(np.median(v1)), v2 * float(np.median(f2_obs)) / float(np.median(v2))


def driver_length(visit_mjd, dt: float = GRID_STEP_DAYS) -> int:
    """Number of grid points of the driver for an object (burn-in before the first visit, tail after the last)."""
    return int(math.ceil((visit_mjd[-1] - visit_mjd[0] + BURN_IN_DAYS + TAIL_DAYS) / dt)) + 1


def sample_statistics(rows: list[dict]) -> dict:
    """Sample-level statistics of section 4 of the plan, from per-object rows
    (the output of ``fit_object`` merged with ``offsets``). The same function is
    applied to the real sample and to every synthetic sample."""
    fit = [r for r in rows if r.get("usable")]
    out = {"n_fit": len(fit)}
    if not fit:
        return out
    tau = np.array([r["tau_t1"] for r in fit])
    out.update(median_tau=float(np.median(tau)), frac_tau_positive=float(np.mean(tau > 0)),
               frac_tau_negative=float(np.mean(tau < 0)),
               mean_ds_01=float(np.mean([r["ds_01"] for r in fit])), mean_ds_02=float(np.mean([r["ds_02"] for r in fit])),
               mean_ds_12=float(np.mean([r["ds_12"] for r in fit])),
               n_well=int(sum(r["well_determined"] for r in fit)),
               median_width_t2=float(np.median([r["width_t2"] for r in fit])))
    for name in ("t0", "t1", "t2"):
        v = np.array([r[f"ds_eta_{name}"] for r in fit], dtype=float)
        e = np.array([r[f"eta_colour_{name}"] for r in fit], dtype=float)
        out[f"mean_ds_eta_{name}"] = float(np.nanmean(v)) if np.isfinite(v).any() else float("nan")
        out[f"median_eta_colour_{name}"] = float(np.nanmedian(e)) if np.isfinite(e).any() else float("nan")
    off = [r for r in fit if r.get("usable_obs") and r.get("comparable_obs")
           and all(r.get(f"usable_{k}") for k in ("t0", "t1", "t2"))]
    out["n_offset"] = len(off)
    if off:
        obs = np.array([r["d_colour_obs"] for r in off])
        out.update(median_d_colour_obs=float(np.median(obs)), frac_obs_bluer=float(np.mean(obs < 0)))
        for name in ("t0", "t1", "t2"):
            pred = np.array([r[f"d_colour_{name}"] for r in off])
            res = obs - pred
            out[f"median_d_colour_pred_{name}"] = float(np.median(pred))
            out[f"median_d_colour_res_{name}"] = float(np.median(res))
            out[f"frac_res_bluer_{name}"] = float(np.mean(res < 0))
            sig = np.array([r.get("sigma_colour_obs", np.nan) for r in off], dtype=float)
            out[f"n_res_beyond_3sigma_{name}"] = int(np.sum(np.abs(res) > 3.0 * sig))
            out[f"ratio_pred_obs_{name}"] = float(np.median(pred) / np.median(obs)) if np.median(obs) != 0 else float("nan")
    return out


def lyu_lags_days(log_lbol: float) -> tuple[float, float]:
    """Rest-frame W1 and W2 lags (Lyu, Rieke & Smith 2019, Eqs. 15-16)."""
    x = 10.0 ** log_lbol / 3.828e33 / 1e11
    return 10.0 ** 2.10 * x ** 0.47, 10.0 ** 2.20 * x ** 0.45
