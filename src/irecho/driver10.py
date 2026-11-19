"""Stage 10: an optical driver (ZTF) and fixed kernels for W1 and W2 (docs/PREDICTIONS_10.md).

Model for one object and one infrared band, at the visits whose kernel window is
covered by the optical light curve:

    F_W(t_i) = a + b * (K * D)(t_i)  [+ alpha' * D(t_i)]  [+ eta * rising_i]

D is the optical flux. Between the optical nights it is the conditional mean of a
Gaussian process (damped random walk plus white noise, constant mean with a flat
prior) given the nightly fluxes; the conditional covariance P of D on the grid goes
into the covariance of the residuals:

    Cov = diag(sigma_W^2) + b^2 c P c^T

where c holds the weights of the kernel on the grid. Because the conditional mean is
linear in the nightly fluxes, (K * D)(t_i) = W y with weights W that depend only on
the observing times, the errors and the hyper-parameters: they are computed once per
object (``prepare``) and reused for every synthetic realisation.

The second half of the module makes the synthetic light curves of section 5 of the plan.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.optimize import minimize
from scipy.signal import fftconvolve

from . import breath, hyst8
from . import transfer9 as t9

# --- numbers fixed in docs/PREDICTIONS_10.md before the ZTF light curves were opened for this stage ---
SUPPORT_PAST_DAYS = 700.0       # the kernel may reach this far into the past of a visit
SUPPORT_FUTURE_DAYS = 100.0     # and this far into its future (negative lags are allowed)
MAX_GAP_DAYS = 250.0            # no optical gap longer than this inside the window of a visit
TAU_STEP_DAYS = 5.0
SHAPES = ("f0", "f05", "f1", "f2", "gauss")
SHAPE_WIDTH = {"f0": 0.0, "f05": 0.5, "f1": 1.0, "f2": 2.0}     # full width of the top-hat / |tau|
BASE = "f1"
GAUSS_SIGMA = 1.0 / math.sqrt(12.0)   # sigma / |tau|: the variance of the base top-hat
GAUSS_CUT = 3.0
FIT_GRID_STEP_DAYS = 2.0
MIN_VISITS = 8
N_ITERATIONS = 3
JOINT_STEP_DAYS = 10.0
MA_MIN_TAU_DAYS = 10.0          # the delayed term of the model with an undelayed term
TRUTH_STEP_DAYS = 1.0
BURN_IN_DAYS = 3000.0
MIN_NIGHTS = 100
DRW_TAU_BOUNDS = (10.0, 5000.0)


# --- optical light curve: nights, hyper-parameters, conditional mean and covariance -----------------
def nightly(mjd, mag, magerr) -> dict:
    """Inverse-variance mean flux per night (relative to the flux of the median magnitude)."""
    mjd, mag, magerr = (np.asarray(v, dtype=float) for v in (mjd, mag, magerr))
    order = np.argsort(mjd)
    mjd, mag, magerr = mjd[order], mag[order], magerr[order]
    ref = float(np.median(mag)) if mag.size else float("nan")
    flux = 10.0 ** (-0.4 * (mag - ref))
    err = flux * magerr * math.log(10.0) / 2.5
    # a night: consecutive points less than 0.5 d apart
    night = np.concatenate([[0], np.cumsum(np.diff(mjd) > 0.5)]) if mjd.size else np.zeros(0, int)
    t, y, s = [], [], []
    for k in np.unique(night):
        sel = night == k
        w = 1.0 / err[sel] ** 2
        t.append(float(np.sum(w * mjd[sel]) / w.sum()))
        y.append(float(np.sum(w * flux[sel]) / w.sum()))
        s.append(float(1.0 / math.sqrt(w.sum())))
    return {"t": np.array(t), "y": np.array(y), "s": np.array(s), "mag_ref": ref}


def clip_nights(t, y, half_window: float = 15.0, n_sigma: float = 5.0, min_neighbours: int = 3) -> np.ndarray:
    """True for the nights kept: a night is dropped when it differs from the median of the other
    nights within ``half_window`` days by more than ``n_sigma`` robust standard deviations of all
    such differences (one pass)."""
    t, y = np.asarray(t, dtype=float), np.asarray(y, dtype=float)
    diff = np.full(t.size, np.nan)
    for i in range(t.size):
        sel = (np.abs(t - t[i]) <= half_window)
        sel[i] = False
        if sel.sum() >= min_neighbours:
            diff[i] = y[i] - np.median(y[sel])
    ok = np.isfinite(diff)
    if ok.sum() < 10:
        return np.ones(t.size, dtype=bool)
    scale = 1.4826 * float(np.median(np.abs(diff[ok] - np.median(diff[ok]))))
    return ~(ok & (np.abs(diff) > n_sigma * scale)) if scale > 0 else np.ones(t.size, dtype=bool)


def _drw_nll(theta, dt_abs, y, s2):
    var, tau, jit = math.exp(theta[0]), math.exp(theta[1]), math.exp(theta[2])
    cov = var * np.exp(-dt_abs / tau)
    cov[np.diag_indices_from(cov)] += s2 + jit
    try:
        chol = np.linalg.cholesky(cov)
    except np.linalg.LinAlgError:
        return 1e30
    one = np.ones(y.size)
    a = np.linalg.solve(chol, np.column_stack([y, one]))
    mean = float(a[:, 0] @ a[:, 1] / (a[:, 1] @ a[:, 1]))
    res = a[:, 0] - mean * a[:, 1]
    # flat prior on the mean: the extra term ln(1^T C^-1 1)
    return float(res @ res + 2.0 * np.log(np.diag(chol)).sum() + math.log(a[:, 1] @ a[:, 1]))


def fit_drw(t, y, s) -> dict:
    """Maximum likelihood damped random walk plus white noise ('jitter') for the nightly fluxes.
    The time scale is kept inside ``DRW_TAU_BOUNDS``."""
    t, y, s = (np.asarray(v, dtype=float) for v in (t, y, s))
    dt_abs = np.abs(t[:, None] - t[None, :])
    s2 = s ** 2
    v0 = max(float(np.var(y) - np.mean(s2)), 1e-3 * float(np.var(y)) + 1e-12)
    lo, hi = math.log(DRW_TAU_BOUNDS[0]), math.log(DRW_TAU_BOUNDS[1])

    def nll(theta):
        th = np.array([theta[0], min(max(theta[1], lo), hi), theta[2]])
        return _drw_nll(th, dt_abs, y, s2) + 1e3 * (theta[1] - th[1]) ** 2

    best = None
    for tau0 in (100.0, 600.0):
        res = minimize(nll, [math.log(v0), math.log(tau0), math.log(max(float(np.median(s2)), 1e-12))],
                       method="Nelder-Mead", options={"xatol": 1e-3, "fatol": 1e-4, "maxiter": 600})
        if best is None or res.fun < best.fun:
            best = res
    return {"var": float(math.exp(best.x[0])), "tau": float(math.exp(min(max(best.x[1], lo), hi))),
            "jitter_var": float(math.exp(best.x[2])), "nll": float(best.fun)}


def driver_grid(t_first: float, t_last: float, step: float = FIT_GRID_STEP_DAYS) -> np.ndarray:
    n = int(math.ceil((t_last - t_first) / step)) + 1
    return t_first + step * np.arange(n)


def gp_driver(t, s, hyper: dict, grid) -> tuple[np.ndarray, np.ndarray]:
    """Gain (G x N) and conditional covariance P (G x G) of the driver on ``grid``:
    mean = gain @ y. Constant mean with a flat prior (universal kriging)."""
    t, s, grid = (np.asarray(v, dtype=float) for v in (t, s, grid))
    var, tau = hyper["var"], hyper["tau"]
    koo = var * np.exp(-np.abs(t[:, None] - t[None, :]) / tau)
    koo[np.diag_indices_from(koo)] += s ** 2 + hyper["jitter_var"]
    kgo = var * np.exp(-np.abs(grid[:, None] - t[None, :]) / tau)
    kgg = var * np.exp(-np.abs(grid[:, None] - grid[None, :]) / tau)
    ainv = np.linalg.inv(koo)
    one = np.ones(t.size)
    a1 = ainv @ one
    denom = float(one @ a1)
    wm = a1 / denom                                # mean estimate = wm @ y
    ka = kgo @ ainv
    u = 1.0 - ka @ one
    gain = ka + np.outer(u, wm)
    p = kgg - ka @ kgo.T + np.outer(u, u) / denom
    return gain, 0.5 * (p + p.T)


def linear_driver(t, s, grid) -> tuple[np.ndarray, np.ndarray]:
    """The same pair for straight lines between the nights (no model of the gaps)."""
    gain = t9.interp_weights(t, np.clip(grid, t[0], t[-1]))
    return gain, (gain * np.asarray(s, dtype=float) ** 2) @ gain.T


# --- kernels on the even grid -------------------------------------------------------------------------
def _interp_rows(g0: float, h: float, n: int, u) -> np.ndarray:
    u = np.atleast_1d(np.asarray(u, dtype=float))
    x = (u - g0) / h
    if np.any(x < -1e-9) or np.any(x > n - 1 + 1e-9):
        raise ValueError("extrapolation is not allowed")
    x = np.clip(x, 0.0, n - 1.0)
    k = np.minimum(np.floor(x).astype(int), n - 2)
    fr = x - k
    w = np.zeros((u.size, n))
    rows = np.arange(u.size)
    w[rows, k] = 1.0 - fr
    w[rows, k + 1] += fr
    return w


def _antiderivative_rows(g0: float, h: float, n: int, u) -> np.ndarray:
    """Weights of integral_{g0}^{u} of the linear interpolant on the grid (exact)."""
    u = np.atleast_1d(np.asarray(u, dtype=float))
    x = (u - g0) / h
    if np.any(x < -1e-9) or np.any(x > n - 1 + 1e-9):
        raise ValueError("extrapolation is not allowed")
    x = np.clip(x, 0.0, n - 1.0)
    k = np.minimum(np.floor(x).astype(int), n - 2)
    fr = x - k
    idx = np.arange(n)[None, :]
    w = np.where(idx < k[:, None], h, 0.0)
    w[:, 0] = np.where(k >= 1, 0.5 * h, 0.0)
    rows = np.arange(u.size)
    w[rows, k] = np.where(k >= 1, 0.5 * h, 0.0) + h * (fr - 0.5 * fr * fr)
    w[rows, k + 1] += h * 0.5 * fr * fr
    return w


def support(shape: str, tau: float) -> tuple[float, float]:
    """Lags covered by the kernel (days)."""
    half = GAUSS_CUT * GAUSS_SIGMA * abs(tau) if shape == "gauss" else 0.5 * SHAPE_WIDTH[shape] * abs(tau)
    return tau - half, tau + half


def kernel_rows(g0: float, h: float, n: int, t_eval, tau: float, shape: str = BASE) -> np.ndarray:
    """c (len(t_eval) x n): (K * D)(t_eval) = c @ D_grid for a kernel of mean lag ``tau``."""
    t_eval = np.atleast_1d(np.asarray(t_eval, dtype=float))
    lo, hi = support(shape, tau)
    if hi - lo < 1e-9:
        return _interp_rows(g0, h, n, t_eval - tau)
    if shape == "gauss":
        z = np.linspace(-GAUSS_CUT, GAUSS_CUT, 41)
        w = np.exp(-0.5 * z * z)
        w /= w.sum()
        sig = GAUSS_SIGMA * abs(tau)
        out = np.zeros((t_eval.size, n))
        for zq, wq in zip(z, w):
            out += wq * _interp_rows(g0, h, n, t_eval - tau - zq * sig)
        return out
    return (_antiderivative_rows(g0, h, n, t_eval - lo) - _antiderivative_rows(g0, h, n, t_eval - hi)) / (hi - lo)


def grid_points() -> list[tuple[str, float]]:
    """(shape, tau) pairs whose support lies inside [-SUPPORT_FUTURE_DAYS, SUPPORT_PAST_DAYS]."""
    out = []
    n = int(round((SUPPORT_PAST_DAYS + SUPPORT_FUTURE_DAYS) / TAU_STEP_DAYS))
    for shape in SHAPES:
        for i in range(n + 1):
            tau = -SUPPORT_FUTURE_DAYS + i * TAU_STEP_DAYS
            lo, hi = support(shape, tau)
            if lo >= -SUPPORT_FUTURE_DAYS - 1e-9 and hi <= SUPPORT_PAST_DAYS + 1e-9:
                out.append((shape, float(tau)))
    return out


def usable_visits(visit_t, opt_t) -> np.ndarray:
    """Visits whose whole window [t - past, t + future] lies inside the optical light curve and
    holds no optical gap longer than ``MAX_GAP_DAYS``. The same set is used by every model."""
    visit_t, opt_t = np.asarray(visit_t, dtype=float), np.asarray(opt_t, dtype=float)
    lo, hi = visit_t - SUPPORT_PAST_DAYS, visit_t + SUPPORT_FUTURE_DAYS
    ok = (lo >= opt_t[0]) & (hi <= opt_t[-1])
    gaps = np.diff(opt_t)
    for k in np.nonzero(gaps > MAX_GAP_DAYS)[0]:
        ok &= ~((lo < opt_t[k + 1]) & (hi > opt_t[k]))
    return ok


def gp_state(opt_t, opt_s, hyper: dict | None, method: str = "gp") -> dict:
    """The driver on the 2-day grid: gain and conditional covariance (see ``gp_driver``)."""
    opt_t, opt_s = np.asarray(opt_t, dtype=float), np.asarray(opt_s, dtype=float)
    grid = driver_grid(opt_t[0], opt_t[-1])
    gain, p = gp_driver(opt_t, opt_s, hyper, grid) if method == "gp" else linear_driver(opt_t, opt_s, grid)
    in_gap = np.zeros(grid.size)        # grid points inside optical gaps longer than 30 days
    for k in np.nonzero(np.diff(opt_t) > 30.0)[0]:
        in_gap[(grid > opt_t[k]) & (grid < opt_t[k + 1])] = 1.0
    return {"opt_t": opt_t, "g0": float(grid[0]), "n": int(grid.size), "gain": gain, "p": p, "in_gap": in_gap}


def prepare(opt_t, opt_s, hyper: dict | None, visit_t, method: str = "gp", joint: bool = False,
            shapes=SHAPES, state: dict | None = None) -> dict:
    """Everything that depends on the observing times only.

    Returns the usable visits ``idx``, the grid points, and for every grid point the weights
    W (m x N) of the nightly fluxes and the covariance V (m x m) of the prediction; the same
    for the undelayed term (W0, V0) and its cross term with the base shape (VX).
    With ``joint`` the cross terms between base-shape lags on the coarse grid are kept (VJ).
    ``state`` (from ``gp_state``) can be passed to reuse the driver of one object for the visits
    of another."""
    opt_t, visit_t = np.asarray(opt_t, dtype=float), np.asarray(visit_t, dtype=float)
    idx = np.nonzero(usable_visits(visit_t, opt_t))[0]
    pts = [p for p in grid_points() if p[0] in shapes]
    out = {"idx": idx, "shapes": np.array([p[0] for p in pts]), "taus": np.array([p[1] for p in pts]), "W": None,
           "n_nights": int(opt_t.size)}
    if idx.size < MIN_VISITS:
        return out
    state = gp_state(opt_t, opt_s, hyper, method) if state is None else state
    g0, h, n, gain, p = state["g0"], FIT_GRID_STEP_DAYS, state["n"], state["gain"], state["p"]
    tv = visit_t[idx]
    m = idx.size
    c = np.stack([kernel_rows(g0, h, n, tv, tau, shape) for shape, tau in pts])        # n_pts x m x G
    c0 = _interp_rows(g0, h, n, tv)
    flat = c.reshape(-1, n)
    cp = (flat @ p).reshape(len(pts), m, n)
    out["W"] = (flat @ gain).reshape(len(pts), m, -1)
    out["V"] = np.einsum("pmg,png->pmn", cp, c)
    out["W0"] = c0 @ gain
    out["V0"] = c0 @ p @ c0.T
    base = out["shapes"] == BASE
    out["VX"] = cp[base] @ c0.T                                                           # n_base x m x m
    # share of the kernel window of each grid point that falls in optical gaps longer than 30 days
    out["gap_share"] = (c @ state["in_gap"]).mean(axis=1)
    if joint:
        sub = np.nonzero(base & (np.abs(np.round(out["taus"] / JOINT_STEP_DAYS) * JOINT_STEP_DAYS - out["taus"]) < 1e-9))[0]
        vj = cp[sub].reshape(-1, n) @ c[sub].reshape(-1, n).T
        out["joint_index"] = sub
        out["VJ"] = vj.reshape(sub.size, m, sub.size, m).transpose(0, 2, 1, 3)           # i, j, m, m
    return out


# --- generalised least squares on the stack of grid points --------------------------------------------
def _solve(cov, a_mat, y):
    cinv = np.linalg.inv(cov)
    ca = cinv @ a_mat
    normal = np.swapaxes(a_mat, -1, -2) @ ca
    rhs = np.einsum("gmp,m->gp", ca, y)
    theta = np.linalg.solve(normal, rhs[..., None])[..., 0]
    chi2 = np.einsum("m,gmk,k->g", y, cinv, y) - np.einsum("gp,gp->g", theta, rhs)
    return theta, chi2, normal


def gls_stack(pred, v, y, s, extra=None, n_iter: int = N_ITERATIONS) -> dict:
    """Fit y = a + b pred for a stack of predictions (G x m) with Cov = diag(s^2) + b^2 V.
    b starts at 0 and is updated ``n_iter`` times. S = chi2 + ln det Cov."""
    g, m = pred.shape
    eye = np.eye(m) * s ** 2
    a_mat = np.stack([np.ones((g, m)), pred], axis=-1)
    b = np.zeros(g)
    for _ in range(n_iter + 1):
        cov = eye + b[:, None, None] ** 2 * v
        theta, chi2, _ = _solve(cov, a_mat, y)
        b = theta[:, 1]
    cov = eye + b[:, None, None] ** 2 * v
    theta, chi2, normal = _solve(cov, a_mat, y)
    dvar = b[:, None] ** 2 * np.einsum("gmm->gm", v)          # variance of the prediction that comes from the driver
    out = {"s": chi2 + np.linalg.slogdet(cov)[1], "chi2": chi2, "theta": theta,
           "pred": theta[:, :1] + theta[:, 1:2] * pred, "var_share": dvar / (s ** 2 + dvar)}
    if extra is not None:
        a3 = np.concatenate([a_mat, np.broadcast_to(np.asarray(extra, dtype=float), (g, m))[..., None]], axis=-1)
        theta3, chi3, normal3 = _solve(cov, a3, y)
        out.update(theta_extra=theta3, chi2_extra=chi3, theta_extra_cov=np.linalg.inv(normal3))
    return out


def gls_two(pred, pred0, v, vx, v0, y, s, n_iter: int = N_ITERATIONS) -> dict:
    """Fit y = a + b pred + c pred0 (the undelayed term) with
    Cov = diag(s^2) + b^2 V + b c (VX + VX^T) + c^2 V0."""
    g, m = pred.shape
    eye = np.eye(m) * s ** 2
    a_mat = np.stack([np.ones((g, m)), pred, np.broadcast_to(pred0, (g, m))], axis=-1)
    sym = vx + np.swapaxes(vx, -1, -2)
    b, c = np.zeros(g), np.zeros(g)
    for _ in range(n_iter + 1):
        cov = eye + b[:, None, None] ** 2 * v + (b * c)[:, None, None] * sym + c[:, None, None] ** 2 * v0
        theta, chi2, normal = _solve(cov, a_mat, y)
        b, c = theta[:, 1], theta[:, 2]
    cov = eye + b[:, None, None] ** 2 * v + (b * c)[:, None, None] * sym + c[:, None, None] ** 2 * v0
    theta, chi2, normal = _solve(cov, a_mat, y)
    return {"s": chi2 + np.linalg.slogdet(cov)[1], "chi2": chi2, "theta": theta, "theta_cov": np.linalg.inv(normal),
            "pred": theta[:, :1] + theta[:, 1:2] * pred + theta[:, 2:3] * pred0}


def _quantiles(taus, s) -> tuple[float, float, float]:
    """16, 50 and 84 per cent points of the distribution exp(-S / 2) on the lag grid."""
    w = np.exp(-0.5 * (s - s.min()))
    cdf = np.cumsum(w) / w.sum()
    return tuple(float(taus[min(int(np.searchsorted(cdf, q)), taus.size - 1)]) for q in (0.16, 0.5, 0.84))


def weighted_corr(x, y, s) -> float:
    w = 1.0 / s ** 2
    mx, my = np.sum(w * x) / w.sum(), np.sum(w * y) / w.sum()
    sxy = np.sum(w * (x - mx) * (y - my))
    den = math.sqrt(float(np.sum(w * (x - mx) ** 2) * np.sum(w * (y - my) ** 2)))
    return float(sxy / den) if den > 0 else 0.0


def fit_band(prep: dict, y_opt, f, s, rising=None, full: bool = True) -> dict:
    """All models for one infrared band. ``f`` and ``s`` are the fluxes and errors at every
    visit of the object (the usable ones are selected with prep['idx']).

    C   constant                 (a)
    M0  no lag                   (a, b)
    M1  base kernel              (a, b, tau)        -- the reference model
    M2  shape free               (a, b, tau, f)     f in 0, 0.5, 1, 2
    MG  Gaussian kernel          (a, b, tau)
    MA  M1 + undelayed term      (a, b, alpha', tau)
    ME  M1 + branch offset eta   (a, b, eta, tau)
    """
    out = {"usable": False, "n_used": int(prep["idx"].size)}
    if prep["W"] is None:
        return out
    idx = prep["idx"]
    y, sy = np.asarray(f, dtype=float)[idx], np.asarray(s, dtype=float)[idx]
    m = idx.size
    shapes, taus = prep["shapes"], prep["taus"]
    p_all = prep["W"] @ y_opt
    p0 = prep["W0"] @ y_opt
    w = 1.0 / sy ** 2
    mean = float(np.sum(w * y) / w.sum())
    chi2_c = float(np.sum(w * (y - mean) ** 2))
    s_c = chi2_c + float(np.sum(np.log(sy ** 2)))
    lab = None
    if rising is not None:
        lab = np.asarray(rising, dtype=bool)[idx].astype(float)
        if not 0 < lab.sum() < m:
            lab = None
    scale = float(np.median(y))
    fit = gls_stack(p_all, prep["V"], y, sy, extra=None if lab is None else lab * scale)
    base = np.nonzero(shapes == BASE)[0]
    i1 = int(base[np.argmin(fit["s"][base])])
    tb, sb = taus[base], fit["s"][base]
    q16, q50, q84 = _quantiles(tb, sb)
    level = max(1.0, float(fit["chi2"][i1]) / max(m - 3, 1))
    inside = sb - sb.min() <= level
    zero = gls_stack(p0[None, :], prep["V0"][None], y, sy)
    out.update(usable=True, idx=idx, s_c=s_c, chi2_c=chi2_c, s_m0=float(zero["s"][0]), s_m1=float(fit["s"][i1]),
               chi2_m1=float(fit["chi2"][i1]), tau=float(taus[i1]), a=float(fit["theta"][i1, 0]),
               b=float(fit["theta"][i1, 1]), pred=fit["pred"][i1], tau_q16=q16, tau_q50=q50, tau_q84=q84,
               tau_lo=float(tb[inside].min()), tau_hi=float(tb[inside].max()),
               tau_at_edge=bool(taus[i1] <= tb.min() + 1e-9 or taus[i1] >= tb.max() - 1e-9),
               r_w=weighted_corr(fit["pred"][i1], y, sy), var_share=float(np.median(fit["var_share"][i1])),
               gap_share=float(prep["gap_share"][i1]),
               ds_c1=float(s_c - fit["s"][i1]), ds_01=float(zero["s"][0] - fit["s"][i1]))
    if not full:
        return out
    tophat = np.nonzero(shapes != "gauss")[0]
    i2 = int(tophat[np.argmin(fit["s"][tophat])])
    gauss = np.nonzero(shapes == "gauss")[0]
    out.update(s_m2=float(fit["s"][i2]), tau_m2=float(taus[i2]), shape_m2=str(shapes[i2]),
               ds_12=float(fit["s"][i1] - fit["s"][i2]))
    for name in SHAPE_WIDTH:        # the best lag of every top-hat width (sensitivity to the shape)
        sel = np.nonzero(shapes == name)[0]
        j = int(sel[np.argmin(fit["s"][sel])])
        out[f"tau_{name}"], out[f"s_{name}"] = float(taus[j]), float(fit["s"][j])
    if gauss.size:
        ig = int(gauss[np.argmin(fit["s"][gauss])])
        out.update(s_mg=float(fit["s"][ig]), tau_mg=float(taus[ig]))
    # MA: the undelayed term next to the base kernel
    # (the delayed term needs tau >= MA_MIN_TAU_DAYS: at tau = 0 the two terms are the same column)
    pos = tb >= MA_MIN_TAU_DAYS - 1e-9
    two = gls_two(p_all[base][pos], p0, prep["V"][base][pos], prep["VX"][pos], prep["V0"], y, sy)
    ja = int(np.argmin(two["s"]))
    b_a, c_a = float(two["theta"][ja, 1]), float(two["theta"][ja, 2])
    sig_c = float(math.sqrt(max(two["theta_cov"][ja, 2, 2], 0.0)))
    out.update(s_ma=float(two["s"][ja]), tau_ma=float(tb[pos][ja]), b_ma=b_a, c_ma=c_a, sigma_c_ma=sig_c,
               z_alpha=float(c_a / sig_c) if sig_c > 0 else float("nan"),
               alpha=float(c_a / (b_a + c_a)) if (b_a + c_a) != 0 else float("nan"),
               ds_1a=float(fit["s"][i1] - two["s"][ja]), pred_ma=two["pred"][ja])
    # ME: the branch offset, with the covariance of the fit without it (same ln det)
    out.update(eta=float("nan"), sigma_eta=float("nan"), ds_1e=float("nan"), tau_me=float("nan"))
    if lab is not None:
        s_eta = fit["s"] - fit["chi2"] + fit["chi2_extra"]
        je = int(base[np.argmin(s_eta[base])])
        out.update(eta=float(fit["theta_extra"][je, 2]),                      # in units of the median flux
                   sigma_eta=float(math.sqrt(fit["theta_extra_cov"][je, 2, 2])),
                   ds_1e=float(fit["s"][i1] - s_eta[je]), tau_me=float(taus[je]))
    return out


def fit_joint(prep: dict, y_opt, f1, s1, f2, s2, r: float, fit1: dict, fit2: dict) -> dict:
    """Both bands at once on the coarse grid of base-shape lags. The two bands share the
    uncertainty of the driver (cross term b1 b2 c1 P c2^T) and the errors of one visit are
    correlated (r). Gives the distribution of tau2 - tau1."""
    idx = prep["idx"]
    y1, y2 = np.asarray(f1, dtype=float)[idx], np.asarray(f2, dtype=float)[idx]
    e1, e2 = np.asarray(s1, dtype=float)[idx], np.asarray(s2, dtype=float)[idx]
    m = idx.size
    sub = prep["joint_index"]
    taus = prep["taus"][sub]
    n = sub.size
    p = prep["W"][sub] @ y_opt                                           # n x m
    v = prep["V"][sub]
    vj = prep["VJ"]                                                      # n x n x m x m
    y = np.concatenate([y1, y2])
    ii, jj = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    ii, jj = ii.ravel(), jj.ravel()
    g = ii.size
    a_mat = np.zeros((g, 2 * m, 4))
    a_mat[:, :m, 0] = 1.0
    a_mat[:, :m, 1] = p[ii]
    a_mat[:, m:, 2] = 1.0
    a_mat[:, m:, 3] = p[jj]
    # start from the slopes of the separate fits at every lag
    sep1 = gls_stack(p, v, y1, e1)["theta"][:, 1]
    sep2 = gls_stack(p, v, y2, e2)["theta"][:, 1]
    b1, b2 = sep1[ii], sep2[jj]
    cross_err = np.diag(r * e1 * e2)

    def cov_of(b1, b2):
        cov = np.zeros((g, 2 * m, 2 * m))
        cov[:, :m, :m] = np.diag(e1 ** 2) + b1[:, None, None] ** 2 * v[ii]
        cov[:, m:, m:] = np.diag(e2 ** 2) + b2[:, None, None] ** 2 * v[jj]
        off = cross_err + (b1 * b2)[:, None, None] * vj[ii, jj]
        cov[:, :m, m:] = off
        cov[:, m:, :m] = np.swapaxes(off, -1, -2)
        return cov

    for _ in range(2):
        cov = cov_of(b1, b2)
        theta, chi2, _ = _solve(cov, a_mat, y)
        b1, b2 = theta[:, 1], theta[:, 3]
    s = (chi2 + np.linalg.slogdet(cov)[1]).reshape(n, n)
    w = np.exp(-0.5 * (s - s.min()))
    w /= w.sum()
    d = (taus[None, :] - taus[:, None]).ravel()
    order = np.argsort(d, kind="stable")
    cdf = np.cumsum(w.ravel()[order])
    q = [float(d[order][min(int(np.searchsorted(cdf, x)), d.size - 1)]) for x in (0.025, 0.16, 0.5, 0.84, 0.975)]
    k = int(np.argmin(s))
    return {"tau1_joint": float(taus[k // n]), "tau2_joint": float(taus[k % n]), "dtau_q025": q[0], "dtau_q16": q[1],
            "dtau_q50": q[2], "dtau_q84": q[3], "dtau_q975": q[4], "s_joint": float(s.min()),
            "p_dtau_positive": float(w.ravel()[d > 0].sum())}


def colour_offsets(f1, s1, f2, s2, r: float, fit1: dict, fit2: dict, rising=None, key: str = "pred") -> dict:
    """Direct offset of stage 8 on the visits used by the fit: measured, and predicted with both
    bands replaced by the fitted models. d_colour < 0: the rising branch is bluer."""
    out = {f"{k}_{n}": (False if k in ("usable", "comparable") else float("nan"))
           for n in ("obs", "pred") for k in ("delta", "sigma_delta", "d_colour", "sigma_colour", "usable", "comparable")}
    if not (fit1.get("usable") and fit2.get("usable")):
        return out
    idx = fit1["idx"]
    f1, s1, f2, s2 = (np.asarray(v, dtype=float) for v in (f1, s1, f2, s2))
    rising = hyst8.branch_rising(f1, f2) if rising is None else np.asarray(rising, dtype=bool)
    lab = rising[idx]
    obs = hyst8.direct_offset(f1[idx], s1[idx], f2[idx], s2[idx], rising=lab, r=r)
    for k in ("delta", "sigma_delta", "d_colour", "sigma_colour", "usable", "comparable"):
        out[f"{k}_obs"] = obs[k]
    if not obs["usable"]:
        return out
    inside = (f2[idx] >= obs["lo"]) & (f2[idx] <= obs["hi"])
    # the model W2 need not span the same range: the visits are those of the measured comparison
    d = hyst8.direct_offset(fit1[key], s1[idx], fit2[key], s2[idx], rising=lab, r=r, inside=inside)
    for k in ("delta", "sigma_delta", "d_colour", "sigma_colour", "usable", "comparable"):
        out[f"{k}_pred"] = d[k]
    return out


# --- literature relations (rest frame, days) --------------------------------------------------------
MANDAL = {"w1": (0.39, 2.66), "w2": (0.33, 2.76), "log_l0": 45.7}     # Mandal et al. 2024, Table 5 (BCES, ICCF)


def lyu_lags_days(log_lbol: float) -> tuple[float, float]:
    return t9.lyu_lags_days(log_lbol)


def mandal_lags_days(log_lbol: float) -> tuple[float, float]:
    x = log_lbol - MANDAL["log_l0"]
    return 10.0 ** (MANDAL["w1"][1] + MANDAL["w1"][0] * x), 10.0 ** (MANDAL["w2"][1] + MANDAL["w2"][0] * x)


# --- synthetic light curves (section 5 of the plan) ---------------------------------------------------
def kernel_array(shape: str, tau: float, dt: float = TRUTH_STEP_DAYS) -> np.ndarray:
    """Causal kernel on the truth grid (index = lag / dt), unit sum.
    'f0', 'f05', 'f1', 'f2': top-hats of full width f * tau centred on tau;
    'gwide': Gaussian of sigma = tau / 2 cut at lag 0;
    'exp': half the response is the 'f1' top-hat of mean lag tau / 2 and half an exponential of
    e-folding time 1.5 tau that starts at lag 0 (mean lag of the sum = tau)."""
    tau = max(float(tau), 0.0)
    if tau < dt:
        return np.array([1.0])
    if shape in SHAPE_WIDTH:
        half = 0.5 * SHAPE_WIDTH[shape] * tau
        lo, hi = tau - half, tau + half
        if hi - lo < dt:
            k = np.zeros(int(math.floor(tau / dt)) + 2)
            x = tau / dt
            i = int(math.floor(x))
            k[i], k[i + 1] = 1.0 - (x - i), x - i
            return k
        edges = np.arange(int(math.ceil(hi / dt)) + 2) * dt - 0.5 * dt
        k = np.clip(np.minimum(edges[1:], hi) - np.maximum(edges[:-1], lo), 0.0, None)
        return k / k.sum()
    lags = np.arange(int(math.ceil(6.0 * tau / dt)) + 1) * dt
    if shape == "gwide":
        k = np.exp(-0.5 * ((lags - tau) / (0.5 * tau)) ** 2)
        return k / k.sum()
    if shape == "exp":
        lags = np.arange(int(math.ceil(12.0 * tau / dt)) + 1) * dt
        box = kernel_array("f1", 0.5 * tau, dt)
        tail = np.exp(-lags / (1.5 * tau))
        k = 0.5 * tail / tail.sum()
        k[: box.size] += 0.5 * box
        return k
    raise ValueError(shape)


def kernel_mean_lag(k, dt: float = TRUTH_STEP_DAYS) -> float:
    return float(np.sum(np.arange(k.size) * dt * k) / k.sum())


def respond(ell, k) -> np.ndarray:
    """ell convolved with a causal kernel; the start is filled with the first value."""
    if k.size == 1:
        return ell.copy()
    pad = np.concatenate([np.full(k.size - 1, ell[0]), ell])
    return fftconvolve(pad, k, mode="valid")


def truth_grid(opt_t, visit_t, dt: float = TRUTH_STEP_DAYS) -> np.ndarray:
    t_end = max(float(opt_t[-1]), float(visit_t[-1])) + 5.0
    t_start = min(float(opt_t[0]), float(visit_t[0])) - BURN_IN_DAYS
    return t_start + dt * np.arange(int(math.ceil((t_end - t_start) / dt)) + 1)


def lognormal_amp(frac_rms: float) -> float:
    """Amplitude (mag) of a unit DRW in magnitude that gives this fractional rms of the flux."""
    return 2.5 / math.log(10.0) * math.sqrt(math.log(1.0 + frac_rms ** 2))


def scale_band(v, target_frac_rms: float, median_flux: float) -> np.ndarray:
    """Values at the visits with the variable part rescaled so that the fractional rms equals the
    target; the median is set to ``median_flux``."""
    v = np.asarray(v, dtype=float)
    d = v / float(np.mean(v)) - 1.0
    rms = float(np.std(d))
    g = float(np.clip(target_frac_rms / rms, 0.02, 10.0)) if rms > 0 else 0.0
    out = 1.0 + g * d
    return median_flux * out / float(np.median(out))


def linear_truth(ell, lag1: float, lag2: float, shape: str = "f1", alpha1: float = 0.0, alpha2: float = 0.0):
    """Noise-free relative W1 and W2 on the truth grid for fixed kernels (+ undelayed shares)."""
    d1 = respond(ell, kernel_array(shape, lag1))
    d2 = respond(ell, kernel_array(shape, lag2))
    if alpha1:
        d1 = (1.0 - alpha1) * d1 + alpha1 * ell
    if alpha2:
        d2 = (1.0 - alpha2) * d2 + alpha2 * ell
    return d1, d2


def boundary_truth(mags, grid, visit_t, frac_rms_w2: float, lag1: float, tau_re_yr: float):
    """Stage 9's history-dependent response (breath.ToyModel 'boundary'), with the amplitude of the
    driver chosen so that the fractional rms of W2 at the visits equals the target.
    Returns ell, W1, W2 on the truth grid."""
    toy = breath.ToyModel("boundary", float(tau_re_yr), host_w1=t9.BOUNDARY_HOST_W1, lag_days=lag1)
    tg = grid - grid[0]
    amp, best = 0.3, None
    for _ in range(6):
        ell = 10.0 ** (-0.4 * amp * mags)
        flux = breath.toy_fluxes(toy, tg, ell)
        v2 = np.interp(visit_t, grid, flux["w2"])
        got = float(np.std(v2) / np.mean(v2))
        best = (ell, flux["w1"], flux["w2"])
        if got <= 0 or abs(got / frac_rms_w2 - 1.0) < 0.02:
            break
        amp = float(np.clip(amp * frac_rms_w2 / got, 1e-3, 3.0))
    return best


def excess_frac_rms(y, s, floor_frac: float = 0.5) -> float:
    """Fractional rms after subtracting the measurement variance (not below ``floor_frac`` of the
    median error)."""
    y, s = np.asarray(y, dtype=float), np.asarray(s, dtype=float)
    return math.sqrt(max(float(np.var(y) - np.mean(s ** 2)), (floor_frac * float(np.median(s))) ** 2)) / float(np.mean(y))
