"""Stage 5b: slope of the variable component (d log T / d log L), toy models
and the hysteresis statistic. Definitions: docs/PREDICTIONS_5b.md sections 6-7,
derivations: docs/THEORY.md sections 11-13.

The estimator never regresses a per-visit temperature on a per-visit
luminosity (both would come from the same two fluxes). It fits a curve in the
(W2 excess flux, W1 excess flux) plane; the temperature only appears as a
parameter of that curve.

    excess       e_b,k = F_b,k - A_b          A = reference (faint) state
    curve        T_c(e2) = T_p (e2 / e2_p)^beta,   e1 = e2 * rho(T_c)
    rho(T)       blackbody F_nu(W1) / F_nu(W2)

A receding sublimation boundary (constant temperature, changing area) gives
beta = 0 whatever the host contributes. Dust heated in place gives beta > 0.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq, least_squares

from .config import LAMBDA_ISO_UM
from .photometry import C_LIGHT, H_PLANCK, K_BOLTZ

DAYS_PER_YEAR = 365.25
_NU = {b: C_LIGHT / (LAMBDA_ISO_UM[b] * 1e-6) for b in ("w1", "w2")}
_X1 = H_PLANCK * _NU["w1"] / K_BOLTZ  # h nu / k in kelvin
_X2 = H_PLANCK * _NU["w2"] / K_BOLTZ
_NU3 = (_NU["w1"] / _NU["w2"]) ** 3
T_MIN, T_MAX = 150.0, 30000.0


# --- blackbody pieces (arrays allowed) ------------------------------------------
def planck_rel(band: str, temperature):
    """B_nu(T) in arbitrary units common to both bands (nu^3 / (e^x - 1))."""
    x = (_X1 if band == "w1" else _X2) / np.asarray(temperature, dtype=float)
    return (_NU[band] / _NU["w2"]) ** 3 / np.expm1(x)


def rho(temperature):
    """F_nu(W1) / F_nu(W2) of a blackbody."""
    t = np.asarray(temperature, dtype=float)
    return _NU3 * np.expm1(_X2 / t) / np.expm1(_X1 / t)


def gamma(band: str, temperature):
    """d ln B_nu / d ln T = x e^x / (e^x - 1)."""
    x = (_X1 if band == "w1" else _X2) / np.asarray(temperature, dtype=float)
    return x / -np.expm1(-x)


def beta_to_s(beta: float, t_pivot: float) -> float:
    """Convert d log T / d log e2 into d log T / d log L (L = Omega sigma T^4,
    Omega from the W2 flux): s = beta / (1 + beta (4 - gamma_2))."""
    den = 1.0 + beta * (4.0 - float(gamma("w2", t_pivot)))
    return float(beta / den) if den > 0 else float("nan")


# --- estimator --------------------------------------------------------------------
@dataclass(frozen=True)
class SlopeSettings:
    """Fixed before the frozen light curves were opened (PREDICTIONS_5b.md 8)."""

    n_ref: int = 3  # visits that define the reference (faint) state
    min_snr: float = 3.0  # combined excess S/N of a used visit
    min_used: int = 6  # visits needed for a slope
    sigma_scale: float = 1.0  # multiplies the fit uncertainty (from the toy tests)
    separation: float = 3.0  # 'measurable': beta_A,min >= separation * sigma_beta
    hyst_min_each: int = 4  # used visits needed on each branch
    hyst_max_sigma: float = 0.10  # 'measurable' hysteresis: sigma_D below this
    beta_bound: float = 3.0


SETTINGS = SlopeSettings()


def _prepare(f1, s1, f2, s2, settings: SlopeSettings) -> dict | None:
    f1, s1, f2, s2 = (np.asarray(x, dtype=float) for x in (f1, s1, f2, s2))
    n = f1.size
    if n < settings.n_ref + settings.min_used:
        return None
    z = f1 / np.median(f1) + f2 / np.median(f2)
    ref = np.argsort(z, kind="stable")[: settings.n_ref]
    is_ref = np.zeros(n, dtype=bool)
    is_ref[ref] = True
    w1, w2 = 1.0 / s1[ref] ** 2, 1.0 / s2[ref] ** 2
    a1, a2 = float(np.sum(w1 * f1[ref]) / w1.sum()), float(np.sum(w2 * f2[ref]) / w2.sum())
    sa1, sa2 = float(w1.sum() ** -0.5), float(w2.sum() ** -0.5)
    e1, e2 = f1 - a1, f2 - a2
    snr = (e1 / np.hypot(s1, sa1) + e2 / np.hypot(s2, sa2)) / np.sqrt(2.0)
    used = ~is_ref & (snr >= settings.min_snr) & (e1 > 0) & (e2 > 0)
    # central difference of the normalised sum over all visits: rising or declining
    dz = np.empty(n)
    dz[1:-1] = z[2:] - z[:-2]
    dz[0], dz[-1] = z[1] - z[0], z[-1] - z[-2]
    return {"f1": f1, "s1": s1, "f2": f2, "s2": s2, "a": (a1, a2), "sa": (sa1, sa2),
            "used": used, "is_ref": is_ref, "rising": dz > 0, "z": z}


def _curve_residuals(theta, f1, s1, f2, s2, a, sa, x_pivot, branch, n_beta):
    """Residuals of the excess-flux curve; the last two entries are the priors
    on the shifts of the reference state."""
    ln_t = theta[0]
    betas = theta[1:1 + n_beta]
    d1, d2 = theta[1 + n_beta], theta[2 + n_beta]
    e1 = f1 - (a[0] + d1)
    e2 = np.maximum(f2 - (a[1] + d2), 1e-12 * a[1])
    beta = betas[branch]
    t_c = np.clip(np.exp(ln_t + beta * (np.log(e2) - x_pivot)), T_MIN, T_MAX)
    r = rho(t_c)
    slope = r * (1.0 + beta * (gamma("w1", t_c) - gamma("w2", t_c)))  # d e1 / d e2
    res = (e1 - e2 * r) / np.sqrt(s1**2 + (slope * s2) ** 2)
    return np.concatenate([res, [d1 / sa[0], d2 / sa[1]]])


def _fit(prep: dict, branch: np.ndarray, n_beta: int, x_pivot: float,
         settings: SlopeSettings) -> dict | None:
    u = prep["used"]
    args = (prep["f1"][u], prep["s1"][u], prep["f2"][u], prep["s2"][u], prep["a"], prep["sa"],
            x_pivot, branch, n_beta)
    e1, e2 = prep["f1"][u] - prep["a"][0], prep["f2"][u] - prep["a"][1]
    ratio = float(np.sum(e1) / np.sum(e2))
    lo, hi = float(rho(T_MIN)), float(rho(T_MAX))
    t0 = brentq(lambda t: rho(t) - ratio, T_MIN, T_MAX) if lo < ratio < hi else 1000.0
    theta0 = np.array([np.log(t0)] + [0.0] * n_beta + [0.0, 0.0])
    bound = settings.beta_bound
    lower = [np.log(T_MIN)] + [-bound] * n_beta + [-np.inf, -np.inf]
    upper = [np.log(T_MAX)] + [bound] * n_beta + [np.inf, np.inf]
    try:
        sol = least_squares(_curve_residuals, theta0, args=args, bounds=(lower, upper),
                            x_scale=[0.1] + [0.1] * n_beta + [prep["sa"][0], prep["sa"][1]])
    except (ValueError, FloatingPointError):
        return None
    jtj = sol.jac.T @ sol.jac
    try:
        cov = np.linalg.inv(jtj)
    except np.linalg.LinAlgError:
        return None
    if not np.all(np.isfinite(cov)) or np.any(np.diag(cov) <= 0):
        return None
    return {"theta": sol.x, "cov": cov, "chi2": float(2.0 * sol.cost),
            "dof": int(u.sum() - 1 - n_beta)}


def heating_in_place_excess(e2, a2: float, t_ref: float):
    """Noise-free W1 excess of dust heated in place (fixed area) whose W2 flux
    in the reference state is ``a2`` at temperature ``t_ref``."""
    e2 = np.asarray(e2, dtype=float)
    x_ref = _X2 / t_ref
    x_k = np.log1p(np.expm1(x_ref) / (1.0 + e2 / a2))
    t_k = _X2 / x_k
    omega = a2 / planck_rel("w2", t_ref)
    return omega * (planck_rel("w1", t_k) - planck_rel("w1", t_ref))


def _beta_in_place_min(prep: dict, t_pivot: float, x_pivot: float,
                       settings: SlopeSettings) -> float:
    """Smallest beta that heating in place can give for this object: all of the
    reference flux is hot dust (no host), and the reference temperature is set
    so that the model has the fitted pivot temperature."""
    u = prep["used"]
    e2 = prep["f2"][u] - prep["a"][1]
    a1, a2 = prep["a"]
    branch = np.zeros(int(u.sum()), dtype=int)

    def fit_model(t_ref: float):
        e1 = heating_in_place_excess(e2, a2, t_ref)
        fake = dict(prep, f1=prep["f1"].copy())
        fake["f1"][u] = a1 + e1
        return _fit(fake, branch, 1, x_pivot, settings)

    def miss(t_ref: float) -> float:
        fit = fit_model(t_ref)
        return (fit["theta"][0] if fit else np.log(T_MAX)) - np.log(t_pivot)

    lo, hi = 200.0, min(t_pivot, 5000.0)
    try:
        if miss(lo) >= 0 or miss(hi) <= 0:
            return float("nan")
        t_ref = brentq(miss, lo, hi, xtol=1.0)
    except (ValueError, RuntimeError):
        return float("nan")
    fit = fit_model(t_ref)
    return float(fit["theta"][1]) if fit else float("nan")


def measure(f1, s1, f2, s2, settings: SlopeSettings = SETTINGS,
            with_in_place: bool = True) -> dict:
    """Slope and hysteresis statistics of one visit-level light curve (fluxes
    in any common unit, time ordered)."""
    out = {"n_visits": int(np.asarray(f1).size), "n_used": 0, "status": "too few visits",
           "beta": np.nan, "sigma_beta": np.nan, "t_pivot_k": np.nan, "s": np.nan,
           "beta_in_place_min": np.nan, "r": np.nan, "measurable": False, "class": "",
           "chi2": np.nan, "dof": 0, "n_rising": 0, "n_declining": 0,
           "beta_rise": np.nan, "beta_decline": np.nan, "d_hyst": np.nan,
           "sigma_d": np.nan, "hyst_measurable": False, "max_excess_ratio": np.nan}
    prep = _prepare(f1, s1, f2, s2, settings)
    if prep is None:
        return out
    u = prep["used"]
    out["n_used"] = int(u.sum())
    if u.sum() < settings.min_used:
        out["status"] = "too few used visits"
        return out
    e2 = prep["f2"][u] - prep["a"][1]
    out["max_excess_ratio"] = float(e2.max() / prep["a"][1])
    x_pivot = float(np.mean(np.log(e2)))
    fit = _fit(prep, np.zeros(int(u.sum()), dtype=int), 1, x_pivot, settings)
    if fit is None:
        out["status"] = "fit failed"
        return out
    out.update(status="ok", beta=float(fit["theta"][1]),
               sigma_beta=float(settings.sigma_scale * np.sqrt(fit["cov"][1, 1])),
               t_pivot_k=float(np.exp(fit["theta"][0])), chi2=fit["chi2"], dof=fit["dof"])
    out["s"] = beta_to_s(out["beta"], out["t_pivot_k"])
    if with_in_place:
        b_min = _beta_in_place_min(prep, out["t_pivot_k"], x_pivot, settings)
        out["beta_in_place_min"] = b_min
        if np.isfinite(b_min) and b_min > 0:
            out["r"] = out["beta"] / b_min
            out["measurable"] = bool(out["sigma_beta"] * settings.separation <= b_min)
    out["class"] = classify(out["beta"], out["sigma_beta"], out["beta_in_place_min"])
    # hysteresis: shared temperature at the largest excess, one exponent per branch
    rising = prep["rising"][u]
    out["n_rising"], out["n_declining"] = int(rising.sum()), int((~rising).sum())
    if min(out["n_rising"], out["n_declining"]) >= settings.hyst_min_each:
        fit2 = _fit(prep, np.where(rising, 0, 1), 2, float(np.log(e2.max())), settings)
        if fit2 is not None:
            c = fit2["cov"]
            out["beta_rise"], out["beta_decline"] = float(fit2["theta"][1]), float(fit2["theta"][2])
            out["d_hyst"] = out["beta_decline"] - out["beta_rise"]
            out["sigma_d"] = float(settings.sigma_scale
                                   * np.sqrt(max(c[1, 1] + c[2, 2] - 2.0 * c[1, 2], 0.0)))
            out["hyst_measurable"] = bool(0 < out["sigma_d"] <= settings.hyst_max_sigma)
    return out


def classify(beta: float, sigma: float, beta_in_place_min: float) -> str:
    """Per-object class of PREDICTIONS_5b.md section 6.2."""
    if not (np.isfinite(beta) and np.isfinite(sigma) and sigma > 0):
        return ""
    if beta <= -2.0 * sigma:
        return "negative"
    if beta >= 2.0 * sigma:
        return "in place"  # the '0.25 side': not consistent with a pure boundary shift
    if np.isfinite(beta_in_place_min) and beta <= beta_in_place_min - 2.0 * sigma:
        return "boundary"  # the '0 side'
    return "unclear"


# --- toy models --------------------------------------------------------------------
@dataclass(frozen=True)
class ToyModel:
    """One dust geometry. All fluxes are relative to the total W2 flux at
    driver luminosity 1.

    kind 'boundary': the inner edge follows the sublimation radius when the
    luminosity rises and returns to it with e-folding time ``tau_re_yr`` when
    the luminosity falls (0 = instantly, inf = never). Dust at the edge has
    T = T_sub (R_sub / R_in)^(1/2); the area scales as R_in^2.
    kind 'fixed': constant area, T = t_fixed * L^(1/4).
    ``f_fixed`` mixes a fixed-area layer (W2 flux fraction at L = 1) into a
    'boundary' model. ``host_w1`` is the constant, colour-zero host as a
    fraction of the total W1 flux at L = 1.
    """

    kind: str = "boundary"
    tau_re_yr: float = 0.0
    t_sub: float = 1500.0
    t_fixed: float = 1200.0
    f_fixed: float = 0.0
    host_w1: float = 0.0
    lag_days: float = 60.0  # delay and width of the dust response
    lag_w2_extra_days: float = 0.0  # W2 responds this much later than W1
    lag_fixed_days: float = 180.0  # the fixed layer of a mix sits further out


HOST_F1_OVER_F2 = 309.540 / 171.787  # W1 - W2 = 0 (Vega)


def drw(rng, t_days, sigma_mag: float, tau_days: float = 1000.0):
    """Damped random walk in magnitude; returns L / L_mean-ish = 10^(-0.4 m)."""
    t = np.asarray(t_days, dtype=float)
    m = np.empty(t.size)
    m[0] = rng.normal(0.0, sigma_mag)
    for i in range(1, t.size):
        a = np.exp(-(t[i] - t[i - 1]) / tau_days)
        m[i] = a * m[i - 1] + rng.normal(0.0, sigma_mag * np.sqrt(1.0 - a * a))
    return 10.0 ** (-0.4 * m)


def _respond(t, ell, lag_days: float):
    """Driver seen by the dust: boxcar of width ``lag_days`` ending ``lag_days/2``
    ago (mean delay = lag_days)."""
    if lag_days <= 0:
        return ell
    dt = t[1] - t[0]
    n = max(int(round(lag_days / dt)), 1)
    kernel = np.ones(n) / n
    sm = np.convolve(ell, kernel, mode="full")[: ell.size]
    sm[: n] = ell[0]
    shift = int(round(0.5 * lag_days / dt))
    return np.concatenate([np.full(shift, sm[0]), sm[: sm.size - shift]]) if shift else sm


def boundary_radius(t_days, ell, tau_re_yr: float):
    """R_in / R_sub(L = 1) for a driver ``ell`` on an even time grid."""
    rsub = np.sqrt(ell)
    r = np.empty(rsub.size)
    r[0] = rsub[0]
    dt = t_days[1] - t_days[0]
    decay = 0.0 if tau_re_yr == 0 else (1.0 if np.isinf(tau_re_yr)
                                        else np.exp(-dt / (tau_re_yr * DAYS_PER_YEAR)))
    for i in range(1, rsub.size):
        r[i] = rsub[i] if rsub[i] >= r[i - 1] else rsub[i] + (r[i - 1] - rsub[i]) * decay
    return r


def _layer(model: ToyModel, kind: str, t, ell, lag: float):
    seen = {b: _respond(t, ell, lag + (model.lag_w2_extra_days if b == "w2" else 0.0))
            for b in ("w1", "w2")}
    out = {}
    for b in ("w1", "w2"):
        if kind == "fixed":
            temp = model.t_fixed * seen[b] ** 0.25
            out[b] = planck_rel(b, temp) / planck_rel("w2", model.t_fixed)
        else:
            r = boundary_radius(t, seen[b], model.tau_re_yr)
            temp = model.t_sub * seen[b] ** 0.25 / np.sqrt(r)
            out[b] = r**2 * planck_rel(b, temp) / planck_rel("w2", model.t_sub)
    return out


def toy_fluxes(model: ToyModel, t_days, ell) -> dict:
    """Total W1 and W2 flux on the grid, total W2 = 1 at L = 1 before the host
    is added (the host adds its own share on top)."""
    if model.kind == "fixed":
        dust = _layer(model, "fixed", t_days, ell, model.lag_days)
    else:
        dust = _layer(model, "boundary", t_days, ell, model.lag_days)
        if model.f_fixed > 0:
            far = _layer(model, "fixed", t_days, ell, model.lag_fixed_days)
            dust = {b: (1.0 - model.f_fixed) * dust[b] + model.f_fixed * far[b] for b in dust}
    if model.host_w1 <= 0:
        return dust
    t_one = model.t_fixed if model.kind == "fixed" else model.t_sub
    dust_w1_at_one = float(rho(t_one))
    if model.kind != "fixed" and model.f_fixed > 0:
        dust_w1_at_one = ((1.0 - model.f_fixed) * float(rho(model.t_sub))
                          + model.f_fixed * float(rho(model.t_fixed)))
    host1 = model.host_w1 / (1.0 - model.host_w1) * dust_w1_at_one
    return {"w1": dust["w1"] + host1, "w2": dust["w2"] + host1 / HOST_F1_OVER_F2}


def synthetic_visits(rng, model: ToyModel, visit_mjd, err_mag: dict, sigma_mag: float,
                     ramp_mag: float = 0.0, ramp_years: float = 3.0, tau_days: float = 1000.0,
                     burn_in_days: float = 2500.0, dt_days: float = 10.0,
                     noise_correlation: float = 0.0) -> dict:
    """One synthetic object observed at ``visit_mjd`` with Gaussian visit errors
    ``err_mag[b]`` (mag). The driver is a damped random walk of asymptotic rms
    ``sigma_mag`` plus an optional one-sided ramp of ``ramp_mag`` that starts at
    a random time in the observed span."""
    visit_mjd = np.asarray(visit_mjd, dtype=float)
    t = np.arange(visit_mjd[0] - burn_in_days, visit_mjd[-1] + dt_days, dt_days)
    ell = drw(rng, t, sigma_mag, tau_days)
    if ramp_mag:
        start = rng.uniform(visit_mjd[0], visit_mjd[-1] - ramp_years * DAYS_PER_YEAR)
        frac = np.clip((t - start) / (ramp_years * DAYS_PER_YEAR), 0.0, 1.0)
        ell = ell * 10.0 ** (0.4 * ramp_mag * frac)
    flux = toy_fluxes(model, t, ell)
    common = rng.normal(size=visit_mjd.size)
    out = {"mjd": visit_mjd}
    for b in ("w1", "w2"):
        true = np.interp(visit_mjd, t, flux[b])
        frac_err = np.log(10.0) / 2.5 * np.broadcast_to(err_mag[b], visit_mjd.shape)
        own = rng.normal(size=visit_mjd.size)
        noise = (np.sqrt(noise_correlation) * common + np.sqrt(1.0 - noise_correlation) * own)
        out[f"{b}_true"] = true
        out[f"{b}_flux"] = true * (1.0 + frac_err * noise)
        out[f"{b}_flux_err"] = out[f"{b}_flux"] * frac_err
    return out
