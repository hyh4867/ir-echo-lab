"""Slow-rise detector: detector (b) of stage 4 (docs/PIPELINE.md section 8).

Detector (a) (flare.py) compares one visit with the earlier ones, so a rise
that never jumps by 5 sigma between consecutive states is absorbed into its
baseline. This detector instead fits the whole visit-level light curve with
"constant, then rising" models and asks whether the rise is significant and
large in both bands.

Models, per band (t in years, m in mag):

    constant            m = a
    linear / hinge      m = a + s * max(0, t - t_on)         t_on scanned over visits
    exponential         f = a + c * exp((t - t_end) / tau)   in relative flux, tau on a grid

Each rising model is compared with the constant by an F-test,

    F = (chi2_const - chi2_model) / (chi2_model / (n - 2)),   p = P(F_{1, n-2} > F),

which uses the scatter around the model rather than the absolute size of the
error bars. The p-value ignores the scan over t_on and tau, so the false-alarm
rate is measured empirically (scripts/run_stage34.py), not read off p.

Decision (fixed in config.SlowRiseCriteria before the injection tests were
run): take the candidate model whose worse band has the smallest p. The light
curve is a slow-rise candidate when, in both W1 and W2,

    p <= p_max,
    fitted cumulative brightening >= min_cum_mag,
    median(visits before the onset) - median(last 3 visits) >= min_cum_mag.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.stats import f as f_dist
from scipy.stats import norm

from .config import SLOW_RISE, SlowRiseCriteria

DAYS_PER_YEAR = 365.25


def mann_kendall_z(y) -> float:
    """Mann-Kendall trend statistic Z (positive = increasing y), no tie correction."""
    y = np.asarray(y, dtype=float)
    n = y.size
    if n < 4:
        return float("nan")
    s = 0.0
    for i in range(n - 1):
        s += np.sign(y[i + 1:] - y[i]).sum()
    var = n * (n - 1) * (2 * n + 5) / 18.0
    if s > 0:
        return float((s - 1) / np.sqrt(var))
    if s < 0:
        return float((s + 1) / np.sqrt(var))
    return 0.0


def theil_sen_slope(t, y) -> float:
    """Median of all pairwise slopes."""
    t, y = np.asarray(t, dtype=float), np.asarray(y, dtype=float)
    i, j = np.triu_indices(t.size, k=1)
    dt = t[j] - t[i]
    ok = dt > 0
    return float(np.median((y[j] - y[i])[ok] / dt[ok])) if ok.any() else float("nan")


def _wls(y: np.ndarray, w: np.ndarray, g: np.ndarray) -> tuple[float, float, float]:
    """Weighted fit y = a + s*g. Returns (a, s, chi2)."""
    sw = w.sum()
    gm, ym = (w * g).sum() / sw, (w * y).sum() / sw
    sgg = (w * (g - gm) ** 2).sum()
    if sgg <= 0:
        return ym, 0.0, float((w * (y - ym) ** 2).sum())
    s = (w * (g - gm) * (y - ym)).sum() / sgg
    a = ym - s * gm
    return float(a), float(s), float((w * (y - a - s * g) ** 2).sum())


def _f_p(chi2_const: float, chi2_model: float, n: int) -> tuple[float, float]:
    if chi2_model <= 0 or n <= 2:
        return float("inf"), 0.0
    f_stat = max(chi2_const - chi2_model, 0.0) / (chi2_model / (n - 2))
    return float(f_stat), float(f_dist.sf(f_stat, 1, n - 2))


def candidate_models(n: int, criteria: SlowRiseCriteria) -> list[tuple[str, float]]:
    """(kind, parameter): ('linear', 0), ('hinge', k) with k the first rising
    visit, ('exp', tau_years)."""
    cands: list[tuple[str, float]] = [("linear", 0.0)]
    for k in range(criteria.min_baseline_visits, n - criteria.min_ramp_visits + 1):
        cands.append(("hinge", float(k)))
    cands += [("exp", float(tau)) for tau in criteria.tau_grid_yr]
    return cands


def _fit_band(t: np.ndarray, m: np.ndarray, err: np.ndarray,
              cands: list[tuple[str, float]]) -> list[dict]:
    n = t.size
    w = 1.0 / err**2
    m0 = (w * m).sum() / w.sum()
    chi2_const_mag = float((w * (m - m0) ** 2).sum())
    # relative flux for the exponential model
    f = 10.0 ** (-0.4 * (m - np.median(m)))
    ferr = f * np.log(10.0) / 2.5 * err
    wf = 1.0 / ferr**2
    f0 = (wf * f).sum() / wf.sum()
    chi2_const_flux = float((wf * (f - f0) ** 2).sum())
    out = []
    for kind, par in cands:
        if kind in ("linear", "hinge"):
            k = int(par)
            t_on = t[0] if k == 0 else t[k - 1]
            g = np.clip(t - t_on, 0.0, None)
            a, s, chi2 = _wls(m, w, g)
            dm = -s * (g[-1] - g[0])  # positive = brightening
            f_stat, p = _f_p(chi2_const_mag, chi2, n)
        else:
            g = np.exp((t - t[-1]) / par)
            a, c, chi2 = _wls(f, wf, g)
            f_start, f_end = a + c * g[0], a + c * g[-1]
            dm = 2.5 * np.log10(f_end / f_start) if f_start > 0 and f_end > 0 else np.nan
            f_stat, p = _f_p(chi2_const_flux, chi2, n)
        if not np.isfinite(dm) or dm <= 0:
            p, f_stat = 1.0, 0.0  # not a brightening
        out.append({"f": f_stat, "p": p, "dm_fit": float(dm) if np.isfinite(dm) else float("nan")})
    return out


@dataclass
class SlowRiseResult:
    flagged: bool
    n_visits: int
    model: str | None = None  # 'linear', 'hinge' or 'exp'
    onset_index: int | None = None  # first rising visit (hinge), 0 (linear), None (exp)
    onset_mjd: float = float("nan")  # last visit before the rise
    tau_yr: float = float("nan")
    p: dict[str, float] = field(default_factory=dict)
    f_stat: dict[str, float] = field(default_factory=dict)
    dm_fit: dict[str, float] = field(default_factory=dict)
    dm_cum: dict[str, float] = field(default_factory=dict)
    mk_z: dict[str, float] = field(default_factory=dict)  # positive = brightening
    theil_sen_mag_per_yr: dict[str, float] = field(default_factory=dict)  # positive = brightening
    criteria: SlowRiseCriteria = field(default_factory=lambda: SLOW_RISE)

    def as_dict(self) -> dict:
        d = {k: getattr(self, k) for k in ("flagged", "n_visits", "model", "onset_index",
                                           "onset_mjd", "tau_yr")}
        for k in ("p", "f_stat", "dm_fit", "dm_cum", "mk_z", "theil_sen_mag_per_yr"):
            for b, val in getattr(self, k).items():
                d[f"{b}_{k}"] = val
        return d


def detect_slow_rise(mjd, w1mag, w1err, w2mag, w2err,
                     criteria: SlowRiseCriteria = SLOW_RISE) -> SlowRiseResult:
    """Run detector (b) on a visit-level light curve (time ordered)."""
    mjd = np.asarray(mjd, dtype=float)
    mags = {"w1": np.asarray(w1mag, dtype=float), "w2": np.asarray(w2mag, dtype=float)}
    errs = {"w1": np.asarray(w1err, dtype=float), "w2": np.asarray(w2err, dtype=float)}
    ok = np.isfinite(mjd)
    for b in mags:
        ok &= np.isfinite(mags[b]) & np.isfinite(errs[b]) & (errs[b] > 0)
    mjd = mjd[ok]
    mags = {b: v[ok] for b, v in mags.items()}
    errs = {b: v[ok] for b, v in errs.items()}
    n = mjd.size
    if n < criteria.min_visits:
        return SlowRiseResult(flagged=False, n_visits=n, criteria=criteria)
    if np.any(np.diff(mjd) < 0):
        raise ValueError("visits must be in time order")
    t = (mjd - mjd[0]) / DAYS_PER_YEAR
    cands = candidate_models(n, criteria)
    fits = {b: _fit_band(t, mags[b], errs[b], cands) for b in mags}
    worst_p = np.array([max(fits["w1"][i]["p"], fits["w2"][i]["p"]) for i in range(len(cands))])
    best = int(np.argmin(worst_p))
    kind, par = cands[best]
    k = int(par) if kind in ("linear", "hinge") else None
    n_before = max(k, 3) if k is not None else 3
    res = SlowRiseResult(flagged=False, n_visits=n, model=kind, onset_index=k,
                         onset_mjd=float(mjd[max(k - 1, 0)]) if k is not None else float("nan"),
                         tau_yr=par if kind == "exp" else float("nan"), criteria=criteria)
    passes = True
    for b in mags:
        fit = fits[b][best]
        cum = float(np.median(mags[b][:n_before]) - np.median(mags[b][-3:]))
        res.p[b], res.f_stat[b], res.dm_fit[b], res.dm_cum[b] = fit["p"], fit["f"], fit["dm_fit"], cum
        res.mk_z[b] = -mann_kendall_z(mags[b])
        res.theil_sen_mag_per_yr[b] = -theil_sen_slope(t, mags[b])
        passes &= (fit["p"] <= criteria.p_max and fit["dm_fit"] >= criteria.min_cum_mag
                   and cum >= criteria.min_cum_mag)
    res.flagged = bool(passes)
    return res


def model_curve(mjd, mag, err, result: SlowRiseResult) -> np.ndarray:
    """Fitted magnitudes of the selected model for one band (for plotting)."""
    mjd, mag, err = (np.asarray(x, dtype=float) for x in (mjd, mag, err))
    t = (mjd - mjd[0]) / DAYS_PER_YEAR
    w = 1.0 / err**2
    if result.model in ("linear", "hinge"):
        k = int(result.onset_index)
        g = np.clip(t - (t[0] if k == 0 else t[k - 1]), 0.0, None)
        a, s, _ = _wls(mag, w, g)
        return a + s * g
    if result.model == "exp":
        ref = np.median(mag)
        f = 10.0 ** (-0.4 * (mag - ref))
        ferr = f * np.log(10.0) / 2.5 * err
        g = np.exp((t - t[-1]) / result.tau_yr)
        a, c, _ = _wls(f, 1.0 / ferr**2, g)
        return ref - 2.5 * np.log10(a + c * g)
    return np.full(mjd.size, np.nan)
