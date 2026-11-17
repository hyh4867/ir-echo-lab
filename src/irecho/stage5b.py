"""Stage 5b: pre-registered analysis of the frozen sample_v1 light curves.

Every definition here is the code form of docs/PREDICTIONS_5b.md (committed
before the frozen files were opened). Section numbers in comments refer to
that document.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy.stats import beta as beta_dist
from scipy.stats import fisher_exact, mannwhitneyu, norm, wilcoxon

from . import breath, sample_io, ztf
from .config import BANDS, FLARE, PROJECT_ROOT, SEED, SLOW_RISE
from .epochs import MEDIAN_SE_FACTOR
from .flare import detect_flares
from .injection import Q1_MIN_MAG
from .slowrise import DAYS_PER_YEAR, detect_slow_rise
from .stats import robust_std

FROZEN_DIR = PROJECT_ROOT / "data" / "frozen" / "sample_v1"
LN10_OVER_2P5 = np.log(10.0) / 2.5

# section 5 (ZTF)
ZTF_MIN_BINS = 6
ZTF_MIN_IR_CHANGE_MAG = 0.1
ZTF_BIN_ERR_FLOOR_MAG = 0.02
ZTF_MIN_CHANGE_MAG = 0.05
ZTF_MIN_SIGMA = 3.0
ZTF_DEC_LIMIT_DEG = -28.0  # only used to label why an object has no ZTF data
# section 6.1
LAMBDA_EDD_THRESHOLD = 0.02
# sections 3-7: minimum sample sizes
MIN_BRIGHTENED = 10
MIN_PER_GROUP = 5
MIN_HYSTERESIS = 8
N_BOOTSTRAP = 10000
COMPLETENESS_FLOOR = 0.2


# --- frozen light curves -------------------------------------------------------------
def read_frozen(name: str, with_floor: bool = True,
                frozen_dir: Path | None = None) -> dict[str, np.ndarray]:
    """Visit table of one target. ``with_floor=False`` rebuilds the errors
    without the galaxy error floor (reference results only). ``frozen_dir``
    selects another frozen sample (stage 6B: data/frozen/sample_v2)."""
    with open((frozen_dir or FROZEN_DIR) / f"{name}.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    v = {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}
    if not with_floor:
        for b in BANDS:
            v[f"{b}_err"] = np.hypot(v[f"{b}_err_stat"], v[f"{b}_offset_err"])
            v[f"{b}_flux_err_jy"] = v[f"{b}_flux_jy"] * LN10_OVER_2P5 * v[f"{b}_err"]
        v["color_err"] = np.hypot(v["w1_err"], v["w2_err"])
    return v


# --- detectors (sections 2-3) -----------------------------------------------------------
def run_detectors(v: dict) -> dict:
    a = detect_flares(v["mjd"], v["w1_mag"], v["w1_err"], v["w2_mag"], v["w2_err"], FLARE)
    b = detect_slow_rise(v["mjd"], v["w1_mag"], v["w1_err"], v["w2_mag"], v["w2_err"], SLOW_RISE)
    out = {"n_visits": int(v["mjd"].size), "a_flagged": a.detected,
           "a_first_index": a.onset_index if a.detected else None,
           "a_n_flagged": int(a.flagged.sum()),
           "a_peak_w1": a.peak("w1")[0], "a_peak_w2": a.peak("w2")[0],
           "b_flagged": b.flagged, "b_model": b.model or "", "b_onset_index": b.onset_index,
           "b_tau_yr": b.tau_yr}
    for band in BANDS:
        out[f"b_{band}_dm_fit"] = b.dm_fit.get(band, np.nan)
        out[f"b_{band}_dm_cum"] = b.dm_cum.get(band, np.nan)
        out[f"b_{band}_p"] = b.p.get(band, np.nan)
    dm = [out[f"b_{band}_dm_fit"] for band in BANDS]
    out["q1"] = bool(b.flagged and min(dm) >= Q1_MIN_MAG)
    out["q1_w2_only"] = bool(b.flagged and out["b_w2_dm_fit"] >= Q1_MIN_MAG)
    span = (v["mjd"][-1] - v["mjd"][0]) / DAYS_PER_YEAR
    if b.model in ("linear", "hinge"):
        k = b.onset_index or 0
        out["b_rise_years"] = float((v["mjd"][-1] - v["mjd"][max(k - 1, 0)]) / DAYS_PER_YEAR)
    elif b.model == "exp":
        out["b_rise_years"] = float(min(2.0 * b.tau_yr, span))
    else:
        out["b_rise_years"] = np.nan
    out["brightened"] = bool(a.detected or b.flagged)
    starts = []
    if a.detected:
        starts.append(a.onset_index)
    if b.flagged:
        starts.append(b.onset_index if b.onset_index is not None else 0)
    out["start_index"] = min(starts) if starts else None
    return out


# --- colour change (section 4) -------------------------------------------------------------
def colour_change(v: dict, start_index: int) -> dict:
    n = v["mjd"].size
    n_base = max(int(start_index), 3)
    out = {"n_base": n_base, "delta_color": np.nan, "delta_color_err": np.nan,
           "dm_w1": np.nan, "dm_w2": np.nan}
    if n - n_base < 3:
        return out
    mean_mag = 0.5 * (v["w1_mag"] + v["w2_mag"])
    rest = np.arange(n_base, n)
    bright = rest[np.argsort(mean_mag[rest], kind="stable")[:3]]
    base = np.arange(n_base)
    w = 1.0 / v["color_err"] ** 2

    def wmean(idx):
        return float(np.sum(w[idx] * v["color"][idx]) / w[idx].sum()), float(w[idx].sum() ** -0.5)

    c_bright, e_bright = wmean(bright)
    c_base, e_base = wmean(base)
    out.update(delta_color=c_bright - c_base, delta_color_err=float(np.hypot(e_bright, e_base)),
               dm_w1=float(np.median(v["w1_mag"][base]) - np.median(v["w1_mag"][bright])),
               dm_w2=float(np.median(v["w2_mag"][base]) - np.median(v["w2_mag"][bright])))
    return out


# --- optical (section 5) -------------------------------------------------------------------
def _group_se(mag, err) -> float:
    err = np.hypot(err, ZTF_BIN_ERR_FLOOR_MAG)
    scatter = MEDIAN_SE_FACTOR * robust_std(mag) / np.sqrt(mag.size) if mag.size > 1 else 0.0
    return float(max(scatter, np.sqrt(np.sum(err**2)) / mag.size))


def ztf_band(v: dict, path: Path) -> dict:
    """Optical change between the IR-faint and IR-bright thirds of the WISE
    visits that have a ZTF bin, for one ZTF band."""
    out = {"status": "", "n_rows_used": 0, "n_bins": 0, "ir_change_w2": np.nan,
           "dm_opt": np.nan, "dm_opt_err": np.nan, "verdict": ""}
    meta = path.with_suffix(".csv.meta.json")
    if not meta.exists():
        out["status"] = "not fetched"
        return out
    if not path.exists():
        out["status"] = "no rows returned"
        return out
    lc = ztf.read_clean(path)
    out["n_rows_used"] = int(lc["n_rows_used"])
    if lc["n_rows_used"] == 0:
        out["status"] = "no clean rows"
        return out
    bins = ztf.bin_to_visits(lc["mjd"], lc["mag"], v["mjd"])
    ok = np.flatnonzero(np.isfinite(bins["mag"]))
    out["n_bins"] = int(ok.size)
    if ok.size < ZTF_MIN_BINS:
        out["status"] = "fewer than 6 visits with a ZTF bin"
        return out
    mean_mag = 0.5 * (v["w1_mag"][ok] + v["w2_mag"][ok])
    order = np.argsort(mean_mag, kind="stable")
    n3 = max(2, ok.size // 3)
    bright, faint = ok[order[:n3]], ok[order[-n3:]]
    out["ir_change_w2"] = float(np.median(v["w2_mag"][faint]) - np.median(v["w2_mag"][bright]))
    if out["ir_change_w2"] < ZTF_MIN_IR_CHANGE_MAG:
        out["status"] = "IR change below 0.1 mag in the ZTF era"
        return out
    dm = float(np.median(bins["mag"][faint]) - np.median(bins["mag"][bright]))
    err = float(np.hypot(_group_se(bins["mag"][faint], bins["err"][faint]),
                         _group_se(bins["mag"][bright], bins["err"][bright])))
    out.update(status="usable", dm_opt=dm, dm_opt_err=err)
    if dm >= ZTF_MIN_CHANGE_MAG and dm / err >= ZTF_MIN_SIGMA:
        out["verdict"] = "accompanies"
    elif dm <= -ZTF_MIN_CHANGE_MAG and dm / err <= -ZTF_MIN_SIGMA:
        out["verdict"] = "opposite"
    else:
        out["verdict"] = "no change"
    return out


def ztf_compare(v: dict, name: str, has_object: dict[str, bool], dec_deg: float) -> dict:
    out = {}
    verdicts = []
    for band in ("g", "r"):
        if not has_object.get(band, False):
            res = {"status": ("outside the ZTF declination range"
                              if dec_deg < ZTF_DEC_LIMIT_DEG else "no ZTF object within 3 arcsec"),
                   "n_rows_used": 0, "n_bins": 0, "ir_change_w2": np.nan, "dm_opt": np.nan,
                   "dm_opt_err": np.nan, "verdict": ""}
        else:
            res = ztf_band(v, sample_io.S5_DIR / "ztf" / f"{name}_{band}.csv")
        for k, val in res.items():
            out[f"ztf_{band}_{k}"] = val
        if res["status"] == "usable":
            verdicts.append(res["verdict"])
    out["ztf_usable"] = bool(verdicts)
    out["ztf_verdict"] = ("accompanies" if "accompanies" in verdicts
                          else "opposite" if "opposite" in verdicts
                          else "no change" if verdicts else "")
    return out


# --- intervals and tests -------------------------------------------------------------------
def wilson(k: int, n: int, level: float = 0.95) -> tuple[float, float]:
    if n == 0:
        return float("nan"), float("nan")
    z = float(norm.ppf(0.5 + level / 2.0))
    p = k / n
    den = 1.0 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return float(max(0.0, centre - half)), float(min(1.0, centre + half))


def clopper_pearson(k: int, n: int, level: float = 0.95) -> tuple[float, float]:
    if n == 0:
        return float("nan"), float("nan")
    a = (1.0 - level) / 2.0
    lo = 0.0 if k == 0 else float(beta_dist.ppf(a, k, n - k + 1))
    hi = 1.0 if k == n else float(beta_dist.ppf(1.0 - a, k + 1, n - k))
    return lo, hi


def verdict_q1(lo: float, hi: float) -> str:
    """Section 3: the prediction is 5-20%."""
    if lo >= 0.05 and hi <= 0.20:
        return "지지"
    if hi < 0.05 or lo > 0.20:
        return "기각"
    return "판정불가"


def verdict_majority(k: int, n: int, n_min: int = MIN_BRIGHTENED) -> tuple[str, float, float]:
    """Sections 4-5: 'most' = more than half; Clopper-Pearson interval."""
    lo, hi = clopper_pearson(k, n)
    if n < n_min:
        return "판정불가 (표본 부족)", lo, hi
    if lo > 0.5:
        return "지지", lo, hi
    if hi < 0.5:
        return "기각", lo, hi
    return "판정불가", lo, hi


def nearest(grid, value: float, log: bool = False) -> float:
    g = np.asarray(grid, dtype=float)
    d = np.abs(np.log(g) - np.log(value)) if log else np.abs(g - value)
    return float(g[int(np.argmin(d))])


def corrected_fraction(weights, n_total: int, seed_tag: int) -> dict:
    """Section 3: f = sum(1/C_i) / N with an object-level bootstrap."""
    w = np.zeros(n_total)
    w[: len(weights)] = weights
    rng = np.random.default_rng([SEED, 53, seed_tag])
    boot = rng.choice(w, size=(N_BOOTSTRAP, n_total), replace=True).mean(axis=1)
    lo, hi = np.percentile(boot, [2.5, 97.5])
    return {"fraction": float(w.mean()), "ci95": [float(lo), float(hi)],
            "verdict_reference": verdict_q1(float(lo), float(hi))}


def compare_groups(r_above, r_below, seed_tag: int) -> dict:
    """Section 6.3: delta_r = median(below) - median(above), one-sided
    Mann-Whitney (above < below), bootstrap interval, verdict."""
    a, b = np.asarray(r_above, dtype=float), np.asarray(r_below, dtype=float)
    out = {"n_above": int(a.size), "n_below": int(b.size), "median_above": None,
           "median_below": None, "delta_r": None, "delta_r_ci95": None, "p_one_sided": None}
    if a.size:
        out["median_above"] = float(np.median(a))
    if b.size:
        out["median_below"] = float(np.median(b))
    if min(a.size, b.size) < MIN_PER_GROUP:
        out["verdict"] = "판정불가 (표본 부족)"
        if a.size and b.size:
            out["delta_r"] = out["median_below"] - out["median_above"]
        return out
    rng = np.random.default_rng([SEED, 54, seed_tag])
    boot = (np.median(rng.choice(b, size=(N_BOOTSTRAP, b.size)), axis=1)
            - np.median(rng.choice(a, size=(N_BOOTSTRAP, a.size)), axis=1))
    lo, hi = (float(x) for x in np.percentile(boot, [2.5, 97.5]))
    p = float(mannwhitneyu(a, b, alternative="less").pvalue)
    delta = out["median_below"] - out["median_above"]
    out.update(delta_r=float(delta), delta_r_ci95=[lo, hi], p_one_sided=p)
    if p < 0.05 and delta >= 0.5:
        out["verdict"] = "지지"
    elif hi < 0.5:
        out["verdict"] = "기각"
    else:
        out["verdict"] = "판정불가"
    return out


def hysteresis_test(d, d_min: float, seed_tag: int) -> dict:
    """Section 7.2: one-sided Wilcoxon signed-rank test for D > 0."""
    d = np.asarray(d, dtype=float)
    out = {"n": int(d.size), "median_d": float(np.median(d)) if d.size else None,
           "median_d_ci95": None, "p_one_sided": None, "d_min": d_min}
    if d.size < MIN_HYSTERESIS:
        out["verdict"] = "판정불가 (표본 부족)"
        return out
    rng = np.random.default_rng([SEED, 55, seed_tag])
    boot = np.median(rng.choice(d, size=(N_BOOTSTRAP, d.size)), axis=1)
    lo, hi = (float(x) for x in np.percentile(boot, [2.5, 97.5]))
    p = float(wilcoxon(d, alternative="greater").pvalue)
    out.update(median_d_ci95=[lo, hi], p_one_sided=p)
    if p < 0.05 and out["median_d"] > 0:
        out["verdict"] = "지지"
    elif hi < d_min:
        out["verdict"] = "기각"
    else:
        out["verdict"] = "판정불가"
    return out


def fisher_2x2(a: int, b: int, c: int, d: int) -> float | None:
    """Two-sided Fisher exact p for [[a, b], [c, d]]; None when a margin is 0."""
    if min(a + b, c + d, a + c, b + d) == 0:
        return None
    return float(fisher_exact([[a, b], [c, d]])[1])


def weighted_mean(values, sigmas) -> tuple[float, float]:
    v, s = np.asarray(values, dtype=float), np.asarray(sigmas, dtype=float)
    ok = np.isfinite(v) & np.isfinite(s) & (s > 0)
    if not ok.any():
        return float("nan"), float("nan")
    w = 1.0 / s[ok] ** 2
    return float(np.sum(w * v[ok]) / w.sum()), float(w.sum() ** -0.5)
