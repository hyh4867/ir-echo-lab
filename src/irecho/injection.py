"""Stage-4 injection tests (docs/PIPELINE.md section 9).

A host light curve (a control galaxy or a comparison star, differentially
corrected) is cut into segments of ``L`` baseline visits followed by
``post_visits`` visits. A brightening of amplitude A (mag, the same in W1 and
W2) is added from the onset on, either at once or as a linear ramp that
reaches A after ``rise`` years. The same segments without injection give the
false-alarm rate.
"""

from __future__ import annotations

import numpy as np

from .config import BANDS, FLARE, INJECTION, SLOW_RISE, InjectionGrid
from .flare import detect_flares
from .slowrise import detect_slow_rise, mann_kendall_z, theil_sen_slope
from .stats import robust_std

DAYS_PER_YEAR = 365.25
# Stage 5b, Q1: a "slow rise" is a detector-(b) flag whose fitted brightening is
# at least this large in both bands (docs/PREDICTIONS_5b.md section 3).
Q1_MIN_MAG = 0.2
SEGMENT_KEYS = ("mjd", "w1_mag", "w1_err", "w2_mag", "w2_err")

FEATURE_NAMES = tuple(
    f"{b}_{name}" for name in (
        "max_brightening_mag", "max_brightening_sigma", "theil_sen_change_mag", "mk_z",
        "cum_change_mag", "max_step_mag", "excess_scatter", "von_neumann",
    ) for b in BANDS
) + ("w1_w2_correlation", "color_change_mag")


def segment_onsets(n_visits: int, baseline_visits: int, post_visits: int) -> list[int]:
    """Onset indices k for which visits [k - L, k + P) exist."""
    return list(range(baseline_visits, n_visits - post_visits + 1))


def cut_segment(v: dict[str, np.ndarray], onset: int, baseline_visits: int,
                post_visits: int) -> dict[str, np.ndarray]:
    sl = slice(onset - baseline_visits, onset + post_visits)
    return {k: np.array(v[k][sl], dtype=float) for k in SEGMENT_KEYS}


def inject(seg: dict[str, np.ndarray], onset_in_segment: int, amp_mag: float,
           rise_years: float) -> dict[str, np.ndarray]:
    """Brighten both bands by ``amp_mag`` from ``onset_in_segment`` on. The
    ramp starts at the last baseline visit, so with ``rise_years`` > 0 the first
    post-onset visit has only reached a fraction of the amplitude."""
    out = {k: v.copy() for k, v in seg.items()}
    mjd = seg["mjd"]
    if rise_years <= 0:
        frac = np.ones(mjd.size)
    else:
        frac = np.clip((mjd - mjd[onset_in_segment - 1]) / (rise_years * DAYS_PER_YEAR), 0.0, 1.0)
    frac[:onset_in_segment] = 0.0
    for b in BANDS:
        out[f"{b}_mag"] = seg[f"{b}_mag"] - amp_mag * frac
    return out


def run_detectors(seg: dict[str, np.ndarray], onset_in_segment: int | None) -> dict[str, bool]:
    """Detector (a) and (b) on one segment. For an injected segment, (a) only
    counts when a visit at or after the onset is flagged."""
    res_a = detect_flares(seg["mjd"], seg["w1_mag"], seg["w1_err"], seg["w2_mag"],
                          seg["w2_err"], FLARE)
    if onset_in_segment is None:
        hit_a = bool(res_a.flagged.any())
    else:
        hit_a = bool(res_a.flagged[onset_in_segment:].any())
    res_b = detect_slow_rise(seg["mjd"], seg["w1_mag"], seg["w1_err"], seg["w2_mag"],
                             seg["w2_err"], SLOW_RISE)
    q1 = bool(res_b.flagged) and min(res_b.dm_fit.values()) >= Q1_MIN_MAG
    return {"a": hit_a, "b": bool(res_b.flagged), "q1": bool(q1)}


def features(seg: dict[str, np.ndarray]) -> np.ndarray:
    """Light-curve features for the anomaly-detection comparison. The reference
    level is the median of the first four visits (the shortest baseline)."""
    t = (seg["mjd"] - seg["mjd"][0]) / DAYS_PER_YEAR
    per_band: dict[str, list[float]] = {}
    resid = {}
    for b in BANDS:
        m, e = seg[f"{b}_mag"], seg[f"{b}_err"]
        ref = float(np.median(m[:4]))
        scat = max(robust_std(m[:4]), float(np.median(e)))
        bright = ref - m  # positive = brighter
        d = np.diff(m)
        var = float(np.var(m, ddof=1))
        per_band[b] = [
            float(bright.max()),
            float(bright.max() / scat),
            float(-theil_sen_slope(t, m) * (t[-1] - t[0])),
            float(-mann_kendall_z(m)),
            float(ref - np.median(m[-3:])),
            float((-d).max()),
            float(np.std(m, ddof=1) / np.median(e)),
            float(np.mean(d**2) / var) if var > 0 else 2.0,
        ]
        resid[b] = m - m.mean()
    feats = [per_band[b][i] for i in range(8) for b in BANDS]
    denom = np.sqrt((resid["w1"] ** 2).sum() * (resid["w2"] ** 2).sum())
    feats.append(float((resid["w1"] * resid["w2"]).sum() / denom) if denom > 0 else 0.0)
    color = seg["w1_mag"] - seg["w2_mag"]
    feats.append(float(np.median(color[-3:]) - np.median(color[:4])))
    return np.array(feats)


def run_grid(hosts: list[dict], grid: InjectionGrid = INJECTION,
             with_features: bool = True) -> list[dict]:
    """Every host x baseline length x onset x (null + amplitude x rise).

    ``hosts``: dicts with ``name``, ``kind`` and ``visits`` (a visit table).
    Returns one row per trial with the detector outcomes and the feature vector
    (``with_features=False`` leaves the features out; they are only needed for
    the anomaly-model comparison).
    """
    feats = features if with_features else (lambda seg: None)
    rows = []
    for host in hosts:
        v = host["visits"]
        n = v["mjd"].size
        for n_base in grid.baseline_visits:
            for onset in segment_onsets(n, n_base, grid.post_visits):
                seg = cut_segment(v, onset, n_base, grid.post_visits)
                common = {"host": host["name"], "kind": host["kind"], "baseline_visits": n_base,
                          "onset": onset}
                rows.append({**common, "amp": 0.0, "rise": -1.0, **run_detectors(seg, None),
                             "features": feats(seg)})
                for rise in grid.rise_years:
                    for amp in grid.amplitudes_mag:
                        inj = inject(seg, n_base, amp, rise)
                        rows.append({**common, "amp": amp, "rise": rise,
                                     **run_detectors(inj, n_base), "features": feats(inj)})
    return rows
