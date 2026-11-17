"""Stage 5a: a per-galaxy error floor from two independent photometries.

Problem (stage 3-4 report, remaining problem 2): the visit errors of a galaxy
nucleus are underestimated. The pipeline error is the scatter of the single
exposures inside a visit; it does not contain errors that are common to all
exposures of a visit but change from visit to visit (how much host light the
profile fit picks up at a given scan angle, deblending, residual zero point).

Estimate. The same visit is measured twice from the same exposures: by the
pipeline (profile-fit on single exposures, visit median) and by unTimely
(PSF fit on the visit coadd). The visit-to-visit scatter of the difference
d_k = m_pipe,k - m_unTimely,k does not contain the source's variability, only
measurement error:

    s_d^2 = sigma_pipe^2 + sigma_unT^2            (errors of the two measurements)

For a point source this scatter is what the two photometries give on a clean
target; it is measured on the comparison stars (``reference_scatter``). A
galaxy's excess over that is attributed to the pipeline measurement of the
galaxy:

    floor^2 = max(0, s_d(galaxy)^2 - s_ref^2)

This is conservative in one direction: part of the excess belongs to unTimely
(which also models an extended nucleus as a point source), so the pipeline
floor is, if anything, overestimated. The floor is one number per object and
band, added in quadrature to every visit error (``apply_floor``).

Nothing here looks at how a target changes with time: the input is the
difference of two measurements of the same visit.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import BANDS
from .photometry import LN10_OVER_2P5, color_w1_w2
from .stats import robust_std

MIN_PAIRS = 8  # visits measured by both photometries needed for a per-object estimate
STAR_MAG_EDGES = (10.5,)  # reference bins: W1 < 10.5 and W1 >= 10.5 (as the star draw)


@dataclass(frozen=True)
class PairScatter:
    n_pairs: int
    scatter_mag: float  # robust std of the differences (NaN below 2 pairs)
    median_diff_mag: float  # pipeline - unTimely, constant offset between the two
    usable: bool  # n_pairs >= MIN_PAIRS


def pair_scatter(mag_pipeline, mag_untimely, min_pairs: int = MIN_PAIRS) -> PairScatter:
    """Robust visit-to-visit scatter of (pipeline - unTimely) for matched visits.

    The robust estimator (1.4826 MAD) is used so that one visit with a bad
    coadd does not set the floor.
    """
    a = np.asarray(mag_pipeline, dtype=float)
    b = np.asarray(mag_untimely, dtype=float)
    d = (a - b)[np.isfinite(a) & np.isfinite(b)]
    if d.size < 2:
        return PairScatter(int(d.size), float("nan"),
                           float(np.median(d)) if d.size else float("nan"), False)
    return PairScatter(int(d.size), robust_std(d), float(np.median(d)), d.size >= min_pairs)


def reference_scatter(star_scatters, star_w1_mags, edges=STAR_MAG_EDGES) -> dict:
    """Median pair scatter of the comparison stars in W1 bins.

    Returns {'edges': [...], 'scatter_mag': [per bin], 'n_stars': [per bin]}.
    A bin without stars gets the median over all stars.
    """
    s = np.asarray(star_scatters, dtype=float)
    m = np.asarray(star_w1_mags, dtype=float)
    ok = np.isfinite(s) & np.isfinite(m)
    s, m = s[ok], m[ok]
    which = np.digitize(m, edges)
    overall = float(np.median(s)) if s.size else float("nan")
    scatter, counts = [], []
    for k in range(len(edges) + 1):
        sel = which == k
        counts.append(int(sel.sum()))
        scatter.append(float(np.median(s[sel])) if sel.any() else overall)
    return {"edges": [float(e) for e in edges], "scatter_mag": scatter, "n_stars": counts}


def reference_for(w1_mag: float, reference: dict) -> float:
    """Reference scatter for an object of magnitude ``w1_mag`` (objects outside
    the stars' range use the nearest bin)."""
    return float(reference["scatter_mag"][int(np.digitize(w1_mag, reference["edges"]))])


def floor_from_scatter(scatter_mag: float, reference_mag: float) -> float:
    """floor = sqrt(max(0, s^2 - s_ref^2)); NaN in, NaN out."""
    if not (np.isfinite(scatter_mag) and np.isfinite(reference_mag)):
        return float("nan")
    return float(np.sqrt(max(0.0, scatter_mag**2 - reference_mag**2)))


def adopt_floors(per_object: dict[str, dict[str, float]],
                 usable: dict[str, dict[str, bool]]) -> tuple[dict, dict]:
    """Final floor per object and band.

    An object with at least MIN_PAIRS matched visits keeps its own estimate.
    Any other object (no unTimely detection, too few pairs) gets the median of
    the usable objects in that band. Returns (floors, fallback) where
    ``fallback[band]`` is that median.
    """
    fallback = {}
    for b in BANDS:
        vals = [per_object[name][b] for name in per_object
                if usable[name][b] and np.isfinite(per_object[name][b])]
        fallback[b] = float(np.median(vals)) if vals else float("nan")
    floors = {}
    for name in per_object:
        floors[name] = {}
        for b in BANDS:
            own = per_object[name][b]
            good = usable[name][b] and np.isfinite(own)
            floors[name][b] = float(own) if good else fallback[b]
            floors[name][f"{b}_source"] = "own" if good else "sample median"
    return floors, fallback


def apply_floor(v: dict[str, np.ndarray], floor_mag: dict[str, float]) -> dict[str, np.ndarray]:
    """Add the floor in quadrature to the visit errors of a visit table.

    Magnitudes and fluxes are unchanged; ``w?_err``, ``w?_flux_err_jy`` and the
    colour error are recomputed. The input table is not modified. A NaN floor
    is refused: an object must not silently keep its too-small errors.
    """
    out = {k: np.array(val, copy=True) for k, val in v.items()}
    for b in BANDS:
        f = float(floor_mag[b])
        if not np.isfinite(f) or f < 0:
            raise ValueError(f"invalid error floor for {b}: {floor_mag[b]!r}")
        out[f"{b}_err"] = np.hypot(v[f"{b}_err"], f)
        out[f"{b}_flux_err_jy"] = v[f"{b}_flux_jy"] * LN10_OVER_2P5 * out[f"{b}_err"]
        out[f"{b}_floor"] = np.full(len(v[f"{b}_err"]), f)
    out["color"], out["color_err"] = color_w1_w2(out["w1_mag"], out["w1_err"],
                                                 out["w2_mag"], out["w2_err"])
    return out
