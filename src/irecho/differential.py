"""Differential photometry: visit-to-visit zero-point offsets from comparison
stars in the same field.

Stars within ~1.5 deg of a target are observed in the same WISE visits. If
the photometric zero point of a visit is off by ``o_k`` magnitudes, every
constant star is shifted by the same ``o_k``. The median residual of the stars
therefore estimates ``o_k``, and subtracting it removes the part of the
visit-to-visit scatter that is shared by the field (docs/THEORY.md section 8).

Sign convention: residual = magnitude - (the star's own median magnitude), so a
positive offset means "everything looks fainter in this visit". The corrected
magnitude is ``m - o_k``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import BANDS, BINNING
from .epochs import MEDIAN_SE_FACTOR, assign_epochs
from .photometry import LN10_OVER_2P5, color_w1_w2
from .stats import robust_std

MATCH_DAYS = 15.0  # a target visit uses the field visit nearest in time within this window
MIN_STARS = 3  # fewer stars than this in a visit -> no offset for that visit


@dataclass
class FieldOffsets:
    mjd: np.ndarray  # mean epoch of each field visit
    offset: dict[str, np.ndarray]  # per band, mag, positive = fainter
    offset_err: dict[str, np.ndarray]
    n_stars: dict[str, np.ndarray]

    def nearest(self, mjd: float, match_days: float = MATCH_DAYS) -> int | None:
        if self.mjd.size == 0:
            return None
        k = int(np.argmin(np.abs(self.mjd - mjd)))
        return k if abs(self.mjd[k] - mjd) <= match_days else None


def field_offsets(star_visits: list[dict[str, np.ndarray]],
                  gap_days: float = BINNING.gap_days,
                  min_stars: int = MIN_STARS) -> FieldOffsets:
    """Median residual of the comparison stars in every field visit.

    ``star_visits`` holds one visit-level table per star (keys ``mjd``,
    ``w1_mag``, ``w2_mag``). The uncertainty of the offset is the standard
    error of the median, sqrt(pi/2) * robust_std / sqrt(n).
    """
    star_visits = [v for v in star_visits if v["mjd"].size]
    if not star_visits:
        empty = np.zeros(0)
        return FieldOffsets(empty, {b: empty for b in BANDS}, {b: empty for b in BANDS},
                            {b: empty.astype(int) for b in BANDS})
    all_mjd = np.concatenate([v["mjd"] for v in star_visits])
    visit_of = assign_epochs(all_mjd, gap_days)
    n_visits = int(visit_of.max()) + 1
    mjd = np.array([all_mjd[visit_of == k].mean() for k in range(n_visits)])
    offset, offset_err, n_stars = {}, {}, {}
    for b in BANDS:
        resid = np.concatenate([v[f"{b}_mag"] - np.nanmedian(v[f"{b}_mag"]) for v in star_visits])
        off = np.full(n_visits, np.nan)
        err = np.full(n_visits, np.nan)
        cnt = np.zeros(n_visits, dtype=int)
        for k in range(n_visits):
            r = resid[(visit_of == k) & np.isfinite(resid)]
            cnt[k] = r.size
            if r.size >= min_stars:
                off[k] = np.median(r)
                err[k] = MEDIAN_SE_FACTOR * robust_std(r) / np.sqrt(r.size)
        offset[b], offset_err[b], n_stars[b] = off, err, cnt
    return FieldOffsets(mjd, offset, offset_err, n_stars)


def apply_offsets(v: dict[str, np.ndarray], offsets: FieldOffsets, *,
                  err_key: str = "err", floor_mag: dict[str, float] | None = None,
                  match_days: float = MATCH_DAYS) -> dict[str, np.ndarray]:
    """Subtract the field offsets from a visit-level table.

    Returns a new table with the same keys. ``w?_mag`` and ``w?_flux_jy`` are
    corrected; ``w?_err`` becomes the quadrature sum of ``w?_<err_key>``, the
    offset uncertainty and an optional systematic floor. Visits without a
    usable field visit are dropped.
    """
    floor_mag = floor_mag or {}
    n = v["mjd"].size
    idx = [offsets.nearest(t, match_days) for t in v["mjd"]]
    keep = np.array([
        k is not None and all(np.isfinite(offsets.offset[b][k]) for b in BANDS) for k in idx
    ], dtype=bool) if n else np.zeros(0, dtype=bool)
    out = {key: np.array(val)[keep] for key, val in v.items()}
    sel = np.array([k for k, ok in zip(idx, keep) if ok], dtype=int)
    for b in BANDS:
        base = f"{b}_{err_key}" if err_key != "err" else f"{b}_err"
        off = offsets.offset[b][sel]
        off_err = offsets.offset_err[b][sel]
        mag = out[f"{b}_mag"] - off
        err = np.sqrt(out[base] ** 2 + off_err**2 + floor_mag.get(b, 0.0) ** 2)
        flux = out[f"{b}_flux_jy"] * 10.0 ** (0.4 * off)
        out[f"{b}_mag"], out[f"{b}_err"] = mag, err
        out[f"{b}_flux_jy"] = flux
        out[f"{b}_flux_err_jy"] = flux * LN10_OVER_2P5 * err
        out[f"{b}_offset"] = off
        out[f"{b}_offset_err"] = off_err
    out["color"], out["color_err"] = color_w1_w2(out["w1_mag"], out["w1_err"],
                                                 out["w2_mag"], out["w2_err"])
    return out


def leave_one_out(star_visits: list[dict[str, np.ndarray]], **kwargs) -> list[dict]:
    """Correct every star with the offsets measured from all the *other* stars,
    so that a star never removes its own noise."""
    out = []
    for i, v in enumerate(star_visits):
        others = star_visits[:i] + star_visits[i + 1:]
        out.append(apply_offsets(v, field_offsets(others), **kwargs))
    return out
