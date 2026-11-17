"""ZTF public light curves (g, r) from the IRSA light-curve service.

The catalogue photometry is PSF photometry for point sources. For the nucleus
of an extended galaxy it mixes host light and depends on seeing, so only
medians over many epochs are used (docs/PREDICTIONS.md).
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import numpy as np

from .config import CACHE_DIR
from .epochs import MEDIAN_SE_FACTOR
from .stats import robust_std

BIN_HALF_WIDTH_DAYS = 91.3  # half of the WISE cadence
MIN_POINTS = 5


def cache_path(name: str, band: str) -> Path:
    return CACHE_DIR / "ztf" / f"{name}_{band}.csv"


def read_clean(path: Path) -> dict:
    """Rows with catflags == 0 of the object ID that has the most such rows
    (different IDs are different field/CCD combinations with their own zero
    points)."""
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    good = [r for r in rows if r["catflags"] == "0"]
    counts = Counter(r["oid"] for r in good)
    main = counts.most_common(1)[0][0] if counts else None
    keep = sorted((r for r in good if r["oid"] == main), key=lambda r: float(r["mjd"]))
    return {
        "n_rows_total": len(rows), "n_rows_clean": len(good), "oid": main,
        "n_rows_used": len(keep), "oid_counts": dict(counts),
        "mjd": np.array([float(r["mjd"]) for r in keep]),
        "mag": np.array([float(r["mag"]) for r in keep]),
        "magerr": np.array([float(r["magerr"]) for r in keep]),
    }


def bin_to_visits(mjd, mag, visit_mjd, half_width: float = BIN_HALF_WIDTH_DAYS,
                  min_points: int = MIN_POINTS) -> dict:
    """Median magnitude of the optical points nearest to each WISE visit
    (within ``half_width`` days). Bins with fewer than ``min_points`` are NaN."""
    mjd, mag, visit_mjd = (np.asarray(x, dtype=float) for x in (mjd, mag, visit_mjd))
    nearest = np.array([int(np.argmin(np.abs(visit_mjd - t))) for t in mjd], dtype=int)
    inside = np.abs(visit_mjd[nearest] - mjd) <= half_width if mjd.size else np.zeros(0, bool)
    n = visit_mjd.size
    out = {"n": np.zeros(n, dtype=int), "mag": np.full(n, np.nan), "err": np.full(n, np.nan),
           "mjd": np.full(n, np.nan)}
    for k in range(n):
        sel = inside & (nearest == k)
        out["n"][k] = int(sel.sum())
        if sel.sum() >= min_points:
            out["mag"][k] = float(np.median(mag[sel]))
            out["err"][k] = float(MEDIAN_SE_FACTOR * robust_std(mag[sel]) / np.sqrt(sel.sum()))
            out["mjd"][k] = float(np.median(mjd[sel]))
    return out
