"""Stage 6B (B4): a direct estimate of the mean W2-minus-W1 lag from sparsely
sampled light curves (docs/PREDICTIONS_6.md section 9).

For a lag delta much shorter than the variability time scale,

    F2(t) = a + b F1(t - delta) ~ a + b [F1(t) - delta dF1/dt],

so the lag is a regression coefficient on the time derivative of W1. With one
visit every six months the derivative is a central difference over a year, so
the estimate is biased when the light curve changes faster than that; the
synthetic tests decide whether the method is usable at all.
"""

from __future__ import annotations

import numpy as np


def derivative(mjd, flux) -> np.ndarray:
    """Central-difference time derivative (per day); NaN at the two ends."""
    t, f = np.asarray(mjd, dtype=float), np.asarray(flux, dtype=float)
    g = np.full(t.size, np.nan)
    g[1:-1] = (f[2:] - f[:-2]) / (t[2:] - t[:-2])
    return g


def object_terms(mjd, f1, f2, s2) -> dict | None:
    """Sufficient statistics of one object for the profile fit: the weighted
    least-squares problem F2 = a + b (F1 - delta G) for any delta."""
    mjd, f1, f2, s2 = (np.asarray(x, dtype=float) for x in (mjd, f1, f2, s2))
    g = derivative(mjd, f1)
    ok = np.isfinite(g) & np.isfinite(f1) & np.isfinite(f2) & (s2 > 0)
    if ok.sum() < 6:
        return None
    scale1, scale2 = np.median(f1[ok]), np.median(f2[ok])
    return {"f1": f1[ok] / scale1, "g": g[ok] / scale1, "f2": f2[ok] / scale2,
            "w": (scale2 / s2[ok]) ** 2, "n": int(ok.sum())}


def chi2_object(terms: dict, delta: float) -> float:
    """Minimum chi-square over (a, b) of one object at lag ``delta`` (days)."""
    x = terms["f1"] - delta * terms["g"]
    w, y = terms["w"], terms["f2"]
    sw, sx, sy = w.sum(), (w * x).sum(), (w * y).sum()
    sxx, sxy, syy = (w * x * x).sum(), (w * x * y).sum(), (w * y * y).sum()
    den = sw * sxx - sx * sx
    if den <= 0:
        return float(syy - sy * sy / sw)
    b = (sw * sxy - sx * sy) / den
    a = (sy - b * sx) / sw
    return float(syy - a * sy - b * sxy)


def stacked_lag(objects: list[dict], scale=None, grid=None) -> dict:
    """Common lag of a set of objects by profile chi-square. ``scale`` (one
    number per object) lets the lag of object i be delta * scale_i, e.g. the
    (1+z) L^(1/2) scaling of model 1 normalised to its median. The quoted error
    is the half-width of the delta-chi-square = 1 interval multiplied by
    sqrt(chi2_min / dof) (the visit errors do not include the W1 errors)."""
    grid = np.arange(-300.0, 600.5, 2.0) if grid is None else np.asarray(grid, dtype=float)
    scale = np.ones(len(objects)) if scale is None else np.asarray(scale, dtype=float)
    chi2 = np.array([sum(chi2_object(o, d * s) for o, s in zip(objects, scale)) for d in grid])
    k = int(np.argmin(chi2))
    dof = sum(o["n"] for o in objects) - 2 * len(objects) - 1
    inside = np.flatnonzero(chi2 <= chi2[k] + 1.0)
    half = 0.5 * (grid[inside[-1]] - grid[inside[0]]) if inside.size > 1 else float(grid[1] - grid[0])
    inflate = float(np.sqrt(max(chi2[k] / dof, 1.0))) if dof > 0 else float("nan")
    return {"lag_days": float(grid[k]), "lag_err_days": float(max(half, grid[1] - grid[0]) * inflate),
            "chi2_min": float(chi2[k]), "dof": int(dof), "n_objects": len(objects),
            "at_grid_edge": bool(k in (0, grid.size - 1))}
