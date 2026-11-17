"""Group single exposures into visits (~6 month cadence) and summarise them."""

from __future__ import annotations

import numpy as np

from .config import BINNING, EpochBinning
from .photometry import (flux_err_to_mag_err, flux_jy_to_mag, mag_err_to_flux_err,
                         mag_to_flux_jy)
from .stats import MAD_TO_SIGMA

# Asymptotic ratio of the standard error of the median to that of the mean for
# Gaussian noise, sqrt(pi/2). For n ~ 5-15 the true ratio is a few per cent
# smaller, so this is slightly conservative.
MEDIAN_SE_FACTOR = float(np.sqrt(np.pi / 2.0))


def assign_epochs(mjd, gap_days: float = BINNING.gap_days) -> np.ndarray:
    """Visit index for each time stamp: a gap longer than ``gap_days`` starts a
    new visit. Input order is preserved; indices increase with time."""
    mjd = np.asarray(mjd, dtype=float)
    if mjd.size == 0:
        return np.zeros(0, dtype=int)
    order = np.argsort(mjd, kind="stable")
    steps = np.diff(mjd[order]) > gap_days
    ids_sorted = np.concatenate([[0], np.cumsum(steps)])
    ids = np.empty(mjd.size, dtype=int)
    ids[order] = ids_sorted
    return ids


def summarise_visit(mjd, mag, mag_err, band: str,
                    binning: EpochBinning = BINNING) -> dict[str, float]:
    """Median flux and uncertainty of one visit in one band.

    Statistics are computed in flux space. Two uncertainties are returned:

    * ``err_formal``  : standard error of the median predicted by the catalogue
      errors, sqrt(pi/2) * sqrt(mean(sigma_i^2)) / sqrt(n)
    * ``err_scatter`` : the same with the observed in-visit standard deviation,
      sqrt(pi/2) * s / sqrt(n)

    The adopted error is the larger of the two.
    """
    nan = float("nan")
    empty = {
        "n": 0, "n_clipped": 0, "mjd": nan, "flux": nan, "flux_err": nan,
        "mag": nan, "err": nan, "err_formal": nan, "err_scatter": nan,
        "scatter": nan, "mean_cat_err": nan,
    }
    mjd = np.asarray(mjd, dtype=float)
    mag = np.asarray(mag, dtype=float)
    mag_err = np.asarray(mag_err, dtype=float)
    if mag.size < binning.min_frames:
        return {**empty, "n": int(mag.size)}

    flux = mag_to_flux_jy(mag, band)
    flux_err = mag_err_to_flux_err(flux, mag_err)

    # One-pass outlier rejection around the median (cosmic rays, residual artefacts).
    med = np.median(flux)
    width = max(MAD_TO_SIGMA * np.median(np.abs(flux - med)), np.median(flux_err))
    keep = np.abs(flux - med) <= binning.clip_sigma * width
    n_clipped = int((~keep).sum())
    mjd, flux, flux_err, mag = mjd[keep], flux[keep], flux_err[keep], mag[keep]
    n = int(flux.size)
    if n < binning.min_frames:
        return {**empty, "n": n, "n_clipped": n_clipped}

    med = float(np.median(flux))
    scatter_flux = float(np.std(flux, ddof=1))
    err_formal = MEDIAN_SE_FACTOR * float(np.sqrt(np.mean(flux_err**2))) / np.sqrt(n)
    err_scatter = MEDIAN_SE_FACTOR * scatter_flux / np.sqrt(n)
    err = max(err_formal, err_scatter)
    return {
        "n": n,
        "n_clipped": n_clipped,
        "mjd": float(np.mean(mjd)),
        "flux": med,
        "flux_err": float(err),
        "mag": float(flux_jy_to_mag(med, band)),
        "err": float(flux_err_to_mag_err(med, err)),
        "err_formal": float(flux_err_to_mag_err(med, err_formal)),
        "err_scatter": float(flux_err_to_mag_err(med, err_scatter)),
        "scatter": float(np.std(mag, ddof=1)),  # in-visit std of single exposures, mag
        "mean_cat_err": float(np.sqrt(np.mean(mag_err[keep] ** 2))),
    }


def bin_band(mjd, mag, mag_err, epoch_id, band: str, n_epochs: int,
             binning: EpochBinning = BINNING) -> dict[str, np.ndarray]:
    """Apply :func:`summarise_visit` to every visit. Arrays have length
    ``n_epochs``; visits with too few frames are NaN."""
    mjd = np.asarray(mjd, dtype=float)
    mag = np.asarray(mag, dtype=float)
    mag_err = np.asarray(mag_err, dtype=float)
    epoch_id = np.asarray(epoch_id, dtype=int)
    rows = [
        summarise_visit(mjd[epoch_id == k], mag[epoch_id == k], mag_err[epoch_id == k],
                        band, binning)
        for k in range(n_epochs)
    ]
    keys = rows[0].keys() if rows else ()
    return {key: np.array([r[key] for r in rows], dtype=float) for key in keys}
