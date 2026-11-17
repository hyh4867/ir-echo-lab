"""Brightening-event ("flare candidate") detector on visit-binned light curves.

Definition (fixed before the known-answer tests, see docs/PIPELINE.md section 5):

For visit k, the baseline is every earlier visit that was not itself flagged.
With m_base the baseline median and s_base its robust scatter (never smaller
than the median visit error),

    dm_k = m_base - m_k                      (positive = brighter)
    S_k  = dm_k / sqrt(s_base^2 + err_k^2)

Visit k is flagged when, in both W1 and W2, S_k >= n_sigma and dm_k >= min_amp.
The onset is the first flagged visit.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .config import FLARE, FlareCriteria
from .stats import robust_std


@dataclass
class FlareResult:
    mjd: np.ndarray
    dmag: dict[str, np.ndarray]  # brightening relative to the running baseline
    significance: dict[str, np.ndarray]
    baseline_mag: dict[str, np.ndarray]
    baseline_scatter: dict[str, np.ndarray]
    flagged: np.ndarray
    criteria: FlareCriteria = field(default_factory=lambda: FLARE)

    @property
    def detected(self) -> bool:
        return bool(self.flagged.any())

    @property
    def onset_index(self) -> int | None:
        idx = np.flatnonzero(self.flagged)
        return int(idx[0]) if idx.size else None

    @property
    def onset_mjd(self) -> float:
        i = self.onset_index
        return float(self.mjd[i]) if i is not None else float("nan")

    def peak(self, band: str) -> tuple[float, float]:
        """(max brightening in mag among flagged visits, its MJD)."""
        if not self.detected:
            return float("nan"), float("nan")
        idx = np.flatnonzero(self.flagged)
        best = idx[np.nanargmax(self.dmag[band][idx])]
        return float(self.dmag[band][best]), float(self.mjd[best])

    def max_significance(self, band: str) -> float:
        s = self.significance[band]
        return float(np.nanmax(s)) if np.isfinite(s).any() else float("nan")


def detect_flares(mjd, w1mag, w1err, w2mag, w2err,
                  criteria: FlareCriteria = FLARE) -> FlareResult:
    """Run the detector. Inputs are visit-level arrays in time order; visits
    with NaN in either band are never flagged and never enter a baseline."""
    mjd = np.asarray(mjd, dtype=float)
    mags = {"w1": np.asarray(w1mag, dtype=float), "w2": np.asarray(w2mag, dtype=float)}
    errs = {"w1": np.asarray(w1err, dtype=float), "w2": np.asarray(w2err, dtype=float)}
    n = mjd.size
    if n > 1 and np.any(np.diff(mjd[np.isfinite(mjd)]) < 0):
        raise ValueError("visits must be in time order")

    valid = np.ones(n, dtype=bool)
    for b in mags:
        valid &= np.isfinite(mags[b]) & np.isfinite(errs[b]) & (errs[b] > 0)

    nan_arr = lambda: np.full(n, np.nan)  # noqa: E731
    dmag = {b: nan_arr() for b in mags}
    signif = {b: nan_arr() for b in mags}
    base_mag = {b: nan_arr() for b in mags}
    base_sc = {b: nan_arr() for b in mags}
    flagged = np.zeros(n, dtype=bool)

    for k in range(n):
        if not valid[k]:
            continue
        base = valid[:k] & ~flagged[:k]
        if base.sum() < criteria.min_baseline_epochs:
            continue
        passes = True
        for b in mags:
            m0 = float(np.median(mags[b][:k][base]))
            s0 = max(robust_std(mags[b][:k][base]), float(np.median(errs[b][:k][base])))
            dm = m0 - mags[b][k]
            s = dm / np.hypot(s0, errs[b][k])
            base_mag[b][k], base_sc[b][k] = m0, s0
            dmag[b][k], signif[b][k] = dm, s
            passes &= (s >= criteria.n_sigma) and (dm >= criteria.min_amp_mag)
        flagged[k] = passes

    return FlareResult(mjd=mjd, dmag=dmag, significance=signif, baseline_mag=base_mag,
                       baseline_scatter=base_sc, flagged=flagged, criteria=criteria)
