"""unTimely (time-domain unWISE; Meisner et al. 2023, AJ 165, 36) light curves.

Input: data/cache/untimely/<name>.csv written by scripts/fetch_untimely.py.
Fluxes are PSF-fit fluxes on the per-visit unWISE coadds in Vega nanomaggies,
mag = 22.5 - 2.5 log10(FLUX) (Meisner et al. 2023, Table 1). DFLUX is
"statistical only".
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from .config import BANDS, CACHE_DIR, F0_JY

MATCH_RADIUS_ARCSEC = 5.0
NANOMAGGY_ZP = 22.5
BAND_CODE = {"w1": 1, "w2": 2}


def cache_path(name: str) -> Path:
    return CACHE_DIR / "untimely" / f"{name}.csv"


def nanomaggy_to_jy(flux_nmgy, band: str):
    """Vega nanomaggies -> Jy with the same zero point as the pipeline."""
    return F0_JY[band] * 10.0 ** (-0.4 * NANOMAGGY_ZP) * np.asarray(flux_nmgy, dtype=float)


def read_rows(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        for k in ("sep_arcsec", "MJDMEAN", "MJDMIN", "MJDMAX", "flux", "dflux", "fracflux",
                  "rchi2", "qf"):
            r[k] = float(r[k])
        for k in ("band", "EPOCH", "FORWARD", "primary", "flags_unwise", "flags_info", "nm"):
            r[k] = int(float(r[k]))
    return rows


def read_field(path: Path) -> dict[int, list[dict]]:
    """Stage 5a field file (scripts/fetch_untimely.py --sample): rows grouped by
    position id (0 = target, k = k-th comparison star)."""
    out: dict[int, list[dict]] = {}
    for r in read_rows(path):
        out.setdefault(int(r["pid"]), []).append(r)
    return out


def light_curve(rows: list[dict], radius_arcsec: float = MATCH_RADIUS_ARCSEC) -> dict:
    """One point per unTimely epoch: the primary detection nearest to the
    target within ``radius_arcsec``. ``w?_n_det`` counts all detections inside
    the radius (more than one would mean the galaxy was split)."""
    epochs = sorted({r["EPOCH"] for r in rows})
    out: dict[str, list] = {"epoch": [], "mjd": []}
    for b in BANDS:
        for k in ("mag", "err", "flux_jy", "flux_err_jy", "n_det", "fracflux", "flags_info"):
            out[f"{b}_{k}"] = []
    for e in epochs:
        mjds = []
        for b in BANDS:
            cand = [r for r in rows if r["EPOCH"] == e and r["band"] == BAND_CODE[b]
                    and r["sep_arcsec"] <= radius_arcsec and r["primary"] == 1 and r["flux"] > 0]
            if not cand:
                for k, val in (("mag", np.nan), ("err", np.nan), ("flux_jy", np.nan),
                               ("flux_err_jy", np.nan), ("n_det", 0), ("fracflux", np.nan),
                               ("flags_info", -1)):
                    out[f"{b}_{k}"].append(val)
                continue
            best = min(cand, key=lambda r: r["sep_arcsec"])
            mjds.append(best["MJDMEAN"])
            out[f"{b}_mag"].append(NANOMAGGY_ZP - 2.5 * np.log10(best["flux"]))
            out[f"{b}_err"].append(2.5 / np.log(10.0) * best["dflux"] / best["flux"])
            out[f"{b}_flux_jy"].append(float(nanomaggy_to_jy(best["flux"], b)))
            out[f"{b}_flux_err_jy"].append(float(nanomaggy_to_jy(best["dflux"], b)))
            out[f"{b}_n_det"].append(len(cand))
            out[f"{b}_fracflux"].append(best["fracflux"])
            out[f"{b}_flags_info"].append(best["flags_info"])
        out["epoch"].append(e)
        out["mjd"].append(float(np.mean(mjds)) if mjds else np.nan)
    return {k: np.array(v) for k, v in out.items()}


def load(name: str) -> dict:
    return light_curve(read_rows(cache_path(name)))


def match_visits(mjd_a, mjd_b, max_days: float = 15.0) -> list[tuple[int, int]]:
    """Pairs (i, j) with |mjd_a[i] - mjd_b[j]| <= max_days, nearest neighbour."""
    mjd_a, mjd_b = np.asarray(mjd_a, float), np.asarray(mjd_b, float)
    pairs = []
    for i, t in enumerate(mjd_a):
        if not np.isfinite(t) or mjd_b.size == 0:
            continue
        j = int(np.nanargmin(np.abs(mjd_b - t)))
        if abs(mjd_b[j] - t) <= max_days:
            pairs.append((i, j))
    return pairs
