"""From cached IRSA tables to a visit-binned two-band light curve."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import irsa, quality
from .config import BANDS, BINNING, QUALITY, EpochBinning, QualityCuts
from .epochs import assign_epochs, bin_band
from .photometry import color_w1_w2

EPOCH_COLUMNS = (
    "epoch", "survey", "mjd",
    "w1_n", "w1_mag", "w1_err", "w1_err_formal", "w1_err_scatter", "w1_scatter",
    "w1_flux_jy", "w1_flux_err_jy",
    "w2_n", "w2_mag", "w2_err", "w2_err_formal", "w2_err_scatter", "w2_scatter",
    "w2_flux_jy", "w2_flux_err_jy",
    "color", "color_err",
)


@dataclass
class LightCurve:
    name: str
    ra: float
    dec: float
    obs: dict[str, np.ndarray]  # single exposures, both surveys, time ordered
    masks: dict[str, np.ndarray]  # per-band boolean masks into ``obs``
    cut_report: dict[str, dict[str, dict[str, int]]]  # survey -> band -> counts
    epoch_id: np.ndarray  # visit index of every row of ``obs``
    epochs: dict[str, np.ndarray]  # visit-level table (see EPOCH_COLUMNS)

    @property
    def n_epochs(self) -> int:
        return int(len(self.epochs["mjd"]))

    def select(self, survey: str | None = None) -> np.ndarray:
        """Mask of visits from one survey ('allwise' / 'neowise') or all."""
        if survey is None:
            return np.ones(self.n_epochs, dtype=bool)
        return self.epochs["survey"] == survey


def concat_obs(parts: list[dict[str, np.ndarray]]) -> dict[str, np.ndarray]:
    parts = [p for p in parts if len(p["mjd"])]
    if not parts:
        raise ValueError("no single-exposure rows")
    out = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    order = np.argsort(out["mjd"], kind="stable")
    return {k: v[order] for k, v in out.items()}


def build_lightcurve(name: str, ra: float, dec: float,
                     tables: dict[str, dict[str, np.ndarray]],
                     cuts: QualityCuts = QUALITY,
                     binning: EpochBinning = BINNING) -> LightCurve:
    """Pure function: raw tables ({'neowise': ..., 'allwise': ...}) -> LightCurve."""
    parts, report = [], {}
    for survey, table in tables.items():
        if table is None or len(table["mjd"]) == 0:
            continue
        obs = quality.normalise(table, survey)
        one = quality.nearest_per_frame(obs, ra, dec)
        obs = {k: v[one] for k, v in obs.items()}
        report[survey] = {b: quality.band_mask(obs, b, cuts)[1] for b in BANDS}
        parts.append(obs)
    obs = concat_obs(parts)
    masks = {b: quality.band_mask(obs, b, cuts)[0] for b in BANDS}

    epoch_id = assign_epochs(obs["mjd"], binning.gap_days)
    n_epochs = int(epoch_id.max()) + 1
    epochs: dict[str, np.ndarray] = {"epoch": np.arange(n_epochs)}
    epochs["survey"] = np.array(
        [obs["survey"][epoch_id == k][0] for k in range(n_epochs)], dtype=str
    )
    per_band = {}
    for b in BANDS:
        m = masks[b]
        per_band[b] = bin_band(obs["mjd"][m], obs[f"{b}mag"][m], obs[f"{b}err"][m],
                               epoch_id[m], b, n_epochs, binning)
    # Visit time: mean over the bands that have data, else mean of all frames.
    band_mjd = np.vstack([per_band[b]["mjd"] for b in BANDS])
    all_mjd = np.array([obs["mjd"][epoch_id == k].mean() for k in range(n_epochs)])
    have = np.isfinite(band_mjd).any(axis=0)
    mjd = all_mjd.copy()
    mjd[have] = np.nanmean(band_mjd[:, have], axis=0)
    epochs["mjd"] = mjd
    for b in BANDS:
        s = per_band[b]
        epochs[f"{b}_n"] = s["n"]
        epochs[f"{b}_mag"] = s["mag"]
        epochs[f"{b}_err"] = s["err"]
        epochs[f"{b}_err_formal"] = s["err_formal"]
        epochs[f"{b}_err_scatter"] = s["err_scatter"]
        epochs[f"{b}_scatter"] = s["scatter"]
        epochs[f"{b}_flux_jy"] = s["flux"]
        epochs[f"{b}_flux_err_jy"] = s["flux_err"]
    epochs["color"], epochs["color_err"] = color_w1_w2(
        epochs["w1_mag"], epochs["w1_err"], epochs["w2_mag"], epochs["w2_err"]
    )
    return LightCurve(name=name, ra=ra, dec=dec, obs=obs, masks=masks,
                      cut_report=report, epoch_id=epoch_id, epochs=epochs)


def load_target(name: str, ra: float, dec: float, *, cache_only: bool = True,
                cuts: QualityCuts = QUALITY,
                binning: EpochBinning = BINNING,
                surveys: tuple[str, ...] = ("neowise", "allwise")) -> LightCurve:
    """Build a light curve from the local cache. With ``cache_only=False`` the
    tables are downloaded first when missing. ``surveys`` restricts the input
    tables (the field stars of the control galaxies have NEOWISE data only)."""
    paths = {s: p for s, p in irsa.cache_paths(name).items() if s in surveys}
    if not cache_only:
        irsa.fetch_target(name, ra, dec, surveys=surveys)
    missing = [str(p) for p in paths.values() if not p.exists()]
    if missing:
        raise FileNotFoundError(
            "IRSA cache missing: " + ", ".join(missing) + " (run scripts/fetch_data.py)"
        )
    tables = {s: irsa.read_csv_table(p) for s, p in paths.items()}
    return build_lightcurve(name, ra, dec, tables, cuts, binning)


def write_epochs_csv(lc: LightCurve, path: Path) -> None:
    """Visit-level table as CSV (derived product, small enough to commit)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(EPOCH_COLUMNS)
        for k in range(lc.n_epochs):
            row = []
            for col in EPOCH_COLUMNS:
                v = lc.epochs[col][k]
                if col in ("epoch",) or col.endswith("_n"):
                    row.append(int(v))
                elif col == "survey":
                    row.append(str(v))
                elif col == "mjd":
                    row.append(f"{v:.4f}")
                elif "flux" in col:
                    row.append("" if not np.isfinite(v) else f"{v:.6e}")
                else:
                    row.append("" if not np.isfinite(v) else f"{v:.4f}")
            writer.writerow(row)


def mjd_to_year(mjd):
    """Decimal year (Julian-year approximation; for plotting and tables)."""
    return 2000.0 + (np.asarray(mjd, dtype=float) - 51544.5) / 365.25
