"""Stage 5a: blind quality checks of the sample_v1 fields.

What this module may compute for a *target*: counts (exposures, visits,
rejected rows), saturation, time coverage, error sizes, and the scatter of the
difference between two photometries of the same visit. What it must not
compute or return for a target: magnitudes as a function of time, amplitudes,
colours, trends, detector outputs. The comparison stars are not blind: their
light curves are used freely (chi-square, offsets, injection tests).
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from . import differential, errorfloor, untimely
from .analysis import mjd_to_date, visit_arrays
from .config import BANDS, SATURATION_WARN_MAG
from .lightcurve import LightCurve
from .stats import constant_chi2

ERR_KEY = "err_scatter"  # adopted visit error, as in stage 3-4 (stage34.ERR_KEY)
FROZEN_COLUMNS = (
    "mjd", "w1_n", "w1_mag", "w1_err", "w1_err_stat", "w1_offset", "w1_offset_err", "w1_floor",
    "w1_flux_jy", "w1_flux_err_jy",
    "w2_n", "w2_mag", "w2_err", "w2_err_stat", "w2_offset", "w2_offset_err", "w2_floor",
    "w2_flux_jy", "w2_flux_err_jy", "color", "color_err",
)


# --- comparison stars ---------------------------------------------------------
def star_tables(stars: list[LightCurve]) -> list[dict]:
    return [visit_arrays(lc, "neowise") for lc in stars]


def star_qa(names: list[str], w1_cat: list[float], sv: list[dict], loo: list[dict]) -> list[dict]:
    """Constant-source chi-square of every comparison star before and after the
    leave-one-out differential correction."""
    rows = []
    for name, w1, before, after in zip(names, w1_cat, sv, loo):
        row = {"star": name, "w1_allwise": w1, "n_visits": int(before["mjd"].size),
               "n_visits_corrected": int(after["mjd"].size)}
        for b in BANDS:
            c0 = constant_chi2(before[f"{b}_mag"], before[f"{b}_{ERR_KEY}"])
            c1 = constant_chi2(after[f"{b}_mag"], after[f"{b}_err"])
            row.update({
                f"{b}_chi2_red_before": c0.chi2_red, f"{b}_chi2_red_after": c1.chi2_red,
                f"{b}_p_before": c0.p_value, f"{b}_p_after": c1.p_value,
                f"{b}_rms_before_mmag": 1e3 * c0.rms, f"{b}_rms_after_mmag": 1e3 * c1.rms,
                f"{b}_median_err_mmag": 1e3 * c0.median_err,
            })
        rows.append(row)
    return rows


def offset_rows(field_name: str, offsets: differential.FieldOffsets) -> list[dict]:
    """Per-visit zero-point offsets of one field (mmag; positive = fainter)."""
    rows = []
    for k in range(offsets.mjd.size):
        row = {"field": field_name, "visit_date": mjd_to_date(offsets.mjd[k]),
               "mjd": round(float(offsets.mjd[k]), 3)}
        for b in BANDS:
            row[f"{b}_offset_mmag"] = _round(1e3 * offsets.offset[b][k], 2)
            row[f"{b}_offset_err_mmag"] = _round(1e3 * offsets.offset_err[b][k], 2)
            row[f"{b}_n_stars"] = int(offsets.n_stars[b][k])
        rows.append(row)
    return rows


def _round(x: float, digits: int):
    return round(float(x), digits) if np.isfinite(x) else None


def field_summary(star_rows: list[dict], offsets: differential.FieldOffsets) -> dict:
    out: dict = {"n_stars": len(star_rows), "n_field_visits": int(offsets.mjd.size)}
    for b in BANDS:
        before = np.array([r[f"{b}_chi2_red_before"] for r in star_rows], dtype=float)
        after = np.array([r[f"{b}_chi2_red_after"] for r in star_rows], dtype=float)
        off = offsets.offset[b]
        usable = np.isfinite(off)
        out[b] = {
            "chi2_red_median_before": _med(before), "chi2_red_median_after": _med(after),
            "n_p_below_0.01_before": int(sum(r[f"{b}_p_before"] < 0.01 for r in star_rows)),
            "n_p_below_0.01_after": int(sum(r[f"{b}_p_after"] < 0.01 for r in star_rows)),
            "n_visits_with_offset": int(usable.sum()),
            "offset_std_mmag": float(1e3 * np.std(off[usable], ddof=1)) if usable.sum() > 1
            else None,
            "offset_max_abs_mmag": float(1e3 * np.max(np.abs(off[usable]))) if usable.any()
            else None,
            "offset_err_median_mmag": float(1e3 * np.median(offsets.offset_err[b][usable]))
            if usable.any() else None,
            "min_stars_per_visit": int(offsets.n_stars[b].min()) if offsets.mjd.size else 0,
        }
    return out


def _med(x) -> float | None:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return float(np.median(x)) if x.size else None


# --- targets (counts only) ------------------------------------------------------
def target_counts(lc: LightCurve, corrected: dict[str, np.ndarray]) -> dict:
    """Blind-safe bookkeeping for one target. ``corrected`` is the visit table
    after the differential correction; only its size, dates and errors are used."""
    e = lc.epochs
    neo, allw = lc.select("neowise"), lc.select("allwise")
    both = np.isfinite(e["w1_mag"]) & np.isfinite(e["w2_mag"])
    out: dict = {
        "n_frames": {s: int((lc.obs["survey"] == s).sum()) for s in ("neowise", "allwise")},
        "cut_report": lc.cut_report,
        "n_visits_neowise": int(neo.sum()),
        "n_visits_neowise_w1": int((neo & np.isfinite(e["w1_mag"])).sum()),
        "n_visits_neowise_w2": int((neo & np.isfinite(e["w2_mag"])).sum()),
        "n_visits_neowise_both": int((neo & both).sum()),
        "n_visits_allwise_both": int((allw & both).sum()),
        "n_visits_corrected": int(corrected["mjd"].size),
    }
    mjd = corrected["mjd"]
    out["first_visit_date"] = mjd_to_date(mjd.min()) if mjd.size else None
    out["last_visit_date"] = mjd_to_date(mjd.max()) if mjd.size else None
    out["max_gap_days"] = float(np.diff(np.sort(mjd)).max()) if mjd.size > 1 else None
    durations = [float(np.ptp(lc.obs["mjd"][lc.epoch_id == k])) for k in range(lc.n_epochs)]
    out["max_visit_duration_days"] = max(durations) if durations else None
    for b in BANDS:
        report = lc.cut_report.get("neowise", {}).get(b, {})
        total = report.get("total_rows", 0)
        out[b] = {
            "median_frames_per_visit": _med(e[f"{b}_n"][neo & both]),
            "saturated_frame_fraction": (report.get("saturation", 0) / total) if total else None,
            # A count, not a light curve: how many visits sit above the saturation margin.
            "n_visits_brighter_than_saturation_margin": int(
                (corrected[f"{b}_mag"] < SATURATION_WARN_MAG[b]).sum()),
            "median_stat_err_mmag": _scaled_med(corrected.get(f"{b}_{ERR_KEY}"), 1e3),
            "median_offset_err_mmag": _scaled_med(corrected.get(f"{b}_offset_err"), 1e3),
            "single_exposure_rchi2_median": _med(lc.obs[f"{b}rchi2"][lc.masks[b]]),
        }
    return out


def _scaled_med(x, scale: float) -> float | None:
    m = _med(x) if x is not None else None
    return None if m is None else scale * m


# --- pipeline vs unTimely ---------------------------------------------------------
def untimely_pairs(v: dict[str, np.ndarray], rows: list[dict] | None) -> dict:
    """Scatter of (pipeline - unTimely) over the visits measured by both.

    ``v`` is a differentially corrected visit table. The returned numbers are
    difference statistics; no magnitude of either photometry is returned.
    """
    if not rows:
        return {"status": "no unTimely detection"}
    u = untimely.light_curve(rows)
    pairs = untimely.match_visits(v["mjd"], u["mjd"])
    i = np.array([p[0] for p in pairs], dtype=int)
    j = np.array([p[1] for p in pairs], dtype=int)
    out: dict = {"status": "ok", "n_untimely_epochs": int(u["mjd"].size),
                 "untimely_last_date": mjd_to_date(np.nanmax(u["mjd"])) if u["mjd"].size else None}
    for b in BANDS:
        ps = errorfloor.pair_scatter(v[f"{b}_mag"][i], u[f"{b}_mag"][j])
        ok = np.isfinite(v[f"{b}_mag"][i]) & np.isfinite(u[f"{b}_mag"][j])
        out[b] = {
            "n_pairs": ps.n_pairs, "usable": ps.usable,
            "pair_scatter_mmag": _round(1e3 * ps.scatter_mag, 2),
            "median_pipeline_minus_untimely_mag": _round(ps.median_diff_mag, 4),
            "median_pipeline_err_mmag": _scaled_med(v[f"{b}_err"][i][ok], 1e3),
            "median_untimely_err_mmag": _scaled_med(u[f"{b}_err"][j][ok], 1e3),
            "max_detections_within_5arcsec": int(u[f"{b}_n_det"][j].max()) if j.size else 0,
            "median_fracflux": _med(u[f"{b}_fracflux"][j]),
        }
    return out


# --- frozen light curves ----------------------------------------------------------
def write_frozen(path: Path, v: dict[str, np.ndarray]) -> None:
    """Final visit table of a target (differentially corrected, floor applied).
    ``w?_err`` is the total error; its three parts are kept in separate columns."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(FROZEN_COLUMNS)
        for k in range(v["mjd"].size):
            row = []
            for col in FROZEN_COLUMNS:
                key = col.replace("_err_stat", f"_{ERR_KEY}")
                val = v[key][k]
                if col == "mjd":
                    row.append(f"{val:.4f}")
                elif col.endswith("_n"):
                    row.append(int(val))
                elif "flux" in col:
                    row.append(f"{val:.6e}")
                else:
                    row.append(f"{val:.5f}")
            writer.writerow(row)
