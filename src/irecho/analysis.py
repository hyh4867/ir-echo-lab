"""Known-answer tests (stage 2): reference stars, NGC 3786, quiet-AGN controls.

Every function here works on :class:`irecho.lightcurve.LightCurve` objects and
returns plain dicts, so scripts/run_analysis.py only has to write them out.
"""

from __future__ import annotations

import csv
import json
from datetime import datetime, timedelta

import numpy as np

from . import config
from .config import BANDS, FLARE, FlareCriteria
from .epochs import MEDIAN_SE_FACTOR, assign_epochs
from .flare import FlareResult, detect_flares
from .lightcurve import LightCurve, mjd_to_year
from .photometry import blackbody_temperature
from .stats import constant_chi2, error_floor, pooled_error_floor, robust_std

# The last WISE point in Fig. 1 of Son et al. (2022) sits just before the 2021
# tick. Visits after this date were not available to the paper.
PAPER_DATA_END_MJD = 59215.0  # 2021-01-01

ERROR_MODELS = ("err_formal", "err_scatter", "err")  # 'err' = adopted (larger of the two)
INJECTION_AMPLITUDES_MAG = (0.1, 0.2, 0.3, 0.5, 0.8)


def mjd_to_date(mjd: float) -> str:
    if not np.isfinite(mjd):
        return ""
    return (datetime(1858, 11, 17) + timedelta(days=float(mjd))).strftime("%Y-%m-%d")


# --- targets ----------------------------------------------------------------
def load_ngc3786() -> dict:
    return json.loads((config.TARGETS_DIR / "ngc3786.json").read_text())


def load_reference_stars() -> list[dict]:
    with open(config.TARGETS_DIR / "reference_stars.csv") as fh:
        return [{**r, "ra": float(r["ra"]), "dec": float(r["dec"])} for r in csv.DictReader(fh)]


def load_controls() -> list[dict]:
    with open(config.TARGETS_DIR / "control_agn.csv") as fh:
        rows = [r for r in csv.DictReader(fh) if r["selected"] == "yes"]
    return [{**r, "ra": float(r["ra"]), "dec": float(r["dec"])} for r in rows]


def load_field_stars() -> list[dict]:
    """Comparison stars around the control galaxies (NEOWISE only), with the
    name of the field they belong to. Empty when the list has not been made."""
    path = config.TARGETS_DIR / "field_stars.csv"
    if not path.exists():
        return []
    with open(path, newline="") as fh:
        return [{"name": r["name"], "field": r["field"], "ra": float(r["ra"]),
                 "dec": float(r["dec"])} for r in csv.DictReader(fh)]


def all_targets() -> list[dict]:
    t = load_ngc3786()
    return [{"name": t["name"], "ra": t["ra"], "dec": t["dec"]}] + [
        {"name": r["name"], "ra": r["ra"], "dec": r["dec"]}
        for r in load_reference_stars() + load_controls()
    ]


# --- helpers ----------------------------------------------------------------
def visit_arrays(lc: LightCurve, survey: str | None = "neowise",
                 max_mjd: float | None = None) -> dict[str, np.ndarray]:
    """Visit-level arrays for visits with both bands measured."""
    e = lc.epochs
    ok = lc.select(survey) & np.isfinite(e["w1_mag"]) & np.isfinite(e["w2_mag"])
    if max_mjd is not None:
        ok &= e["mjd"] < max_mjd
    return {k: v[ok] for k, v in e.items()}


def run_detector(v: dict[str, np.ndarray], criteria: FlareCriteria = FLARE) -> FlareResult:
    return detect_flares(v["mjd"], v["w1_mag"], v["w1_err"], v["w2_mag"], v["w2_err"], criteria)


def variability_stats(v: dict[str, np.ndarray]) -> dict:
    out = {}
    for b in BANDS:
        r = constant_chi2(v[f"{b}_mag"], v[f"{b}_err"])
        out[b] = {
            "n_visits": r.n, "median_mag": float(np.median(v[f"{b}_mag"])),
            "rms_mag": r.rms, "median_err_mag": r.median_err,
            "chi2_red": r.chi2_red, "p_value": r.p_value,
            "floor_mag": error_floor(v[f"{b}_mag"], v[f"{b}_err"]),
        }
    return out


def detector_summary(res: FlareResult) -> dict:
    out = {
        "n_visits": int(res.mjd.size),
        "n_flagged": int(res.flagged.sum()),
        "flagged_dates": [mjd_to_date(m) for m in res.mjd[res.flagged]],
        "onset_mjd": res.onset_mjd if res.detected else None,
        "onset_date": mjd_to_date(res.onset_mjd) if res.detected else None,
    }
    for b in BANDS:
        s, d = res.significance[b], res.dmag[b]
        if np.isfinite(s).any():
            k = int(np.nanargmax(s))
            out[f"{b}_max_significance"] = float(s[k])
            out[f"{b}_dmag_at_max_significance"] = float(d[k])
            out[f"{b}_max_brightening_mag"] = float(np.nanmax(d))
        else:
            out[f"{b}_max_significance"] = out[f"{b}_dmag_at_max_significance"] = None
            out[f"{b}_max_brightening_mag"] = None
    return out


# --- (a) reference stars ------------------------------------------------------
def reference_star_analysis(lcs: list[LightCurve]) -> dict:
    """Is the visit-to-visit scatter of constant stars explained by the errors?"""
    per_star = []
    values = {(b, m, s): [] for b in BANDS for m in ERROR_MODELS for s in ("neowise", "all")}
    errors = {key: [] for key in values}
    n_flagged = 0
    for lc in lcs:
        row = {"name": lc.name}
        for survey in ("neowise", "all"):
            v = visit_arrays(lc, None if survey == "all" else survey)
            for b in BANDS:
                for model in ERROR_MODELS:
                    r = constant_chi2(v[f"{b}_mag"], v[f"{b}_{model}"])
                    values[(b, model, survey)].append(v[f"{b}_mag"])
                    errors[(b, model, survey)].append(v[f"{b}_{model}"])
                    if survey == "neowise":
                        row[f"{b}_chi2_red_{model}"] = r.chi2_red
                        if model == "err":
                            row[f"{b}_median_mag"] = float(np.median(v[f"{b}_mag"]))
                            row[f"{b}_n_visits"] = r.n
                            row[f"{b}_rms_mag"] = r.rms
                            row[f"{b}_median_err_mag"] = r.median_err
                            row[f"{b}_p_value"] = r.p_value
                            row[f"{b}_floor_mag"] = error_floor(v[f"{b}_mag"], v[f"{b}_err"])
                    elif model == "err":
                        row[f"{b}_chi2_red_err_all_visits"] = r.chi2_red
        res = run_detector(visit_arrays(lc, "neowise"))
        row["n_flagged"] = int(res.flagged.sum())
        row["w1_max_significance"] = res.max_significance("w1")
        row["w2_max_significance"] = res.max_significance("w2")
        n_flagged += row["n_flagged"]
        per_star.append(row)

    summary: dict = {"n_stars": len(lcs), "n_stars_flagged_by_detector": int(
        sum(r["n_flagged"] > 0 for r in per_star)), "n_visits_flagged": int(n_flagged)}
    for b in BANDS:
        band = {}
        for model in ERROR_MODELS:
            chi = np.array([r[f"{b}_chi2_red_{model}"] for r in per_star])
            band[model] = {
                "chi2_red_median": float(np.median(chi)),
                "chi2_red_mean": float(np.mean(chi)),
                "chi2_red_min": float(np.min(chi)),
                "chi2_red_max": float(np.max(chi)),
                "pooled_floor_mag_neowise": pooled_error_floor(
                    values[(b, model, "neowise")], errors[(b, model, "neowise")]),
                "pooled_floor_mag_all_visits": pooled_error_floor(
                    values[(b, model, "all")], errors[(b, model, "all")]),
            }
        p = np.array([r[f"{b}_p_value"] for r in per_star])
        band["n_stars_p_below_0.01_adopted"] = int((p < 0.01).sum())
        band["median_rms_mag"] = float(np.median([r[f"{b}_rms_mag"] for r in per_star]))
        band["median_visit_err_mag"] = float(
            np.median([r[f"{b}_median_err_mag"] for r in per_star]))
        summary[b] = band
    return {"per_star": per_star, "summary": summary}


def common_mode_analysis(lcs: list[LightCurve]) -> dict:
    """Visit-to-visit zero-point offsets shared by all stars of the field.

    All reference stars lie within 1.5 deg, so they are observed in the same
    visits. For each visit the median residual over stars estimates a common
    offset; subtracting it shows how much of the excess scatter is shared
    (calibration-like) rather than specific to each star.
    """
    per = [visit_arrays(lc, "neowise") for lc in lcs]
    all_mjd = np.concatenate([v["mjd"] for v in per])
    star_of = np.concatenate([np.full(v["mjd"].size, i) for i, v in enumerate(per)])
    visit_of = assign_epochs(all_mjd, config.BINNING.gap_days)
    n_visits = int(visit_of.max()) + 1
    out: dict = {"visit_mjd": [float(all_mjd[visit_of == k].mean()) for k in range(n_visits)]}
    out["visit_year"] = [float(mjd_to_year(m)) for m in out["visit_mjd"]]
    out["n_stars_per_visit"] = [int(np.unique(star_of[visit_of == k]).size)
                                for k in range(n_visits)]
    for b in BANDS:
        resid = np.concatenate([v[f"{b}_mag"] - np.median(v[f"{b}_mag"]) for v in per])
        err = np.concatenate([v[f"{b}_err"] for v in per])
        offset = np.array([np.median(resid[visit_of == k]) for k in range(n_visits)])
        n_k = np.array([(visit_of == k).sum() for k in range(n_visits)])
        offset_err = np.array([
            MEDIAN_SE_FACTOR * robust_std(resid[visit_of == k]) / np.sqrt(n_k[k])
            for k in range(n_visits)
        ])
        corrected = resid - offset[visit_of]
        corrected_err = np.hypot(err, offset_err[visit_of])
        sets_v = [corrected[star_of == i] for i in range(len(per))]
        sets_e = [corrected_err[star_of == i] for i in range(len(per))]
        chi = np.array([constant_chi2(v, e).chi2_red for v, e in zip(sets_v, sets_e)])
        out[b] = {
            "offset_mag": offset.tolist(),
            "offset_err_mag": offset_err.tolist(),
            "offset_std_mag": float(np.std(offset, ddof=1)),
            "offset_max_abs_mag": float(np.max(np.abs(offset))),
            "chi2_red_median_after_subtraction": float(np.median(chi)),
            "pooled_floor_mag_after_subtraction": pooled_error_floor(sets_v, sets_e),
        }
    return out


# --- (b) NGC 3786 -------------------------------------------------------------
def excess_temperature(v: dict[str, np.ndarray], baseline: np.ndarray, k: int,
                       n_draws: int = 10000, seed: int = config.SEED) -> dict:
    """Blackbody colour temperature of the flux added on top of the baseline at
    visit ``k`` (observed frame, monochromatic; see photometry.blackbody_temperature)."""
    rng = np.random.default_rng(seed)
    draws = {}
    point = {}
    for b in BANDS:
        base = v[f"{b}_flux_jy"][baseline]
        base_med = float(np.median(base))
        base_err = MEDIAN_SE_FACTOR * float(np.std(base, ddof=1)) / np.sqrt(base.size)
        point[b] = float(v[f"{b}_flux_jy"][k]) - base_med
        draws[b] = (rng.normal(v[f"{b}_flux_jy"][k], v[f"{b}_flux_err_jy"][k], n_draws)
                    - rng.normal(base_med, base_err, n_draws))
    ratio = point["w1"] / point["w2"]
    with np.errstate(invalid="ignore", divide="ignore"):
        temps = np.array([blackbody_temperature(r) for r in draws["w1"] / draws["w2"]])
    good = np.isfinite(temps)
    lo, med, hi = (np.percentile(temps[good], [16, 50, 84]) if good.any()
                   else (np.nan, np.nan, np.nan))
    return {
        "date": mjd_to_date(v["mjd"][k]),
        "excess_flux_mjy": {b: 1e3 * point[b] for b in BANDS},
        "excess_flux_ratio_w1_w2": float(ratio),
        "temperature_k": blackbody_temperature(ratio),
        "temperature_k_p16": float(lo), "temperature_k_p50": float(med),
        "temperature_k_p84": float(hi),
        "n_draws": n_draws, "seed": seed, "valid_draw_fraction": float(good.mean()),
    }


def flare_measurements(v: dict[str, np.ndarray], res: FlareResult) -> dict:
    """Onset, amplitude and colour change measured from one detector run."""
    out = detector_summary(res)
    if not res.detected:
        return out
    base = ~res.flagged & (np.arange(res.mjd.size) < res.onset_index)
    last = int(np.flatnonzero(res.flagged)[-1])
    peak_w2 = int(np.flatnonzero(res.flagged)[np.nanargmax(res.dmag["w2"][res.flagged])])
    out["n_baseline_visits"] = int(base.sum())
    out["baseline_span"] = [mjd_to_date(v["mjd"][base].min()), mjd_to_date(v["mjd"][base].max())]
    color_base = float(np.median(v["color"][base]))
    out["baseline"] = {
        "w1_mag": float(np.median(v["w1_mag"][base])),
        "w2_mag": float(np.median(v["w2_mag"][base])),
        "w1_std_mag": float(np.std(v["w1_mag"][base], ddof=1)),
        "w2_std_mag": float(np.std(v["w2_mag"][base], ddof=1)),
        "color_mag": color_base,
        "color_std_mag": float(np.std(v["color"][base], ddof=1)),
    }
    # First visit that is individually >= 3 sigma brighter in both bands
    # (weaker than the flare criterion; a hint of when the rise began).
    early = np.flatnonzero((res.significance["w1"] >= 3) & (res.significance["w2"] >= 3))
    out["first_3sigma_both_bands_date"] = mjd_to_date(res.mjd[early[0]]) if early.size else None

    def at(k: int) -> dict:
        d = {"date": mjd_to_date(v["mjd"][k]), "mjd": float(v["mjd"][k]),
             "year": float(mjd_to_year(v["mjd"][k]))}
        for b in BANDS:
            d[f"{b}_mag"] = float(v[f"{b}_mag"][k])
            d[f"{b}_err_mag"] = float(v[f"{b}_err"][k])
            d[f"{b}_dmag"] = float(res.dmag[b][k])
            d[f"{b}_significance"] = float(res.significance[b][k])
            d[f"{b}_dmag_over_baseline_std"] = float(
                res.dmag[b][k] / out["baseline"][f"{b}_std_mag"])
        d["color_mag"] = float(v["color"][k])
        d["color_err_mag"] = float(v["color_err"][k])
        d["delta_color_mag"] = float(v["color"][k] - color_base)
        d["excess_blackbody"] = excess_temperature(v, base, k)
        return d

    out["onset_visit"] = at(res.onset_index)
    out["last_flagged_visit"] = at(last)
    out["peak_w2_visit"] = at(peak_w2)
    return out


def ngc3786_analysis(lc: LightCurve) -> dict:
    full = visit_arrays(lc, "neowise")
    paper = visit_arrays(lc, "neowise", max_mjd=PAPER_DATA_END_MJD)
    res_full, res_paper = run_detector(full), run_detector(paper)
    table = []
    for k in range(res_full.mjd.size):
        table.append({
            "date": mjd_to_date(full["mjd"][k]), "year": round(float(mjd_to_year(full["mjd"][k])), 3),
            "w1_mag": float(full["w1_mag"][k]), "w1_err": float(full["w1_err"][k]),
            "w2_mag": float(full["w2_mag"][k]), "w2_err": float(full["w2_err"][k]),
            "color": float(full["color"][k]), "color_err": float(full["color_err"][k]),
            "w1_dmag": float(res_full.dmag["w1"][k]), "w2_dmag": float(res_full.dmag["w2"][k]),
            "w1_significance": float(res_full.significance["w1"][k]),
            "w2_significance": float(res_full.significance["w2"][k]),
            "flagged": bool(res_full.flagged[k]),
        })
    return {
        "paper_data_end_mjd": PAPER_DATA_END_MJD,
        "criteria": {"n_sigma": FLARE.n_sigma, "min_amp_mag": FLARE.min_amp_mag,
                     "min_baseline_epochs": FLARE.min_baseline_epochs},
        "paper_range": flare_measurements(paper, res_paper),
        "full_range": flare_measurements(full, res_full),
        "variability_neowise": variability_stats(full),
        "cut_report": lc.cut_report,
        "n_visits": {"allwise": int(lc.select("allwise").sum()),
                     "neowise": int(lc.select("neowise").sum())},
        "single_exposure_rchi2_median": {
            b: float(np.nanmedian(lc.obs[f"{b}rchi2"][lc.masks[b]])) for b in BANDS},
        "visit_table": table,
    }


# --- (c) controls -------------------------------------------------------------
def injection_recovery(v: dict[str, np.ndarray],
                       amplitudes=INJECTION_AMPLITUDES_MAG) -> dict:
    """Inject a step brightening of equal size in both bands at every admissible
    visit and count how often the detector puts the onset exactly there."""
    n = v["mjd"].size
    out = {}
    for amp in amplitudes:
        hits = trials = 0
        for k in range(FLARE.min_baseline_epochs, n):
            w1, w2 = v["w1_mag"].copy(), v["w2_mag"].copy()
            w1[k:] -= amp
            w2[k:] -= amp
            res = detect_flares(v["mjd"], w1, v["w1_err"], w2, v["w2_err"])
            trials += 1
            hits += int(res.onset_index == k)
        out[f"{amp:.1f}"] = {"recovered": hits, "trials": trials}
    return out


def control_analysis(lcs: list[LightCurve]) -> dict:
    rows = []
    for lc in lcs:
        v = visit_arrays(lc, "neowise")
        res = run_detector(v)
        row = {"name": lc.name, **detector_summary(res),
               "variability": variability_stats(v),
               "injection_recovery": injection_recovery(v),
               "single_exposure_rchi2_median": {
                   b: float(np.nanmedian(lc.obs[f"{b}rchi2"][lc.masks[b]])) for b in BANDS}}
        rows.append(row)
    pooled = {}
    for amp in INJECTION_AMPLITUDES_MAG:
        key = f"{amp:.1f}"
        pooled[key] = {
            "recovered": int(sum(r["injection_recovery"][key]["recovered"] for r in rows)),
            "trials": int(sum(r["injection_recovery"][key]["trials"] for r in rows)),
        }
    return {
        "per_object": rows,
        "n_objects": len(rows),
        "n_objects_flagged": int(sum(r["n_flagged"] > 0 for r in rows)),
        "n_visits_tested": int(sum(r["n_visits"] for r in rows)),
        "injection_recovery_pooled": pooled,
    }
