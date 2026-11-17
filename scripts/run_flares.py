"""Stage 6A: measure the flare candidates and fill the consistency table.

    python scripts/run_flares.py

Rules: docs/FLARES_PLAN.md (committed first). Output: results/s6a_flares.json.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from irecho import config, flares, sample_io, untimely  # noqa: E402
from irecho import stage5b as s5  # noqa: E402
from irecho.config import BANDS  # noqa: E402

TARGETS = [("bat0463", "main"), ("bat0582", "main"), ("bat0074", "also flagged by (a)"),
           ("bat0686", "also flagged by (a)"), ("bat0560", "reference (NGC 3786)")]
OUT = config.RESULTS_DIR / "s6a_flares.json"


def _clean(x):
    if isinstance(x, dict):
        return {k: _clean(v) for k, v in x.items() if not k.startswith("_")}
    if isinstance(x, (list, tuple)):
        return [_clean(v) for v in x]
    if isinstance(x, (np.floating, float)):
        return None if not np.isfinite(x) else float(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x


def main() -> None:
    sample = {t["name"]: t for t in sample_io.load_sample()}
    stars = sample_io.load_stars()
    with open(config.RESULTS_DIR / "s5b_objects.csv", newline="") as fh:
        obj = {r["name"]: r for r in csv.DictReader(fh)}
    with open(config.TARGETS_DIR / "sample_v1_props.csv", newline="") as fh:
        props = {r["name"]: r for r in csv.DictReader(fh)}
    offsets: dict[str, list[dict]] = {}
    with open(config.RESULTS_DIR / "qa_offsets.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            offsets.setdefault(r["field"], []).append(r)
    out = {}
    for tag, (name, role) in enumerate(TARGETS):
        t, o, p = sample[name], obj[name], props[name]
        z = float(t["z"])
        v = s5.read_frozen(name)
        n_base = max(int(float(o["a_first_index"])), 4)
        base, base_err = flares.baseline(v, n_base)
        ex = flares.excess(v, base, base_err)
        runs = flares.episodes(ex, n_base)
        ep = flares.main_episode(ex, runs)
        start, end, peak = ep
        phys = flares.visit_physics(v, n_base, z, tag)
        en = flares.energy(v, phys["_l_draws"], start, end, z)
        ts = flares.time_scales(v, ex, n_base, ep)
        event_idx = np.arange(start, end + 1)
        base_idx = np.arange(n_base)
        field = sample_io.load_field(t, stars[name])
        base_mag = {b: float(np.median(v[f"{b}_mag"][:n_base] + v[f"{b}_offset"][:n_base]))
                    for b in BANDS}  # single exposures are not offset-corrected
        expo = flares.exposure_checks(field["target"], v, base_idx, event_idx, base_mag)
        # zero-point offsets of the field at the event visits
        off = offsets[name]
        off_mjd = np.array([float(r["mjd"]) for r in off])
        zp = {}
        for b in BANDS:
            vals = np.array([float(r[f"{b}_offset_mmag"]) if r[f"{b}_offset_mmag"] else np.nan
                             for r in off])
            std = float(np.nanstd(vals, ddof=1))
            at_event = [float(vals[int(np.argmin(np.abs(off_mjd - v["mjd"][k])))]) for k in event_idx]
            zp[b] = {"field_offset_std_mmag": std,
                     "max_abs_offset_at_event_mmag": float(np.nanmax(np.abs(at_event))),
                     "within_3_std": bool(np.nanmax(np.abs(at_event)) <= 3.0 * std)}
        rows = untimely.read_field(sample_io.field_paths(name)["untimely"]).get(0)
        if rows:
            u = untimely.light_curve(rows)
            unt = flares.untimely_check(v, u, n_base, event_idx,
                                        untimely.match_visits(v["mjd"], u["mjd"]))
        else:
            unt = {"status": "not covered", "n_event_pairs": 0}
        # colour at the peak (O7): brightening in each band and its error
        mag_err_base = {b: 2.5 / np.log(10.0) * base_err[b] / base[b] for b in BANDS}
        dm = {b: float(-ex[f"{b}_dmag"][peak]) for b in BANDS}
        dm_err = {b: float(np.hypot(v[f"{b}_err"][peak], mag_err_base[b])) for b in BANDS}
        peak_row = next(r for r in phys["visits"] if r["index"] == peak)
        l_bol = 10.0 ** float(p["dr2_loglbol"]) if p["dr2_loglbol"] else float("nan")
        chi2_base = {b: float(np.sum(((v[f"{b}_mag"][:n_base] - np.mean(v[f"{b}_mag"][:n_base]))
                                      / v[f"{b}_err"][:n_base]) ** 2) / (n_base - 1)) for b in BANDS}
        below_half = np.flatnonzero(ex["w2"][:peak] < 0.5 * ex["w2"][peak])
        half_rise = (float((v["mjd"][peak] - v["mjd"][below_half[-1]]) / flares.DAYS_PER_YEAR)
                     if below_half.size else None)
        o1 = {
            "both_bands_in_same_exposures":
                expo["fraction_event_exposures_brighter_in_both_bands"] >= flares.BOTH_BANDS_MIN_FRACTION,
            "position_within_3_arcsec": expo["position_shift_arcsec"] < flares.ARTEFACT_MAX_ARCSEC,
            "field_zero_point_normal": zp["w1"]["within_3_std"] and zp["w2"]["within_3_std"],
            "untimely_agrees": (None if unt["status"] != "covered"
                                else unt["max_abs_difference_mag"] <= flares.UNTIMELY_MAX_DIFF_MAG),
            "kept_fraction_not_halved":
                expo["kept_fraction_event"] >= 0.5 * expo["kept_fraction_baseline"],
        }
        measured = {
            "o1_checks": o1, "nu_l_nu_w2_peak": peak_row["nu_l_nu_w2_erg_s"],
            "e_dust": en["e_dust_erg"], "n_episodes": len(runs),
            "dm_w2_minus_dm_w1": dm["w2"] - dm["w1"], "dm_diff_err": float(np.hypot(dm_err["w1"], dm_err["w2"])),
            "position_shift_arcsec": expo["position_shift_arcsec"],
            "position_shift_err_arcsec": expo["position_shift_err_arcsec"],
            "rebrightening_after_peak": ts["rebrightening_after_peak"],
            "t_rise_yr": ts["t_rise_yr"], "t_decay_efold_yr": ts["t_decay_efold_yr"],
        }
        out[name] = {
            "counterpart": t["counterpart"], "role": role, "type": p["type"], "z": z,
            "log_lbol_bass": float(p["dr2_loglbol"]) if p["dr2_loglbol"] else None,
            "log_edd_bass": float(p["dr2_logedd"]) if p["dr2_logedd"] else None,
            "log_nh_bass": float(p["dr1_lognh"]) if p["dr1_lognh"] else None,
            "n_visits": int(v["mjd"].size), "n_baseline_visits": n_base,
            "baseline_chi2_red": chi2_base,
            "episodes": [{"start_mjd": float(v["mjd"][a]), "end_mjd": float(v["mjd"][b]),
                          "n_visits": int(b - a + 1),
                          "peak_w2_brightening_mag": float(np.max(-ex["w2_dmag"][a: b + 1]))}
                         for a, b in runs],
            "main_episode": {"start_mjd": float(v["mjd"][start]), "end_mjd": float(v["mjd"][end]),
                             "peak_mjd": float(v["mjd"][peak]), "peak_index": int(peak),
                             "brightening_at_peak_mag": dm, "brightening_err_mag": dm_err},
            "time_scales": ts, "physics": phys, "energy": en,
            "peak": {**{k: peak_row[k] for k in peak_row if k != "index"},
                     "l_dust_over_lbol": peak_row["l_dust_erg_s"] / l_bol if np.isfinite(l_bol) else None},
            "exposure_checks": expo, "field_zero_point": zp, "untimely": unt,
            "ztf_from_stage5b": {k: o[k] for k in ("ztf_g_status", "ztf_g_dm_opt", "ztf_g_dm_opt_err",
                                                  "ztf_r_status", "ztf_r_dm_opt", "ztf_r_dm_opt_err",
                                                  "ztf_verdict")},
            "descriptive_not_in_plan": {
                "half_rise_yr": half_rise,
                "note": "time from the last visit below half of the peak W2 excess to the peak"},
            "measured_for_table": measured, "table": flares.classify(measured),
        }
        tb = out[name]["table"]
        print(name, t["counterpart"], "| episodes", len(runs), "| peak dm", {b: round(dm[b], 2) for b in BANDS},
              "| " + " ".join(f"{k}:{tb[k]['verdict']}" for k in tb))
    OUT.write_text(json.dumps(_clean(out), indent=1, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
