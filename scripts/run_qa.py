"""Stage 5a: blind quality checks of the sample_v1 fields, error floors, and
the freeze of the target light curves.

    python scripts/run_qa.py

Reads the cached field tables (scripts/fetch_sample.py, fetch_untimely.py
--sample). For the targets it writes counts and difference statistics only;
no target magnitude, amplitude, colour or trend is printed, plotted or stored
in results/. The final target light curves go to data/frozen/sample_v1/
(git-ignored) and only their SHA-256 is committed (data/freeze_v1.json).

Outputs
    results/qa_fields.json       per field: comparison-star chi-square, offset statistics
    results/qa_offsets.csv       per field and visit: zero-point offsets
    results/qa_stars.csv         per comparison star
    results/qa_targets.json      per target: counts, saturation, unTimely pair scatter, flags
    results/qa_error_floor.json  reference scatter of the stars and the adopted floors
    results/qa_summary.json      sample-level numbers quoted in the report
    results/fig_qa_*.png         QA figures (no target light curve)
"""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from irecho import config, differential, errorfloor, irsa, qa, sample_io, untimely  # noqa: E402
from irecho.analysis import visit_arrays  # noqa: E402
from irecho.config import BANDS  # noqa: E402

FROZEN_DIR = sample_io.FROZEN_DIR
FREEZE_PATH = sample_io.FREEZE_PATH
QA = sample_io.QA_PREFIX  # "qa_" for sample_v1, "qa2_" for sample_v2
MIN_STARS_FOR_FIELD = differential.MIN_STARS + 1  # leave-one-out needs one more

# QA flags. They mark objects for attention; nothing is removed from the sample here.
FLAG_FEW_VISITS = 15  # fewer corrected two-band NEOWISE visits than this
FLAG_SATURATED_FRACTION = 0.10  # more than this fraction of frames with saturated pixels
FLAG_STAR_CHI2 = 2.0  # field median chi2_red of the stars after correction above this
FLAG_LONG_VISIT_DAYS = 10.0


def quantiles(values, qs=(0, 5, 25, 50, 75, 95, 100)) -> dict:
    v = np.asarray([x for x in values if x is not None], dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return {"n": 0}
    return {"n": int(v.size), **{f"p{q}": float(np.percentile(v, q)) for q in qs}}


def write_csv(path: Path, rows: list[dict]) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        for r in rows:
            writer.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()})


def main() -> None:
    targets = sample_io.load_sample()
    stars = sample_io.load_stars()
    fields_out, targets_out, star_rows, offset_rows = {}, {}, [], []
    corrected: dict[str, dict] = {}
    star_pairs: list[dict] = []

    for n, t in enumerate(targets, start=1):
        name = t["name"]
        trow: dict = {"role": t["role"], "seen_before": t["seen_before"], "status": "ok"}
        targets_out[name] = trow
        if not sample_io.field_is_cached(name):
            trow["status"] = "not downloaded"
            continue
        field = sample_io.load_field(t, stars.get(name, []))
        if len(field["stars"]) < MIN_STARS_FOR_FIELD:
            trow["status"] = f"only {len(field['stars'])} comparison stars with data"
            continue
        sv = qa.star_tables(field["stars"])
        loo = differential.leave_one_out(sv, err_key=qa.ERR_KEY)
        offsets = differential.field_offsets(sv)
        rows = qa.star_qa(field["star_names"], field["star_w1"], sv, loo)
        for r in rows:
            r["field"] = name
        fields_out[name] = qa.field_summary(rows, offsets)
        offset_rows.extend(qa.offset_rows(name, offsets))

        by_pid = {}
        upath = sample_io.field_paths(name)["untimely"]
        if upath.exists():
            by_pid = untimely.read_field(upath)
        for r, pid, v in zip(rows, field["star_pids"], loo):
            pair = qa.untimely_pairs(v, by_pid.get(pid))
            for b in BANDS:
                ok = pair["status"] == "ok"
                r[f"{b}_untimely_n_pairs"] = pair[b]["n_pairs"] if ok else 0
                r[f"{b}_untimely_pair_scatter_mmag"] = (pair[b]["pair_scatter_mmag"]
                                                        if ok else None)
                r[f"{b}_untimely_usable"] = bool(ok and pair[b]["usable"])
            star_pairs.append(r)
        star_rows.extend(rows)

        if field["target"] is None:
            trow["status"] = "no single-exposure detection of the target"
            continue
        v_t = differential.apply_offsets(visit_arrays(field["target"], "neowise"), offsets,
                                         err_key=qa.ERR_KEY)
        trow.update(qa.target_counts(field["target"], v_t))
        trow["untimely"] = (qa.untimely_pairs(v_t, by_pid.get(0)) if upath.exists()
                            else {"status": "no unTimely file"})
        corrected[name] = v_t
        if n % 25 == 0:
            print(f"[{n}/{len(targets)}] fields processed", flush=True)

    # --- error floor -----------------------------------------------------------
    reference, per_object, usable = {}, {}, {}
    for b in BANDS:
        good = [r for r in star_pairs if r[f"{b}_untimely_usable"]]
        reference[b] = errorfloor.reference_scatter(
            [1e-3 * r[f"{b}_untimely_pair_scatter_mmag"] for r in good],
            [r["w1_allwise"] for r in good])
    w1_cat = {t["name"]: float(t["w1mpro"]) for t in targets}
    for name in corrected:
        u = targets_out[name]["untimely"]
        per_object[name], usable[name] = {}, {}
        for b in BANDS:
            ok = u["status"] == "ok" and u[b]["usable"]
            usable[name][b] = bool(ok)
            per_object[name][b] = errorfloor.floor_from_scatter(
                1e-3 * u[b]["pair_scatter_mmag"], errorfloor.reference_for(w1_cat[name],
                                                                             reference[b])
            ) if ok else float("nan")
    floors, fallback = errorfloor.adopt_floors(per_object, usable)
    # The adopted floor is the value rounded to 0.01 mmag, i.e. exactly the number
    # written to qa_error_floor.json: the frozen files can be rebuilt from it.
    for name in floors:
        for b in BANDS:
            floors[name][b] = 1e-3 * round(1e3 * floors[name][b], 2)
    for name in corrected:
        targets_out[name]["error_floor_mmag"] = {
            b: round(1e3 * floors[name][b], 2) for b in BANDS}
        targets_out[name]["error_floor_source"] = {
            b: floors[name][f"{b}_source"] for b in BANDS}

    # --- freeze -------------------------------------------------------------------
    frozen = {}
    for name, v_t in corrected.items():
        final = errorfloor.apply_floor(v_t, floors[name])
        path = FROZEN_DIR / f"{name}.csv"
        qa.write_frozen(path, final)
        frozen[name] = {"sha256": irsa.sha256_file(path), "n_visits": int(final["mjd"].size)}

    # --- flags ----------------------------------------------------------------------
    for name, trow in targets_out.items():
        flags = []
        if trow["status"] != "ok":
            flags.append("no_light_curve")
        else:
            f = fields_out[name]
            if trow["n_visits_corrected"] < FLAG_FEW_VISITS:
                flags.append("few_visits")
            if any((trow[b]["saturated_frame_fraction"] or 0) > FLAG_SATURATED_FRACTION
                   or trow[b]["n_visits_brighter_than_saturation_margin"] > 0 for b in BANDS):
                flags.append("saturation")
            if any((f[b]["chi2_red_median_after"] or 0) > FLAG_STAR_CHI2 for b in BANDS):
                flags.append("star_chi2_high")
            if (trow["max_visit_duration_days"] or 0) > FLAG_LONG_VISIT_DAYS:
                flags.append("long_visit")
            u = trow["untimely"]
            if u["status"] != "ok":
                flags.append("no_untimely")
            else:
                if not all(u[b]["usable"] for b in BANDS):
                    flags.append("untimely_few_pairs")
                if any(u[b]["max_detections_within_5arcsec"] > 1 for b in BANDS):
                    flags.append("untimely_split")
        trow["flags"] = flags

    # --- write ---------------------------------------------------------------------
    res = config.RESULTS_DIR
    (res / f"{QA}fields.json").write_text(json.dumps(fields_out, indent=1) + "\n")
    (res / f"{QA}targets.json").write_text(json.dumps(targets_out, indent=1) + "\n")
    write_csv(res / f"{QA}offsets.csv", offset_rows)
    star_fields = ["field", "star"] + [k for k in star_rows[0] if k not in ("field", "star")]
    write_csv(res / f"{QA}stars.csv", [{k: r[k] for k in star_fields} for r in star_rows])
    floor_doc = {
        "method": "errorfloor.py: floor^2 = max(0, s_pair(target)^2 - s_pair(stars)^2); "
                  "own value with >= %d matched visits, else the sample median" % errorfloor.MIN_PAIRS,
        "reference_scatter_of_stars": reference,
        "fallback_floor_mmag": {b: round(1e3 * fallback[b], 2) for b in BANDS},
        "floors_mmag": {name: {b: round(1e3 * floors[name][b], 2) for b in BANDS}
                        for name in floors},
        "source": {name: {b: floors[name][f"{b}_source"] for b in BANDS} for name in floors},
    }
    (res / f"{QA}error_floor.json").write_text(json.dumps(floor_doc, indent=1) + "\n")

    ok_targets = [t for t in targets_out.values() if t["status"] == "ok"]
    summary: dict = {
        "n_sample": len(targets),
        "n_with_light_curve": len(ok_targets),
        "status_counts": _count(t["status"] for t in targets_out.values()),
        "n_comparison_stars_with_data": len(star_rows),
        "flag_counts": _count(f for t in targets_out.values() for f in t["flags"]),
        "n_targets_without_flags": int(sum(1 for t in targets_out.values() if not t["flags"])),
        "flag_thresholds": {"few_visits": FLAG_FEW_VISITS,
                            "saturated_frame_fraction": FLAG_SATURATED_FRACTION,
                            "star_chi2_red_after": FLAG_STAR_CHI2,
                            "long_visit_days": FLAG_LONG_VISIT_DAYS},
        "targets": {
            "n_visits_corrected": quantiles([t["n_visits_corrected"] for t in ok_targets]),
            "n_visits_allwise_both": quantiles([t["n_visits_allwise_both"] for t in ok_targets]),
            "n_frames_neowise": quantiles([t["n_frames"]["neowise"] for t in ok_targets]),
            "max_gap_days": quantiles([t["max_gap_days"] for t in ok_targets]),
            "max_visit_duration_days": quantiles([t["max_visit_duration_days"]
                                                  for t in ok_targets]),
        },
        "stars": {}, "fields": {}, "untimely": {}, "error_floor": {},
    }
    for b in BANDS:
        kept = sum(t["cut_report"]["neowise"][b]["kept"] for t in ok_targets)
        total = sum(t["cut_report"]["neowise"][b]["total_rows"] for t in ok_targets)
        crit = {k: int(sum(t["cut_report"]["neowise"][b].get(k, 0) for t in ok_targets))
                for k in ok_targets[0]["cut_report"]["neowise"][b] if k not in ("kept", "total_rows")}
        summary["targets"][b] = {
            "neowise_frames_total": int(total), "neowise_frames_kept": int(kept),
            "neowise_kept_fraction": kept / total if total else None,
            "neowise_frames_failing_each_cut": crit,
            "saturated_frame_fraction": quantiles([t[b]["saturated_frame_fraction"]
                                                   for t in ok_targets]),
            "n_targets_saturated_fraction_above_threshold": int(sum(
                (t[b]["saturated_frame_fraction"] or 0) > FLAG_SATURATED_FRACTION
                for t in ok_targets)),
            "n_targets_with_visits_above_saturation_margin": int(sum(
                t[b]["n_visits_brighter_than_saturation_margin"] > 0 for t in ok_targets)),
            "median_stat_err_mmag": quantiles([t[b]["median_stat_err_mmag"] for t in ok_targets]),
            "median_offset_err_mmag": quantiles([t[b]["median_offset_err_mmag"]
                                                 for t in ok_targets]),
        }
        before = [r[f"{b}_chi2_red_before"] for r in star_rows]
        after = [r[f"{b}_chi2_red_after"] for r in star_rows]
        summary["stars"][b] = {
            "chi2_red_before": quantiles(before), "chi2_red_after": quantiles(after),
            "n_p_below_0.01_before": int(sum(r[f"{b}_p_before"] < 0.01 for r in star_rows)),
            "n_p_below_0.01_after": int(sum(r[f"{b}_p_after"] < 0.01 for r in star_rows)),
            "rms_before_mmag": quantiles([r[f"{b}_rms_before_mmag"] for r in star_rows]),
            "rms_after_mmag": quantiles([r[f"{b}_rms_after_mmag"] for r in star_rows]),
        }
        summary["fields"][b] = {
            "chi2_red_median_before": quantiles([f[b]["chi2_red_median_before"]
                                                 for f in fields_out.values()]),
            "chi2_red_median_after": quantiles([f[b]["chi2_red_median_after"]
                                                for f in fields_out.values()]),
            "offset_std_mmag": quantiles([f[b]["offset_std_mmag"] for f in fields_out.values()]),
            "offset_max_abs_mmag": quantiles([f[b]["offset_max_abs_mmag"]
                                              for f in fields_out.values()]),
            "offset_err_median_mmag": quantiles([f[b]["offset_err_median_mmag"]
                                                 for f in fields_out.values()]),
            "min_stars_per_visit": quantiles([f[b]["min_stars_per_visit"]
                                              for f in fields_out.values()]),
        }
        tu = [t["untimely"][b] for t in ok_targets if t["untimely"]["status"] == "ok"]
        summary["untimely"][b] = {
            "n_targets_with_detection": len(tu),
            "n_targets_usable": int(sum(u["usable"] for u in tu)),
            "target_n_pairs": quantiles([u["n_pairs"] for u in tu]),
            "target_pair_scatter_mmag": quantiles([u["pair_scatter_mmag"] for u in tu
                                                   if u["usable"]]),
            "target_median_pipeline_minus_untimely_mag": quantiles(
                [u["median_pipeline_minus_untimely_mag"] for u in tu if u["usable"]]),
            "n_stars_usable": int(sum(r[f"{b}_untimely_usable"] for r in star_pairs)),
            "star_pair_scatter_mmag": quantiles([r[f"{b}_untimely_pair_scatter_mmag"]
                                                 for r in star_pairs
                                                 if r[f"{b}_untimely_usable"]]),
        }
        own = [1e3 * floors[n][b] for n in floors if floors[n][f"{b}_source"] == "own"]
        summary["error_floor"][b] = {
            "reference_scatter_mmag": [round(1e3 * x, 2) for x in reference[b]["scatter_mag"]],
            "reference_bin_edges_w1": reference[b]["edges"],
            "reference_n_stars": reference[b]["n_stars"],
            "n_own": len(own), "n_sample_median": len(floors) - len(own),
            "own_floor_mmag": quantiles(own),
            "n_own_floor_zero": int(sum(x == 0 for x in own)),
            "fallback_floor_mmag": round(1e3 * fallback[b], 2),
        }
    (res / f"{QA}summary.json").write_text(json.dumps(summary, indent=1) + "\n")

    commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True,
                            cwd=config.PROJECT_ROOT).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "--", "src", "scripts"],
                           capture_output=True, text=True, cwd=config.PROJECT_ROOT).stdout.strip()
    FREEZE_PATH.write_text(json.dumps({
        "what": "SHA-256 of the frozen target light curves (data/frozen/sample_" + sample_io.SAMPLE_TAG + "/<name>.csv, "
                "not in git). Columns: " + ", ".join(qa.FROZEN_COLUMNS),
        "frozen_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "code_commit": commit, "code_tree_clean": not dirty,
        "n_light_curves": len(frozen), "light_curves": frozen,
    }, indent=1, sort_keys=True) + "\n")
    print(json.dumps({k: summary[k] for k in ("n_sample", "n_with_light_curve", "status_counts",
                                               "n_comparison_stars_with_data", "flag_counts",
                                               "n_targets_without_flags")}, indent=1))


def _count(items) -> dict:
    out: dict = {}
    for x in items:
        out[x] = out.get(x, 0) + 1
    return dict(sorted(out.items()))


if __name__ == "__main__":
    main()
