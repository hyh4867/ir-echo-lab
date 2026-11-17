"""Run stage 0 (corrections), stage 3 (physical quantities) and stage 4
(detectors, injection tests, ML comparison) from the local cache.

    python scripts/run_stage34.py

Needs data/cache/ complete (scripts/fetch_data.py, scripts/fetch_untimely.py).
Writes results/stage34_*.json, results/injection_recovery.csv and figures.
"""

from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from irecho import analysis, config, differential, plots34, stage34, untimely  # noqa: E402
from irecho.analysis import visit_arrays  # noqa: E402
from irecho.flare import detect_flares  # noqa: E402
from irecho.slowrise import detect_slow_rise  # noqa: E402

R = config.RESULTS_DIR


def write_visits_csv(path: Path, v: dict) -> None:
    cols = ("mjd", "w1_mag", "w1_err", "w1_offset", "w1_offset_err", "w1_flux_jy",
            "w1_flux_err_jy", "w2_mag", "w2_err", "w2_offset", "w2_offset_err", "w2_flux_jy",
            "w2_flux_err_jy", "color", "color_err")
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for i in range(v["mjd"].size):
            w.writerow([f"{v[c][i]:.6e}" if "flux" in c else f"{v[c][i]:.4f}" for c in cols])


def main() -> None:
    start = time.time()
    partial = "--partial" in sys.argv  # development only: use whatever is cached
    fields = stage34.load_fields(allow_partial=partial)
    print("fields:", {k: len(f["stars"]) for k, f in fields.items()}, flush=True)
    ngc = fields["ngc3786"]
    z = json.loads((config.TARGETS_DIR / "ngc3786.json").read_text())["ned_z"]

    # --- stage 0a
    diff = stage34.differential_validation(ngc["stars"])
    diff["control_fields"] = stage34.field_star_validation(fields)
    stage34.write_json("stage34_differential.json", diff)
    sv = stage34.star_visits(ngc["stars"])
    plots34.plot_differential(sv, differential.leave_one_out(sv, err_key=stage34.ERR_KEY),
                              R / "fig_differential_chi2.png")
    (R / "corrected").mkdir(exist_ok=True)
    corrected = {name: stage34.corrected_target(f) for name, f in fields.items() if f["stars"]}
    for name, v in corrected.items():
        write_visits_csv(R / "corrected" / f"{name}.csv", v)

    # --- stage 0b
    v_ngc = corrected["ngc3786"]
    raw = visit_arrays(ngc["target"], "neowise")
    res_raw = detect_flares(raw["mjd"], raw["w1_mag"], raw["w1_err"], raw["w2_mag"], raw["w2_err"])
    unt = {"per_object": stage34.untimely_comparison(fields),
           "ngc3786_flare": stage34.untimely_flare(ngc, res_raw),
           "ngc3786_followup": stage34.followup_comparisons(ngc)}
    stage34.write_json("stage34_untimely.json", unt)
    plots34.plot_untimely(raw, untimely.load("ngc3786"), R / "fig_untimely_ngc3786.png")

    # --- stage 3
    phys = stage34.physics_analysis(v_ngc, z)
    stage34.write_json("stage34_physics_ngc3786.json", phys)
    plots34.plot_physics(phys, R / "fig_physics_ngc3786.png")
    hyp = stage34.hypothesis_tests(phys)
    base_w1 = phys["baseline"]["w1_flux_mjy"]
    hyp["alternative_central_source"] = stage34.optical_comparison(v_ngc, base_w1)
    hyp["exploratory"] = stage34.exploratory_temperature_checks(
        phys, hyp["alternative_central_source"])
    stage34.write_json("stage34_hypotheses.json", hyp)
    plots34.plot_optical("NGC 3786", v_ngc, base_w1, hyp["alternative_central_source"],
                         R / "fig_optical_ngc3786.png")

    # --- stage 4
    det = stage34.detector_analysis(fields)
    stage34.write_json("stage34_detectors.json", det)
    for name in ("ngc2273", "ngc3786"):
        if name not in corrected:
            continue
        v = corrected[name]
        res_b = detect_slow_rise(v["mjd"], v["w1_mag"], v["w1_err"], v["w2_mag"], v["w2_err"])
        plots34.plot_slow_rise(name, v, res_b, R / f"fig_slowrise_{name}.png")
    if "ngc2273" in corrected:  # extra check, not in docs/PREDICTIONS.md
        v2 = corrected["ngc2273"]
        base2 = float(1e3 * np.median(v2["w1_flux_jy"][v2["mjd"] < 59000.0]))
        opt2 = stage34.optical_comparison(v2, base2, name="ngc2273")
        opt2["note"] = ("same rule as for NGC 3786, applied after the fact to the slow-rise "
                        "candidate; baseline = median W1 flux before MJD 59000")
        stage34.write_json("stage34_ngc2273_optical.json", opt2)
        if opt2["bands"].get("g", {}).get("status") == "ok":
            plots34.plot_optical("NGC 2273", v2, base2, opt2, R / "fig_optical_ngc2273.png")
    hosts = stage34.injection_hosts(fields)
    inj = stage34.injection_analysis(hosts)
    inj["dimming_check"] = stage34.dimming_check(hosts)
    stage34.write_json("stage34_injection.json", inj)
    with open(R / "injection_recovery.csv", "w", newline="") as fh:
        cols = ["kind", "baseline_visits", "rise_years", "amp_mag", "n", *stage34.METHODS]
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(inj["recovery"])
    plots34.plot_injection(inj, R / "fig_injection_recovery.png")

    print(f"done in {time.time() - start:.0f}s")
    m = diff["models"][stage34.ERR_KEY]
    print("chi2_red median before -> after:",
          {b: (round(m[b]["chi2_red_median_before"], 2), round(m[b]["chi2_red_median_after"], 2))
           for b in ("w1", "w2")})
    print("hypotheses:", {k: hyp[k]["verdict"] for k in hyp if k.startswith("H")},
          "| alternative:", hyp["alternative_central_source"]["verdict"])
    print("galaxies flagged by (a):", [k for k, g in det["galaxies"].items() if g["a_n_flagged"]],
          "| by (b):", [k for k, g in det["galaxies"].items() if g["b"]["flagged"]])
    print("stars flagged (a)/(b):", det["stars"]["n_flagged_a"], det["stars"]["n_flagged_b"],
          "of", det["stars"]["n"])


if __name__ == "__main__":
    main()
