"""Run the stage-2 known-answer tests from the local cache and write results/.

    python scripts/run_analysis.py

Needs data/cache/ (scripts/fetch_data.py). No network access.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from irecho import analysis, config, plots  # noqa: E402
from irecho.lightcurve import load_target, write_epochs_csv  # noqa: E402


def to_jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_jsonable(v) for v in obj]
    if isinstance(obj, (np.floating, float)):
        return None if not np.isfinite(obj) else round(float(obj), 6)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    return obj


def write_json(name: str, data: dict) -> None:
    path = config.RESULTS_DIR / name
    path.write_text(json.dumps(to_jsonable(data), indent=1) + "\n")
    print(f"wrote {path.relative_to(config.PROJECT_ROOT)}")


def write_star_table(rows: list[dict]) -> None:
    path = config.RESULTS_DIR / "reference_stars_chi2.csv"
    fields = list(rows[0].keys())
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: (f"{v:.5g}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"wrote {path.relative_to(config.PROJECT_ROOT)}")


def main() -> None:
    config.RESULTS_DIR.mkdir(exist_ok=True)
    (config.RESULTS_DIR / "epochs").mkdir(exist_ok=True)

    # (a) reference stars
    stars = [load_target(t["name"], t["ra"], t["dec"]) for t in analysis.load_reference_stars()]
    star_result = analysis.reference_star_analysis(stars)
    common = analysis.common_mode_analysis(stars)
    write_star_table(star_result["per_star"])
    write_json("reference_stars_summary.json",
               {"summary": star_result["summary"], "common_mode": common})
    plots.plot_reference_stars(star_result, config.RESULTS_DIR / "fig_refstars_chi2.png")
    plots.plot_common_mode(common, config.RESULTS_DIR / "fig_refstars_common_mode.png")

    # (b) NGC 3786
    t = analysis.load_ngc3786()
    target = load_target(t["name"], t["ra"], t["dec"])
    write_epochs_csv(target, config.RESULTS_DIR / "epochs" / "ngc3786.csv")
    ngc = analysis.ngc3786_analysis(target)
    write_json("ngc3786_summary.json", ngc)
    plots.plot_ngc3786(target, ngc, config.RESULTS_DIR / "fig_ngc3786_lightcurve.png")

    # (c) controls
    controls = [load_target(c["name"], c["ra"], c["dec"]) for c in analysis.load_controls()]
    for lc in controls:
        write_epochs_csv(lc, config.RESULTS_DIR / "epochs" / f"{lc.name}.csv")
    ctrl = analysis.control_analysis(controls)
    write_json("controls_summary.json", ctrl)
    plots.plot_controls(target, controls, config.RESULTS_DIR / "fig_controls.png")

    s = star_result["summary"]
    print("\n(a) reference stars:", s["n_stars"], "stars,",
          "median chi2_red W1/W2 =",
          f"{s['w1']['err']['chi2_red_median']:.2f}/{s['w2']['err']['chi2_red_median']:.2f},",
          "flagged by detector:", s["n_stars_flagged_by_detector"])
    p = ngc["paper_range"]
    print("(b) NGC 3786 (paper range): onset", p.get("onset_date"),
          "| dW1 =", round(p["last_flagged_visit"]["w1_dmag"], 3) if p.get("onset_date") else None,
          "| dW2 =", round(p["last_flagged_visit"]["w2_dmag"], 3) if p.get("onset_date") else None)
    print("(c) controls:", ctrl["n_objects_flagged"], "of", ctrl["n_objects"], "flagged")


if __name__ == "__main__":
    main()
