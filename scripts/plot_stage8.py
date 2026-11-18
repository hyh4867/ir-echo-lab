"""Stage 8 figures.

    python scripts/plot_stage8.py

fig_s8_null.png (B1), fig_s8_direct.png (B2), fig_s8_examples.png (B4). The
examples are the sample_v2 objects with the largest D and with D closest to
the sample median, among the objects that are also directly comparable
(docs/VALIDATION_8.md section 3.4). Needs the frozen light curves.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import run_stage8 as r8  # noqa: E402
from irecho import config, plots8  # noqa: E402

RES = config.RESULTS_DIR


def main() -> None:
    made = [plots8.plot_null(RES / "s8_null.json", config.CACHE_DIR / "s8", RES / "fig_s8_null.png"),
            plots8.plot_direct(RES / "s8_direct_objects.csv", RES / "s8_direct.json", RES / "s8_null.json", RES / "fig_s8_direct.png")]
    with open(RES / "s8_direct_objects.csv", newline="") as fh:
        rows = [r for r in csv.DictReader(fh) if r["sample"] == "v2" and r["hyst_measurable"] == "1" and r["comparable"] == "1"]
    d = np.array([float(r["d_hyst"]) for r in rows])
    all_d = np.array([o["d"] for o in r8.sample_names("v2") if o["d"] is not None])
    pick = [rows[int(np.argmax(d))], rows[int(np.argmin(np.abs(d - np.median(all_d))))]]
    objs = {o["name"]: o for o in r8.load_sample("v2")}
    cp = {}
    with open(config.TARGETS_DIR / "sample_v2.csv", newline="") as fh:
        for t in csv.DictReader(fh):
            cp[t["name"]] = t["counterpart"]
    examples = []
    for r, what in zip(pick, ("largest D", "D closest to the median")):
        o = objs[r["name"]]
        examples.append({"name": r["name"], "label": f"{what}\n{cp[r['name']]} ({r['name']})",
                         "mjd": o["mjd"], "f1": o["f1"], "s1": o["s1"], "f2": o["f2"], "s2": o["s2"], "lo": float(r["lo"]),
                         "hi": float(r["hi"]), "d": float(r["d_hyst"]), "d_colour": float(r["d_colour"])})
    made.append(plots8.plot_examples(examples, RES / "fig_s8_examples.png"))
    (RES / "s8_examples.json").write_text(json.dumps(
        {"rule": "sample_v2, measurable D and directly comparable: largest D; D closest to the median of the 94",
         "objects": [{"name": e["name"], "label": e["label"], "d": e["d"], "d_colour_mag": e["d_colour"]} for e in examples]},
        indent=1, ensure_ascii=False) + "\n")
    for m in made:
        print(m.relative_to(config.PROJECT_ROOT))


if __name__ == "__main__":
    main()
