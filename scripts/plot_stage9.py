"""Stage 9 figures.

    python scripts/plot_stage9.py ident        fig_s9_identifiability.png (synthetic samples only)
    python scripts/plot_stage9.py              all three figures (needs the fit, the verdict and the frozen light curves)

Example objects (sample_v2, fit usable and branch offset comparable), chosen by a fixed rule:
the largest improvement from T0 to T1 among the 'well determined' lags, the 'well determined'
lag closest to the sample median, and the largest residual offset in units of its error.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import run_stage9 as r9  # noqa: E402
import stage9_verdict as sv  # noqa: E402
from irecho import config, plots9  # noqa: E402
from irecho import transfer9 as t9  # noqa: E402

RES = config.RESULTS_DIR


def main() -> None:
    ident = json.loads((RES / "s9_identifiability.json").read_text())
    only_ident = len(sys.argv) > 1 and sys.argv[1] == "ident"
    observed = None
    if not only_ident:
        tr = json.loads((RES / "s9_transfer.json").read_text())
        observed = {tag: tr["samples"][tag]["statistics"]["median_d_colour_obs"] for tag in r9.TAGS}
    made = [plots9.plot_identifiability(ident, observed, RES / "fig_s9_identifiability.png")]
    if not only_ident:
        verdict = json.loads((RES / "s9_verdict.json").read_text())
        rows = sv.real_rows()
        objs = {tag: r9.load(tag) for tag in r9.TAGS}
        ok = {tag: [(o, r) for o, r in zip(objs[tag], rows[tag]) if r["usable"] and r["usable_obs"] and r["comparable_obs"]
                    and all(r[f"usable_{k}"] for k in ("t0", "t1", "t2"))] for tag in r9.TAGS}
        # --- examples
        cand = ok["v2"]
        well = [(o, r) for o, r in cand if r["well_determined"]]
        med = float(np.median([r["tau_t1"] for r in rows["v2"] if r["usable"]]))
        picks = [(max(well, key=lambda x: x[1]["ds_01"]), "largest improvement from no lag to a lag"),
                 (min(well, key=lambda x: (abs(x[1]["tau_t1"] - med), -x[1]["ds_01"])), "lag closest to the sample median"),
                 (max(cand, key=lambda x: abs(x[1]["d_colour_obs"] - x[1]["d_colour_t2"]) / x[1]["sigma_colour_obs"]),
                  "largest residual offset in units of its error")]
        cp = {}
        with open(config.TARGETS_DIR / "sample_v2.csv", newline="") as fh:
            for t in csv.DictReader(fh):
                cp[t["name"]] = t["counterpart"]
        examples, record = [], []
        for (o, r), what in picks:
            fit = t9.fit_object(o["mjd"], o["f1"], o["s1"], o["f2"], o["s2"], o["r"], prep=o["prep"])
            examples.append({"label": f"{what}\n{cp[o['name']]} ({o['name']})", "mjd": o["mjd"], "f1": o["f1"], "s1": o["s1"], "f2": o["f2"],
                             "s2": o["s2"], "idx": fit["idx"], "pred_t0": fit["pred_t0"], "pred_t1": fit["pred_t1"], "pred_t2": fit["pred_t2"],
                             "tau_t1": fit["tau_t1"], "tau_lo": fit["tau_lo"], "tau_hi": fit["tau_hi"], "tau_t2": fit["tau_t2"],
                             "width_t2": fit["width_t2"], "ds_01": fit["ds_01"], "d_colour_obs": r["d_colour_obs"], "d_colour_t2": r["d_colour_t2"]})
            record.append({"name": o["name"], "counterpart": cp[o["name"]], "rule": what, "tau_t1_days": fit["tau_t1"],
                           "tau_interval_days": [fit["tau_lo"], fit["tau_hi"]], "tau_t2_days": fit["tau_t2"], "width_t2_days": fit["width_t2"],
                           "ds_01": fit["ds_01"], "d_colour_obs_mag": r["d_colour_obs"], "d_colour_t2_mag": r["d_colour_t2"]})
        made.append(plots9.plot_examples(examples, RES / "fig_s9_examples.png"))
        (RES / "s9_examples.json").write_text(json.dumps({"rule": __doc__.split("Example objects")[1].strip(), "objects": record},
                                                         indent=1, ensure_ascii=False) + "\n")
        # --- offsets and lag against luminosity
        lum, tau, wl = [], [], []
        for tag in r9.TAGS:
            for o, r in zip(objs[tag], rows[tag]):
                if r["usable"] and o["log_lbol"] is not None:
                    lum.append(o["log_lbol"])
                    tau.append(r["tau_t1"] / (1.0 + o["z"]))
                    wl.append(bool(r["well_determined"]))
        x = np.linspace(min(lum) - 0.1, max(lum) + 0.1, 60)
        lagl = {"log_lbol": lum, "tau_rest": tau, "well": wl, "curve_x": x.tolist(),
                "lyu": [float(np.diff(t9.lyu_lags_days(v))[0]) for v in x], "mandal": [sv.mandal_delay_rest(v) for v in x]}
        made.append(plots9.plot_offsets({tag: [r for _, r in ok[tag]] for tag in r9.TAGS}, verdict, lagl, RES / "fig_s9_offsets.png"))
    for m in made:
        print(m.relative_to(config.PROJECT_ROOT))


if __name__ == "__main__":
    main()
