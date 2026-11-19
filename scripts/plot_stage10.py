"""Stage 10 figures.

    python scripts/plot_stage10.py ident    fig_s10_identifiability.png (synthetic samples only)
    python scripts/plot_stage10.py          all four figures (needs the fit, the verdict and the frozen light curves)

Example objects, chosen by a fixed rule among the objects detected in both bands: the largest, the
median and the smallest value of min(r_w of W1, r_w of W2). When fewer than three objects are detected
the rule is applied to all fitted objects.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import run_stage10 as r10  # noqa: E402
import stage10_verdict as sv  # noqa: E402
from irecho import driver10 as d10  # noqa: E402
from irecho import plots10  # noqa: E402

RES = r10.RES


def table() -> list[dict]:
    with open(RES / "s10_objects.csv", newline="") as fh:
        return list(csv.DictReader(fh))


def num(v) -> float:
    return float(v) if v not in ("", None) else float("nan")


def main() -> None:
    ident = json.loads((RES / "s10_identifiability.json").read_text())
    if len(sys.argv) > 1 and sys.argv[1] == "ident":
        print(plots10.plot_identifiability(ident, RES / "fig_s10_identifiability.png"))
        return
    fit = json.loads((RES / "s10_fit.json").read_text())
    observed = {"frac_detected": fit["statistics"]["all"]["frac_detected"], "med_alpha_w1": fit["statistics"]["d"]["med_alpha_w1"]}
    made = [plots10.plot_identifiability(ident, RES / "fig_s10_identifiability.png", observed)]
    verdict = json.loads((RES / "s10_verdict.json").read_text())
    objs = {o["name"]: o for o in r10.load()}
    rows = table()
    examples = json.loads((r10.CKPT / "fit_examples.json").read_text())
    # --- examples
    cand = [r for r in rows if r["in_d"] == "1"] or rows
    cand.sort(key=lambda r: min(num(r["r_w_w1"]), num(r["r_w_w2"])))
    picks = [(cand[-1], "largest"), (cand[len(cand) // 2], "median"), (cand[0], "smallest")] if len(cand) >= 3 else [(c, "") for c in cand]
    shown = []
    for r, why in picks:
        o = objs[r["name"]]
        ex = examples[r["name"]]
        state = d10.gp_state(o["opt_t"], o["opt_s"], o["hyper"])
        grid_t = state["g0"] + d10.FIT_GRID_STEP_DAYS * np.arange(state["n"])
        shown.append(dict(o, idx=ex["idx"], pred_w1=ex["pred_w1"], pred_w2=ex["pred_w2"], grid_t=grid_t,
                          grid_mean=state["gain"] @ o["opt_y"], grid_sd=np.sqrt(np.clip(np.diag(state["p"]), 0, None)),
                          tau_w1=num(r["tau_w1"]), tau_w2=num(r["tau_w2"]),
                          q_w1=(num(r["tau_q16_w1"]), num(r["tau_q84_w1"])), q_w2=(num(r["tau_q16_w2"]), num(r["tau_q84_w2"])),
                          title=f"{r['name']} ({why} min r_w among the detected objects): r_w = {num(r['r_w_w1']):.2f} (W1), "
                                f"{num(r['r_w_w2']):.2f} (W2); {r['n_used']} visits used, {r['n_nights']} ZTF {o['band']} nights",
                          kernel_title=f"tau2 - tau1 = {num(r['dtau']):.0f} d (joint 16-84%: {num(r['dtau_q16']):.0f} to {num(r['dtau_q84']):.0f})"))
    made.append(plots10.plot_examples(shown, RES / "fig_s10_examples.png"))
    # --- Q1
    q_rows = [{"in_d": r["in_d"] == "1", "tau9": objs[r["name"]]["tau9"] if objs[r["name"]]["tau9"] is not None else float("nan"),
               "dtau": num(r["dtau"]), "dtau_q16": num(r["dtau_q16"]), "dtau_q84": num(r["dtau_q84"])} for r in rows]
    made.append(plots10.plot_dtau(q_rows, verdict, RES / "fig_s10_dtau_vs_stage9.png"))
    # --- Q2
    l_rows = []
    for r in rows:
        o = objs[r["name"]]
        l_rows.append({"in_d": r["in_d"] == "1", "logl": o["log_lbol"] if o["log_lbol"] is not None else float("nan"),
                       "rest_w1": num(r["tau_w1"]) / (1.0 + o["z"]), "rest_w2": num(r["tau_w2"]) / (1.0 + o["z"])})

    def rel(fn, i):
        return lambda x: np.array([fn(v)[i] for v in x])

    def non_adc(b):
        a, c = sv.MANDAL_NON_ADC[b]
        return lambda x: 10.0 ** (c + a * (x - d10.MANDAL["log_l0"]))

    relations = {"lyu": {"w1": rel(d10.lyu_lags_days, 0), "w2": rel(d10.lyu_lags_days, 1)},
                 "mandal": {"w1": rel(d10.mandal_lags_days, 0), "w2": rel(d10.mandal_lags_days, 1)},
                 "mandal_non_adc": {"w1": non_adc("w1"), "w2": non_adc("w2")},
                 "edge": lambda x: np.full(np.asarray(x).shape, 465.0 / 1.05)}
    made.append(plots10.plot_lag_luminosity(l_rows, relations, verdict, RES / "fig_s10_lag_luminosity.png"))
    print("\n".join(str(p) for p in made))


if __name__ == "__main__":
    main()
