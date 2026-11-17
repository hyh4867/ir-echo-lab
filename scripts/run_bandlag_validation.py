"""Stage 6B (B4): can the mean W2-minus-W1 lag be estimated from six-monthly
visits? Synthetic test with the pass conditions of PREDICTIONS_6.md section 9.

    python scripts/run_bandlag_validation.py [prelim|final]

Synthetic samples as in run_d_lum_power.py (one toy light curve per sample_v2
object with a BASS L_bol; real visit times and errors), with a constant
inter-band lag of 0, 20, 40 or 80 days, or with no lag and a re-formation time
of 1, 2 or 5 yr. Output: results/s6b_bandlag_validation_<stage>.json.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import run_d_lum_power as power  # noqa: E402
from irecho import bandlag, breath, config  # noqa: E402
from irecho.breath import ToyModel  # noqa: E402

N_SAMPLES = 40
CASES = [("lag 0 d", 0.0, 0.0), ("lag 20 d", 20.0, 0.0), ("lag 40 d", 40.0, 0.0),
         ("lag 80 d", 80.0, 0.0), ("re-formation 1 yr, no lag", 0.0, 1.0),
         ("re-formation 2 yr, no lag", 0.0, 2.0), ("re-formation 5 yr, no lag", 0.0, 5.0)]


def synthetic_sample(rng, objects, cond, lag: float, tau: float) -> list[dict]:
    out = []
    for o in objects:
        if o["name"] not in cond:
            continue
        mjd, err = cond[o["name"]]
        dt1, _ = power.lyu_lags(o["log_lbol"])
        toy = ToyModel("boundary", tau, host_w1=float(rng.uniform(0.0, 0.8)),
                       lag_days=(1.0 + o["z"]) * dt1, lag_w2_extra_days=lag)
        ramp = float(rng.uniform(0.5, 1.5)) if rng.random() < 0.3 else 0.0
        v = breath.synthetic_visits(rng, toy, mjd, err, float(rng.uniform(0.15, 0.6)),
                                    ramp_mag=ramp, burn_in_days=3500.0)
        terms = bandlag.object_terms(v["mjd"], v["w1_flux"], v["w2_flux"], v["w2_flux_err"])
        if terms is not None:
            out.append(terms)
    return out


def main() -> None:
    stage = sys.argv[1] if len(sys.argv) > 1 else "prelim"
    objects = power.load_objects()
    cond = power.load_conditions(stage, [o["name"] for o in objects])
    rows = []
    for k, (label, lag, tau) in enumerate(CASES):
        rng = np.random.default_rng([config.SEED, 63, k, 0 if stage == "prelim" else 1])
        est, err = [], []
        for _ in range(N_SAMPLES):
            res = bandlag.stacked_lag(synthetic_sample(rng, objects, cond, lag, tau))
            est.append(res["lag_days"])
            err.append(res["lag_err_days"])
        est = np.array(est)
        rows.append({"case": label, "true_lag_days": lag, "tau_re_yr": tau, "n_samples": N_SAMPLES,
                     "mean_estimate_days": float(est.mean()), "std_of_estimates_days": float(est.std(ddof=1)),
                     "median_quoted_error_days": float(np.median(err)),
                     "bias_days": float(est.mean() - lag)})
        print(rows[-1], flush=True)
    lag_cases = [r for r in rows if r["tau_re_yr"] == 0.0]
    re_cases = [r for r in rows if r["tau_re_yr"] > 0.0]
    cond1 = all(abs(r["bias_days"]) <= max(10.0, 0.3 * r["true_lag_days"]) for r in lag_cases)
    cond2 = all(r["std_of_estimates_days"] <= 15.0 for r in lag_cases)
    cond3 = all(abs(r["mean_estimate_days"]) <= 10.0 for r in re_cases)
    out = {"stage": stage, "n_objects": len([o for o in objects if o["name"] in cond]),
           "cases": rows,
           "pass": {"1_bias": bool(cond1), "2_precision": bool(cond2),
                    "3_not_fooled_by_reformation": bool(cond3),
                    "all": bool(cond1 and cond2 and cond3)}}
    (config.RESULTS_DIR / f"s6b_bandlag_validation_{stage}.json").write_text(json.dumps(out, indent=1) + "\n")
    print(out["pass"])


if __name__ == "__main__":
    main()
