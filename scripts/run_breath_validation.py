"""Stage 5b, step 1: test the slope and hysteresis estimators on toy models
before any frozen light curve is opened.

    python scripts/run_breath_validation.py            # all cells, resumable
    python scripts/run_breath_validation.py --aggregate-only

Synthetic objects take the visit times of a random sample_v1 field
(results/qa_offsets.csv) and the visit errors of a random main-sample object
(results/qa_targets.json: median statistical and offset errors;
results/qa_error_floor.json: error floor). Those are QA values of stage 5a; no
target magnitude is read.

One checkpoint per cell in data/cache/s5/breath_validation/. Output:
results/s5b_validation.json.
"""

from __future__ import annotations

import csv
import json
import sys
import time
import zlib
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.stats import wilcoxon

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from irecho import breath, config, sample_io  # noqa: E402
from irecho.breath import ToyModel  # noqa: E402

CHECKPOINT_DIR = sample_io.S5_DIR / "breath_validation"
OUT_PATH = config.RESULTS_DIR / "s5b_validation.json"
N_PER_CELL = 400
SEED_TAG = 52

# driver classes: (rms of the damped random walk in mag, one-sided ramp in mag over 3 yr)
DRIVERS = {"drw0.15": (0.15, 0.0), "drw0.3": (0.3, 0.0), "drw0.6": (0.6, 0.0),
           "ramp1.0": (0.2, 1.0), "ramp2.0": (0.2, 2.0)}
STRONG = ("drw0.3", "drw0.6", "ramp1.0", "ramp2.0")
TAUS = (0.25, 0.5, 1.0, 2.0, 5.0, 10.0, float("inf"))


def cells() -> list[dict]:
    out = []
    for host in (0.0, 0.5, 0.8):
        for d in DRIVERS:
            out.append({"group": "slope", "model": "boundary", "host": host, "driver": d,
                        "toy": ToyModel("boundary", 0.0, host_w1=host)})
            out.append({"group": "slope", "model": "in place", "host": host, "driver": d,
                        "toy": ToyModel("fixed", host_w1=host)})
            out.append({"group": "slope", "model": "mix 50%", "host": host, "driver": d,
                        "toy": ToyModel("boundary", 0.0, f_fixed=0.5, host_w1=host)})
    for lag in (10.0, 30.0, 90.0):
        for d in STRONG:
            out.append({"group": "band lag", "model": f"boundary, W2 lag {lag:.0f} d",
                        "host": 0.0, "driver": d, "lag": lag,
                        "toy": ToyModel("boundary", 0.0, lag_w2_extra_days=lag)})
    for d in STRONG:
        out.append({"group": "correlated noise", "model": "boundary, noise corr 0.5",
                    "host": 0.0, "driver": d, "noise_correlation": 0.5,
                    "toy": ToyModel("boundary", 0.0)})
    for tau in TAUS:
        for host in (0.0, 0.5):
            for d in STRONG:
                out.append({"group": "reformation", "model": f"reformation {tau:g} yr",
                            "tau_re_yr": tau, "host": host, "driver": d,
                            "toy": ToyModel("boundary", tau, host_w1=host)})
    for c in out:
        c["key"] = f"{c['group']}|{c['model']}|host{c['host']:.1f}|{c['driver']}"
    return out


def load_real_conditions() -> tuple[list[np.ndarray], list[dict]]:
    """Visit times per field and visit error (mag) per main-sample object."""
    cadence: dict[str, list[float]] = {}
    with open(config.RESULTS_DIR / "qa_offsets.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            if r["w1_offset_mmag"] and r["w2_offset_mmag"]:
                cadence.setdefault(r["field"], []).append(float(r["mjd"]))
    qa = json.loads((config.RESULTS_DIR / "qa_targets.json").read_text())
    floors = json.loads((config.RESULTS_DIR / "qa_error_floor.json").read_text())["floors_mmag"]
    errors = []
    for t in sample_io.load_sample():
        if t["role"] != "science" or t["seen_before"]:
            continue
        q = qa[t["name"]]
        errors.append({b: 1e-3 * float(np.sqrt(q[b]["median_stat_err_mmag"] ** 2
                                               + q[b]["median_offset_err_mmag"] ** 2
                                               + floors[t["name"]][b] ** 2))
                       for b in ("w1", "w2")})
    return [np.array(sorted(v)) for v in cadence.values()], errors


KEEP = ("status", "n_used", "beta", "sigma_beta", "beta_in_place_min", "r", "measurable",
        "class", "d_hyst", "sigma_d", "hyst_measurable", "max_excess_ratio", "t_pivot_k")


def run_cell(cell: dict, cadences, errors) -> list[dict]:
    rng = np.random.default_rng([config.SEED, SEED_TAG, zlib.crc32(cell["key"].encode())])
    sigma, ramp = DRIVERS[cell["driver"]]
    rows = []
    for _ in range(N_PER_CELL):
        mjd = cadences[rng.integers(len(cadences))]
        err = errors[rng.integers(len(errors))]
        v = breath.synthetic_visits(rng, cell["toy"], mjd, err, sigma, ramp_mag=ramp,
                                    noise_correlation=cell.get("noise_correlation", 0.0))
        m = breath.measure(v["w1_flux"], v["w1_flux_err"], v["w2_flux"], v["w2_flux_err"])
        rows.append({k: m[k] for k in KEEP})
    return rows


def _pct(x, q):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return [round(float(v), 4) for v in np.percentile(x, q)] if x.size else None


def summarise(cell: dict, rows: list[dict]) -> dict:
    ok = [r for r in rows if r["status"] == "ok"]
    meas = [r for r in ok if r["measurable"]]
    hm = [r for r in ok if r["hyst_measurable"]]
    beta = np.array([r["beta"] for r in ok])
    sig = np.array([r["sigma_beta"] for r in ok])
    w = 1.0 / sig**2 if ok else np.zeros(0)
    bmin = np.array([r["beta_in_place_min"] for r in ok])
    fin = np.isfinite(bmin)
    out = {k: cell[k] for k in ("group", "model", "host", "driver")}
    out.update({
        "n": len(rows), "n_ok": len(ok), "n_measurable": len(meas),
        "beta_p16_p50_p84": _pct(beta, [16, 50, 84]),
        "sigma_beta_median": float(np.median(sig)) if ok else None,
        "pull_std_about_zero": float(np.std(beta / sig)) if ok else None,
        "stack_beta_weighted_mean": float(np.sum(w * beta) / w.sum()) if ok else None,
        "stack_beta_weighted_mean_err": float(w.sum() ** -0.5) if ok else None,
        "stack_beta_in_place_min_weighted_mean":
            float(np.sum(w[fin] * bmin[fin]) / w[fin].sum()) if fin.any() else None,
        "measurable_class_counts": dict(Counter(r["class"] for r in meas)),
        "measurable_r_p16_p50_p84": _pct([r["r"] for r in meas], [16, 50, 84]),
        "measurable_beta_p16_p50_p84": _pct([r["beta"] for r in meas], [16, 50, 84]),
        "measurable_max_excess_ratio_median":
            float(np.median([r["max_excess_ratio"] for r in meas])) if meas else None,
        "hyst_n_measurable": len(hm),
        "hyst_d_p16_p50_p84": _pct([r["d_hyst"] for r in hm], [16, 50, 84]),
        "hyst_sigma_d_median": float(np.median([r["sigma_d"] for r in hm])) if hm else None,
        "hyst_pull_std": float(np.std([r["d_hyst"] / r["sigma_d"] for r in hm])) if hm else None,
        "hyst_n_above_2sigma": int(sum(r["d_hyst"] > 2 * r["sigma_d"] for r in hm)),
        "hyst_n_below_minus_2sigma": int(sum(r["d_hyst"] < -2 * r["sigma_d"] for r in hm)),
    })
    for key in ("tau_re_yr", "lag", "noise_correlation"):
        if key in cell:
            out[key] = cell[key]
    return out


def power(d_values, n_objects: int, rng, n_trials: int = 2000) -> float:
    """Fraction of trials in which a one-sided Wilcoxon signed-rank test on
    ``n_objects`` values drawn from ``d_values`` gives p < 0.05 for D > 0."""
    d = np.asarray(d_values, dtype=float)
    if d.size < n_objects:
        return float("nan")
    hits = 0
    for _ in range(n_trials):
        draw = rng.choice(d, size=n_objects, replace=True)
        hits += wilcoxon(draw, alternative="greater").pvalue < 0.05
    return hits / n_trials


def aggregate() -> None:
    summaries, pooled = [], {}
    for cell in cells():
        path = CHECKPOINT_DIR / (cell["key"].replace("|", "__").replace(" ", "_") + ".json")
        if not path.exists():
            continue
        rows = json.loads(path.read_text())
        summaries.append(summarise(cell, rows))
        tag = None
        if cell["group"] == "reformation":
            tag = f"tau {cell['tau_re_yr']:g}"
        elif cell["group"] == "slope" and cell["model"] == "boundary" and cell["host"] < 0.6 \
                and cell["driver"] in STRONG:
            tag = "tau 0"
        elif cell["group"] == "band lag":
            tag = f"lag {cell['lag']:.0f}"
        if tag:
            pooled.setdefault((tag, cell["driver"]), []).extend(
                r["d_hyst"] for r in rows if r["status"] == "ok" and r["hyst_measurable"])
    rng = np.random.default_rng([config.SEED, SEED_TAG, 1])
    hyst = []
    tags = sorted({t for t, _ in pooled}, key=lambda s: (s.split()[0], float(s.split()[1])))
    for tag in tags:
        for driver in list(STRONG) + ["all"]:
            d = (np.concatenate([np.array(pooled.get((tag, x), [])) for x in STRONG])
                 if driver == "all" else np.array(pooled.get((tag, driver), [])))
            if not d.size:
                continue
            hyst.append({"model": tag, "driver": driver, "n_measurable": int(d.size),
                         "d_p16_p50_p84": _pct(d, [16, 50, 84]),
                         "power_n8": power(d, 8, rng), "power_n15": power(d, 15, rng),
                         "power_n30": power(d, 30, rng)})
    inf_all = next((h for h in hyst if h["model"] == "tau inf" and h["driver"] == "all"), None)
    out = {
        "what": "toy-model test of breath.measure; no target light curve was read",
        "settings": {k: getattr(breath.SETTINGS, k) for k in breath.SETTINGS.__dataclass_fields__},
        "n_per_cell": N_PER_CELL, "drivers": DRIVERS,
        "cells": summaries, "hysteresis": hyst,
        "d_min": round(0.5 * inf_all["d_p16_p50_p84"][1], 4) if inf_all else None,
    }
    OUT_PATH.write_text(json.dumps(out, indent=1) + "\n")
    print(f"{len(summaries)} cells -> {OUT_PATH.name}; d_min = {out['d_min']}")


def main() -> None:
    if "--aggregate-only" not in sys.argv:
        cadences, errors = load_real_conditions()
        CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)
        todo = cells()
        for n, cell in enumerate(todo, 1):
            path = CHECKPOINT_DIR / (cell["key"].replace("|", "__").replace(" ", "_") + ".json")
            if path.exists():
                continue
            start = time.time()
            rows = run_cell(cell, cadences, errors)
            path.write_text(json.dumps(rows) + "\n")
            print(f"[{n}/{len(todo)}] {cell['key']}: {time.time() - start:.0f}s", flush=True)
    aggregate()


if __name__ == "__main__":
    main()
