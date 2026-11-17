"""Stage 6B (B0): what the two models predict for the relation between the
hysteresis statistic D and luminosity, and the power to tell them apart.

    python scripts/run_d_lum_power.py prelim   # before the download: v1 cadence and errors
    python scripts/run_d_lum_power.py final    # after the freeze: sample_v2 cadence and errors
    python scripts/run_d_lum_power.py exact    # check added after the result: visit times of the
                                               # frozen light curves themselves (models 1 and 2)

Synthetic samples: one toy light curve for every sample_v2 object that has a
BASS L_bol, with that object's redshift and luminosity.

    model 1    W2 lags W1 by (1+z)[dt_W2(L) - dt_W1(L)] (Lyu et al. 2019, Eqs. 15-16),
               also with 0.5x and 2x that lag; no re-formation delay
    model 2    no inter-band lag; re-formation time 1, 2 or 5 yr, independent of L

In every model the dust response itself is delayed and smoothed by (1+z) dt_W1(L).
No target light curve is read: 'final' uses only the visit times of the fields
(qa2_offsets.csv) and the QA error sizes (qa2_targets.json, qa2_error_floor.json).
Output: results/s6b_power_<stage>.json. Checkpoints: data/cache/s6b/power_<stage>/.
"""

from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from irecho import breath, config, sample  # noqa: E402
from irecho.breath import ToyModel  # noqa: E402

L_SUN = 3.828e33  # erg/s (IAU 2015 nominal)
N_REALISATIONS = 200
SEED_TAG = 62
MODELS = {"model1": 1.0, "model1_half": 0.5, "model1_double": 2.0, "model2": None}
TAU_RE_CHOICES = (1.0, 2.0, 5.0)
N_GRID = (20, 25, 30, 35, 40, 50, 60, 70, 80, 100)


def lyu_lags(log_lbol: float) -> tuple[float, float]:
    """Rest-frame W1 and W2 lags in days (Lyu, Rieke & Smith 2019, Eqs. 15-16)."""
    x = 10.0**log_lbol / L_SUN / 1e11
    return 10.0**2.10 * x**0.47, 10.0**2.20 * x**0.45


def load_objects() -> list[dict]:
    dr2 = sample.parse_vizier_tsv((config.CACHE_DIR / "catalogs" / "bass_dr2"
                                   / "J_ApJS_261_2_table9.tsv").read_text(errors="replace"))
    by_id: dict[str, list[dict]] = {}
    for r in dr2:
        by_id.setdefault(r["ID"].strip(), []).append(r)
    out = []
    with open(config.TARGETS_DIR / "sample_v2.csv", newline="") as fh:
        for t in csv.DictReader(fh):
            cands = by_id.get(t["dr2_id"].strip().rstrip("AB"), []) if t["dr2_id"] else []
            vals = [float(c["logLbol"]) for c in cands if c["logLbol"].strip()]
            if vals:
                z = float(t["z"])
                out.append({"name": t["name"], "z": z, "log_lbol": vals[0],
                            "x": float(np.log10(1.0 + z) + 0.5 * vals[0])})
    return out


def load_conditions(stage: str, names: list[str]) -> dict[str, tuple[np.ndarray, dict]]:
    """Visit times and visit errors (mag) for every object. 'prelim': a random
    sample_v1 field and error set; 'final': the object's own field and QA errors."""
    if stage == "exact":
        # 'final' takes the visits of the field (those with a zero-point offset); for some
        # targets the frozen light curve has fewer. Here the frozen mjd column is used.
        cond = load_conditions("final", names)
        frozen = config.PROJECT_ROOT / "data" / "frozen" / "sample_v2"
        out = {}
        for n, (_, e) in cond.items():
            with open(frozen / f"{n}.csv", newline="") as fh:
                out[n] = (np.array(sorted(float(r["mjd"]) for r in csv.DictReader(fh))), e)
        return out
    prefix = "qa_" if stage == "prelim" else "qa2_"
    cadence: dict[str, list[float]] = {}
    with open(config.RESULTS_DIR / f"{prefix}offsets.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            if r["w1_offset_mmag"] and r["w2_offset_mmag"]:
                cadence.setdefault(r["field"], []).append(float(r["mjd"]))
    qa = json.loads((config.RESULTS_DIR / f"{prefix}targets.json").read_text())
    floors = json.loads((config.RESULTS_DIR / f"{prefix}error_floor.json").read_text())["floors_mmag"]

    def err(name: str) -> dict:
        q = qa[name]
        return {b: 1e-3 * float(np.sqrt(q[b]["median_stat_err_mmag"] ** 2
                                        + q[b]["median_offset_err_mmag"] ** 2
                                        + floors[name][b] ** 2)) for b in ("w1", "w2")}

    if stage == "final":
        return {n: (np.array(sorted(cadence[n])), err(n)) for n in names if n in cadence}
    rng = np.random.default_rng([config.SEED, SEED_TAG, 0])
    pool = sorted(cadence)
    return {n: (np.array(sorted(cadence[k])), err(k))
            for n, k in zip(names, rng.choice(pool, size=len(names)))}


def run_model(model: str, objects: list[dict], cond: dict, stage: str) -> list[list[dict]]:
    out_dir = config.CACHE_DIR / "s6b" / f"power_{stage}"
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{model}.json"
    if path.exists():
        return json.loads(path.read_text())
    rng = np.random.default_rng([config.SEED, SEED_TAG, list(MODELS).index(model) + 1,
                                 {"prelim": 0, "final": 1, "exact": 2}[stage]])
    scale = MODELS[model]
    reals = []
    for _ in range(N_REALISATIONS):
        rows = []
        for o in objects:
            if o["name"] not in cond:
                continue
            mjd, err = cond[o["name"]]
            dt1, dt2 = lyu_lags(o["log_lbol"])
            lag = (1.0 + o["z"]) * dt1
            extra = 0.0 if scale is None else scale * (1.0 + o["z"]) * (dt2 - dt1)
            tau = float(rng.choice(TAU_RE_CHOICES)) if scale is None else 0.0
            toy = ToyModel("boundary", tau, host_w1=float(rng.uniform(0.0, 0.8)), lag_days=lag,
                           lag_w2_extra_days=extra)
            ramp = float(rng.uniform(0.5, 1.5)) if rng.random() < 0.3 else 0.0
            v = breath.synthetic_visits(rng, toy, mjd, err, float(rng.uniform(0.15, 0.6)),
                                        ramp_mag=ramp, burn_in_days=3500.0)
            m = breath.measure(v["w1_flux"], v["w1_flux_err"], v["w2_flux"], v["w2_flux_err"],
                               with_in_place=False)
            rows.append({"x": o["x"], "log_lbol": o["log_lbol"], "d": m["d_hyst"],
                         "ok": bool(m["hyst_measurable"])})
        reals.append(rows)
    path.write_text(json.dumps(reals) + "\n")
    return reals


def rho_of(rows: list[dict]) -> float:
    use = [r for r in rows if r["ok"]]
    if len(use) < 8:
        return float("nan")
    return float(spearmanr([r["x"] for r in use], [r["d"] for r in use]).statistic)


def summarise(reals: list[list[dict]], rng) -> dict:
    n_meas = np.array([sum(r["ok"] for r in rows) for rows in reals])
    rho = np.array([rho_of(rows) for rows in reals])
    pooled = [r for rows in reals for r in rows if r["ok"]]
    d = np.array([r["d"] for r in pooled])
    x = np.array([r["x"] for r in pooled])
    hi = x >= np.median(x)
    out = {"n_objects_per_realisation": len(reals[0]),
           "n_measurable_p5_p50_p95": [float(v) for v in np.percentile(n_meas, [5, 50, 95])],
           "rho_p5_p50_p95_at_own_n": [float(v) for v in np.nanpercentile(rho, [5, 50, 95])],
           "d_median": float(np.median(d)),
           "d_median_high_minus_low_luminosity": float(np.median(d[hi]) - np.median(d[~hi])),
           "by_n": {}}
    edges = np.quantile(x, np.linspace(0.0, 1.0, 7))
    out["d_by_x_bin"] = []
    for lo, hi_edge in zip(edges[:-1], edges[1:]):
        sel = (x >= lo) & (x <= hi_edge)
        p16, p50, p84 = np.percentile(d[sel], [16, 50, 84])
        out["d_by_x_bin"].append({"x_mid": float(np.median(x[sel])), "d_p16": float(p16),
                                  "d_p50": float(p50), "d_p84": float(p84)})
    for n in N_GRID:
        if len(pooled) < n:
            continue
        rhos = []
        for _ in range(3000):
            idx = rng.choice(len(pooled), size=n, replace=False)
            rhos.append(spearmanr(x[idx], d[idx]).statistic)
        out["by_n"][str(n)] = {"rho_p5": float(np.percentile(rhos, 5)),
                               "rho_p50": float(np.percentile(rhos, 50)),
                               "rho_p95": float(np.percentile(rhos, 95)), "_rhos": rhos}
    return out


def main() -> None:
    stage = sys.argv[1] if len(sys.argv) > 1 else "prelim"
    if stage not in ("prelim", "final", "exact"):
        raise SystemExit(__doc__)
    objects = load_objects()
    cond = load_conditions(stage, [o["name"] for o in objects])
    rng = np.random.default_rng([config.SEED, SEED_TAG, 99])
    res = {}
    for model in (("model1", "model2") if stage == "exact" else MODELS):
        start = time.time()
        res[model] = summarise(run_model(model, objects, cond, stage), rng)
        print(f"{model}: {time.time() - start:.0f}s", flush=True)
    # null threshold q2(n): 95th percentile of rho under no association (permutation of
    # ranks), and the power of each model-1 variant to exceed it
    null = {}
    for n in N_GRID:
        perm = [spearmanr(np.arange(n), rng.permutation(n)).statistic for _ in range(10000)]
        null[str(n)] = float(np.percentile(perm, 95))
    for model, s in res.items():
        for n, b in s["by_n"].items():
            rhos = np.array(b.pop("_rhos"))
            b["q2_null_p95"] = null[n]
            b["power_rho_above_q2"] = float(np.mean(rhos > null[n]))
            b["models_separable"] = bool(b["rho_p5"] > null[n])
    out = {"stage": stage, "n_objects_with_lbol": len(objects), "n_realisations": N_REALISATIONS,
           "assumptions": {"driver_rms_mag": "uniform 0.15-0.6", "ramp": "30% of objects, 0.5-1.5 mag over 3 yr",
                           "host_w1_fraction": "uniform 0-0.8", "drw_time_scale_days": 1000,
                           "tau_re_model2_yr": list(TAU_RE_CHOICES)},
           "models": res}
    (config.RESULTS_DIR / f"s6b_power_{stage}.json").write_text(json.dumps(out, indent=1) + "\n")
    for model, s in res.items():
        print(model, "n_meas", s["n_measurable_p5_p50_p95"], "rho", [round(v, 2) for v in s["rho_p5_p50_p95_at_own_n"]],
              "D med", round(s["d_median"], 3), "dD hi-lo", round(s["d_median_high_minus_low_luminosity"], 3))
        print("   ", {n: (round(b["rho_p5"], 2), round(b["rho_p50"], 2), "q2", round(b["q2_null_p95"], 2),
                          "power", round(b["power_rho_above_q2"], 2)) for n, b in s["by_n"].items()})


if __name__ == "__main__":
    main()
