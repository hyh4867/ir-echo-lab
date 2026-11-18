"""Stage 9: is the colour asymmetry explained by an inter-band lag alone? (docs/PREDICTIONS_9.md)

    python scripts/run_stage9.py ident [workers]   B1: identifiability on synthetic samples (real cadence and errors)
    python scripts/run_stage9.py ident_summary     B1: criteria of section 3.2
    python scripts/run_stage9.py fit               B2: the ladder T0-T3 on the frozen light curves
    python scripts/run_stage9.py null [workers]    B2: synthetic samples with the fitted lags as truth ('fitted')
    python scripts/run_stage9.py verdict           B2-B4: model comparison, verdicts, lag against luminosity
    python scripts/run_stage9.py mock5 [workers]   B5: lag-luminosity mock with a driver matched in amplitude

Checkpoints: data/cache/s9/ (one file per synthetic sample; a rerun skips what exists).
Needs the frozen light curves and results/s8_error_correlation_*.json.
"""

from __future__ import annotations

import csv
import json
import multiprocessing
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import run_stage8 as r8  # noqa: E402
from irecho import breath, config, hyst8  # noqa: E402
from irecho import transfer9 as t9  # noqa: E402

RES = config.RESULTS_DIR
CKPT = config.CACHE_DIR / "s9"
SEED_TAG = 9
TAGS = ("v1", "v2")
DRIVERS = {"fast": 200.0, "slow": 1000.0}        # DRW time scale in days
N_IDENT = 50
N_MOCK5 = 100
SCATTERS = (0.0, 0.15, 0.30)
RHO_OBSERVED = 0.129
OBSERVED_V2 = {"n_measurable_with_lbol": 70, "median_d": 0.268}
ROW_KEYS = ("tau_t1", "tau_lo", "tau_hi", "tau_half_width", "tau_at_edge", "well_determined", "tau_t2", "width_t2",
            "ds_01", "ds_02", "ds_12", "chi2_t0", "chi2_t1", "chi2_t2", "a_t1", "b_t1",
            "eta_t0", "eta_t1", "eta_t2", "sigma_eta_t1", "sigma_eta_t2", "ds_eta_t0", "ds_eta_t1", "ds_eta_t2",
            "eta_colour_t0", "eta_colour_t1", "eta_colour_t2")
OFF_KEYS = tuple(f"{k}_{n}" for n in ("obs", "t0", "t1", "t2")
                 for k in ("delta", "sigma_delta", "d_colour", "sigma_colour", "usable", "comparable"))
_OBJ: dict = {}
_DELAYS: dict = {}


# ----------------------------------------------------------------------------- objects
def load(tag: str) -> list[dict]:
    """Frozen light curves with what the synthetic samples need: the W1 lag and the Lyu delay (observed frame)."""
    if tag in _OBJ:
        return _OBJ[tag]
    try:
        objs = r8.load_sample(tag)
    except ValueError as exc:
        raise SystemExit("the object tables have blank catalogue columns (the public copy leaves them out). "
                         "Run scripts/restore_catalog_columns.py first.") from exc
    lag1, delay = [], []
    for o in objs:
        if o["log_lbol"] is not None:
            w1, w2 = t9.lyu_lags_days(o["log_lbol"])
            o["lyu_w1_rest"], o["lyu_delay_rest"] = w1, w2 - w1
            lag1.append((1.0 + o["z"]) * w1)
            delay.append((1.0 + o["z"]) * (w2 - w1))
    med_lag, med_delay = float(np.median(lag1)), float(np.median(delay))
    for o in objs:
        has = o["log_lbol"] is not None
        o["lag1"] = (1.0 + o["z"]) * o["lyu_w1_rest"] if has else med_lag
        o["lyu_delay"] = (1.0 + o["z"]) * o["lyu_delay_rest"] if has else med_delay
        o["prep"] = t9.prepare(o["mjd"])
    _OBJ[tag] = objs
    return objs


def slim(fit: dict, off: dict) -> dict:
    row = {"usable": bool(fit.get("usable")), "n_used": int(fit.get("n_used", 0))}
    if row["usable"]:
        for k in ROW_KEYS:
            v = fit[k]
            row[k] = bool(v) if isinstance(v, (bool, np.bool_)) else float(v)
    for k in OFF_KEYS:
        v = off[k]
        row[k] = bool(v) if isinstance(v, (bool, np.bool_)) else float(v)
    return row


def lag_luminosity(objs: list[dict], rows: list[dict]) -> dict:
    """Section 5: k = median of tau_rest / delay_Lyu and Spearman rho of tau_rest against log L_bol."""
    tau, lum, ratio = [], [], []
    for o, r in zip(objs, rows):
        if r is None or not r.get("usable") or o["log_lbol"] is None:
            continue
        rest = r["tau_t1"] / (1.0 + o["z"])
        tau.append(rest)
        lum.append(o["log_lbol"])
        ratio.append(rest / o["lyu_delay_rest"])
    if len(tau) < 10:
        return {"n": len(tau), "k": None, "rho": None}
    return {"n": len(tau), "k": float(np.median(ratio)), "rho": float(spearmanr(lum, tau).statistic),
            "median_tau_rest": float(np.median(tau))}


# ----------------------------------------------------------------------------- synthetic samples
def one_sample(args) -> dict:
    """One synthetic sample of every requested truth model for (sample tag, driver, index)."""
    tag, driver, index, truths = args
    objs = load(tag)
    seed = [config.SEED, SEED_TAG, TAGS.index(tag), list(DRIVERS).index(driver), index]
    rows = {name: [] for name in truths}
    for k, o in enumerate(objs):
        mags = t9.drw_magnitudes(np.random.default_rng(seed + [k, 0]), t9.driver_length(o["mjd"]), t9.GRID_STEP_DAYS,
                                 DRIVERS[driver])
        for name in truths:
            if name == "fitted":
                kind, delay, tau_re = "linear", _DELAYS[tag].get(o["name"]), None
            else:
                kind, delay, tau_re = t9.TRUTHS[name]
                delay = o["lyu_delay"] if delay == "lyu" else delay
            row = None
            if delay is not None and o["prep"]["c"] is not None:
                truth = t9.mock_truth(mags, o["mjd"], o["f1"], o["s1"], o["f2"], o["s2"], kind, o["lag1"], delay, tau_re)
                if truth is not None:
                    f1, f2 = hyst8.add_noise(np.random.default_rng(seed + [k, 1]), truth[0], o["s1"], truth[1], o["s2"], o["r"])
                    if np.all(f1 > 0) and np.all(f2 > 0):
                        fit = t9.fit_object(o["mjd"], f1, o["s1"], f2, o["s2"], o["r"], prep=o["prep"])
                        row = slim(fit, t9.offsets(f1, o["s1"], f2, o["s2"], o["r"], fit))
                        row["delay"] = float(delay)
            rows[name].append(row)
    out = {}
    for name in truths:
        good = [r for r in rows[name] if r is not None]
        out[name] = {"stats": t9.sample_statistics(good), "lagl": lag_luminosity(objs, rows[name]),
                     "tau": [None if r is None or not r["usable"] else r["tau_t1"] for r in rows[name]],
                     "well": [None if r is None or not r["usable"] else bool(r["well_determined"]) for r in rows[name]],
                     "delay": [None if r is None else r["delay"] for r in rows[name]]}
    return {"tag": tag, "driver": driver, "index": index, "truths": out}


def run_samples(kind: str, truths: list[str], workers: int, n: int = N_IDENT) -> None:
    CKPT.mkdir(parents=True, exist_ok=True)
    for tag in TAGS:
        load(tag)
    todo = [(tag, driver, i, truths) for tag in TAGS for driver in DRIVERS for i in range(n)
            if not (CKPT / f"{kind}_{tag}_{driver}_{i:03d}.json").exists()]
    print(f"{kind}: {len(todo)} synthetic samples to make", flush=True)
    t0 = time.time()
    with multiprocessing.get_context("fork").Pool(workers) as pool:
        for done, res in enumerate(pool.imap_unordered(one_sample, todo), start=1):
            path = CKPT / f"{kind}_{res['tag']}_{res['driver']}_{res['index']:03d}.json"
            path.write_text(json.dumps(res["truths"]) + "\n")
            if done % 10 == 0 or done == len(todo):
                print(f"{time.time() - t0:7.0f} s  {done}/{len(todo)}", flush=True)


def read_samples(kind: str, tag: str, driver: str, n: int = N_IDENT) -> list[dict]:
    out = []
    for i in range(n):
        path = CKPT / f"{kind}_{tag}_{driver}_{i:03d}.json"
        if path.exists():
            out.append(json.loads(path.read_text()))
    return out


def band(values) -> list[float] | None:
    """'Null range' of the plan: the sorted values without the smallest and the largest one."""
    v = np.sort(np.array([x for x in values if x is not None and np.isfinite(x)], dtype=float))
    return [float(v[1]), float(v[-2])] if v.size >= 10 else None


def stat(samples: list[dict], truth: str, key: str) -> np.ndarray:
    return np.array([s[truth]["stats"].get(key, np.nan) for s in samples], dtype=float)


def frac_outside(values, lo: float, hi: float, side: str = "both") -> float:
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if side == "above":
        return float(np.mean(v > hi))
    return float(np.mean((v < lo) | (v > hi)))


LINEAR_TRUTHS = ("T0", "T1_20", "T1_40", "T1_80", "T1_160", "T1_lyu", "T2_40", "T2_80")


def ident_summary() -> None:
    out = {"note": "stage 9 B1: identifiability with the real visit times, errors and error correlations",
           "n_samples_per_configuration": N_IDENT, "drivers_tau_days": DRIVERS,
           "null_range": "sorted values of the synthetic samples without the smallest and the largest one", "configs": {}}
    flags = {"lag_recovered": [], "t0_vs_t1": [], "t1_vs_t2": [], "fixed_vs_t3": []}
    for tag in TAGS:
        objs = load(tag)
        for driver in DRIVERS:
            samples = read_samples("ident", tag, driver)
            rec = {"n_samples": len(samples), "truths": {}}
            for name in t9.TRUTHS:
                per = {}
                for key in ("n_fit", "n_offset", "n_well", "median_tau", "frac_tau_positive", "mean_ds_01", "mean_ds_02", "mean_ds_12",
                            "mean_ds_eta_t0", "mean_ds_eta_t1", "mean_ds_eta_t2", "median_d_colour_obs", "median_d_colour_pred_t1",
                            "median_d_colour_pred_t2", "median_d_colour_res_t1", "median_d_colour_res_t2", "frac_res_bluer_t1",
                            "frac_res_bluer_t2", "ratio_pred_obs_t1", "ratio_pred_obs_t2", "median_width_t2"):
                    v = stat(samples, name, key)
                    per[key] = {"median": float(np.nanmedian(v)), "null_range": band(v),
                                "p16_p84": [float(x) for x in np.nanpercentile(v, [16, 84])]}
                # per-object accuracy of the lag
                err_all, err_well, n_well, n_fit = [], [], 0, 0
                for s in samples:
                    for tau, well, delay in zip(s[name]["tau"], s[name]["well"], s[name]["delay"]):
                        if tau is None:
                            continue
                        n_fit += 1
                        err_all.append(tau - delay)
                        if well:
                            n_well += 1
                            err_well.append(tau - delay)
                per["per_object"] = {"rms_error_days": float(np.sqrt(np.mean(np.square(err_all)))),
                                     "median_error_days": float(np.median(err_all)),
                                     "frac_well_determined": n_well / max(n_fit, 1),
                                     "rms_error_well_days": float(np.sqrt(np.mean(np.square(err_well)))) if err_well else None}
                per["lag_luminosity"] = {"k_median": float(np.nanmedian([s[name]["lagl"]["k"] for s in samples])),
                                         "k_null_range": band([s[name]["lagl"]["k"] for s in samples]),
                                         "rho_null_range": band([s[name]["lagl"]["rho"] for s in samples])}
                rec["truths"][name] = per
            # --- criteria of section 3.2
            crit = {}
            ok = []
            for name, delay in (("T1_40", 40.0), ("T1_80", 80.0)):
                v = stat(samples, name, "median_tau") - delay
                p16, p84 = np.percentile(v, [16, 84])
                good = bool(abs(np.median(v)) <= max(10.0, 0.25 * delay) and 0.5 * (p84 - p16) <= 15.0)
                crit[f"lag_bias_{name}"] = {"median_error_days": float(np.median(v)), "half_width_16_84_days": float(0.5 * (p84 - p16)),
                                            "pass": good}
                ok.append(good)
            crit["lag_recovered"] = bool(all(ok))
            pw = {}
            for key in ("median_tau", "mean_ds_01"):
                lo, hi = band(stat(samples, "T0", key))
                pw[key] = {name: frac_outside(stat(samples, name, key), lo, hi, "above")
                           for name in ("T1_20", "T1_40", "T1_80", "T1_160", "T1_lyu", "T2_40", "T2_80", "T3_1yr", "T3_3yr", "T0b")}
            crit["t0_vs_t1_power"] = pw
            crit["t0_vs_t1"] = bool(max(pw["median_tau"]["T1_lyu"], pw["mean_ds_01"]["T1_lyu"]) >= 0.8)
            p12 = {}
            for a, b in (("T1_40", "T2_40"), ("T1_80", "T2_80")):
                lo, hi = band(stat(samples, a, "mean_ds_12"))
                p12[b] = frac_outside(stat(samples, b, "mean_ds_12"), lo, hi, "above")
            crit["t1_vs_t2_power"] = p12
            crit["t1_vs_t2"] = bool(min(p12.values()) >= 0.8)
            p3 = {}
            for key in ("median_d_colour_res_t2", "mean_ds_eta_t2"):
                bands = [band(stat(samples, name, key)) for name in LINEAR_TRUTHS]
                lo, hi = min(b[0] for b in bands), max(b[1] for b in bands)
                p3[key] = {"envelope_of_fixed_kernel_truths": [lo, hi],
                           **{name: frac_outside(stat(samples, name, key), lo, hi) for name in ("T3_1yr", "T3_3yr", "T0b")}}
            crit["fixed_vs_t3_power"] = p3
            crit["fixed_vs_t3"] = bool(all(max(p3[k][name] for k in p3) >= 0.8 for name in ("T3_1yr", "T3_3yr")))
            for k in flags:
                flags[k].append(crit[k])
            rec["criteria"] = crit
            rec["n_objects"] = len(objs)
            out["configs"][f"{tag}_{driver}"] = rec
    # a question counts as answered 'yes' only if every sample and driver agrees (section 3.2)
    out["overall"] = {k: bool(all(v)) for k, v in flags.items()}
    out["overall_by_config"] = {k: dict(zip(out["configs"], v)) for k, v in flags.items()}
    (RES / "s9_identifiability.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out["overall"], indent=1))
    for name, rec in out["configs"].items():
        c = rec["criteria"]
        print(name, "lag bias", {k: (round(c[f"lag_bias_{k}"]["median_error_days"], 1), round(c[f"lag_bias_{k}"]["half_width_16_84_days"], 1))
                                 for k in ("T1_40", "T1_80")}, "| T0 vs T1_lyu", c["t0_vs_t1_power"]["median_tau"]["T1_lyu"],
              c["t0_vs_t1_power"]["mean_ds_01"]["T1_lyu"], "| T1 vs T2", c["t1_vs_t2_power"],
              "| T3", {k: (v["T3_1yr"], v["T3_3yr"]) for k, v in c["fixed_vs_t3_power"].items()})


# ----------------------------------------------------------------------------- real data
def fit_real(r_override: float | None = None) -> dict:
    out = {}
    for tag in TAGS:
        rows = []
        for o in load(tag):
            r = o["r"] if r_override is None else r_override
            fit = t9.fit_object(o["mjd"], o["f1"], o["s1"], o["f2"], o["s2"], r, prep=o["prep"])
            row = slim(fit, t9.offsets(o["f1"], o["s1"], o["f2"], o["s2"], r, fit))
            row["name"] = o["name"]
            rows.append(row)
        out[tag] = rows
    return out


def run_fit() -> None:
    real = fit_real()
    zero = fit_real(0.0)
    s8 = {}
    with open(RES / "s8_direct_objects.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            s8[(r["sample"], r["name"])] = r
    fields = ["sample", "name", "usable", "n_used", *ROW_KEYS, *OFF_KEYS, "d_colour_res_t1", "d_colour_res_t2", "stage8_d_colour"]
    with open(RES / "s9_objects.csv", "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(fields)
        for tag in TAGS:
            for row in real[tag]:
                rec = dict(row, sample=tag)
                for k in ("t1", "t2"):
                    rec[f"d_colour_res_{k}"] = row["d_colour_obs"] - row[f"d_colour_{k}"]
                rec["stage8_d_colour"] = s8[(tag, row["name"])]["d_colour"]
                w.writerow(["" if rec.get(k) is None or (isinstance(rec.get(k), float) and not np.isfinite(rec[k]))
                            else (f"{rec[k]:.6g}" if isinstance(rec[k], float) else (int(rec[k]) if isinstance(rec[k], bool) else rec[k]))
                            for k in fields])
    out = {"note": "stage 9 B2: transfer W1 -> W2 with a top-hat kernel on linearly interpolated W1; statistics of section 3.1",
           "samples": {}}
    for tag in TAGS:
        rows = real[tag]
        st = t9.sample_statistics(rows)
        off = [r for r in rows if r["usable"] and r["usable_obs"] and r["comparable_obs"] and all(r[f"usable_{k}"] for k in ("t0", "t1", "t2"))]
        rec = {"n_objects": len(rows), "statistics": st, "statistics_r_zero": t9.sample_statistics(zero[tag])}
        obs = np.array([r["d_colour_obs"] for r in off])
        for k in ("t1", "t2"):
            pred = np.array([r[f"d_colour_{k}"] for r in off])
            res = obs - pred
            sig = np.array([r["sigma_colour_obs"] for r in off])
            big = np.abs(res) > 3.0 * sig
            rec[f"offsets_{k}"] = {
                "spearman_pred_obs": float(spearmanr(pred, obs).statistic),
                "n_res_beyond_3sigma": int(big.sum()),
                "objects_res_beyond_3sigma": [{"name": r["name"], "res_mag": float(x), "sigma_mag": float(s), "obs_mag": float(o_)}
                                              for r, x, s, o_, b in zip(off, res, sig, obs, big) if b],
                "res_p16_p50_p84_mag": [float(x) for x in np.percentile(res, [16, 50, 84])],
                "frac_objects_abs_res_smaller_than_abs_obs": float(np.mean(np.abs(res) < np.abs(obs)))}
        both = [(float(s8[(tag, r["name"])]["d_colour"]), r["d_colour_obs"]) for r in off if s8[(tag, r["name"])]["d_colour"]]
        rec["stage8_against_stage9_visits"] = {
            "n": len(both), "median_stage8_all_visits_mag": float(np.median([b[0] for b in both])),
            "median_stage9_used_visits_mag": float(np.median([b[1] for b in both])),
            "spearman": float(spearmanr([b[0] for b in both], [b[1] for b in both]).statistic)}
        fit = [r for r in rows if r["usable"]]
        rec["n_used_visits_p5_p50_p95"] = [float(x) for x in np.percentile([r["n_used"] for r in fit], [5, 50, 95])]
        rec["chi2_per_dof_t1_median"] = float(np.median([r["chi2_t1"] / max(r["n_used"] - 3, 1) for r in fit]))
        rec["chi2_per_dof_t0_median"] = float(np.median([r["chi2_t0"] / max(r["n_used"] - 2, 1) for r in fit]))
        rec["tau_t1_p5_p16_p50_p84_p95"] = [float(x) for x in np.percentile([r["tau_t1"] for r in fit], [5, 16, 50, 84, 95])]
        rec["frac_tau_at_edge"] = float(np.mean([r["tau_at_edge"] for r in fit]))
        rec["width_t2_counts"] = {str(int(w)): int(sum(r["width_t2"] == w for r in fit)) for w in t9.WIDTHS_DAYS}
        out["samples"][tag] = rec
        print(tag, json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in st.items()}))
    (RES / "s9_transfer.json").write_text(json.dumps(out, indent=1) + "\n")


def load_fitted_delays() -> None:
    with open(RES / "s9_objects.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            if r["usable"] == "1":
                _DELAYS.setdefault(r["sample"], {})[r["name"]] = float(r["tau_t1"])


# ----------------------------------------------------------------------------- B5
def mock5_sample(args) -> dict:
    driver, scatter, index = args
    objs = [o for o in load("v2") if o["log_lbol"] is not None]
    seed = [config.SEED, SEED_TAG, 50, list(DRIVERS).index(driver), SCATTERS.index(scatter), index]
    rng = np.random.default_rng(seed)
    xs, ds = [], []
    for k, o in enumerate(objs):
        mags = t9.drw_magnitudes(np.random.default_rng(seed + [k, 0]), t9.driver_length(o["mjd"]), t9.GRID_STEP_DAYS, DRIVERS[driver])
        f = 10.0 ** rng.normal(0.0, scatter) if scatter > 0 else 1.0
        truth = t9.mock_truth(mags, o["mjd"], o["f1"], o["s1"], o["f2"], o["s2"], "linear", o["lag1"] * f, o["lyu_delay"] * f)
        if truth is None:
            continue
        f1, f2 = hyst8.add_noise(np.random.default_rng(seed + [k, 1]), truth[0], o["s1"], truth[1], o["s2"], o["r"])
        if np.any(f1 <= 0) or np.any(f2 <= 0):
            continue
        m = breath.measure(f1, o["s1"], f2, o["s2"], with_in_place=False)
        if m["hyst_measurable"]:
            xs.append(np.log10(1.0 + o["z"]) + 0.5 * o["log_lbol"])
            ds.append(m["d_hyst"])
    return {"driver": driver, "scatter": scatter, "index": index, "n_measurable": len(ds),
            "median_d": float(np.median(ds)) if ds else None,
            "rho": float(spearmanr(xs, ds).statistic) if len(ds) >= 10 else None}


def run_mock5(workers: int) -> None:
    CKPT.mkdir(parents=True, exist_ok=True)
    load("v2")
    out = {"note": "stage 9 B5: inter-band lag from Lyu et al. 2019 (linear response, pure shift) with a driver whose amplitude "
                   "is matched to the observed W2 rms of every object; sample_v2 objects with a bolometric luminosity",
           "n_objects": int(sum(o["log_lbol"] is not None for o in load("v2"))), "observed": dict(OBSERVED_V2, rho=RHO_OBSERVED),
           "configs": {}}
    t0 = time.time()
    with multiprocessing.get_context("fork").Pool(workers) as pool:
        for driver in DRIVERS:
            for scatter in SCATTERS:
                key = f"{driver}_scatter{scatter:.2f}"
                path = CKPT / f"mock5_{key}.json"
                if path.exists():
                    rows = json.loads(path.read_text())
                else:
                    rows = pool.map(mock5_sample, [(driver, scatter, i) for i in range(N_MOCK5)], chunksize=2)
                    path.write_text(json.dumps(rows) + "\n")
                    print(f"{time.time() - t0:6.0f} s  {key}", flush=True)
                n = np.array([r["n_measurable"] for r in rows], dtype=float)
                d = np.array([r["median_d"] for r in rows if r["median_d"] is not None], dtype=float)
                rho = np.array([r["rho"] for r in rows if r["rho"] is not None], dtype=float)
                n_rng, d_rng = np.percentile(n, [5, 95]), np.percentile(d, [5, 95])
                rep = bool(n_rng[0] <= OBSERVED_V2["n_measurable_with_lbol"] <= n_rng[1] and d_rng[0] <= OBSERVED_V2["median_d"] <= d_rng[1])
                rec = {"driver_tau_days": DRIVERS[driver], "scatter_dex": scatter, "n_samples": len(rows),
                       "n_measurable_p5_p50_p95": [float(x) for x in np.percentile(n, [5, 50, 95])],
                       "median_d_p5_p50_p95": [float(x) for x in np.percentile(d, [5, 50, 95])], "reproduces_the_sample": rep}
                if rho.size >= 20:
                    rec.update(rho_p5_p50_p95=[float(x) for x in np.percentile(rho, [5, 50, 95])],
                               fraction_at_or_below_observed=float(np.mean(rho <= RHO_OBSERVED)),
                               outside=bool(np.mean(rho <= RHO_OBSERVED) < 0.05))
                out["configs"][key] = rec
    (RES / "s9_mock5.json").write_text(json.dumps(out, indent=1) + "\n")
    for k, v in out["configs"].items():
        print(k, v)


def main() -> None:
    what = sys.argv[1] if len(sys.argv) > 1 else ""
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    if what == "ident":
        run_samples("ident", list(t9.TRUTHS), workers)
        ident_summary()
    elif what == "ident_summary":
        ident_summary()
    elif what == "fit":
        run_fit()
    elif what == "null":
        load_fitted_delays()
        run_samples("fitted", ["fitted"], workers)
    elif what == "verdict":
        import stage9_verdict
        stage9_verdict.main()
    elif what == "mock5":
        run_mock5(workers)
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
