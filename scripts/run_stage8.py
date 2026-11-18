"""Stage 8: validation of the hysteresis statistic (docs/VALIDATION_8.md).

    python scripts/run_stage8.py null [workers]   B1 + the null test of the direct estimator -> results/s8_null.json
    python scripts/run_stage8.py direct           B2 on the frozen data, unTimely, scan parity -> results/s8_direct*.{json,csv}
    python scripts/run_stage8.py lag [workers]    B3 -> results/s8_lag_model.json
    python scripts/run_stage8.py select           B4 -> results/s8_selection.json

Needs the frozen light curves (data/frozen/) and results/s8_error_correlation_*.json
(scripts/run_stage8_noise.py). Long stages keep checkpoints in data/cache/s8/.
"""

from __future__ import annotations

import csv
import json
import sys
import time
import multiprocessing
from pathlib import Path

import numpy as np
from scipy.stats import binomtest, mannwhitneyu, spearmanr, wilcoxon

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from irecho import breath, config, hyst8  # noqa: E402
from irecho import stage5b as s5  # noqa: E402

RES = config.RESULTS_DIR
CKPT = config.CACHE_DIR / "s8"
FROZEN = {"v1": config.PROJECT_ROOT / "data" / "frozen" / "sample_v1",
          "v2": config.PROJECT_ROOT / "data" / "frozen" / "sample_v2"}
SEED_TAG = 8
N_NULL, N_VARIANT = 200, 100
NULLS = ("linear", "single_t")
OBSERVED = {"v1": {"n": 71, "median_d": 0.228}, "v2": {"n": 94, "median_d": 0.268}}


# ----------------------------------------------------------------------------- samples
def sample_names(tag: str) -> list[dict]:
    """Objects of the registered analyses: sample_v1 main sample (5b) and sample_v2 (6B)."""
    if tag == "v1":
        with open(RES / "s5b_objects.csv", newline="") as fh:
            return [{"name": r["name"], "type_group": r["type_group"], "z": float(r["z"]), "qa_flags": r["qa_flags"],
                     "d": float(r["sl_d_hyst"]) if r["sl_hyst_measurable"] == "1" else None,
                     "sigma_d": float(r["sl_sigma_d"]) if r["sl_hyst_measurable"] == "1" else None,
                     "log_lbol": float(r["loglbol"]) if r["loglbol"] else None, "type": r["type"]}
                    for r in csv.DictReader(fh) if r["sample"] == "main"]
    with open(RES / "s6b_objects.csv", newline="") as fh:
        return [{"name": r["name"], "type": r["type"], "z": float(r["z"]), "qa_flags": r["qa_flags"],
                 "d": float(r["d_hyst"]) if r["hyst_measurable"] == "1" else None,
                 "sigma_d": float(r["sigma_d"]) if r["hyst_measurable"] == "1" else None,
                 "log_lbol": float(r["log_lbol"]) if r["log_lbol"] else None}
                for r in csv.DictReader(fh)]


def load_sample(tag: str) -> list[dict]:
    corr = json.loads((RES / f"s8_error_correlation_{tag}.json").read_text())
    r_med = corr["target_r_median"]
    out = []
    for o in sample_names(tag):
        v = s5.read_frozen(o["name"], frozen_dir=FROZEN[tag])
        ok = np.isfinite(v["w1_flux_jy"]) & np.isfinite(v["w2_flux_jy"])
        rec = corr["objects"].get(o["name"], {})
        r = rec.get("r")
        out.append(dict(o, mjd=v["mjd"][ok], f1=v["w1_flux_jy"][ok], s1=v["w1_flux_err_jy"][ok],
                        f2=v["w2_flux_jy"][ok], s2=v["w2_flux_err_jy"][ok],
                        r=float(np.clip(r if r is not None else r_med, -0.95, 0.95)), r_estimated=r is not None))
    return out


# ----------------------------------------------------------------------------- B1: null test
VARIANTS = {
    "primary": {},
    "r_zero": {"r": 0.0},
    "r_0p8": {"r": 0.8},
    "smooth_truth": {"smooth": True},
    "pivot_mean": {"pivot": "mean"},
    "truth_reference": {"reference": "truth"},
    "truth_branches": {"branches": "truth"},
    "noise_w1_only": {"bands": "w1"},
    "noise_w2_only": {"bands": "w2"},
}
_OBJECTS: dict = {}


def prepare_truths(obj: dict, smooth: bool) -> dict | None:
    t2 = hyst8.smooth3(obj["f2"]) if smooth else obj["f2"]
    lin = hyst8.truth_linear(obj["f1"], obj["s1"], t2, obj["s2"])
    one = hyst8.truth_single_temperature(obj["f1"], obj["s1"], t2, obj["s2"])
    return {"t2": t2, "linear": lin, "single_t": one}


def one_realisation(args) -> dict:
    tag, kind, variant, index = args
    opts = VARIANTS[variant]
    objs = _OBJECTS[(tag, bool(opts.get("smooth")))]
    rng = np.random.default_rng([config.SEED, SEED_TAG, ("v1", "v2").index(tag), NULLS.index(kind),
                                 list(VARIANTS).index(variant), index])
    d, dsig, delta_sig, dcol, per_obj = [], [], [], [], []
    for k, o in enumerate(objs):
        t1 = o["truth"][kind]
        if t1 is None:
            continue
        t2 = o["truth"]["t2"]
        r = opts.get("r", o["r"])
        f1, f2 = hyst8.add_noise(rng, t1, o["s1"], t2, o["s2"], r, opts.get("bands", "both"))
        if np.any(f1 <= 0) or np.any(f2 <= 0):
            continue
        ref = None
        if opts.get("reference") == "truth":
            prep = breath._prepare(t1, o["s1"], t2, o["s2"], breath.SETTINGS)
            ref = prep["a"] if prep else None
        rising = hyst8.branch_rising(t1, t2) if opts.get("branches") == "truth" else None
        if variant in ("primary", "r_zero", "r_0p8", "smooth_truth", "noise_w1_only", "noise_w2_only"):
            m = breath.measure(f1, o["s1"], f2, o["s2"], with_in_place=False)
            res = {"d": m["d_hyst"], "sigma_d": m["sigma_d"], "measurable": m["hyst_measurable"]}
        else:
            res = hyst8.d_statistic(f1, o["s1"], f2, o["s2"], pivot=opts.get("pivot", "peak"), reference=ref, rising=rising)
        if res["measurable"]:
            d.append(res["d"])
            dsig.append(res["sigma_d"])
            per_obj.append((k, res["d"]))
        if variant in ("primary", "r_zero", "r_0p8", "smooth_truth"):
            dr = hyst8.direct_offset(f1, o["s1"], f2, o["s2"], r=r)
            if dr["comparable"]:
                delta_sig.append(dr["delta"] / dr["sigma_delta"])
                dcol.append(dr["d_colour"])
    d, dcol, delta_sig = np.array(d), np.array(dcol), np.array(delta_sig)
    out = {"n_measurable": int(d.size), "median_d": float(np.median(d)) if d.size else None,
           "frac_positive": float(np.mean(d > 0)) if d.size else None,
           "wilcoxon_p": float(wilcoxon(d, alternative="greater").pvalue) if d.size >= 8 else None,
           "median_sigma_d": float(np.median(dsig)) if d.size else None, "per_obj": per_obj}
    if dcol.size:
        out.update(n_comparable=int(dcol.size), frac_delta_positive=float(np.mean(delta_sig > 0)),
                   median_d_colour=float(np.median(dcol)), mean_pull=float(np.mean(delta_sig)), std_pull=float(np.std(delta_sig)),
                   wilcoxon_p_delta=float(wilcoxon(delta_sig, alternative="greater").pvalue) if dcol.size >= 8 else None)
    return out


def summarise_null(rows: list[dict], n_objects: int, observed_median: float | None) -> dict:
    med = np.array([r["median_d"] for r in rows if r["median_d"] is not None])
    pos = np.array([r["frac_positive"] for r in rows if r["frac_positive"] is not None])
    nm = np.array([r["n_measurable"] for r in rows])
    p = np.array([r["wilcoxon_p"] for r in rows if r["wilcoxon_p"] is not None])
    out = {"n_realisations": len(rows), "n_objects": n_objects,
           "n_measurable_p5_p50_p95": [float(x) for x in np.percentile(nm, [5, 50, 95])],
           "median_d_p2p5_p16_p50_p84_p97p5_p99": [float(x) for x in np.percentile(med, [2.5, 16, 50, 84, 97.5, 99])],
           "frac_positive_mean": float(pos.mean()), "frac_positive_p5_p95": [float(x) for x in np.percentile(pos, [5, 95])],
           "wilcoxon_rejection_rate": float(np.mean(p < 0.05)) if p.size else None,
           "median_sigma_d": float(np.median([r["median_sigma_d"] for r in rows if r["median_sigma_d"] is not None]))}
    m50 = out["median_d_p2p5_p16_p50_p84_p97p5_p99"][2]
    rej = out["wilcoxon_rejection_rate"] or 0.0
    if abs(m50) <= 0.03 and 0.45 <= out["frac_positive_mean"] <= 0.55 and rej <= 0.10:
        out["label"] = "치우침 없음"
    elif m50 >= 0.10 or rej >= 0.50:
        out["label"] = "치우침 큼"
    else:
        out["label"] = "부분적"
    if observed_median is not None:
        q = out["median_d_p2p5_p16_p50_p84_p97p5_p99"]
        out["observed_median_d"] = observed_median
        out["observed_percentile_in_null"] = float(np.mean(med < observed_median) * 100.0)
        out["observed_vs_null"] = ("null로 설명되지 않음" if observed_median > q[5] else
                                   "null로 설명됨" if q[0] <= observed_median <= q[4] else "일부 설명됨")
        out["observed_minus_null_median"] = observed_median - m50
    dr = [r for r in rows if "n_comparable" in r]
    if dr:
        fp = np.array([r["frac_delta_positive"] for r in dr])
        pw = np.array([r["wilcoxon_p_delta"] for r in dr if r["wilcoxon_p_delta"] is not None])
        out["direct"] = {"n_comparable_p5_p50_p95": [float(x) for x in np.percentile([r["n_comparable"] for r in dr], [5, 50, 95])],
                         "frac_delta_positive_mean": float(fp.mean()),
                         "frac_delta_positive_p2p5_p50_p97p5_p99": [float(x) for x in np.percentile(fp, [2.5, 50, 97.5, 99])],
                         "median_d_colour_p2p5_p50_p97p5": [float(x) for x in np.percentile([r["median_d_colour"] for r in dr], [2.5, 50, 97.5])],
                         "pull_mean": float(np.mean([r["mean_pull"] for r in dr])), "pull_std": float(np.mean([r["std_pull"] for r in dr])),
                         "wilcoxon_rejection_rate": float(np.mean(pw < 0.05)) if pw.size else None}
        out["direct"]["usable_under_null"] = bool(0.45 <= out["direct"]["frac_delta_positive_mean"] <= 0.55
                                                  and (out["direct"]["wilcoxon_rejection_rate"] or 0.0) <= 0.10)
    return out


def run_null(workers: int) -> None:
    CKPT.mkdir(parents=True, exist_ok=True)
    samples = {tag: load_sample(tag) for tag in ("v1", "v2")}
    for tag, objs in samples.items():
        for smooth in (False, True):
            _OBJECTS[(tag, smooth)] = [dict(o, truth=prepare_truths(o, smooth)) for o in objs]
    out = {"note": "stage 8 B1: null light curves on the real visit times, errors and W2 light curves",
           "settings": {"n_null": N_NULL, "n_variant": N_VARIANT, "dust_fraction": hyst8.DUST_FRACTION, "variants": VARIANTS},
           "samples": {}}
    t0 = time.time()
    # fork: the workers inherit the prepared light curves (_OBJECTS)
    with multiprocessing.get_context("fork").Pool(workers) as pool:
        for tag, objs in samples.items():
            rec = {"n_objects": len(objs), "n_r_estimated": int(sum(o["r_estimated"] for o in objs)),
                   "r_median": float(np.median([o["r"] for o in objs])), "nulls": {}}
            for kind in NULLS:
                rec["nulls"][kind] = {"n_objects_with_truth": int(sum(o["truth"][kind] is not None for o in _OBJECTS[(tag, False)]))}
                for variant in VARIANTS:
                    path = CKPT / f"null_{tag}_{kind}_{variant}.json"
                    if path.exists():
                        rows = json.loads(path.read_text())
                    else:
                        n = N_NULL if variant == "primary" else N_VARIANT
                        rows = pool.map(one_realisation, [(tag, kind, variant, i) for i in range(n)], chunksize=4)
                        path.write_text(json.dumps(rows) + "\n")
                        print(f"{time.time() - t0:6.0f} s  {tag} {kind} {variant}", flush=True)
                    s = summarise_null(rows, len(objs), OBSERVED[tag]["median_d"] if variant in ("primary", "r_zero", "r_0p8", "smooth_truth") else None)
                    if variant == "primary":
                        # (f) objects whose two branches have the same number of visits (within one)
                        # is not available per realisation; per-object means are kept instead
                        acc: dict[int, list[float]] = {}
                        for r in rows:
                            for k, dval in r["per_obj"]:
                                acc.setdefault(int(k), []).append(dval)
                        names = [o["name"] for o in objs]
                        s["per_object"] = {names[k]: {"n_measurable": len(v), "mean_d": float(np.mean(v)), "median_d": float(np.median(v))}
                                           for k, v in acc.items() if len(v) >= 20}
                    rec["nulls"][kind][variant] = s
            out["samples"][tag] = rec
    (RES / "s8_null.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    for tag in out["samples"]:
        for kind in NULLS:
            for variant in VARIANTS:
                s = out["samples"][tag]["nulls"][kind][variant]
                print(tag, kind, variant, "n", s["n_measurable_p5_p50_p95"][1], "med D %.3f" % s["median_d_p2p5_p16_p50_p84_p97p5_p99"][2],
                      "pos %.2f" % s["frac_positive_mean"], "rej", s["wilcoxon_rejection_rate"], s["label"], s.get("observed_vs_null", ""))


# ----------------------------------------------------------------------------- B2: direct comparison
HALF_YEAR = 182.62


def visit_parity(mjd) -> np.ndarray:
    """0/1 label that alternates every half year (WISE scans the same sky in opposite directions)."""
    mjd = np.asarray(mjd, dtype=float)
    return (np.round((mjd - mjd[0]) / HALF_YEAR).astype(int) % 2).astype(float)


def boot_median(x, rng, n: int = 10000) -> list[float]:
    x = np.asarray(x, dtype=float)
    b = np.median(rng.choice(x, size=(n, x.size)), axis=1)
    return [float(v) for v in np.percentile(b, [2.5, 97.5])]


def sample_stats(rows: list[dict], key_ok: str, prefix: str, rng) -> dict:
    use = [r for r in rows if r.get(key_ok)]
    out = {"n": len(use)}
    if len(use) < 8:
        return out
    delta = np.array([r[f"{prefix}delta"] for r in use])
    sig = np.array([r[f"{prefix}sigma_delta"] for r in use])
    col = np.array([r[f"{prefix}d_colour"] for r in use])
    k = int(np.sum(delta > 0))
    out.update(n_delta_positive=k, frac_delta_positive=k / len(use),
               binomial_p_one_sided=float(binomtest(k, len(use), 0.5, alternative="greater").pvalue),
               median_d_colour_mag=float(np.median(col)), median_d_colour_ci95=boot_median(col, rng),
               median_frac_delta=float(np.median([r[f"{prefix}frac_delta"] for r in use])),
               wilcoxon_p_one_sided=float(wilcoxon(delta / sig, alternative="greater").pvalue),
               n_above_2sigma=int(np.sum(delta > 2 * sig)), n_below_minus_2sigma=int(np.sum(delta < -2 * sig)),
               median_sigma_colour_mag=float(np.median([r[f"{prefix}sigma_colour"] for r in use])))
    return out


def run_direct() -> None:
    rng = np.random.default_rng([config.SEED, SEED_TAG, 20])
    null = json.loads((RES / "s8_null.json").read_text()) if (RES / "s8_null.json").exists() else None
    out = {"note": "stage 8 B2: offset of W1 between the rising and the declining branch at the same W2 flux; "
                   "no reference visits, no shared peak, no blackbody", "samples": {}}
    all_rows = []
    for tag in ("v1", "v2"):
        objs = load_sample(tag)
        ut = json.loads((CKPT / f"untimely_lc_{tag}.json").read_text())
        rows = []
        parity_agree = []
        for o in objs:
            dr = hyst8.direct_offset(o["f1"], o["s1"], o["f2"], o["s2"], r=o["r"])
            par = visit_parity(o["mjd"])
            dp = hyst8.direct_offset(o["f1"], o["s1"], o["f2"], o["s2"], r=o["r"], parity=par)
            d0 = hyst8.direct_offset(o["f1"], o["s1"], o["f2"], o["s2"], r=0.0)
            dt = hyst8.direct_offset(o["f1"], o["s1"], o["f2"], o["s2"], r=o["r"], trend=(o["mjd"] - o["mjd"].mean()) / 365.25)
            row = {"sample": tag, "name": o["name"], "type": o.get("type", ""), "z": o["z"], "qa_flags": o["qa_flags"],
                   "n_visits": int(o["f1"].size), "r_error": o["r"], "d_hyst": o["d"], "sigma_d": o["sigma_d"],
                   "hyst_measurable": o["d"] is not None}
            for k in ("usable", "comparable", "n_rise", "n_fade", "n_rise_in", "n_fade_in", "lo", "hi", "delta", "sigma_delta",
                      "frac_delta", "d_colour", "sigma_colour", "slope", "delta_two_lines", "sigma_two_lines"):
                row[k] = dr.get(k)
            row.update(parity_delta=dp["delta"] if dp["usable"] else None, parity_sigma_delta=dp["sigma_delta"] if dp["usable"] else None,
                       parity_d_colour=dp["d_colour"] if dp["usable"] else None, parity_frac_delta=dp.get("frac_delta"),
                       parity_sigma_colour=dp["sigma_colour"] if dp["usable"] else None,
                       parity_gamma=dp["gamma"] if dp["usable"] else None, parity_sigma_gamma=dp["sigma_gamma"] if dp["usable"] else None,
                       parity_ok=bool(dr["comparable"] and dp["usable"] and np.isfinite(dp["gamma"])),
                       r0_delta=d0["delta"] if d0["usable"] else None,
                       trend_ok=bool(dr["comparable"] and dt["usable"] and np.isfinite(dt["gamma"])),
                       trend_delta=dt["delta"] if dt["usable"] else None, trend_sigma_delta=dt["sigma_delta"] if dt["usable"] else None,
                       trend_d_colour=dt["d_colour"] if dt["usable"] else None, trend_frac_delta=dt.get("frac_delta"),
                       trend_sigma_colour=dt["sigma_colour"] if dt["usable"] else None,
                       trend_w1_fraction_per_year=dt["gamma"] if dt["usable"] else None,
                       trend_sigma=dt["sigma_gamma"] if dt["usable"] else None,
                       years_rise_minus_fade=dt.get("mean_trend_rise_minus_fade"))
            # unTimely light curve of the same object
            u = ut.get(o["name"])
            row.update(ut_ok=False)
            if u:
                arr = {k: np.array([np.nan if x is None else x for x in u[k]], dtype=float)
                       for k in ("mjd", "w1_flux_jy", "w1_flux_err_jy", "w2_flux_jy", "w2_flux_err_jy")}
                ok = np.all([np.isfinite(arr[k]) for k in arr], axis=0) & (arr["w1_flux_jy"] > 0) & (arr["w2_flux_jy"] > 0)
                fwd = np.array(u["forward"], dtype=float)[ok]
                if ok.sum() >= 10:
                    order = np.argsort(arr["mjd"][ok])
                    a = {k: arr[k][ok][order] for k in arr}
                    fwd = fwd[order]
                    du = hyst8.direct_offset(a["w1_flux_jy"], a["w1_flux_err_jy"], a["w2_flux_jy"], a["w2_flux_err_jy"], r=0.0)
                    dup = hyst8.direct_offset(a["w1_flux_jy"], a["w1_flux_err_jy"], a["w2_flux_jy"], a["w2_flux_err_jy"], r=0.0,
                                              parity=fwd if np.all(fwd >= 0) else None)
                    row.update(ut_n_epochs=int(ok.sum()), ut_ok=bool(du["usable"]), ut_delta=du["delta"], ut_sigma_delta=du["sigma_delta"],
                               ut_d_colour=du["d_colour"], ut_frac_delta=du.get("frac_delta"), ut_sigma_colour=du["sigma_colour"],
                               ut_n_rise_in=du["n_rise_in"], ut_n_fade_in=du["n_fade_in"],
                               ut_forward_delta=dup["delta"] if dup["usable"] else None, ut_forward_gamma=dup["gamma"] if dup["usable"] else None)
                    # do the half-year parity labels of the pipeline visits follow the unTimely scan flag?
                    if np.all(fwd >= 0):
                        pv = visit_parity(a["mjd"])
                        parity_agree.append(max(np.mean(pv == fwd), np.mean(pv != fwd)))
            rows.append(row)
        all_rows += rows
        st = {"n_objects": len(rows), "n_usable": int(sum(bool(r["usable"]) for r in rows)),
              "direct": sample_stats(rows, "comparable", "", rng),
              "direct_with_scan_parity_term": sample_stats(rows, "parity_ok", "parity_", rng),
              "untimely": sample_stats(rows, "ut_ok", "ut_", rng),
              "not_in_plan_with_linear_time_term": sample_stats(rows, "trend_ok", "trend_", rng)}
        tr = [r for r in rows if r["trend_ok"]]
        st["not_in_plan_with_linear_time_term"].update(
            note="added after the plan: F1 = c + m F2 + delta [rise] + kappa (t - mean t); guards against a slow drift of the colour",
            median_w1_fraction_per_year=float(np.median([r["trend_w1_fraction_per_year"] for r in tr])) if tr else None,
            median_years_rise_minus_fade=float(np.median([r["years_rise_minus_fade"] for r in tr])) if tr else None,
            frac_rise_earlier_than_fade=float(np.mean([r["years_rise_minus_fade"] < 0 for r in tr])) if tr else None)
        comp = [r for r in rows if r["comparable"]]
        # two separate lines (auxiliary)
        tl = np.array([r["delta_two_lines"] for r in comp if r["delta_two_lines"] is not None and np.isfinite(r["delta_two_lines"])])
        st["two_lines"] = {"n": int(tl.size), "frac_positive": float(np.mean(tl > 0)) if tl.size else None}
        r0 = np.array([r["r0_delta"] for r in comp if r["r0_delta"] is not None])
        st["with_r_zero"] = {"n": int(r0.size), "frac_positive": float(np.mean(r0 > 0)) if r0.size else None}
        # agreement with D where both exist
        both = [r for r in comp if r["hyst_measurable"]]
        st["with_d"] = {"n_both": len(both), "n_same_sign": int(sum((r["delta"] > 0) == (r["d_hyst"] > 0) for r in both)),
                        "frac_delta_positive_among_d_measurable": float(np.mean([r["delta"] > 0 for r in both])) if both else None,
                        "spearman_d_vs_minus_colour": float(spearmanr([r["d_hyst"] for r in both], [-r["d_colour"] for r in both]).statistic) if len(both) >= 8 else None,
                        "n_d_measurable_not_comparable": int(sum(r["hyst_measurable"] and not r["comparable"] for r in rows))}
        # objects seen by both photometries
        bo = [r for r in rows if r["comparable"] and r["ut_ok"]]
        st["untimely_against_pipeline"] = {"n_both": len(bo), "n_same_sign": int(sum((r["delta"] > 0) == (r["ut_delta"] > 0) for r in bo))}
        gam = np.array([r["parity_gamma"] / r["parity_sigma_gamma"] for r in rows if r["parity_ok"]])
        d_plain = st["direct"].get("median_d_colour_mag")
        d_par = st["direct_with_scan_parity_term"].get("median_d_colour_mag")
        st["scan_parity"] = {"n": int(gam.size), "gamma_pull_mean": float(gam.mean()) if gam.size else None,
                             "gamma_pull_std": float(gam.std()) if gam.size else None,
                             "frac_abs_pull_above_2": float(np.mean(np.abs(gam) > 2)) if gam.size else None,
                             "median_d_colour_plain": d_plain, "median_d_colour_with_parity": d_par,
                             "relative_change": abs(d_par - d_plain) / abs(d_plain) if d_plain and d_par is not None else None,
                             "parity_follows_untimely_forward_flag_median_agreement": float(np.median(parity_agree)) if parity_agree else None}
        if st["scan_parity"]["relative_change"] is not None:
            st["scan_parity"]["label"] = "scan 방향과 무관" if st["scan_parity"]["relative_change"] <= 0.30 else "scan 방향에 기댐"
        u = st["untimely"]
        if u.get("n", 0) >= 8:
            u["label"] = "같은 방향" if u["frac_delta_positive"] > 0.5 and u["binomial_p_one_sided"] < 0.05 else "같은 방향 아님"
        # registered labels (section 3.2), using the null distribution of the sign fraction
        d = st["direct"]
        if null and d.get("n", 0) >= 8:
            nd = [null["samples"][tag]["nulls"][k]["primary"]["direct"] for k in NULLS]
            p99 = max(x["frac_delta_positive_p2p5_p50_p97p5_p99"][3] for x in nd)
            lo = min(x["frac_delta_positive_p2p5_p50_p97p5_p99"][0] for x in nd)
            hi = max(x["frac_delta_positive_p2p5_p50_p97p5_p99"][2] for x in nd)
            d["null_usable"] = bool(all(x["usable_under_null"] for x in nd))
            d["null_frac_positive_p99"], d["null_frac_positive_p2p5_p97p5"] = p99, [lo, hi]
            ci = d["median_d_colour_ci95"]
            if d["frac_delta_positive"] >= 0.65 and d["binomial_p_one_sided"] < 0.01 and ci[1] < 0 and d["frac_delta_positive"] > p99:
                d["label"] = "비대칭 지지"
            elif ci[0] <= 0 <= ci[1] and lo <= d["frac_delta_positive"] <= hi:
                d["label"] = "사라짐"
            else:
                d["label"] = "약해짐"
        out["samples"][tag] = st
    labels = [out["samples"][t]["direct"].get("label") for t in ("v1", "v2")]
    out["overall"] = ("지지" if all(x == "비대칭 지지" for x in labels) else "사라짐" if all(x == "사라짐" for x in labels) else "약해짐")
    cols = list(all_rows[0])
    for r in all_rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    with open(RES / "s8_direct_objects.csv", "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(cols)
        for r in all_rows:
            w.writerow(["" if r.get(k) is None or (isinstance(r.get(k), float) and not np.isfinite(r.get(k)))
                        else (f"{r[k]:.6g}" if isinstance(r[k], float) else (int(r[k]) if isinstance(r[k], (bool, np.bool_)) else r[k]))
                        for k in cols])
    (RES / "s8_direct.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    for tag in ("v1", "v2"):
        st = out["samples"][tag]
        print(tag, json.dumps({k: st[k] for k in ("direct", "untimely", "with_d", "two_lines", "scan_parity", "untimely_against_pipeline")}, ensure_ascii=False)[:1500])
    print("overall", out["overall"])


# ----------------------------------------------------------------------------- B3: lag model re-test
LAG_CONFIGS = {f"{drv}_scatter{sc:.2f}": {"driver": drv, "scatter_dex": sc}
               for drv in ("macleod", "macleod_sf_burke_tau") for sc in (0.0, 0.15, 0.30)}
RHO_OBSERVED = 0.129
_LAG: dict = {}


def lag_objects() -> list[dict]:
    import run_d_lum_power as power

    objs = power.load_objects()
    mbh = {}
    with open(config.TARGETS_DIR / "sample_v2.csv", newline="") as fh:
        for t in csv.DictReader(fh):
            mbh[t["name"]] = float(t["dr2_logmbh"]) if t["dr2_logmbh"].strip() else None
    known = [v for v in mbh.values() if v is not None]
    med = float(np.median(known))
    real = {r["name"] for r in sample_names("v2") if r["d"] is not None and r["log_lbol"] is not None}
    cond = power.load_conditions("exact", [o["name"] for o in objs])
    out = []
    for o in objs:
        if o["name"] not in cond:
            continue
        m = mbh.get(o["name"])
        out.append(dict(o, log_mbh=m if m is not None else med, mbh_known=m is not None, mjd=cond[o["name"]][0],
                        err=cond[o["name"]][1], real_measurable=o["name"] in real))
    return out


def lag_realisation(args) -> dict:
    import run_d_lum_power as power

    key, index = args
    cfg = LAG_CONFIGS[key]
    rng = np.random.default_rng([config.SEED, SEED_TAG, 30, list(LAG_CONFIGS).index(key), index])
    xs, ds, real = [], [], []
    for o in _LAG["objects"]:
        dt1, dt2 = power.lyu_lags(o["log_lbol"])
        f = 10.0 ** rng.normal(0.0, cfg["scatter_dex"]) if cfg["scatter_dex"] > 0 else 1.0
        lag = (1.0 + o["z"]) * dt1 * f
        extra = (1.0 + o["z"]) * (dt2 - dt1) * f
        sigma, tau = hyst8.macleod_drw(o["log_lbol"], o["log_mbh"])
        if cfg["driver"] == "macleod_sf_burke_tau":
            tau = hyst8.burke_tau(o["log_mbh"])
        tau_obs = tau * (1.0 + o["z"])
        toy = breath.ToyModel("boundary", 0.0, host_w1=float(rng.uniform(0.0, 0.8)), lag_days=lag, lag_w2_extra_days=extra)
        v = breath.synthetic_visits(rng, toy, o["mjd"], o["err"], sigma, tau_days=tau_obs, burn_in_days=3500.0)
        m = breath.measure(v["w1_flux"], v["w1_flux_err"], v["w2_flux"], v["w2_flux_err"], with_in_place=False)
        if m["hyst_measurable"]:
            xs.append(o["x"])
            ds.append(m["d_hyst"])
            real.append(o["real_measurable"])
    xs, ds, real = np.array(xs), np.array(ds), np.array(real, dtype=bool)
    out = {"n_measurable": int(xs.size), "median_d": float(np.median(ds)) if xs.size else None,
           "rho": float(spearmanr(xs, ds).statistic) if xs.size >= 10 else None, "n_matched": int(real.sum()),
           "rho_matched": float(spearmanr(xs[real], ds[real]).statistic) if real.sum() >= 10 else None}
    return out


def run_lag(workers: int) -> None:
    CKPT.mkdir(parents=True, exist_ok=True)
    _LAG["objects"] = lag_objects()
    objs = _LAG["objects"]
    drw = np.array([hyst8.macleod_drw(o["log_lbol"], o["log_mbh"]) for o in objs])
    out = {"note": "stage 8 B3: model 1 (inter-band lag from Lyu et al. 2019) with a driver that depends on mass and luminosity",
           "n_objects": len(objs), "n_mbh_known": int(sum(o["mbh_known"] for o in objs)),
           "n_real_measurable_with_lbol": int(sum(o["real_measurable"] for o in objs)), "rho_observed": RHO_OBSERVED,
           "driver": {"sigma_mag_p5_p50_p95": [float(x) for x in np.percentile(drw[:, 0], [5, 50, 95])],
                      "tau_rest_days_macleod_p5_p50_p95": [float(x) for x in np.percentile(drw[:, 1], [5, 50, 95])],
                      "tau_rest_days_burke_p5_p50_p95": [float(x) for x in np.percentile([hyst8.burke_tau(o["log_mbh"]) for o in objs], [5, 50, 95])],
                      "stage6b_toy": "sigma uniform 0.15-0.6 mag, tau 1000 d, 30% of objects with a ramp of 0.5-1.5 mag"},
           "configs": {}}
    t0 = time.time()
    with multiprocessing.get_context("fork").Pool(workers) as pool:
        for key in LAG_CONFIGS:
            path = CKPT / f"lag_{key}.json"
            if path.exists():
                rows = json.loads(path.read_text())
            else:
                rows = pool.map(lag_realisation, [(key, i) for i in range(N_NULL)], chunksize=4)
                path.write_text(json.dumps(rows) + "\n")
                print(f"{time.time() - t0:6.0f} s  {key}", flush=True)
            rho = np.array([r["rho"] for r in rows if r["rho"] is not None])
            rm = np.array([r["rho_matched"] for r in rows if r["rho_matched"] is not None])
            rec = {**LAG_CONFIGS[key], "n_realisations": len(rows), "n_with_rho": int(rho.size),
                   "n_measurable_p5_p50_p95": [float(x) for x in np.percentile([r["n_measurable"] for r in rows], [5, 50, 95])],
                   "median_d_p50": float(np.median([r["median_d"] for r in rows if r["median_d"] is not None])) if any(r["median_d"] is not None for r in rows) else None}
            if rho.size >= 20:
                rec.update(rho_p5_p50_p95=[float(x) for x in np.percentile(rho, [5, 50, 95])],
                           fraction_at_or_below_observed=float(np.mean(rho <= RHO_OBSERVED)),
                           outside=bool(RHO_OBSERVED < np.percentile(rho, 5)))
            if rm.size >= 20:
                rec.update(matched_n_with_rho=int(rm.size), matched_rho_p5_p50_p95=[float(x) for x in np.percentile(rm, [5, 50, 95])],
                           matched_fraction_at_or_below_observed=float(np.mean(rm <= RHO_OBSERVED)))
            out["configs"][key] = rec
    flags = [c.get("outside") for c in out["configs"].values()]
    out["verdict"] = ("판정불가 (mock의 측정 가능 수 부족)" if any(f is None for f in flags) else
                      "여전히 바깥" if all(flags) else "바깥이 아님" if not any(flags) else "가정에 달림")
    (RES / "s8_lag_model.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    for k, c in out["configs"].items():
        print(k, {a: c.get(a) for a in ("n_measurable_p5_p50_p95", "median_d_p50", "rho_p5_p50_p95", "fraction_at_or_below_observed", "outside")})
    print(out["verdict"])


# ----------------------------------------------------------------------------- B4: selection and reporting
def run_select() -> None:
    out = {"note": "stage 8 B4: descriptive"}
    # (1) sample_v2: measurable objects with and without L_bol
    rows = sample_names("v2")
    objs = {o["name"]: o for o in load_sample("v2")}
    meas = [r for r in rows if r["d"] is not None]

    def snr(name):
        o = objs[name]
        return float(np.median(np.abs(o["f2"] - np.median(o["f2"])) / o["s2"]))

    def group(g):
        d = np.array([r["d"] for r in g])
        t1 = sum(1 for r in g if r["type"] in ("Sy1", "Sy1.2", "Sy1.5", "Sy1.8"))
        t2 = sum(1 for r in g if r["type"] in ("Sy1.9", "Sy2", "Sy2 candidate"))
        return {"n": len(g), "median_d": float(np.median(d)), "d_p16_p84": [float(x) for x in np.percentile(d, [16, 84])],
                "median_z": float(np.median([r["z"] for r in g])), "n_type1": t1, "n_type2": t2, "n_other_type": len(g) - t1 - t2,
                "median_sigma_d": float(np.median([r["sigma_d"] for r in g])),
                "median_w2_variability_snr": float(np.median([snr(r["name"]) for r in g])),
                "median_n_visits": float(np.median([objs[r["name"]]["f1"].size for r in g]))}

    with_l, without = [r for r in meas if r["log_lbol"] is not None], [r for r in meas if r["log_lbol"] is None]
    out["v2_measurable_with_and_without_lbol"] = {
        "with_lbol": group(with_l), "without_lbol": group(without),
        "mann_whitney_p_two_sided": {
            "d": float(mannwhitneyu([r["d"] for r in with_l], [r["d"] for r in without], alternative="two-sided").pvalue),
            "z": float(mannwhitneyu([r["z"] for r in with_l], [r["z"] for r in without], alternative="two-sided").pvalue),
            "w2_variability_snr": float(mannwhitneyu([snr(r["name"]) for r in with_l], [snr(r["name"]) for r in without], alternative="two-sided").pvalue)},
        "types_of_objects_without_lbol": sorted({r["type"] for r in without})}
    # (2) sample_v1 types and the Q1 / Q2-Q3 object sets
    with open(RES / "s5b_objects.csv", newline="") as fh:
        v1 = [r for r in csv.DictReader(fh) if r["sample"] == "main"]
    groups: dict[str, int] = {}
    for r in v1:
        groups[r["type_group"]] = groups.get(r["type_group"], 0) + 1
    out["v1_types"] = {"n_main": len(v1), "by_type_group": groups,
                       "objects_outside_type_1_and_2": [{"name": r["name"], "counterpart": r["counterpart"], "type": r["type"],
                                                         "type_group": r["type_group"]} for r in v1 if r["type_group"] not in ("1", "2")]}
    q1 = [r for r in v1 if r["q1"] == "1"]
    br = [r for r in v1 if r["brightened"] == "1"]
    out["v1_q1_against_q2_q3"] = {
        "n_q1": len(q1), "n_brightened_for_q2_q3": len(br),
        "n_in_both": int(sum(r["brightened"] == "1" for r in q1)),
        "q1_not_in_q2": [r["name"] for r in q1 if r["brightened"] != "1"],
        "q2_not_in_q1": [r["name"] for r in br if r["q1"] != "1"],
        "n_flagged_a": int(sum(r["a_flagged"] == "1" for r in v1)), "n_flagged_b": int(sum(r["b_flagged"] == "1" for r in v1)),
        "columns": "q1 and brightened in results/s5b_objects.csv"}
    # (3) the quoted p values recomputed with the documented test
    d1 = np.array([float(r["sl_d_hyst"]) for r in v1 if r["sl_hyst_measurable"] == "1"])
    d2 = np.array([r["d"] for r in meas])
    out["p_values"] = {"test": "one-sided Wilcoxon signed-rank test of D against 0 (scipy.stats.wilcoxon, alternative='greater')",
                       "code": {"v1": "src/irecho/stage5b.py hysteresis_test (called by scripts/run_stage5b.py)",
                                "v2": "scripts/run_stage6b.py main (d_replication)"},
                       "v1": {"n": int(d1.size), "p": float(wilcoxon(d1, alternative="greater").pvalue), "n_positive": int(np.sum(d1 > 0))},
                       "v2": {"n": int(d2.size), "p": float(wilcoxon(d2, alternative="greater").pvalue), "n_positive": int(np.sum(d2 > 0))},
                       "smallest_possible_p_when_all_positive": {"v1": float(0.5 ** d1.size), "v2": float(0.5 ** d2.size)}}
    (RES / "s8_selection.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    print(json.dumps(out, indent=1, ensure_ascii=False)[:3500])


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else ""
    if stage == "null":
        run_null(int(sys.argv[2]) if len(sys.argv) > 2 else 2)
    elif stage == "direct":
        run_direct()
    elif stage == "lag":
        run_lag(int(sys.argv[2]) if len(sys.argv) > 2 else 2)
    elif stage == "select":
        run_select()
    elif stage == "time":
        objs = load_sample("v2")
        _OBJECTS[("v2", False)] = [dict(o, truth=prepare_truths(o, False)) for o in objs]
        t0 = time.time()
        r = one_realisation(("v2", "linear", "primary", 0))
        print("one realisation v2:", round(time.time() - t0, 1), "s", {k: v for k, v in r.items() if k != "per_obj"})
        print("truth available:", sum(o["truth"]["linear"] is not None for o in _OBJECTS[("v2", False)]),
              sum(o["truth"]["single_t"] is not None for o in _OBJECTS[("v2", False)]), "of", len(objs))
    else:
        raise SystemExit(__doc__)
