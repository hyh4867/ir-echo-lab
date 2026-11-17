"""Stage 5b: run the pre-registered analysis on the frozen sample_v1 light curves.

    python scripts/run_stage5b.py

Before anything else the frozen files are checked against data/freeze_v1.json.
Rules: docs/PREDICTIONS_5b.md. Outputs:

    results/s5b_objects.csv           one row per object (errors with the floor)
    results/s5b_objects_nofloor.csv   the same without the floor (reference)
    results/s5b_summary.json          verdict tables
    data/cache/s5/ztf_priority.txt    brightened objects first (for fetch_ztf.py --order)
"""

from __future__ import annotations

import csv
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from irecho import breath, config, irsa, sample_io  # noqa: E402
from irecho import stage5b as s5  # noqa: E402

RES = config.RESULTS_DIR
FREEZE_PATH = config.PROJECT_ROOT / "data" / "freeze_v1.json"


def check_freeze() -> dict:
    freeze = json.loads(FREEZE_PATH.read_text())
    bad = [n for n, rec in freeze["light_curves"].items()
           if irsa.sha256_file(s5.FROZEN_DIR / f"{n}.csv") != rec["sha256"]]
    if bad:
        raise SystemExit(f"frozen files differ from data/freeze_v1.json: {bad[:5]}")
    return {"n_files": len(freeze["light_curves"]), "n_mismatch": 0,
            "freeze_code_commit": freeze["code_commit"], "frozen_utc": freeze["frozen_utc"]}


def _num(text) -> float:
    try:
        return float(text)
    except (TypeError, ValueError):
        return float("nan")


def load_context() -> dict:
    with open(config.TARGETS_DIR / "sample_v1_props.csv", newline="") as fh:
        props = {r["name"]: r for r in csv.DictReader(fh)}
    qa = json.loads((RES / "qa_targets.json").read_text())
    floors = json.loads((RES / "qa_error_floor.json").read_text())["floors_mmag"]
    has = {}
    obj_path = sample_io.S5_DIR / "ztf" / "objects_dr24.csv"
    targets = sample_io.load_sample()
    if obj_path.exists():
        with open(obj_path, newline="") as fh:
            for r in csv.DictReader(fh):
                has.setdefault(targets[int(r["pid"])]["name"], {})[r["filtercode"][1]] = True
    completeness = {}
    path = RES / "s5b_q1_completeness.csv"
    if path.exists():
        with open(path, newline="") as fh:
            for r in csv.DictReader(fh):
                completeness[(float(r["floor_mmag"]), float(r["amp_mag"]), float(r["rise_yr"]))] = \
                    float(r["completeness_q1"])
    val_path = RES / "s5b_validation.json"
    d_min = json.loads(val_path.read_text())["d_min"] if val_path.exists() else None
    return {"targets": targets, "props": props, "qa": qa, "floors": floors, "ztf_has": has,
            "completeness": completeness, "d_min": d_min}


def analyse(ctx: dict, with_floor: bool) -> list[dict]:
    rows = []
    for t in ctx["targets"]:
        name = t["name"]
        p = ctx["props"][name]
        v = s5.read_frozen(name, with_floor)
        row = {"name": name, "counterpart": t["counterpart"], "role": t["role"],
               "seen_before": t["seen_before"],
               "sample": "main" if (t["role"] == "science" and not t["seen_before"]) else "separate",
               "qa_flags": ";".join(ctx["qa"][name]["flags"]), "type": p["type"],
               "type_group": p["type_group"], "z": t["z"], "dec": t["dec"],
               "logedd": p["dr2_logedd"], "loglbol": p["dr2_loglbol"], "logmbh": p["dr2_logmbh"],
               "lognh": p["dr1_lognh"], "above_threshold": p["above_threshold"],
               "above_threshold_nh": p["above_threshold_nh"],
               "w1_floor_mmag": ctx["floors"][name]["w1"]}
        det = s5.run_detectors(v)
        row.update(det)
        col = (s5.colour_change(v, det["start_index"]) if det["brightened"]
               else {"n_base": "", "delta_color": np.nan, "delta_color_err": np.nan,
                     "dm_w1": np.nan, "dm_w2": np.nan})
        row.update(col)
        row.update(s5.ztf_compare(v, name, ctx["ztf_has"].get(name, {}), float(t["dec"])))
        m = breath.measure(v["w1_flux_jy"], v["w1_flux_err_jy"], v["w2_flux_jy"],
                           v["w2_flux_err_jy"])
        row.update({f"sl_{k}": val for k, val in m.items() if k != "n_visits"})
        rows.append(row)
    return rows


def _f(x) -> float | None:
    return None if x is None or not np.isfinite(x) else float(x)


def summarise(rows: list[dict], ctx: dict, tag: int) -> dict:
    n = len(rows)
    out: dict = {"n": n}
    # ---- Q1 ------------------------------------------------------------------------
    q1 = [r for r in rows if r["q1"]]
    lo, hi = s5.wilson(len(q1), n)
    q = {"n_slow_rise": len(q1), "fraction": len(q1) / n, "ci95": [lo, hi],
         "verdict": s5.verdict_q1(lo, hi), "objects": [r["name"] for r in q1]}
    for key, cond in (("b_flag_0p1", lambda r: r["b_flagged"]),
                      ("b_flag_w2_0p2", lambda r: r["q1_w2_only"]),
                      ("a_flag", lambda r: r["a_flagged"]),
                      ("a_or_b", lambda r: r["brightened"])):
        k = sum(1 for r in rows if cond(r))
        q[key] = {"n": k, "fraction": k / n, "ci95": list(s5.wilson(k, n))}
    if ctx["completeness"]:
        grid = sorted(ctx["completeness"])
        floors = sorted({g[0] for g in grid})
        amps = sorted({g[1] for g in grid})
        rises = sorted({g[2] for g in grid})
        weights, detail = [], []
        for r in q1:
            key = (s5.nearest(floors, r["w1_floor_mmag"]),
                   s5.nearest(amps, min(r["b_w1_dm_fit"], r["b_w2_dm_fit"]), log=True),
                   s5.nearest(rises, r["b_rise_years"]))
            c = ctx["completeness"][key]
            weights.append(1.0 / max(c, s5.COMPLETENESS_FLOOR))
            detail.append({"name": r["name"], "floor_mmag": key[0], "amp_mag": key[1],
                           "rise_yr": key[2], "completeness": c})
        q["completeness_corrected"] = s5.corrected_fraction(weights, n, tag)
        q["completeness_detail"] = detail
    out["q1"] = q
    # ---- Q2 ------------------------------------------------------------------------
    br = [r for r in rows if r["brightened"] and np.isfinite(r["delta_color"])]
    dc = np.array([r["delta_color"] for r in br])
    red = int((dc > 0).sum())
    verdict, lo, hi = s5.verdict_majority(red, len(br))
    sig = np.array([r["delta_color"] / r["delta_color_err"] for r in br])
    q2 = {"n_brightened": len(br), "n_redder": red,
          "fraction_redder": red / len(br) if br else None, "ci95": [lo, hi], "verdict": verdict,
          "delta_color_p25_p50_p75": [float(x) for x in np.percentile(dc, [25, 50, 75])] if br else None,
          "n_redder_2sigma": int((sig >= 2).sum()), "n_bluer_2sigma": int((sig <= -2).sum()),
          "n_no_change_2sigma": int((np.abs(sig) < 2).sum()), "by_type": {}}
    for g in ("1", "2"):
        sub = [r for r in br if r["type_group"] == g]
        q2["by_type"][g] = {"n": len(sub), "n_redder": sum(1 for r in sub if r["delta_color"] > 0),
                            "median_delta_color":
                                _f(np.median([r["delta_color"] for r in sub])) if sub else None}
    t1, t2 = q2["by_type"]["1"], q2["by_type"]["2"]
    q2["type_fisher_p"] = s5.fisher_2x2(t1["n_redder"], t1["n"] - t1["n_redder"],
                                        t2["n_redder"], t2["n"] - t2["n_redder"])
    out["q2"] = q2
    # ---- Q3 ------------------------------------------------------------------------
    use = [r for r in rows if r["brightened"] and r["ztf_usable"]]
    acc = sum(1 for r in use if r["ztf_verdict"] == "accompanies")
    verdict, lo, hi = s5.verdict_majority(acc, len(use))
    q3 = {"n_brightened": sum(1 for r in rows if r["brightened"]), "n_usable": len(use),
          "n_accompanies": acc, "fraction": acc / len(use) if use else None, "ci95": [lo, hi],
          "verdict": verdict, "verdict_counts": dict(Counter(r["ztf_verdict"] for r in use)),
          "by_type": {}, "all_objects_with_ir_change": {}}
    for g in ("1", "2"):
        sub = [r for r in use if r["type_group"] == g]
        q3["by_type"][g] = {"n": len(sub), **dict(Counter(r["ztf_verdict"] for r in sub))}
    every = [r for r in rows if r["ztf_usable"]]
    q3["all_objects_with_ir_change"] = {
        "n": len(every), **dict(Counter(r["ztf_verdict"] for r in every)),
        "by_type": {g: {"n": sum(1 for r in every if r["type_group"] == g),
                        **dict(Counter(r["ztf_verdict"] for r in every if r["type_group"] == g))}
                    for g in ("1", "2")}}
    out["q3"] = q3
    # ---- H-breath ------------------------------------------------------------------
    ok = [r for r in rows if r["sl_status"] == "ok"]
    meas = [r for r in ok if r["sl_measurable"]]

    def group_stats(key_fn, label_true: str, label_false: str, seed: int) -> dict:
        res = {}
        above = [r for r in meas if key_fn(r) is True]
        below = [r for r in meas if key_fn(r) is False]
        res["test"] = s5.compare_groups([r["sl_r"] for r in above], [r["sl_r"] for r in below],
                                        tag * 10 + seed)
        for label, grp, allg in ((label_true, above, [r for r in ok if key_fn(r) is True]),
                                 (label_false, below, [r for r in ok if key_fn(r) is False])):
            wm, we = s5.weighted_mean([r["sl_beta"] for r in allg], [r["sl_sigma_beta"] for r in allg])
            bm, _ = s5.weighted_mean([r["sl_beta_in_place_min"] for r in allg],
                                     [r["sl_sigma_beta"] for r in allg])
            res[label] = {"n_measurable": len(grp),
                          "classes": dict(Counter(r["sl_class"] for r in grp)),
                          "median_beta_measurable":
                              _f(np.median([r["sl_beta"] for r in grp])) if grp else None,
                          "n_with_slope": len(allg), "stack_beta": _f(wm), "stack_beta_err": _f(we),
                          "stack_beta_in_place_min": _f(bm)}
        ca, cb = res[label_true]["classes"], res[label_false]["classes"]
        res["class_fisher_p"] = s5.fisher_2x2(ca.get("boundary", 0), ca.get("in place", 0),
                                              cb.get("boundary", 0), cb.get("in place", 0))
        return res

    def flag(col):
        return lambda r: (True if r[col] in (1, "1") else False if r[col] in (0, "0") else None)

    def split(col, cut):
        def fn(r):
            x = _num(r[col])
            return None if not np.isfinite(x) else bool(x >= cut)
        return fn

    lbol = [_num(r["loglbol"]) for r in ok if np.isfinite(_num(r["loglbol"]))]
    hb = {"n_with_slope": len(ok), "n_measurable": len(meas),
          "n_measurable_with_logedd": sum(1 for r in meas if r["above_threshold"] in (0, 1, "0", "1")),
          "measurable_objects": [
              {k: (r[k] if isinstance(r[k], str) else _f(r[k])) for k in
               ("name", "counterpart", "type", "logedd", "lognh", "sl_beta", "sl_sigma_beta",
                "sl_beta_in_place_min", "sl_r", "sl_s", "sl_t_pivot_k", "sl_class", "sl_n_used",
                "sl_max_excess_ratio")} for r in meas],
          "all_classes_measurable": dict(Counter(r["sl_class"] for r in meas)),
          "primary": group_stats(flag("above_threshold"), "above", "below", 1),
          "threshold_10^-1.5": group_stats(split("logedd", -1.5), "above", "below", 2),
          "threshold_nh_dependent": group_stats(flag("above_threshold_nh"), "above", "below", 3),
          "type": group_stats(lambda r: {"1": True, "2": False}.get(r["type_group"]),
                              "type 1", "type 2", 4),
          "nh": group_stats(split("lognh", 22.0), "logNH >= 22", "logNH < 22", 5),
          "luminosity": group_stats(split("loglbol", float(np.median(lbol)) if lbol else 0.0),
                                    "logLbol >= median", "logLbol < median", 6),
          "loglbol_median": _f(np.median(lbol)) if lbol else None}
    wm, we = s5.weighted_mean([r["sl_beta"] for r in ok], [r["sl_sigma_beta"] for r in ok])
    bm, _ = s5.weighted_mean([r["sl_beta_in_place_min"] for r in ok],
                             [r["sl_sigma_beta"] for r in ok])
    hb["stack_all"] = {"n": len(ok), "beta": _f(wm), "beta_err": _f(we), "beta_in_place_min": _f(bm)}
    hb["verdict"] = hb["primary"]["test"]["verdict"]
    out["h_breath"] = hb
    # ---- hysteresis ----------------------------------------------------------------
    hm = [r for r in ok if r["sl_hyst_measurable"]]
    hy = s5.hysteresis_test([r["sl_d_hyst"] for r in hm], ctx["d_min"], tag)
    hy.update(n_above_2sigma=sum(1 for r in hm if r["sl_d_hyst"] > 2 * r["sl_sigma_d"]),
              n_below_minus_2sigma=sum(1 for r in hm if r["sl_d_hyst"] < -2 * r["sl_sigma_d"]),
              d_p16_p50_p84=[float(x) for x in np.percentile([r["sl_d_hyst"] for r in hm],
                                                             [16, 50, 84])] if hm else None,
              objects=[{"name": r["name"], "d": _f(r["sl_d_hyst"]), "sigma_d": _f(r["sl_sigma_d"]),
                        "beta_rise": _f(r["sl_beta_rise"]), "beta_decline": _f(r["sl_beta_decline"])}
                       for r in hm])
    out["hysteresis"] = hy
    return out


def final_verdicts(all_s: dict, clean_s: dict) -> dict:
    out = {}
    for key in ("q1", "q2", "q3", "h_breath", "hysteresis"):
        a, c = all_s[key]["verdict"], clean_s[key]["verdict"]
        # "판정불가" and "판정불가 (표본 부족)" are the same verdict with a reason attached.
        base_a, base_c = a.split(" (")[0], c.split(" (")[0]
        out[key] = {"main_all": a, "main_clean": c,
                    "final": (a if a == c else base_a) if base_a == base_c else "갈림"}
    return out


def ztf_coverage(rows: list[dict]) -> dict:
    table = {}
    for band in ("g", "r"):
        table[band] = {g: dict(Counter(r[f"ztf_{band}_status"] for r in rows
                                       if g == "all" or r["type_group"] == g))
                       for g in ("all", "1", "2")}
    return table


def write_rows(path: Path, rows: list[dict]) -> None:
    def fmt(x):
        if isinstance(x, (bool, np.bool_)):
            return int(x)
        if isinstance(x, (float, np.floating)):
            return "" if not np.isfinite(x) else f"{x:.6g}"
        return "" if x is None else x
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(list(rows[0]))
        for r in rows:
            writer.writerow([fmt(r[k]) for k in rows[0]])


def main() -> None:
    freeze = check_freeze()
    ctx = load_context()
    summary = {"freeze_check": freeze, "d_min": ctx["d_min"],
               "settings": {k: getattr(breath.SETTINGS, k)
                            for k in breath.SETTINGS.__dataclass_fields__}}
    priority: list[str] = []
    for mode, with_floor in (("floor", True), ("nofloor", False)):
        rows = analyse(ctx, with_floor)
        write_rows(RES / ("s5b_objects.csv" if with_floor else "s5b_objects_nofloor.csv"), rows)
        main_all = [r for r in rows if r["sample"] == "main"]
        main_clean = [r for r in main_all if not r["qa_flags"]]
        s_all, s_clean = summarise(main_all, ctx, 1), summarise(main_clean, ctx, 2)
        summary[mode] = {"main_all": s_all, "main_clean": s_clean,
                         "verdicts": final_verdicts(s_all, s_clean),
                         "ztf_coverage_all_231": ztf_coverage(rows)}
        priority.extend(r["name"] for r in rows if r["brightened"] and r["name"] not in priority)
    (sample_io.S5_DIR / "ztf_priority.txt").write_text("\n".join(priority) + "\n")
    (RES / "s5b_summary.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False) + "\n")
    print(json.dumps({m: summary[m]["verdicts"] for m in ("floor", "nofloor")}, indent=1,
                     ensure_ascii=False))


if __name__ == "__main__":
    main()
