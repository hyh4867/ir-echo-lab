"""Stage 6B (B3): the hysteresis statistic D in the independent sample and its
relation to luminosity. Rules: docs/PREDICTIONS_6.md (sections 6-8).

    python scripts/run_stage6b.py

The frozen sample_v2 files are checked against data/freeze_v2.json first.
Thresholds: q2 from a permutation of the observed D values; q1 from the
model-1 toy samples made with the real cadence and errors
(results/s6b_power_final.json, committed before D was computed).
Output: results/s6b_objects.csv, results/s6b_summary.json.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr, wilcoxon

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import run_d_lum_power as power  # noqa: E402
from irecho import breath, config, irsa  # noqa: E402
from irecho import stage5b as s5  # noqa: E402

FROZEN_V2 = config.PROJECT_ROOT / "data" / "frozen" / "sample_v2"
FREEZE_V2 = config.PROJECT_ROOT / "data" / "freeze_v2.json"
MIN_OBJECTS = 20
N_PERM = 10000


def check_freeze() -> dict:
    freeze = json.loads(FREEZE_V2.read_text())
    bad = [n for n, rec in freeze["light_curves"].items()
           if irsa.sha256_file(FROZEN_V2 / f"{n}.csv") != rec["sha256"]]
    if bad:
        raise SystemExit(f"frozen files differ from data/freeze_v2.json: {bad[:5]}")
    return {"n_files": len(freeze["light_curves"]), "n_mismatch": 0,
            "freeze_code_commit": freeze["code_commit"], "frozen_utc": freeze["frozen_utc"]}


def interp_threshold(by_n: dict, n: int, key: str) -> float:
    grid = sorted(int(k) for k in by_n)
    return float(np.interp(n, grid, [by_n[str(g)][key] for g in grid]))


def verdict(rho: float, q1: float, q2: float, n: int) -> str:
    """PREDICTIONS_6.md section 7."""
    if n < MIN_OBJECTS:
        return "판정불가 (표본 부족)"
    if rho <= q2 and rho < q1:
        return "사용자 예측 지지"
    if rho > q2 and rho >= q1:
        return "사용자 예측 기각"
    return "판정불가"


def boot_median(x, rng) -> list[float]:
    x = np.asarray(x, dtype=float)
    b = np.median(rng.choice(x, size=(10000, x.size)), axis=1)
    return [float(v) for v in np.percentile(b, [2.5, 97.5])]


def main() -> None:
    freeze = check_freeze()
    objects = {o["name"]: o for o in power.load_objects()}
    with open(config.TARGETS_DIR / "sample_v2.csv", newline="") as fh:
        targets = list(csv.DictReader(fh))
    qa = json.loads((config.RESULTS_DIR / "qa2_targets.json").read_text())
    rows = []
    for t in targets:
        name = t["name"]
        if qa[name]["status"] != "ok":
            continue
        v = s5.read_frozen(name, frozen_dir=FROZEN_V2)
        m = breath.measure(v["w1_flux_jy"], v["w1_flux_err_jy"], v["w2_flux_jy"], v["w2_flux_err_jy"],
                           with_in_place=False)
        o = objects.get(name)
        rows.append({"name": name, "counterpart": t["counterpart"], "type": t["dr2_type"] or t["bat_type"],
                     "z": float(t["z"]), "log_lbol": o["log_lbol"] if o else float("nan"),
                     "x": o["x"] if o else float("nan"), "qa_flags": ";".join(qa[name]["flags"]),
                     "n_visits": m["n_visits"], "n_used": m["n_used"], "status": m["status"],
                     "n_rising": m["n_rising"], "n_declining": m["n_declining"],
                     "beta_rise": m["beta_rise"], "beta_decline": m["beta_decline"],
                     "d_hyst": m["d_hyst"], "sigma_d": m["sigma_d"],
                     "hyst_measurable": bool(m["hyst_measurable"])})
    with open(config.RESULTS_DIR / "s6b_objects.csv", "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(list(rows[0]))
        for r in rows:
            w.writerow(["" if isinstance(x, float) and not np.isfinite(x)
                        else (f"{x:.6g}" if isinstance(x, float) else (int(x) if isinstance(x, bool) else x))
                        for x in r.values()])
    rng = np.random.default_rng([config.SEED, 64])
    meas = [r for r in rows if r["hyst_measurable"]]
    use = [r for r in meas if np.isfinite(r["x"])]
    d_all = np.array([r["d_hyst"] for r in meas])
    out = {"freeze_check": freeze, "n_sample": len(targets), "n_with_light_curve": len(rows),
           "n_measurable": len(meas), "n_measurable_with_lbol": len(use)}
    # D replication (B3, first item)
    rep = {"n": int(d_all.size), "median_d": float(np.median(d_all)) if d_all.size else None}
    if d_all.size >= 8:
        rep.update(median_d_ci95=boot_median(d_all, rng),
                   wilcoxon_p_one_sided=float(wilcoxon(d_all, alternative="greater").pvalue),
                   n_above_2sigma=int(sum(r["d_hyst"] > 2 * r["sigma_d"] for r in meas)),
                   n_below_minus_2sigma=int(sum(r["d_hyst"] < -2 * r["sigma_d"] for r in meas)))
        rep["sign_reproduced"] = bool(rep["wilcoxon_p_one_sided"] < 0.05 and rep["median_d"] > 0)
    rep["stage5b_reference"] = {"median_d": 0.228, "ci95": [0.200, 0.248], "n": 71}
    out["d_replication"] = rep
    # D against luminosity
    final = json.loads((config.RESULTS_DIR / "s6b_power_final.json").read_text())["models"]
    n = len(use)
    rel = {"n": n}
    if n >= 8:
        x = np.array([r["x"] for r in use])
        d = np.array([r["d_hyst"] for r in use])
        lum = np.array([r["log_lbol"] for r in use])
        rho = float(spearmanr(x, d).statistic)
        perm = np.array([spearmanr(x, rng.permutation(d)).statistic for _ in range(N_PERM)])
        q2 = float(np.percentile(perm, 95))
        boot = []
        for _ in range(5000):
            idx = rng.integers(0, n, n)
            if np.unique(x[idx]).size > 2:
                boot.append(spearmanr(x[idx], d[idx]).statistic)
        hi = lum >= np.median(lum)
        rel.update(rho=rho, rho_ci95=[float(v) for v in np.percentile(boot, [2.5, 97.5])],
                   p_one_sided_permutation=float(np.mean(perm >= rho)), q2=q2,
                   rho_with_log_lbol=float(spearmanr(lum, d).statistic),
                   d_median_high_luminosity=float(np.median(d[hi])),
                   d_median_low_luminosity=float(np.median(d[~hi])),
                   d_median_high_minus_low=float(np.median(d[hi]) - np.median(d[~hi])),
                   log_lbol_median=float(np.median(lum)),
                   log_lbol_range=[float(lum.min()), float(lum.max())])
        rel["q1"] = {k: interp_threshold(final[k]["by_n"], n, "rho_p5")
                     for k in ("model1", "model1_half", "model1_double")}
        rel["model1_rho_median"] = interp_threshold(final["model1"]["by_n"], n, "rho_p50")
        rel["power_model1"] = interp_threshold(final["model1"]["by_n"], n, "power_rho_above_q2")
        rel["verdict"] = verdict(rho, rel["q1"]["model1"], q2, n)
        rel["verdict_if_lag_half"] = verdict(rho, rel["q1"]["model1_half"], q2, n)
        rel["verdict_if_lag_double"] = verdict(rho, rel["q1"]["model1_double"], q2, n)
    else:
        rel["verdict"] = "판정불가 (표본 부족)"
    out["d_vs_luminosity"] = rel
    # reference only: both samples together
    with open(config.RESULTS_DIR / "s5b_objects.csv", newline="") as fh:
        v1 = [r for r in csv.DictReader(fh) if r["sample"] == "main" and r["sl_hyst_measurable"] == "1"
              and r["loglbol"]]
    x1 = np.array([np.log10(1 + float(r["z"])) + 0.5 * float(r["loglbol"]) for r in v1])
    d1 = np.array([float(r["sl_d_hyst"]) for r in v1])
    comb = {"n_v1": int(x1.size), "rho_v1_already_seen": float(spearmanr(x1, d1).statistic)}
    if n >= 8:
        xc = np.concatenate([x1, [r["x"] for r in use]])
        dc = np.concatenate([d1, [r["d_hyst"] for r in use]])
        comb.update(n_combined=int(xc.size), rho_combined=float(spearmanr(xc, dc).statistic))
    out["combined_reference_only"] = comb
    (config.RESULTS_DIR / "s6b_summary.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    print(json.dumps({k: out[k] for k in ("n_sample", "n_measurable", "n_measurable_with_lbol")}))
    print(json.dumps(out["d_replication"], ensure_ascii=False))
    print(json.dumps(out["d_vs_luminosity"], ensure_ascii=False))
    print(json.dumps(out["combined_reference_only"]))


if __name__ == "__main__":
    main()
