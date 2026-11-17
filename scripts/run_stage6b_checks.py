"""Stage 6B: checks added AFTER the registered result was seen (not part of
docs/PREDICTIONS_6.md). They do not change the verdict; they show how the
registered numbers move under choices the registration left open.

    python scripts/run_stage6b_checks.py

1. QA-flagged objects excluded; errors without the galaxy floor.
2. Type 1 and type 2 separately (grouping of stage 5b).
3. D against redshift and luminosity separately; bootstrap of the
   high-minus-low luminosity difference.
4. sample_v1 against sample_v2: D distributions, also inside the common
   luminosity range.
5. How much an error on log L_bol lowers the model-1 correlation (toy samples
   of results/s6b_power_final.json, cached realisations).
Output: results/s6b_checks.json.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import mannwhitneyu, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import run_d_lum_power as power  # noqa: E402
import run_stage6b as main6b  # noqa: E402
from irecho import breath, config  # noqa: E402
from irecho import stage5b as s5  # noqa: E402

N_PERM = 10000
LBOL_SCATTER_DEX = (0.0, 0.3, 0.5, 1.0)


def type_group(t: str) -> str:
    """Stage 5b grouping: Sy1-Sy1.8 -> '1'; Sy1.9, Sy2, Sy2 candidate -> '2'."""
    if t in ("Sy1", "Sy1.2", "Sy1.5", "Sy1.8"):
        return "1"
    if t in ("Sy1.9", "Sy2", "Sy2 candidate"):
        return "2"
    return "other"


def relation(rows: list[dict], final: dict, rng) -> dict:
    """Registered statistic and rule on a subset (rows need x and d)."""
    n = len(rows)
    out = {"n": n}
    if n < 8:
        return out
    x = np.array([r["x"] for r in rows])
    d = np.array([r["d"] for r in rows])
    rho = float(spearmanr(x, d).statistic)
    perm = np.array([spearmanr(x, rng.permutation(d)).statistic for _ in range(N_PERM)])
    q2 = float(np.percentile(perm, 95))
    q1 = main6b.interp_threshold(final["model1"]["by_n"], n, "rho_p5")
    out.update(median_d=float(np.median(d)), rho=rho, p_one_sided_permutation=float(np.mean(perm >= rho)),
               q2=q2, q1_model1=q1, verdict_by_the_registered_rule=main6b.verdict(rho, q1, q2, n))
    return out


def main() -> None:
    rng = np.random.default_rng([config.SEED, 65])
    final = json.loads((config.RESULTS_DIR / "s6b_power_final.json").read_text())["models"]
    with open(config.RESULTS_DIR / "s6b_objects.csv", newline="") as fh:
        obj = list(csv.DictReader(fh))
    meas = [dict(r, d=float(r["d_hyst"]), x=float(r["x"]) if r["x"] else float("nan"),
                 lum=float(r["log_lbol"]) if r["log_lbol"] else float("nan"), zz=float(r["z"]),
                 group=type_group(r["type"]))
            for r in obj if r["hyst_measurable"] == "1"]
    use = [r for r in meas if np.isfinite(r["x"])]
    out = {"note": "computed after the registered result was seen; not part of PREDICTIONS_6.md"}

    # 1. flags and floor
    clean_all = [r for r in meas if not r["qa_flags"]]
    clean = [r for r in use if not r["qa_flags"]]
    out["registered"] = relation(use, final, rng)
    out["without_flagged_objects"] = dict(
        relation(clean, final, rng), n_measurable_any_lbol=len(clean_all),
        median_d_any_lbol=float(np.median([r["d"] for r in clean_all])))
    nofloor = []
    for r in obj:
        v = s5.read_frozen(r["name"], with_floor=False, frozen_dir=main6b.FROZEN_V2)
        m = breath.measure(v["w1_flux_jy"], v["w1_flux_err_jy"], v["w2_flux_jy"], v["w2_flux_err_jy"],
                           with_in_place=False)
        if m["hyst_measurable"]:
            nofloor.append({"d": float(m["d_hyst"]), "x": float(r["x"]) if r["x"] else float("nan")})
    out["errors_without_floor"] = dict(
        relation([r for r in nofloor if np.isfinite(r["x"])], final, rng),
        n_measurable_any_lbol=len(nofloor), median_d_any_lbol=float(np.median([r["d"] for r in nofloor])))

    # 2. by type
    out["by_type"] = {}
    for g in ("1", "2"):
        sub_all = [r for r in meas if r["group"] == g]
        out["by_type"][g] = dict(relation([r for r in use if r["group"] == g], final, rng),
                                 n_measurable_any_lbol=len(sub_all),
                                 median_d_any_lbol=float(np.median([r["d"] for r in sub_all])))
    out["n_sample_by_type"] = {g: sum(1 for r in obj if type_group(r["type"]) == g) for g in ("1", "2", "other")}

    # 3. separate predictors, high-minus-low difference
    d = np.array([r["d"] for r in use])
    lum = np.array([r["lum"] for r in use])
    zz = np.array([r["zz"] for r in use])
    hi = lum >= np.median(lum)
    diff = float(np.median(d[hi]) - np.median(d[~hi]))
    boot = []
    for _ in range(10000):
        idx = rng.integers(0, d.size, d.size)
        h = hi[idx]
        if 2 < h.sum() < d.size - 2:
            boot.append(np.median(d[idx][h]) - np.median(d[idx][~h]))
    perm_diff = []
    for _ in range(N_PERM):
        p = rng.permutation(d)
        perm_diff.append(np.median(p[hi]) - np.median(p[~hi]))
    out["separate_predictors"] = {
        "rho_d_log_lbol": float(spearmanr(lum, d).statistic),
        "rho_d_redshift": float(spearmanr(zz, d).statistic),
        "rho_log_lbol_redshift": float(spearmanr(lum, zz).statistic),
        "high_minus_low_median_d": diff,
        "high_minus_low_ci95": [float(v) for v in np.percentile(boot, [2.5, 97.5])],
        "high_minus_low_p_one_sided_permutation": float(np.mean(np.array(perm_diff) >= diff)),
        "model1_expectation": final["model1"]["d_median_high_minus_low_luminosity"],
        "model1_half_expectation": final["model1_half"]["d_median_high_minus_low_luminosity"],
        "model2_expectation": final["model2"]["d_median_high_minus_low_luminosity"]}
    out["measurable_objects_without_lbol"] = {
        "n": len(meas) - len(use),
        "median_d": float(np.median([r["d"] for r in meas if not np.isfinite(r["x"])]))}

    # 4. sample_v1 against sample_v2
    with open(config.RESULTS_DIR / "s5b_objects.csv", newline="") as fh:
        v1 = [r for r in csv.DictReader(fh) if r["sample"] == "main" and r["sl_hyst_measurable"] == "1"]
    d1 = np.array([float(r["sl_d_hyst"]) for r in v1])
    d2 = np.array([r["d"] for r in meas])
    bd = [np.median(rng.choice(d2, d2.size)) - np.median(rng.choice(d1, d1.size)) for _ in range(10000)]
    lum1 = np.array([float(r["loglbol"]) if r["loglbol"] else np.nan for r in v1])
    lo, hi_edge = float(lum.min()), float(np.nanmax(lum1))
    in1 = d1[(lum1 >= lo) & (lum1 <= hi_edge)]
    in2 = d[(lum >= lo) & (lum <= hi_edge)]
    out["v1_against_v2"] = {
        "n_v1": int(d1.size), "n_v2": int(d2.size),
        "measurable_fraction_v1": [int(d1.size), 225], "measurable_fraction_v2": [int(d2.size), len(obj)],
        "median_d_v1": float(np.median(d1)), "median_d_v2": float(np.median(d2)),
        "median_difference_v2_minus_v1": float(np.median(d2) - np.median(d1)),
        "median_difference_ci95": [float(v) for v in np.percentile(bd, [2.5, 97.5])],
        "mann_whitney_p_two_sided": float(mannwhitneyu(d2, d1, alternative="two-sided").pvalue),
        "common_log_lbol_range": [lo, hi_edge],
        "common_range_n_v1": int(in1.size), "common_range_n_v2": int(in2.size),
        "common_range_median_d_v1": float(np.median(in1)) if in1.size else None,
        "common_range_median_d_v2": float(np.median(in2)) if in2.size else None,
        "common_range_mann_whitney_p_two_sided":
            float(mannwhitneyu(in2, in1, alternative="two-sided").pvalue) if min(in1.size, in2.size) >= 5 else None,
        "toy_median_d": {k: final[k]["d_median"] for k in final}}

    # 5. error on log L_bol in the model-1 toy samples
    n = len(use)
    out["lbol_error_in_model1"] = {"n": n, "observed_rho": out["registered"]["rho"], "cases": []}
    for key in ("model1", "model1_half"):
        reals = json.loads((config.CACHE_DIR / "s6b" / "power_final" / f"{key}.json").read_text())
        pooled = [r for rows in reals for r in rows if r["ok"]]
        px = np.array([r["x"] for r in pooled])
        pd_ = np.array([r["d"] for r in pooled])
        for sig in LBOL_SCATTER_DEX:
            rhos = []
            for _ in range(3000):
                idx = rng.choice(px.size, size=n, replace=False)
                rhos.append(spearmanr(px[idx] + 0.5 * sig * rng.standard_normal(n), pd_[idx]).statistic)
            rhos = np.array(rhos)
            out["lbol_error_in_model1"]["cases"].append({
                "model": key, "sigma_log_lbol_dex": sig,
                "rho_p5": float(np.percentile(rhos, 5)), "rho_p50": float(np.percentile(rhos, 50)),
                "fraction_at_or_below_observed": float(np.mean(rhos <= out["registered"]["rho"]))})
    (config.RESULTS_DIR / "s6b_checks.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    print(json.dumps(out, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
