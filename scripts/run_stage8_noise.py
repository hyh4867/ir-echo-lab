"""Stage 8 (B1 input): correlation of the W1 and W2 visit errors, and the unTimely
light curves used in B2. Needs the cached field tables.

    python scripts/run_stage8_noise.py                 (sample_v1)
    IRECHO_SAMPLE=v2 python scripts/run_stage8_noise.py

Per target: the same visit is measured by the pipeline and by unTimely; the
differences d_W1, d_W2 over the visits do not contain the source's variability.
Their Pearson correlation estimates how much of the measurement error is
common to the two bands (docs/VALIDATION_8.md section 3.1).
Output: results/s8_error_correlation_<tag>.json (per object: r, n pairs) and
data/cache/s8/untimely_lc_<tag>.json (not in git).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from irecho import config, differential, qa, sample_io, untimely  # noqa: E402
from irecho.analysis import visit_arrays  # noqa: E402

MIN_PAIRS = 6


def pair_differences(v: dict, u: dict) -> tuple[np.ndarray, np.ndarray]:
    pairs = untimely.match_visits(v["mjd"], u["mjd"])
    i = np.array([p[0] for p in pairs], dtype=int)
    j = np.array([p[1] for p in pairs], dtype=int)
    if i.size == 0:
        return np.array([]), np.array([])
    d1 = v["w1_mag"][i] - u["w1_mag"][j]
    d2 = v["w2_mag"][i] - u["w2_mag"][j]
    ok = np.isfinite(d1) & np.isfinite(d2)
    return d1[ok], d2[ok]


def correlation(d1, d2) -> dict:
    out = {"n_pairs": int(len(d1)), "r": None}
    if len(d1) >= MIN_PAIRS and np.std(d1) > 0 and np.std(d2) > 0:
        out["r"] = float(np.corrcoef(d1, d2)[0, 1])
        out["scatter_w1_mmag"] = float(1e3 * np.std(d1, ddof=1))
        out["scatter_w2_mmag"] = float(1e3 * np.std(d2, ddof=1))
    return out


def forward_flags(rows: list[dict], epochs) -> list[int]:
    flag = {}
    for r in rows:
        flag.setdefault(r["EPOCH"], int(r.get("FORWARD", -1)))
    return [flag.get(e, -1) for e in epochs]


def main() -> None:
    tag = sample_io.SAMPLE_TAG
    targets = sample_io.load_sample()
    stars = sample_io.load_stars()
    out_t, star_r, lcs = {}, [], {}
    for n, t in enumerate(targets, start=1):
        name = t["name"]
        if not sample_io.field_is_cached(name):
            continue
        field = sample_io.load_field(t, stars.get(name, []))
        upath = sample_io.field_paths(name)["untimely"]
        if field["target"] is None or len(field["stars"]) < 3 or not upath.exists():
            continue
        by_pid = untimely.read_field(upath)
        sv = qa.star_tables(field["stars"])
        offsets = differential.field_offsets(sv)
        v_t = differential.apply_offsets(visit_arrays(field["target"], "neowise"), offsets, err_key=qa.ERR_KEY)
        rows = by_pid.get(0)
        if rows:
            u = untimely.light_curve(rows)
            out_t[name] = correlation(*pair_differences(v_t, u))
            lcs[name] = {k: [None if not np.isfinite(x) else float(x) for x in u[k]]
                         for k in ("mjd", "w1_flux_jy", "w1_flux_err_jy", "w2_flux_jy", "w2_flux_err_jy")}
            lcs[name]["forward"] = forward_flags(rows, u["epoch"])
            lcs[name]["w1_n_det"] = [int(x) for x in u["w1_n_det"]]
            lcs[name]["w2_n_det"] = [int(x) for x in u["w2_n_det"]]
        loo = differential.leave_one_out(sv, err_key=qa.ERR_KEY)
        for pid, v in zip(field["star_pids"], loo):
            srows = by_pid.get(pid)
            if srows:
                c = correlation(*pair_differences(v, untimely.light_curve(srows)))
                if c["r"] is not None:
                    star_r.append(c["r"])
        if n % 40 == 0:
            print(f"[{n}/{len(targets)}]", flush=True)
    rs = np.array([v["r"] for v in out_t.values() if v["r"] is not None])
    summary = {"sample": tag, "min_pairs": MIN_PAIRS, "n_targets": len(out_t), "n_targets_with_r": int(rs.size),
               "target_r_percentiles_5_25_50_75_95": [float(x) for x in np.percentile(rs, [5, 25, 50, 75, 95])],
               "target_r_median": float(np.median(rs)), "n_stars_with_r": len(star_r),
               "star_r_percentiles_5_25_50_75_95": [float(x) for x in np.percentile(star_r, [5, 25, 50, 75, 95])],
               "note": "Pearson correlation of (pipeline - unTimely) in W1 and W2 over the visits measured by both; "
                       "it mixes the error correlations of the two photometries",
               "objects": out_t}
    (config.RESULTS_DIR / f"s8_error_correlation_{tag}.json").write_text(json.dumps(summary, indent=1) + "\n")
    out_dir = config.CACHE_DIR / "s8"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f"untimely_lc_{tag}.json").write_text(json.dumps(lcs) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "objects"}, indent=1))


if __name__ == "__main__":
    main()
