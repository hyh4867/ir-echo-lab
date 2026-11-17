"""Stage 5a: injection tests on the comparison stars of the sample_v1 fields.

    python scripts/run_injection_stars.py            # all cached fields
    python scripts/run_injection_stars.py --aggregate-only
    python scripts/run_injection_stars.py --floor-mmag 20   # variant, see below

Detectors (a) and (b) are run on comparison-star light curves only
(leave-one-out differentially corrected), with the stage-4 injection grid
(config.InjectionGrid). No target light curve is read. One checkpoint file per
field is written to data/cache/s5/injection/, so the run can be interrupted.

Variant ``--floor-mmag X``: the stars are made to look like galaxies with an
error floor of X mmag. Independent Gaussian noise of that size is added to
every visit of both bands and the same X is added in quadrature to the visit
errors (as errorfloor.apply_floor does for the targets). Two stars per field
(the first of each brightness bin) are used. This is a model: the extra
visit-level error of a real galaxy need not be Gaussian or independent.

Outputs
    results/qa_injection_stars.json        null rates, completeness by condition
    results/qa_injection_completeness.csv  the same table, one row per condition
    results/qa_injection_stars_floor<X>.json, ..._completeness_floor<X>.csv  (variant)
"""

from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from irecho import config, differential, injection, qa, sample_io  # noqa: E402
from irecho.config import INJECTION  # noqa: E402

CHECKPOINT_DIR = sample_io.S5_DIR / "injection"
W1_BINS = ((9.5, 10.5), (10.5, 11.5))
MIN_STARS_FOR_FIELD = differential.MIN_STARS + 1


FLOOR_STAR_RANKS = (1, 6)  # variant: the first star of each brightness bin
FLOOR_SEED_TAG = 51


def add_floor(v: dict, floor_mag: float, seed) -> dict:
    """Star visit table with extra visit-level noise ``floor_mag`` and the same
    amount added to the errors."""
    rng = np.random.default_rng(seed)
    out = {k: np.array(val, copy=True) for k, val in v.items()}
    for b in ("w1", "w2"):
        out[f"{b}_mag"] = v[f"{b}_mag"] + rng.normal(0.0, floor_mag, v["mjd"].size)
        out[f"{b}_err"] = np.hypot(v[f"{b}_err"], floor_mag)
    return out


def run_field(target: dict, stars: list[dict], floor_mmag: float = 0.0) -> list[dict]:
    """Counts per star and condition: key 'baseline|amp|rise' -> [n, a, b, a_or_b]."""
    field = sample_io.load_field(target, stars)
    if len(field["stars"]) < MIN_STARS_FOR_FIELD:
        return []
    loo = differential.leave_one_out(qa.star_tables(field["stars"]), err_key=qa.ERR_KEY)
    out = []
    for name, pid, w1, v in zip(field["star_names"], field["star_pids"], field["star_w1"], loo):
        if floor_mmag > 0:
            if pid not in FLOOR_STAR_RANKS:
                continue
            v = add_floor(v, 1e-3 * floor_mmag,
                          [config.SEED, FLOOR_SEED_TAG, int(target["bat_index"]), pid,
                           int(round(floor_mmag))])
        counts: dict[str, list[int]] = {}
        for r in injection.run_grid([{"name": name, "kind": "star", "visits": v}],
                                    with_features=False):
            key = f"{r['baseline_visits']}|{r['amp']}|{r['rise']}"
            c = counts.setdefault(key, [0, 0, 0, 0])
            c[0] += 1
            c[1] += int(r["a"])
            c[2] += int(r["b"])
            c[3] += int(r["a"] or r["b"])
        out.append({"star": name, "w1_allwise": w1, "n_visits": int(v["mjd"].size),
                    "counts": counts})
    return out


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion."""
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return float(max(0.0, centre - half)), float(min(1.0, centre + half))


def aggregate(checkpoint_dir: Path = CHECKPOINT_DIR, suffix: str = "",
              floor_mmag: float = 0.0) -> None:
    stars = []
    for path in sorted(checkpoint_dir.glob("*.json")):
        stars.extend(json.loads(path.read_text()))
    groups = {"all": stars}
    for lo, hi in W1_BINS:
        groups[f"{lo}-{hi}"] = [s for s in stars if lo <= s["w1_allwise"] < hi]
    rows = []
    null = {}
    for gname, members in groups.items():
        total: dict[str, np.ndarray] = {}
        for s in members:
            for key, c in s["counts"].items():
                total[key] = total.get(key, np.zeros(4, dtype=int)) + np.array(c)
        any_null_a = sum(1 for s in members
                         if any(c[1] for k, c in s["counts"].items() if k.split("|")[1] == "0.0"))
        any_null_b = sum(1 for s in members
                         if any(c[2] for k, c in s["counts"].items() if k.split("|")[1] == "0.0"))
        n_null = np.zeros(4, dtype=int)
        for key, c in sorted(total.items(), key=lambda kv: [float(x) for x in kv[0].split("|")]):
            base, amp, rise = key.split("|")
            if float(amp) == 0.0:
                n_null += c
                continue
            lo_ci, hi_ci = wilson(int(c[3]), int(c[0]))
            rows.append({
                "star_w1_bin": gname, "baseline_visits": int(base), "amp_mag": float(amp),
                "rise_yr": float(rise), "n_trials": int(c[0]),
                "recovered_a": int(c[1]), "recovered_b": int(c[2]), "recovered_a_or_b": int(c[3]),
                "completeness_a": c[1] / c[0], "completeness_b": c[2] / c[0],
                "completeness_a_or_b": c[3] / c[0],
                "completeness_a_or_b_ci95_low": lo_ci, "completeness_a_or_b_ci95_high": hi_ci,
            })
        upper = 1.0 - 0.05 ** (1.0 / n_null[0]) if n_null[0] and n_null[3] == 0 else None
        null[gname] = {
            "n_stars": len(members), "n_null_segments": int(n_null[0]),
            "flagged_a": int(n_null[1]), "flagged_b": int(n_null[2]),
            "flagged_a_or_b": int(n_null[3]),
            "n_stars_with_a_flag": any_null_a, "n_stars_with_b_flag": any_null_b,
            "false_alarm_rate_a_or_b": n_null[3] / n_null[0] if n_null[0] else None,
            "false_alarm_ci95": list(wilson(int(n_null[3]), int(n_null[0]))),
            "upper_limit_95_if_zero": upper,
        }
    # Table without the baseline split (all baseline lengths pooled).
    pooled = {}
    for r in rows:
        key = (r["star_w1_bin"], r["amp_mag"], r["rise_yr"])
        p = pooled.setdefault(key, np.zeros(4, dtype=int))
        p += np.array([r["n_trials"], r["recovered_a"], r["recovered_b"], r["recovered_a_or_b"]])
    pooled_rows = [{"star_w1_bin": k[0], "amp_mag": k[1], "rise_yr": k[2], "n_trials": int(v[0]),
                    "completeness_a": v[1] / v[0], "completeness_b": v[2] / v[0],
                    "completeness_a_or_b": v[3] / v[0]} for k, v in sorted(pooled.items())]
    with open(config.RESULTS_DIR / f"qa_injection_completeness{suffix}.csv", "w",
              newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        for r in rows:
            writer.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()})
    (config.RESULTS_DIR / f"qa_injection_stars{suffix}.json").write_text(json.dumps({
        "hosts": "comparison stars of the sample_v1 fields, leave-one-out corrected; "
                 "no target light curve is used",
        "added_floor_mmag": floor_mmag,
        "grid": {"amplitudes_mag": INJECTION.amplitudes_mag, "rise_years": INJECTION.rise_years,
                 "baseline_visits": INJECTION.baseline_visits,
                 "post_visits": INJECTION.post_visits},
        "n_fields": len(list(checkpoint_dir.glob("*.json"))),
        "n_stars": len(stars),
        "n_trials_total": int(sum(c[0] for s in stars for c in s["counts"].values())),
        "null": null,
        "completeness_pooled_over_baseline": pooled_rows,
    }, indent=1) + "\n")
    print(json.dumps(null, indent=1))


def main() -> None:
    floor_mmag = 0.0
    if "--floor-mmag" in sys.argv:
        floor_mmag = float(sys.argv[sys.argv.index("--floor-mmag") + 1])
    tag = f"_floor{int(round(floor_mmag))}" if floor_mmag > 0 else ""
    checkpoint_dir = sample_io.S5_DIR / f"injection{tag}"
    if "--aggregate-only" not in sys.argv:
        targets = sample_io.load_sample()
        stars = sample_io.load_stars()
        checkpoint_dir.mkdir(parents=True, exist_ok=True)
        for i, t in enumerate(targets, start=1):
            out = checkpoint_dir / f"{t['name']}.json"
            if out.exists() or not sample_io.field_is_cached(t["name"]):
                continue
            start = time.time()
            rows = run_field(t, stars.get(t["name"], []), floor_mmag)
            tmp = out.with_suffix(".json.part")
            tmp.write_text(json.dumps(rows) + "\n")
            tmp.replace(out)
            print(f"[{i}/{len(targets)}] {t['name']}: {len(rows)} stars "
                  f"({time.time() - start:.0f}s)", flush=True)
    aggregate(checkpoint_dir, tag, floor_mmag)


if __name__ == "__main__":
    main()
