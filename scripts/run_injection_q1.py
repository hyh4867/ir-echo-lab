"""Stage 5b: completeness of the Q1 'slow rise' class from star injections.

    python scripts/run_injection_q1.py            # floors 0, 10, 20, 50 mmag, resumable
    python scripts/run_injection_q1.py --aggregate-only

Same hosts, grid and noise-floor model as scripts/run_injection_stars.py
(two comparison stars per field, leave-one-out corrected), but the count is of
trials that satisfy the Q1 definition: detector (b) flags the segment and its
fitted brightening is >= 0.2 mag in both bands (PREDICTIONS_5b.md section 3).
No target light curve is read.

Output: results/s5b_q1_completeness.csv, results/s5b_q1_completeness.json.
"""

from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import beta as beta_dist

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import run_injection_stars as base  # noqa: E402
from irecho import config, differential, injection, qa, sample_io  # noqa: E402

FLOORS_MMAG = (0, 10, 20, 50)


def checkpoint_dir(floor: int) -> Path:
    return sample_io.S5_DIR / f"injection_q1_floor{floor}"


def run_field(target: dict, stars: list[dict], floor: int) -> list[dict]:
    field = sample_io.load_field(target, stars)
    if len(field["stars"]) < base.MIN_STARS_FOR_FIELD:
        return []
    loo = differential.leave_one_out(qa.star_tables(field["stars"]), err_key=qa.ERR_KEY)
    out = []
    for name, pid, v in zip(field["star_names"], field["star_pids"], loo):
        if pid not in base.FLOOR_STAR_RANKS:
            continue
        if floor > 0:
            v = base.add_floor(v, 1e-3 * floor, [config.SEED, base.FLOOR_SEED_TAG,
                                                 int(target["bat_index"]), pid, floor])
        counts: dict[str, list[int]] = {}
        for r in injection.run_grid([{"name": name, "kind": "star", "visits": v}],
                                    with_features=False):
            c = counts.setdefault(f"{r['amp']}|{r['rise']}", [0, 0, 0])
            c[0] += 1
            c[1] += int(r["q1"])
            c[2] += int(r["b"])
        out.append({"star": name, "counts": counts})
    return out


def clopper_pearson(k: int, n: int, level: float = 0.95) -> tuple[float, float]:
    a = (1.0 - level) / 2.0
    lo = 0.0 if k == 0 else float(beta_dist.ppf(a, k, n - k + 1))
    hi = 1.0 if k == n else float(beta_dist.ppf(1.0 - a, k + 1, n - k))
    return lo, hi


def aggregate() -> None:
    rows, null = [], {}
    for floor in FLOORS_MMAG:
        stars = []
        for path in sorted(checkpoint_dir(floor).glob("*.json")):
            stars.extend(json.loads(path.read_text()))
        total: dict[str, np.ndarray] = {}
        for s in stars:
            for key, c in s["counts"].items():
                total[key] = total.get(key, np.zeros(3, dtype=int)) + np.array(c)
        n_star_q1 = sum(1 for s in stars if s["counts"].get("0.0|-1.0", [0, 0, 0])[1] > 0)
        seg = total.get("0.0|-1.0", np.zeros(3, dtype=int))
        null[str(floor)] = {
            "n_stars": len(stars), "n_stars_with_q1_in_a_null_segment": n_star_q1,
            "star_rate_upper_95": clopper_pearson(n_star_q1, len(stars))[1] if stars else None,
            "n_null_segments": int(seg[0]), "n_null_segments_q1": int(seg[1]),
            "n_null_segments_b": int(seg[2]),
        }
        for key, c in sorted(total.items(), key=lambda kv: [float(x) for x in kv[0].split("|")]):
            amp, rise = (float(x) for x in key.split("|"))
            if amp == 0.0:
                continue
            lo, hi = clopper_pearson(int(c[1]), int(c[0]))
            rows.append({"floor_mmag": floor, "amp_mag": amp, "rise_yr": rise,
                         "n_trials": int(c[0]), "n_q1": int(c[1]), "n_b": int(c[2]),
                         "completeness_q1": f"{c[1] / c[0]:.4f}",
                         "ci95_low": f"{lo:.4f}", "ci95_high": f"{hi:.4f}"})
    with open(config.RESULTS_DIR / "s5b_q1_completeness.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (config.RESULTS_DIR / "s5b_q1_completeness.json").write_text(json.dumps({
        "definition": "detector (b) flag and fitted brightening >= %.1f mag in both bands"
                      % injection.Q1_MIN_MAG,
        "hosts": "two comparison stars per sample_v1 field, leave-one-out corrected; Gaussian "
                 "noise of the floor size added to both bands and to the errors",
        "null": null}, indent=1) + "\n")
    print(json.dumps(null, indent=1))


def main() -> None:
    if "--aggregate-only" not in sys.argv:
        targets, stars = sample_io.load_sample(), sample_io.load_stars()
        for floor in FLOORS_MMAG:
            checkpoint_dir(floor).mkdir(parents=True, exist_ok=True)
            for n, t in enumerate(targets, 1):
                path = checkpoint_dir(floor) / f"{t['name']}.json"
                if path.exists() or not sample_io.field_is_cached(t["name"]):
                    continue
                start = time.time()
                path.write_text(json.dumps(run_field(t, stars[t["name"]], floor)) + "\n")
                print(f"floor {floor} [{n}/{len(targets)}] {t['name']} ({time.time() - start:.0f}s)",
                      flush=True)
    aggregate()


if __name__ == "__main__":
    main()
