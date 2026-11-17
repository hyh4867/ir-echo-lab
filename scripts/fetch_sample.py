"""Stage 5a: download the single-exposure tables of the sample_v1 fields.

    python scripts/fetch_sample.py                 # every field, both surveys
    python scripts/fetch_sample.py --limit 5       # first five fields (timing run)
    python scripts/fetch_sample.py --only bat0560

One TAP query per field and survey returns the target and its comparison
stars together (positions are uploaded as a table). A finished field is never
queried again, so the script can be stopped and restarted. A failed field is
logged and retried once at the end. Timings go to data/fetch_sample_timing.jsonl.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from irecho import config, sample_io  # noqa: E402

TIMING_PATH = config.PROJECT_ROOT / "data" / ("fetch_sample_timing.jsonl" if sample_io.SAMPLE_TAG == "v1"
                                              else f"fetch_sample_timing_{sample_io.SAMPLE_TAG}.jsonl")


def fetch_one(target: dict, stars: list[dict], surveys) -> dict:
    record = {"name": target["name"], "n_positions": 1 + len(stars),
              "utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}
    for survey in surveys:
        path = sample_io.field_paths(target["name"])[survey]
        if path.exists():
            record[survey] = "cached"
            continue
        start = time.time()
        sample_io.fetch_field(target, stars, surveys=(survey,))
        meta = json.loads(path.with_suffix(".csv.meta.json").read_text())
        record[survey] = {"seconds": round(time.time() - start, 1), "n_rows": meta["n_rows"],
                          "bytes": meta["bytes"]}
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, help="only the first N fields")
    parser.add_argument("--only", help="a single field by name")
    parser.add_argument("--survey", choices=("neowise", "allwise", "both"), default="both")
    args = parser.parse_args()
    surveys = ("neowise", "allwise") if args.survey == "both" else (args.survey,)

    targets = sample_io.load_sample()
    stars = sample_io.load_stars()
    if args.only:
        targets = [t for t in targets if t["name"] == args.only]
    if args.limit:
        targets = targets[:args.limit]

    failed: list[dict] = []
    for attempt, todo in ((1, targets), (2, failed)):
        if attempt == 2:
            todo, failed = list(failed), []
            if not todo:
                break
            print(f"retrying {len(todo)} failed fields", flush=True)
        for i, t in enumerate(todo, start=1):
            try:
                record = fetch_one(t, stars.get(t["name"], []), surveys)
            except Exception as exc:  # keep going; the field is retried at the end
                failed.append(t)
                print(f"[{i}/{len(todo)}] {t['name']}: FAILED {type(exc).__name__}: "
                      f"{str(exc)[:200]}", flush=True)
                continue
            if any(isinstance(record.get(s), dict) for s in surveys):
                with open(TIMING_PATH, "a") as fh:
                    fh.write(json.dumps(record) + "\n")
            summary = ", ".join(
                f"{s} {record[s]['seconds']}s/{record[s]['n_rows']} rows"
                if isinstance(record[s], dict) else f"{s} cached" for s in surveys)
            print(f"[{i}/{len(todo)}] {t['name']}: {summary}", flush=True)
    print(f"done; {len(failed)} fields still failing: {[t['name'] for t in failed]}", flush=True)


if __name__ == "__main__":
    main()
