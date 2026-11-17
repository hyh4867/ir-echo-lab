"""Download the single-exposure tables for every target into data/cache/.

    python scripts/fetch_data.py            # everything in targets/
    python scripts/fetch_data.py --only ngc3786

Queries run one at a time; a finished target is never re-downloaded, so the
script can be interrupted and restarted. Each file's ADQL and SHA-256 are
recorded in data/manifest.json.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from irecho import analysis, irsa  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--only", help="fetch a single target by name")
    parser.add_argument("--refresh", action="store_true", help="re-download even if cached")
    args = parser.parse_args()

    # Stage-2 targets need both surveys; the control-field comparison stars
    # are only used on NEOWISE visits, so their AllWISE tables are not fetched.
    targets = [dict(t, surveys=("neowise", "allwise")) for t in analysis.all_targets()]
    targets += [dict(t, surveys=("neowise",)) for t in analysis.load_field_stars()]
    if args.only:
        targets = [t for t in targets if t["name"] == args.only]
    for i, t in enumerate(targets, start=1):
        paths = {s: p for s, p in irsa.cache_paths(t["name"]).items() if s in t["surveys"]}
        cached = all(p.exists() for p in paths.values()) and not args.refresh
        start = time.time()
        irsa.fetch_target(t["name"], t["ra"], t["dec"], refresh=args.refresh,
                          surveys=t["surveys"])
        status = "cached" if cached else f"downloaded in {time.time() - start:.0f}s"
        print(f"[{i}/{len(targets)}] {t['name']}: {status}", flush=True)


if __name__ == "__main__":
    main()
