"""Stage 10: ZTF public light curves (g, r) and exposure metadata for the type 1 objects.

    python scripts/fetch_ztf10.py objects            # which positions have a ZTF object (one TAP upload query)
    python scripts/fetch_ztf10.py curves --workers 3 # light curves that are not on disk yet
    python scripts/fetch_ztf10.py meta               # seeing, airmass per exposure of the fields that are used
    python scripts/fetch_ztf10.py status

Sources (no login): the IRSA ZTF light-curve service (public data releases) and the IRSA TAP
tables ``ztf_objects_dr24`` and ``ztf.ztf_current_meta_sci``. The light curves of the 75 type 1
objects of sample_v1 that were fetched for stage 5b (data/cache/s5/ztf/) are reused, not fetched
again. New files go to data/cache/s10/ztf/ with a sidecar; a file whose sidecar exists is skipped,
so the run can be interrupted. Nothing here is committed to git (docs/DATA.md).
"""

from __future__ import annotations

import csv
import json
import sys
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import fetch_ztf  # noqa: E402
from irecho import config, irsa  # noqa: E402

S10 = config.CACHE_DIR / "s10"
OLD = config.CACHE_DIR / "s5" / "ztf"
NEW = S10 / "ztf"
META = S10 / "ztf_meta"
BANDS = ("g", "r")
TYPE1 = ("Sy1", "Sy1.2", "Sy1.5", "Sy1.8")
META_COLUMNS = ("expid", "obsjd", "seeing", "airmass", "maglimit", "moonillf")


def targets() -> list[dict]:
    """Type 1 objects of the two registered samples with their AllWISE positions."""
    import run_stage8 as r8

    out = []
    for tag in ("v1", "v2"):
        with open(config.TARGETS_DIR / f"sample_{tag}.csv", newline="") as fh:
            pos = {r["name"]: (float(r["ra"]), float(r["dec"])) for r in csv.DictReader(fh)}
        for o in r8.sample_names(tag):
            if o["type"] in TYPE1:
                out.append({"name": o["name"], "sample": tag, "ra": pos[o["name"]][0], "dec": pos[o["name"]][1]})
    if not out:
        raise SystemExit("no type 1 object found: the object tables have blank catalogue columns (the public copy "
                         "leaves them out). Run scripts/restore_catalog_columns.py first.")
    return out


def curve_path(name: str, band: str) -> Path | None:
    """The light-curve file of an object: the stage 5b file when it exists, else the new one."""
    for base in (OLD, NEW):
        path = base / f"{name}_{band}.csv"
        if path.exists():
            return path
    return None


def fetch_objects(tg: list[dict]) -> Path:
    positions = [(i, t["ra"], t["dec"]) for i, t in enumerate(tg)]
    return irsa.fetch_upload_csv(fetch_ztf.OBJECTS_TABLE, fetch_ztf.OBJECTS_COLUMNS, positions, 3.0,
                                 NEW / "objects_dr24_type1.csv", sort_keys=("pid", "oid"), timeout=900)


def have_objects(tg: list[dict]) -> set[tuple[str, str]]:
    with open(fetch_objects(tg), newline="") as fh:
        return {(tg[int(r["pid"])]["name"], r["filtercode"][1]) for r in csv.DictReader(fh)}


def run_curves(workers: int) -> None:
    import requests

    tg = targets()
    have = have_objects(tg)
    jobs = [(t, b) for t in tg for b in BANDS
            if (t["name"], b) in have and curve_path(t["name"], b) is None
            and not (NEW / f"{t['name']}_{b}.csv.meta.json").exists()]
    if "--reverse" in sys.argv:          # a second process can take the list from the other end
        jobs.reverse()
    print(f"{len(jobs)} requests to make, {workers} at a time", flush=True)

    def work(job) -> str | None:
        t, band = job
        if (NEW / f"{t['name']}_{band}.csv.meta.json").exists():
            return None
        start = time.time()
        err = ""
        for _ in range(2):
            try:
                meta = fetch_ztf.fetch_one(requests.Session(), t["ra"], t["dec"], band, NEW / f"{t['name']}_{band}.csv")
            except Exception as exc:  # noqa: BLE001 - one retry, then report
                err = repr(exc)
                time.sleep(5)
                continue
            print(f"{t['name']} {band}: {meta['status']}, {meta['n_rows']} rows in {time.time() - start:.0f}s", flush=True)
            return None
        print(f"{t['name']} {band}: FAILED {err}", flush=True)
        return f"{t['name']}_{band}"

    with ThreadPoolExecutor(max_workers=workers) as pool:
        still = [x for x in pool.map(work, jobs) if x]
    print(f"done; {len(still)} files still failing: {still}", flush=True)


def main_oid_field(path: Path) -> tuple[int, int, int, str] | None:
    """(field, ccdid, qid, filtercode) of the object ID with the most catflags == 0 rows."""
    with open(path, newline="") as fh:
        rows = [r for r in csv.DictReader(fh) if r["catflags"] == "0"]
    if not rows:
        return None
    oid = Counter(r["oid"] for r in rows).most_common(1)[0][0]
    r = next(x for x in rows if x["oid"] == oid)
    return int(r["field"]), int(r["ccdid"], 0), int(r["qid"], 0), r["filtercode"]


def meta_path(field: int, ccdid: int, qid: int, filtercode: str) -> Path:
    return META / f"f{field:06d}_c{ccdid:02d}_q{qid}_{filtercode}.csv"


def run_meta() -> None:
    keys = set()
    for t in targets():
        for b in BANDS:
            path = curve_path(t["name"], b)
            if path is not None:
                key = main_oid_field(path)
                if key is not None:
                    keys.add(key)
    todo = sorted(k for k in keys if not meta_path(*k).exists())
    print(f"{len(keys)} field quadrants, {len(todo)} to fetch", flush=True)
    failed = []
    for n, (field, ccdid, qid, fc) in enumerate(todo, start=1):
        adql = (f"SELECT {', '.join(META_COLUMNS)} FROM ztf.ztf_current_meta_sci WHERE field = {field} "
                f"AND ccdid = {ccdid} AND qid = {qid} AND filtercode = '{fc}'")
        start = time.time()
        try:
            irsa.fetch_csv(adql, meta_path(field, ccdid, qid, fc), timeout=600, retries=2)
        except Exception as exc:  # noqa: BLE001
            failed.append((field, ccdid, qid, fc))
            print(f"{n}/{len(todo)} field {field} ccd {ccdid} q {qid} {fc}: FAILED {exc!r}", flush=True)
            continue
        print(f"{n}/{len(todo)} field {field} ccd {ccdid} q {qid} {fc}: {time.time() - start:.0f}s", flush=True)
    print(f"done; {len(failed)} failed: {failed}", flush=True)


def status() -> None:
    tg = targets()
    have = have_objects(tg)
    n_obj = sum(any((t["name"], b) in have for b in BANDS) for t in tg)
    old = sum(any((OLD / f"{t['name']}_{b}.csv").exists() for b in BANDS) for t in tg)
    got = sum(any(curve_path(t["name"], b) is not None for b in BANDS) for t in tg)
    missing = [f"{t['name']}_{b}" for t in tg for b in BANDS if (t["name"], b) in have and curve_path(t["name"], b) is None]
    print(json.dumps({"type1": len(tg), "with_ztf_object": n_obj, "with_stage5b_file": old, "with_light_curve": got,
                      "requests_missing": len(missing)}))


def main() -> None:
    cmd = sys.argv[1] if len(sys.argv) > 1 else "status"
    if cmd == "objects":
        print(fetch_objects(targets()))
    elif cmd == "curves":
        run_curves(int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 1)
    elif cmd == "meta":
        run_meta()
    else:
        status()


if __name__ == "__main__":
    main()
