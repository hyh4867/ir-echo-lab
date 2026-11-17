"""Fetch the public ZTF light curves (g, r) at the positions of NGC 3786 and NGC 2273.

    python scripts/fetch_ztf.py

Source: IRSA ZTF light-curve service (public data releases; no login),
https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves, cone of 3 arcsec.
Result: data/cache/ztf/<name>_<band>.csv with a .meta.json sidecar; run with
--register-only to copy the sidecars into data/manifest.json.
"""

from __future__ import annotations

import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from irecho import analysis, config, irsa  # noqa: E402

URL = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves"
RADIUS_DEG = 3.0 / 3600.0
BANDS = ("g", "r")


def register_manifest() -> None:
    manifest = json.loads(config.MANIFEST_PATH.read_text())
    for meta_path in sorted((config.CACHE_DIR / "ztf").glob("*.csv.meta.json")):
        manifest["ztf/" + meta_path.name.removesuffix(".meta.json")] = json.loads(
            meta_path.read_text())
    config.MANIFEST_PATH.write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")


def fetch_one(session, ra: float, dec: float, band: str, path: Path) -> dict:
    """One light-curve request; writes the CSV and its sidecar. A response
    without the CSV header is recorded in the sidecar as status 'unexpected'."""
    params = {"POS": f"CIRCLE {ra:.7f} {dec:.7f} {RADIUS_DEG:.7f}",
              "BANDNAME": band, "FORMAT": "csv"}
    resp = session.get(URL, params=params, timeout=600,
                       headers={"User-Agent": "irecho/0.1 (python-requests)"})
    resp.raise_for_status()
    meta = {"source": URL, "params": params,
            "retrieved_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}
    path.parent.mkdir(parents=True, exist_ok=True)
    if resp.text.startswith("oid,"):
        path.write_bytes(resp.content)
        meta.update({"status": "ok", "n_rows": irsa.count_csv_rows(resp.text),
                     "sha256": irsa.sha256_file(path), "bytes": len(resp.content)})
    else:
        meta.update({"status": "unexpected", "n_rows": 0,
                     "response_head": resp.text[:300]})
    path.with_suffix(".csv.meta.json").write_text(json.dumps(meta, indent=1) + "\n")
    return meta


OBJECTS_TABLE = "ztf_objects_dr24"
OBJECTS_COLUMNS = ("oid", "ra", "dec", "filtercode", "nobs", "ngoodobs", "field", "ccdid",
                   "qid", "medianmag")


def fetch_objects() -> Path:
    """One TAP upload query: the ZTF DR24 objects within 3 arcsec of every
    sample_v1 target (pid = row index in targets/sample_v1.csv). Says which
    targets have a ZTF light curve at all, without the slow light-curve service."""
    from irecho import sample_io

    targets = sample_io.load_sample()
    positions = [(i, t["ra"], t["dec"]) for i, t in enumerate(targets)]
    return irsa.fetch_upload_csv(
        OBJECTS_TABLE, OBJECTS_COLUMNS, positions, 3.0,
        sample_io.S5_DIR / "ztf" / "objects_dr24.csv", sort_keys=("pid", "oid"),
        manifest_path=config.MANIFEST_PATH, timeout=900)


def main_sample() -> None:
    """Stage 5b: g and r light curves at the AllWISE positions of the sample_v1
    targets that have an object in the ZTF objects table. One file per target
    and band under data/cache/s5/ztf/; a file whose sidecar exists is skipped,
    so the run can be interrupted. ``--order FILE`` (one target name per line)
    moves those targets to the front; ``--workers N`` runs N requests at a time
    (the service takes 1-4 minutes per request)."""
    import csv
    from concurrent.futures import ThreadPoolExecutor

    import requests

    from irecho import sample_io

    out_dir = sample_io.S5_DIR / "ztf"
    targets = sample_io.load_sample()
    with open(fetch_objects(), newline="") as fh:
        have = {(int(r["pid"]), r["filtercode"][1]) for r in csv.DictReader(fh)}
    first = []
    if "--order" in sys.argv:
        first = Path(sys.argv[sys.argv.index("--order") + 1]).read_text().split()
    rank = {name: i for i, name in enumerate(first)}
    jobs = [(t, band) for i, t in enumerate(targets) for band in BANDS if (i, band) in have]
    jobs.sort(key=lambda j: rank.get(j[0]["name"], len(rank)))
    jobs = [(t, band) for t, band in jobs
            if not (out_dir / f"{t['name']}_{band}.csv.meta.json").exists()]
    if "--limit" in sys.argv:
        jobs = jobs[: int(sys.argv[sys.argv.index("--limit") + 1])]
    workers = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 1
    print(f"{len(jobs)} requests to make, {workers} at a time", flush=True)

    def work(job) -> str | None:
        t, band = job
        start = time.time()
        for attempt in range(2):
            try:
                meta = fetch_one(requests.Session(), float(t["ra"]), float(t["dec"]), band,
                                 out_dir / f"{t['name']}_{band}.csv")
            except Exception as exc:  # noqa: BLE001 - one retry, then report
                err = repr(exc)
                continue
            print(f"{t['name']} {band}: {meta['status']}, {meta['n_rows']} rows in "
                  f"{time.time() - start:.0f}s", flush=True)
            return None
        print(f"{t['name']} {band}: FAILED {err}", flush=True)
        return f"{t['name']}_{band}"

    with ThreadPoolExecutor(max_workers=workers) as pool:
        still = [x for x in pool.map(work, jobs) if x]
    print(f"done; {len(still)} files still failing: {still}", flush=True)


def register_sample_manifest() -> None:
    from irecho import sample_io

    manifest = json.loads(config.MANIFEST_PATH.read_text())
    for meta_path in sorted((sample_io.S5_DIR / "ztf").glob("*.csv.meta.json")):
        manifest["s5/ztf/" + meta_path.name.removesuffix(".meta.json")] = json.loads(
            meta_path.read_text())
    config.MANIFEST_PATH.write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")


def main() -> None:
    import requests

    if "--sample" in sys.argv:
        main_sample()
        return
    if "--register-sample" in sys.argv:
        register_sample_manifest()
        return
    if "--register-only" in sys.argv:
        register_manifest()
        return
    # NGC 3786 is the pre-specified comparison (docs/PREDICTIONS.md). NGC 2273
    # was added afterwards as an extra check of the slow-rise candidate.
    targets = [analysis.load_ngc3786()] + [c for c in analysis.load_controls()
                                           if c["name"] == "ngc2273"]
    for t, band in [(t, band) for t in targets for band in BANDS]:
        path = config.CACHE_DIR / "ztf" / f"{t['name']}_{band}.csv"
        if path.exists():
            print(f"{t['name']} {band}: cached")
            continue
        params = {"POS": f"CIRCLE {t['ra']:.7f} {t['dec']:.7f} {RADIUS_DEG:.7f}",
                  "BANDNAME": band, "FORMAT": "csv"}
        start = time.time()
        resp = requests.get(URL, params=params, timeout=600,
                            headers={"User-Agent": "irecho/0.1 (python-requests)"})
        resp.raise_for_status()
        if not resp.text.startswith("oid,"):
            raise RuntimeError(f"unexpected response: {resp.text[:200]}")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(resp.content)
        meta = {
            "source": URL, "params": params,
            "retrieved_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "n_rows": irsa.count_csv_rows(resp.text),
            "sha256": irsa.sha256_file(path), "bytes": len(resp.content),
        }
        path.with_suffix(".csv.meta.json").write_text(json.dumps(meta, indent=1) + "\n")
        print(f"{t['name']} {band}: {meta['n_rows']} rows in {time.time() - start:.0f}s",
              flush=True)
    if "--no-manifest" not in sys.argv:
        register_manifest()


if __name__ == "__main__":
    main()
