"""Fetch unTimely (time-domain unWISE, Meisner et al. 2023) detections near
NGC 3786 and the control galaxies.

    python scripts/fetch_untimely.py

Source: IRSA's HEALPix-partitioned parquet copy of the unWISE Time-Domain
Catalog (neo7, doi:10.26131/IRSA580) in the public bucket ``nasa-irsa-wise``
(anonymous access). There is no cone-search service for this catalogue, so
the script opens the one partition (HEALPix order 5, NESTED) that contains a
target, reads the coordinate columns, and keeps the rows within
``RADIUS_ARCSEC`` of the target. The other columns are then read one at a
time so that memory stays low. Result: data/cache/untimely/<name>.csv, with
partition path, row count and SHA-256 recorded in data/manifest.json.
"""

from __future__ import annotations

import csv
import io
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from irecho import analysis, config, irsa  # noqa: E402
from irecho.healpix import ang2pix_nest  # noqa: E402

BUCKET_ROOT = ("nasa-irsa-wise/unwise/neo7/catalogs/time_domain/healpix_k5/"
               "unwise-neo7-time_domain-healpix_k5.parquet")
REGION = "us-west-2"
RADIUS_ARCSEC = 12.0  # generous: the analysis applies its own smaller radius
COLUMNS = ("ra", "dec", "band", "EPOCH", "FORWARD", "MJDMIN", "MJDMAX", "MJDMEAN",
           "flux", "dflux", "fluxlbs", "dfluxlbs", "qf", "rchi2", "fracflux", "fwhm",
           "spread_model", "sky", "nm", "primary", "flags_unwise", "flags_info",
           "coadd_id", "unwise_detid")


def partition_dir(ra: float, dec: float) -> str:
    k5 = ang2pix_nest(5, ra, dec)
    return f"{BUCKET_ROOT}/healpix_k0={k5 >> 10}/healpix_k5={k5}"


def fetch_one(name: str, ra: float, dec: float, s3) -> dict:
    import pyarrow.fs as fs
    import pyarrow.parquet as pq

    out_path = config.CACHE_DIR / "untimely" / f"{name}.csv"
    if out_path.exists():
        return {"status": "cached"}
    files = [i.path for i in s3.get_file_info(fs.FileSelector(partition_dir(ra, dec)))
             if i.path.endswith(".parquet")]
    cosd = np.cos(np.deg2rad(dec))
    rows_out: list[dict] = []
    total_rows = 0
    for path in sorted(files):
        pf = pq.ParquetFile(path, filesystem=s3)
        total_rows += pf.metadata.num_rows
        pos = pf.read(columns=["ra", "dec"])
        sep = 3600.0 * np.hypot((pos["ra"].to_numpy() - ra) * cosd,
                                pos["dec"].to_numpy() - dec)
        idx = np.flatnonzero(sep <= RADIUS_ARCSEC)
        del pos
        if idx.size == 0:
            continue
        cols = {"sep_arcsec": np.round(sep[idx], 3)}
        del sep
        for col in COLUMNS:
            arr = pf.read(columns=[col])[col].take(idx).to_pylist()
            cols[col] = [v.decode() if isinstance(v, bytes) else v for v in arr]
        for j in range(idx.size):
            rows_out.append({k: cols[k][j] for k in cols})
    rows_out.sort(key=lambda r: (r["MJDMEAN"], r["band"], r["sep_arcsec"]))
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=["sep_arcsec", *COLUMNS], lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows_out)
    data = buf.getvalue().encode()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_bytes(data)
    meta = {
        "source": "s3://" + partition_dir(ra, dec) + "/ (anonymous, region " + REGION + ")",
        "selection": f"rows within {RADIUS_ARCSEC} arcsec of ra={ra:.7f}, dec={dec:.7f}",
        "partition_files": [p.split("/")[-1] for p in sorted(files)],
        "partition_rows": int(total_rows),
        "retrieved_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "n_rows": len(rows_out),
        "sha256": irsa.sha256_file(out_path),
        "bytes": len(data),
    }
    out_path.with_suffix(".csv.meta.json").write_text(json.dumps(meta, indent=1) + "\n")
    return {"status": "downloaded", "n_rows": len(rows_out), "partition_rows": total_rows}


def fetch_field(name: str, positions, s3, out_path: Path) -> dict:
    """Stage 5a: unTimely rows near the target (pid 0) and near those of its
    comparison stars that fall in the same HEALPix partition as the target.

    Only the target's partition is opened, so the stars cost no extra download;
    stars in a neighbouring partition are skipped and listed in the sidecar.
    """
    import pyarrow.fs as fs
    import pyarrow.parquet as pq

    if out_path.exists():
        return {"status": "cached"}
    _, ra0, dec0 = positions[0]
    part = partition_dir(ra0, dec0)
    inside = [(pid, ra, dec) for pid, ra, dec in positions if partition_dir(ra, dec) == part]
    files = [i.path for i in s3.get_file_info(fs.FileSelector(part)) if i.path.endswith(".parquet")]
    rows_out: list[dict] = []
    total_rows = 0
    for path in sorted(files):
        pf = pq.ParquetFile(path, filesystem=s3)
        total_rows += pf.metadata.num_rows
        pos = pf.read(columns=["ra", "dec"])
        ra_all, dec_all = pos["ra"].to_numpy(), pos["dec"].to_numpy()
        del pos
        idx_parts, pid_parts, sep_parts = [], [], []
        for pid, ra, dec in inside:
            near = np.flatnonzero((np.abs(dec_all - dec) <= RADIUS_ARCSEC / 3600.0))
            sep = 3600.0 * np.hypot((ra_all[near] - ra) * np.cos(np.deg2rad(dec)),
                                    dec_all[near] - dec)
            ok = sep <= RADIUS_ARCSEC
            idx_parts.append(near[ok])
            pid_parts.append(np.full(int(ok.sum()), pid))
            sep_parts.append(np.round(sep[ok], 3))
        del ra_all, dec_all
        idx = np.concatenate(idx_parts)
        if idx.size == 0:
            continue
        order = np.argsort(idx, kind="stable")
        idx, pids, seps = idx[order], np.concatenate(pid_parts)[order], np.concatenate(sep_parts)[order]
        cols = {"pid": pids.tolist(), "sep_arcsec": seps.tolist()}
        for col in COLUMNS:
            arr = pf.read(columns=[col])[col].take(idx).to_pylist()
            cols[col] = [v.decode() if isinstance(v, bytes) else v for v in arr]
        for j in range(idx.size):
            rows_out.append({k: cols[k][j] for k in cols})
    rows_out.sort(key=lambda r: (r["pid"], r["MJDMEAN"], r["band"], r["sep_arcsec"]))
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=["pid", "sep_arcsec", *COLUMNS], lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows_out)
    data = buf.getvalue().encode()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix(".csv.part")
    tmp.write_bytes(data)
    tmp.replace(out_path)
    meta = {
        "source": "s3://" + part + "/ (anonymous, region " + REGION + ")",
        "selection": f"rows within {RADIUS_ARCSEC} arcsec of the uploaded positions",
        "positions": [[int(pid), round(float(ra), 7), round(float(dec), 7)]
                      for pid, ra, dec in inside],
        "positions_in_other_partitions": [int(pid) for pid, ra, dec in positions
                                          if partition_dir(ra, dec) != part],
        "partition_files": [f.split("/")[-1] for f in sorted(files)],
        "partition_rows": int(total_rows),
        "retrieved_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "n_rows": len(rows_out),
        "sha256": irsa.sha256_file(out_path),
        "bytes": len(data),
    }
    out_path.with_suffix(".csv.meta.json").write_text(json.dumps(meta, indent=1) + "\n")
    return {"status": "downloaded", "n_rows": len(rows_out), "n_positions": len(inside),
            "partition_rows": total_rows}


def main_sample() -> None:
    """python scripts/fetch_untimely.py --sample [--limit N]"""
    import pyarrow.fs as fs

    from irecho import sample_io

    targets = sample_io.load_sample()
    stars = sample_io.load_stars()
    if "--limit" in sys.argv:
        targets = targets[:int(sys.argv[sys.argv.index("--limit") + 1])]
    s3 = fs.S3FileSystem(anonymous=True, region=REGION)
    failed = []
    for attempt, todo in ((1, targets), (2, failed)):
        if attempt == 2:
            todo, failed = list(failed), []
        for i, t in enumerate(todo, start=1):
            start = time.time()
            out_path = sample_io.field_paths(t["name"])["untimely"]
            try:
                info = fetch_field(t["name"], sample_io.field_positions(t, stars.get(t["name"], [])),
                                   s3, out_path)
            except Exception as exc:  # report, retry once at the end
                failed.append(t)
                info = {"status": f"FAILED: {type(exc).__name__}: {str(exc)[:200]}"}
            print(f"[{i}/{len(todo)}] {t['name']}: {info} ({time.time() - start:.0f}s)", flush=True)
    print(f"done; {len(failed)} fields still failing: {[t['name'] for t in failed]}", flush=True)


def register_sample_manifest() -> None:
    """Copy the sidecars of the stage-5a unTimely files into data/manifest.json."""
    from irecho import sample_io

    manifest = json.loads(config.MANIFEST_PATH.read_text())
    for meta_path in sorted((sample_io.S5_DIR / "untimely").glob("*.csv.meta.json")):
        key = f"{sample_io.CACHE_KEY}/untimely/" + meta_path.name.removesuffix(".meta.json")
        manifest[key] = json.loads(meta_path.read_text())
    config.MANIFEST_PATH.write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")


def main() -> None:
    import pyarrow.fs as fs

    if "--register-sample" in sys.argv:
        register_sample_manifest()
        return
    if "--sample" in sys.argv:
        main_sample()
        return
    if "--register-only" in sys.argv:
        register_manifest()
        return
    s3 = fs.S3FileSystem(anonymous=True, region=REGION)
    targets = [analysis.load_ngc3786()] + analysis.load_controls()
    for i, t in enumerate(targets, start=1):
        start = time.time()
        try:
            info = fetch_one(t["name"], t["ra"], t["dec"], s3)
        except Exception as exc:  # report and continue with the next target
            info = {"status": f"FAILED: {type(exc).__name__}: {exc}"}
        print(f"[{i}/{len(targets)}] {t['name']}: {info} ({time.time() - start:.0f}s)", flush=True)
    if "--no-manifest" not in sys.argv:
        register_manifest()


def register_manifest() -> None:
    """Copy the sidecar metadata of every cached unTimely file into
    data/manifest.json (done in one step at the end so that this script never
    interleaves manifest writes with a running scripts/fetch_data.py)."""
    manifest = json.loads(config.MANIFEST_PATH.read_text())
    for meta_path in sorted((config.CACHE_DIR / "untimely").glob("*.csv.meta.json")):
        key = "untimely/" + meta_path.name.removesuffix(".meta.json")
        manifest[key] = json.loads(meta_path.read_text())
    config.MANIFEST_PATH.write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
