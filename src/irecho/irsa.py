"""IRSA TAP access with a local cache.

Raw query results are stored as CSV under data/cache/ (not committed). Each
file has a sidecar ``.meta.json`` with the exact ADQL, retrieval time, row
count and SHA-256, and the same record is merged into data/manifest.json
(committed) so that a re-download can be compared against what was analysed.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from . import config

# Columns requested from each table (kept explicit so the query is reproducible).
NEOWISE_COLUMNS = (
    "ra", "dec", "mjd", "scan_id", "frame_num", "source_id",
    "w1mpro", "w1sigmpro", "w1snr", "w1rchi2",
    "w2mpro", "w2sigmpro", "w2snr", "w2rchi2",
    "w1sat", "w2sat", "nb", "na", "cc_flags", "det_bit", "ph_qual", "sso_flg",
    "qual_frame", "qi_fact", "saa_sep", "moon_masked", "w1frtr", "w2frtr",
    "allwise_cntr",
)
ALLWISE_MEP_COLUMNS = (
    "cntr_mf", "source_id_mf", "ra", "dec", "nb", "na", "cc_flags",
    "frame_id", "mjd",
    "w1mpro_ep", "w1sigmpro_ep", "w1rchi2_ep",
    "w2mpro_ep", "w2sigmpro_ep", "w2rchi2_ep",
    "qi_fact", "saa_sep", "moon_masked",
)

# Columns that stay as text; everything else is parsed as float (empty -> NaN).
STRING_COLUMNS = frozenset(
    {
        "scan_id", "source_id", "source_id_mf", "frame_id", "cc_flags",
        "ph_qual", "moon_masked", "w1frtr", "w2frtr", "designation",
        "var_flg", "cntr_mf", "allwise_cntr", "cntr",
    }
)

MIN_SECONDS_BETWEEN_QUERIES = 1.0
_last_query_time = 0.0


def cone_adql(table: str, columns, ra_deg: float, dec_deg: float,
              radius_arcsec: float, order_by: str = "mjd") -> str:
    """ADQL for a cone search. Coordinates are formatted to fixed precision so
    that the query text (and therefore the cache key) is deterministic."""
    cols = ", ".join(columns)
    radius_deg = radius_arcsec / 3600.0
    return (
        f"SELECT {cols} FROM {table} "
        f"WHERE CONTAINS(POINT('ICRS', ra, dec), "
        f"CIRCLE('ICRS', {ra_deg:.7f}, {dec_deg:.7f}, {radius_deg:.7f})) = 1 "
        f"ORDER BY {order_by}"
    )


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(Path(path).read_bytes())


def count_csv_rows(text: str) -> int:
    return max(sum(1 for _ in csv.reader(io.StringIO(text))) - 1, 0)


def _update_manifest(key: str, meta: dict, manifest_path: Path) -> None:
    manifest = {}
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
    manifest[key] = meta
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")


def fetch_csv(adql: str, cache_path: Path, *, refresh: bool = False,
              manifest_path: Path | None = None, timeout: float = 300.0,
              retries: int = 3, session=None) -> Path:
    """Run a synchronous TAP query once and cache the CSV.

    A cached file is reused unless ``refresh`` is set. The stored ADQL is
    compared with the requested one so a changed query never silently reuses
    stale data.
    """
    global _last_query_time
    cache_path = Path(cache_path)
    meta_path = cache_path.with_suffix(cache_path.suffix + ".meta.json")
    if cache_path.exists() and meta_path.exists() and not refresh:
        meta = json.loads(meta_path.read_text())
        if meta.get("adql") != adql:
            raise RuntimeError(
                f"{cache_path} was cached with a different query; "
                "pass refresh=True to re-download"
            )
        return cache_path

    import requests  # imported lazily so offline unit tests never need it

    http = session or requests
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        wait = MIN_SECONDS_BETWEEN_QUERIES - (time.time() - _last_query_time)
        if wait > 0:
            time.sleep(wait)
        try:
            _last_query_time = time.time()
            resp = http.get(
                config.TAP_SYNC_URL,
                params={"QUERY": adql, "FORMAT": "csv"},
                timeout=timeout,
                headers={"User-Agent": "irecho/0.1 (python-requests)"},
            )
            resp.raise_for_status()
            data = resp.content
            head = data[:200].lstrip().lower()
            if head.startswith(b"<") or b"error" in head[:60]:
                raise RuntimeError(f"TAP returned a non-CSV response: {data[:300]!r}")
            break
        except Exception as exc:  # network errors are retried, then re-raised
            last_error = exc
            if attempt == retries:
                raise
            time.sleep(5.0 * attempt)
    else:  # pragma: no cover
        raise RuntimeError(str(last_error))

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = cache_path.with_suffix(cache_path.suffix + ".part")
    tmp.write_bytes(data)
    tmp.replace(cache_path)
    meta = {
        "adql": adql,
        "endpoint": config.TAP_SYNC_URL,
        "format": "csv",
        "retrieved_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "n_rows": count_csv_rows(data.decode("utf-8", errors="replace")),
        "sha256": sha256_bytes(data),
        "bytes": len(data),
    }
    meta_path.write_text(json.dumps(meta, indent=1) + "\n")
    if manifest_path is not None:
        key = str(cache_path.relative_to(config.CACHE_DIR)) if cache_path.is_relative_to(
            config.CACHE_DIR) else cache_path.name
        _update_manifest(key, meta, manifest_path)
    return cache_path


# --- several positions in one query (TAP table upload) -------------------------
UPLOAD_NAME = "pos"


def positions_votable(positions) -> str:
    """VOTable with columns pid, pra, pdec for a TAP upload. ``positions`` is a
    sequence of (pid, ra_deg, dec_deg); coordinates are written with fixed
    precision so that the upload (and its hash) is deterministic."""
    rows = "".join(f"<TR><TD>{int(pid)}</TD><TD>{ra:.7f}</TD><TD>{dec:.7f}</TD></TR>"
                   for pid, ra, dec in positions)
    return (
        '<?xml version="1.0"?><VOTABLE version="1.3" '
        'xmlns="http://www.ivoa.net/xml/VOTable/v1.3"><RESOURCE><TABLE>'
        '<FIELD name="pid" datatype="int"/>'
        '<FIELD name="pra" datatype="double" unit="deg"/>'
        '<FIELD name="pdec" datatype="double" unit="deg"/>'
        f"<DATA><TABLEDATA>{rows}</TABLEDATA></DATA></TABLE></RESOURCE></VOTABLE>"
    )


def upload_cone_adql(table: str, columns, radius_arcsec: float, extra_where: str = "") -> str:
    """ADQL that returns, for every uploaded position, the rows of ``table``
    within ``radius_arcsec``. The first output column ``pid`` says which
    position a row belongs to."""
    cols = ", ".join(f"t.{c}" for c in columns)
    where = (f"CONTAINS(POINT('ICRS', t.ra, t.dec), "
             f"CIRCLE('ICRS', p.pra, p.pdec, {radius_arcsec / 3600.0:.7f})) = 1")
    if extra_where:
        where += f" AND {extra_where}"
    return f"SELECT p.pid, {cols} FROM {table} AS t, TAP_UPLOAD.{UPLOAD_NAME} AS p WHERE {where}"


def sort_csv_rows(data: bytes, keys: tuple[str, ...]) -> bytes:
    """Sort the data rows of a CSV by ``keys`` (numeric where possible). The
    server does not promise an order without ORDER BY, and a sorted file has a
    reproducible hash."""
    reader = csv.reader(io.StringIO(data.decode("utf-8", errors="replace")))
    header = next(reader)
    rows = [r for r in reader if r]
    idx = [header.index(k) for k in keys]

    def sort_key(row):
        out = []
        for i in idx:
            try:
                out.append((0, float(row[i]), ""))
            except ValueError:
                out.append((1, 0.0, row[i]))
        return out

    rows.sort(key=sort_key)
    buf = io.StringIO()
    writer = csv.writer(buf, lineterminator="\n")
    writer.writerow(header)
    writer.writerows(rows)
    return buf.getvalue().encode()


def fetch_upload_csv(table: str, columns, positions, radius_arcsec: float, cache_path: Path, *,
                     extra_where: str = "", sort_keys: tuple[str, ...] = ("pid",),
                     refresh: bool = False, manifest_path: Path | None = None,
                     timeout: float = 900.0, retries: int = 3, session=None) -> Path:
    """One synchronous TAP query for many cones (positions are uploaded as a
    table). Cached like :func:`fetch_csv`; the sidecar also stores the uploaded
    positions, and a cached file is only reused for the same query and positions.
    """
    global _last_query_time
    cache_path = Path(cache_path)
    meta_path = cache_path.with_suffix(cache_path.suffix + ".meta.json")
    adql = upload_cone_adql(table, columns, radius_arcsec, extra_where)
    upload = [[int(pid), round(float(ra), 7), round(float(dec), 7)] for pid, ra, dec in positions]
    if cache_path.exists() and meta_path.exists() and not refresh:
        meta = json.loads(meta_path.read_text())
        if meta.get("adql") != adql or meta.get("upload") != upload:
            raise RuntimeError(
                f"{cache_path} was cached with a different query or position list; "
                "pass refresh=True to re-download"
            )
        return cache_path

    import requests  # imported lazily so offline unit tests never need it

    http = session or requests
    votable = positions_votable(upload)
    for attempt in range(1, retries + 1):
        wait = MIN_SECONDS_BETWEEN_QUERIES - (time.time() - _last_query_time)
        if wait > 0:
            time.sleep(wait)
        try:
            _last_query_time = time.time()
            resp = http.post(
                config.TAP_SYNC_URL,
                data={"QUERY": adql, "FORMAT": "csv",
                      "UPLOAD": f"{UPLOAD_NAME},param:{UPLOAD_NAME}"},
                files={UPLOAD_NAME: ("pos.xml", votable, "application/x-votable+xml")},
                timeout=timeout,
                headers={"User-Agent": "irecho/0.1 (python-requests)"},
            )
            resp.raise_for_status()
            data = resp.content
            head = data[:200].lstrip().lower()
            if head.startswith(b"<") or b"error" in head[:60] or not head.startswith(b"pid"):
                raise RuntimeError(f"TAP returned a non-CSV response: {data[:300]!r}")
            break
        except Exception:  # network errors are retried, then re-raised
            if attempt == retries:
                raise
            time.sleep(10.0 * attempt)
    data = sort_csv_rows(data, sort_keys)

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = cache_path.with_suffix(cache_path.suffix + ".part")
    tmp.write_bytes(data)
    tmp.replace(cache_path)
    meta = {
        "adql": adql,
        "upload": upload,
        "upload_columns": ["pid", "pra", "pdec"],
        "sorted_by": list(sort_keys),
        "endpoint": config.TAP_SYNC_URL,
        "format": "csv",
        "retrieved_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "n_rows": count_csv_rows(data.decode("utf-8", errors="replace")),
        "sha256": sha256_bytes(data),
        "bytes": len(data),
    }
    meta_path.write_text(json.dumps(meta, indent=1) + "\n")
    if manifest_path is not None:
        key = str(cache_path.relative_to(config.CACHE_DIR)) if cache_path.is_relative_to(
            config.CACHE_DIR) else cache_path.name
        _update_manifest(key, meta, manifest_path)
    return cache_path


def split_by_pid(table: dict[str, np.ndarray]) -> dict[int, dict[str, np.ndarray]]:
    """Split an upload-query table into one table per uploaded position."""
    pid = table["pid"].astype(int)
    out = {}
    for value in np.unique(pid):
        sel = pid == value
        out[int(value)] = {k: v[sel] for k, v in table.items() if k != "pid"}
    return out


def read_csv_table(path: Path) -> dict[str, np.ndarray]:
    """Read a TAP CSV into a dict of numpy arrays (float unless listed in
    STRING_COLUMNS; empty or 'null' numeric cells become NaN)."""
    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        rows = [r for r in reader if r]
    table: dict[str, np.ndarray] = {}
    for i, name in enumerate(header):
        raw = [r[i] if i < len(r) else "" for r in rows]
        if name in STRING_COLUMNS:
            table[name] = np.array(raw, dtype=str)
        else:
            table[name] = np.array(
                [float(v) if v not in ("", "null", "NULL", "nan") else np.nan for v in raw],
                dtype=float,
            )
    return table


def cache_paths(name: str) -> dict[str, Path]:
    """Cache file locations for one target."""
    return {
        "neowise": config.CACHE_DIR / "neowise" / f"{name}.csv",
        "allwise": config.CACHE_DIR / "allwise_mep" / f"{name}.csv",
    }


def fetch_target(name: str, ra_deg: float, dec_deg: float, *,
                 radius_arcsec: float = config.CONE_RADIUS_ARCSEC,
                 refresh: bool = False,
                 surveys: tuple[str, ...] = ("neowise", "allwise")) -> dict[str, Path]:
    """Download (or reuse) the single-exposure tables for one target."""
    paths = {s: p for s, p in cache_paths(name).items() if s in surveys}
    spec = {"neowise": (config.TABLE_NEOWISE, NEOWISE_COLUMNS),
            "allwise": (config.TABLE_ALLWISE_MEP, ALLWISE_MEP_COLUMNS)}
    for survey, path in paths.items():
        table, columns = spec[survey]
        fetch_csv(cone_adql(table, columns, ra_deg, dec_deg, radius_arcsec), path,
                  refresh=refresh, manifest_path=config.MANIFEST_PATH)
    return paths
