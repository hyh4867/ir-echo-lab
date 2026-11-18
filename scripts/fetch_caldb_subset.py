"""Stage 7: download only the NuSTAR CALDB files that are valid for one date.

    python scripts/fetch_caldb_subset.py 2022-05-01 [--dry] [--all]

Run inside the HEASoft environment (needs astropy). The full NuSTAR CALDB is
2.3 GB; for every (instrument, detector, code name, boundary) of the index
this takes the good-quality entry with the latest validity start on or
before the date, which is the entry a CALDB query for that date returns.
Without --all the 'nocut' response matrices (0.5 GB) are left out; the
spectral products of this stage do not query them. (The ghost-ray files were
left out at first, but numkarf opens them also for a point source.) Files go to data/cache/s7/caldb/ under their CALDB paths; the list with
sizes and SHA-256 is written to data/cache/s7/caldb/nustar_subset.json.
"""

import hashlib
from concurrent.futures import ThreadPoolExecutor
import json
import sys
import urllib.request
from pathlib import Path

from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
CALDB = ROOT / "data" / "cache" / "s7" / "caldb"
BASE = "https://heasarc.gsfc.nasa.gov/FTP/caldb/"


def select(index_path: Path, date: str) -> list[str]:
    with fits.open(index_path) as h:
        d = h[1].data
        best = {}
        for r in d:
            if r["CAL_QUAL"] != 0:
                continue
            start = f"{r['CAL_VSD']}T{r['CAL_VST']}"
            if start > f"{date}T00:00:00":
                continue
            key = (r["INSTRUME"], r["DETNAM"], r["CAL_CNAM"], r["CAL_CBD"])
            if key not in best or start >= best[key][0]:
                best[key] = (start, f"{r['CAL_DIR']}/{r['CAL_FILE']}")
    return sorted({v[1] for v in best.values()})


def main() -> None:
    date = sys.argv[1]
    dry = "--dry" in sys.argv
    files = select(CALDB / "data" / "nustar" / "fpm" / "caldb.indx", date)
    skipped = []
    if "--all" not in sys.argv:
        skipped = [f for f in files if "nocutdet" in f]
        files = [f for f in files if f not in skipped]

    def get(rel: str) -> dict:
        path = CALDB / rel
        if dry:
            req = urllib.request.Request(BASE + rel, method="HEAD")
            return {"file": rel, "bytes": int(urllib.request.urlopen(req, timeout=60).headers["Content-Length"])}
        if not path.exists():
            path.parent.mkdir(parents=True, exist_ok=True)
            tmp = path.with_suffix(path.suffix + ".part")
            for attempt in range(4):
                try:
                    urllib.request.urlretrieve(BASE + rel, tmp)
                    break
                except Exception as exc:  # noqa: BLE001
                    if attempt == 3:
                        raise
                    print("retry", rel, exc, flush=True)
            tmp.rename(path)
        rec = {"file": rel, "bytes": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        print(f"{rec['bytes'] / 1e6:8.1f} MB  {rel}", flush=True)
        return rec

    with ThreadPoolExecutor(max_workers=6) as pool:
        out = list(pool.map(get, files))
    total = sum(r["bytes"] for r in out)
    print(len(files), "files", f"{total / 1e6:.1f} MB")
    if not dry:
        with fits.open(CALDB / "data" / "nustar" / "fpm" / "caldb.indx") as h:
            version = h[1].header.get("CALDBVER")
        (CALDB / "nustar_subset.json").write_text(json.dumps(
            {"date": date, "caldb_index_version": version, "n_files": len(out), "bytes": total,
             "left_out": skipped, "files": out},
            indent=1) + "\n")


if __name__ == "__main__":
    main()
