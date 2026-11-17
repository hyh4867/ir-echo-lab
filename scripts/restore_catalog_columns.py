"""Fill in the catalogue columns that the public repository leaves blank.

    python scripts/restore_catalog_columns.py

The public copy of this repository keeps object names and IDs but not the
per-object values copied from the Swift-BAT 105-month catalogue (class, type,
redshift) or from the BASS DR1/DR2 tables (type, redshift, distance, black-hole
mass, bolometric luminosity, Eddington ratio, column density), because the terms
for redistributing those tables were not confirmed (see README, "Data").

This script downloads the same tables from their original sources (the Swift
site at NASA GSFC and VizieR), matches them to the objects exactly as
scripts/build_sample.py and scripts/build_props.py did, and writes the blank
columns back into:

    targets/sample_v1.csv, sample_v1_rejected.csv, sample_v2.csv, sample_v2_rejected.csv
    targets/sample_v1_props.csv
    results/s5b_objects.csv, s5b_objects_nofloor.csv, s6b_objects.csv
    results/s5b_summary.json, s6a_flares.json

Nothing else in those files is touched. The downloads are cached under
data/cache/catalogs/ (git-ignored) and recorded in data/manifest.json. Please cite Oh et al. (2018), Koss et al.
(2022) and Ricci et al. (2017) when you use the restored values.
"""

from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import build_props  # noqa: E402
import build_sample  # noqa: E402
from irecho import config, sample  # noqa: E402
from irecho.sample import CRITERIA  # noqa: E402

TARGET_COLUMNS = ("bat_class", "bat_type", "z", "dr2_type", "dr2_z", "dr2_dist_mpc",
                  "dr2_logmbh", "dr2_mbh_method")


def read_csv(path: Path) -> tuple[list[str], list[dict]]:
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        return list(reader.fieldnames), list(reader)


def write_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def restore_targets(bat: list[dict], dr2: list[dict]) -> None:
    by_index = {r["bat_index"]: r for r in bat}
    for name in ("sample_v1", "sample_v1_rejected", "sample_v2", "sample_v2_rejected"):
        path = config.TARGETS_DIR / f"{name}.csv"
        fields, rows = read_csv(path)
        for row in rows:
            src = by_index[int(row["bat_index"])]
            match = sample.match_dr2(src, dr2, CRITERIA.dr2_match_arcsec)
            full = build_sample.base_row(src, match)
            if row["dr2_id"] != full.get("dr2_id", ""):
                raise SystemExit(f"{name}: BASS DR2 match of {row['name']} changed "
                                 f"({row['dr2_id']!r} -> {full.get('dr2_id', '')!r})")
            for col in TARGET_COLUMNS:
                row[col] = full.get(col, "")
        write_csv(path, fields, rows)
        print(f"targets/{name}.csv: {len(rows)} rows")


def restore_props() -> dict[str, dict]:
    path = build_props.OUT_PATH
    _, old = read_csv(path)
    build_props.main()  # rewrites targets/sample_v1_props.csv from the downloaded tables
    _, new = read_csv(path)
    for a, b in zip(old, new):
        if (a["name"], a["dr1_swift"]) != (b["name"], b["dr1_swift"]):
            raise SystemExit(f"props: BASS DR1 match of {a['name']} changed")
    return {r["name"]: r for r in new}


def restore_stage5b(props: dict[str, dict]) -> None:
    _, targets = read_csv(config.TARGETS_DIR / "sample_v1.csv")
    z = {t["name"]: t["z"] for t in targets}
    for name in ("s5b_objects.csv", "s5b_objects_nofloor.csv"):
        path = config.RESULTS_DIR / name
        fields, rows = read_csv(path)
        for row in rows:
            p = props[row["name"]]
            row.update(type=p["type"], type_group=p["type_group"], z=z[row["name"]],
                       logedd=p["dr2_logedd"], loglbol=p["dr2_loglbol"], logmbh=p["dr2_logmbh"],
                       lognh=p["dr1_lognh"], above_threshold=p["above_threshold"],
                       above_threshold_nh=p["above_threshold_nh"])
        write_csv(path, fields, rows)
        print(f"results/{name}: {len(rows)} rows")
    path = config.RESULTS_DIR / "s5b_summary.json"
    summary = json.loads(path.read_text())
    for variant in ("floor", "nofloor"):
        for group in summary[variant].values():
            if not isinstance(group, dict):
                continue
            for obj in group.get("h_breath", {}).get("measurable_objects", []):
                p = props[obj["name"]]
                obj.update(type=p["type"], logedd=p["dr2_logedd"], lognh=p["dr1_lognh"])
    path.write_text(json.dumps(summary, indent=1) + "\n")
    print("results/s5b_summary.json")


def restore_stage6b(dr2: list[dict]) -> None:
    by_id: dict[str, list[dict]] = {}
    for r in dr2:
        by_id.setdefault(r["ID"].strip(), []).append(r)
    _, targets = read_csv(config.TARGETS_DIR / "sample_v2.csv")
    by_name = {t["name"]: t for t in targets}
    path = config.RESULTS_DIR / "s6b_objects.csv"
    fields, rows = read_csv(path)
    for row in rows:
        t = by_name[row["name"]]
        zz = float(t["z"])
        cands = by_id.get(t["dr2_id"].strip().rstrip("AB"), []) if t["dr2_id"] else []
        vals = [float(c["logLbol"]) for c in cands if c["logLbol"].strip()]
        row.update(type=t["dr2_type"] or t["bat_type"], z=f"{zz:.6g}",
                   log_lbol=f"{vals[0]:.6g}" if vals else "",
                   x=f"{float(np.log10(1.0 + zz) + 0.5 * vals[0]):.6g}" if vals else "")
    write_csv(path, fields, rows)
    print(f"results/s6b_objects.csv: {len(rows)} rows")


def restore_flares(props: dict[str, dict]) -> None:
    _, targets = read_csv(config.TARGETS_DIR / "sample_v1.csv")
    z = {t["name"]: float(t["z"]) for t in targets}
    path = config.RESULTS_DIR / "s6a_flares.json"
    flares = json.loads(path.read_text())
    num = lambda text: float(text) if text else None  # noqa: E731
    for name, obj in flares.items():
        p = props[name]
        obj.update(type=p["type"], z=z[name], log_lbol_bass=num(p["dr2_loglbol"]),
                   log_edd_bass=num(p["dr2_logedd"]), log_nh_bass=num(p["dr1_lognh"]))
    path.write_text(json.dumps(flares, indent=1) + "\n")
    print("results/s6a_flares.json")


CATALOG_DIR = config.CACHE_DIR / "catalogs"
DOWNLOADS = (  # url, params, cache path, manifest key, parser, minimum number of rows
    (sample.BAT105_URL, None, CATALOG_DIR / "bat105" / "BAT_105m_catalog_07jul2019.txt",
     "catalogs/bat105/BAT_105m_catalog_07jul2019.txt", sample.parse_bat105, 1000),
    (sample.BASS_DR2_URL, sample.BASS_DR2_PARAMS, CATALOG_DIR / "bass_dr2" / "J_ApJS_261_2_table9.tsv",
     "catalogs/bass_dr2/J_ApJS_261_2_table9.tsv", sample.parse_vizier_tsv, 500),
) + tuple(
    (sample.BASS_DR2_URL, {"-source": f"{build_props.DR1_SOURCE}/{tab}", "-out.max": "unlimited", "-out.all": "1"},
     CATALOG_DIR / "bass_dr1" / f"J_ApJS_233_17_{tab}.tsv", f"catalogs/bass_dr1/J_ApJS_233_17_{tab}.tsv",
     sample.parse_vizier_tsv, 500) for tab in ("agns", "Xraydata"))


def download_checked(tries: int = 6, wait_s: float = 30.0) -> None:
    """Download the four catalogue files, refusing truncated or empty answers.

    VizieR sometimes answers a valid request with a header and no rows. Such a
    file is deleted and requested again instead of being cached.
    """
    for url, params, path, key, parse, min_rows in DOWNLOADS:
        for attempt in range(1, tries + 1):
            try:
                n = len(parse(build_sample.fetch_file(url, params, path, key)))
            except (StopIteration, ValueError, KeyError):
                n = 0  # a header without a table
            if n >= min_rows:
                print(f"{key}: {n} rows")
                break
            path.unlink(missing_ok=True)
            path.with_suffix(path.suffix + ".meta.json").unlink(missing_ok=True)
            print(f"{key}: only {n} rows (attempt {attempt} of {tries}); retrying")
            time.sleep(wait_s)
        else:
            raise SystemExit(f"{key}: the server kept returning fewer than {min_rows} rows; "
                             "try again later. No file in targets/ or results/ was changed.")


def main() -> None:
    download_checked()
    bat, dr2 = build_sample.load_catalogues()
    restore_targets(bat, dr2)
    props = restore_props()
    restore_stage5b(props)
    restore_stage6b(dr2)
    restore_flares(props)
    print("done. The restored values come from Oh et al. (2018), Koss et al. (2022) and "
          "Ricci et al. (2017); do not commit them to a public copy.")


if __name__ == "__main__":
    main()
