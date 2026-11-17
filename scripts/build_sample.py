"""Stage 5a: build the nearby hard-X-ray AGN sample (docs/SAMPLE.md).

    python scripts/build_sample.py targets   # S1-S10 and the redshift rule
    python scripts/build_sample.py stars     # S11: comparison stars per field

Both steps cache every query under data/cache/ and can be interrupted and
restarted: a finished object is not queried again. Only catalogue data are
read here (BAT 105-month, BASS DR2, AllWISE Source Catalog). No single-exposure
or time-domain table is touched.

Outputs
    targets/sample_v1.csv            selected objects (role science / calibration)
    targets/sample_v1_rejected.csv   objects that passed S1 with the first cut they failed
    targets/sample_v1_stars.csv      comparison stars of every field
    targets/sample_v1_funnel.json    counts after every step
"""

from __future__ import annotations

import csv
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import select_targets as st  # noqa: E402
from irecho import analysis, config, irsa, sample, sample_io  # noqa: E402
from irecho.sample import CRITERIA  # noqa: E402

CATALOG_DIR = config.CACHE_DIR / "catalogs"
TARGET_CONE_ARCSEC = CRITERIA.bright_star_radius_arcsec + 15.0
# IRECHO_SAMPLE=v2 builds the independent sample of stage 6B (0.04 < z <= 0.30,
# sample_v1 objects excluded, no redshift rule, no calibration objects). Every
# other criterion and the star selection are those of sample_v1.
TAG = sample_io.SAMPLE_TAG
PREFIX = f"sample_{TAG}"
CAT = sample_io.CACHE_KEY
TARGET_BATCH = 40  # candidate positions per AllWISE query
STAR_FIELD_BATCH = 5  # 1.5 deg fields per AllWISE query
STAR_W1_RANGE = (9.5, 11.5)  # as scripts/select_field_stars.py
STAR_MAG_BINS = ((9.5, 10.5), (10.5, 11.5))
STARS_PER_BIN = 5
SEED_TAG = 5  # rng seed = [config.SEED, SEED_TAG, bat_index]

SAMPLE_FIELDS = [
    "name", "role", "bat_index", "bat_name", "counterpart", "bat_class", "bat_type", "z",
    "gal_lat_deg", "ecl_lat_deg", "dr2_id", "dr2_type", "dr2_z", "dr2_dist_mpc",
    "dr2_logmbh", "dr2_mbh_method", "dr2_sep_arcsec", "allwise_designation", "ra", "dec",
    "sep_ctpt_allwise_arcsec", "w1mpro", "w2mpro", "ext_flg", "cc_flags", "nb", "na",
    "w1rchi2", "w2rchi2", "failed_cut", "n_comparison_stars", "seen_before",
]


def fetch_file(url: str, params: dict | None, path: Path, key: str) -> str:
    """Download a catalogue file once and record it in data/manifest.json."""
    if not path.exists():
        import requests

        resp = requests.get(url, params=params, timeout=300,
                            headers={"User-Agent": "irecho/0.1 (python-requests)"})
        resp.raise_for_status()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(resp.content)
        meta = {
            "url": url, "params": params or {},
            "retrieved_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "n_rows": max(resp.text.count("\n") - 1, 0),
            "sha256": irsa.sha256_bytes(resp.content), "bytes": len(resp.content),
        }
        path.with_suffix(path.suffix + ".meta.json").write_text(json.dumps(meta, indent=1) + "\n")
        irsa._update_manifest(key, meta, config.MANIFEST_PATH)
    return path.read_text(errors="replace")


def load_catalogues() -> tuple[list[dict], list[dict]]:
    bat = sample.parse_bat105(fetch_file(
        sample.BAT105_URL, None, CATALOG_DIR / "bat105" / "BAT_105m_catalog_07jul2019.txt",
        "catalogs/bat105/BAT_105m_catalog_07jul2019.txt"))
    dr2 = sample.parse_vizier_tsv(fetch_file(
        sample.BASS_DR2_URL, sample.BASS_DR2_PARAMS,
        CATALOG_DIR / "bass_dr2" / "J_ApJS_261_2_table9.tsv",
        "catalogs/bass_dr2/J_ApJS_261_2_table9.tsv"))
    return bat, dr2


def base_row(row: dict, match: dict | None) -> dict:
    name = f"bat{row['bat_index']:04d}"
    role = "calibration" if row["counterpart"] in sample.CALIBRATION_COUNTERPARTS else "science"
    out = {
        "name": name, "role": role, "bat_index": row["bat_index"], "bat_name": row["bat_name"],
        "counterpart": row["counterpart"], "bat_class": row["bat_class"],
        "bat_type": row["bat_type"], "z": f"{row['z']:.4f}",
        "gal_lat_deg": f"{sample.galactic_latitude_deg(row['ctpt_ra'], row['ctpt_dec']):.2f}",
        "ecl_lat_deg": f"{sample.ecliptic_latitude_deg(row['ctpt_ra'], row['ctpt_dec']):.2f}",
    }
    if match is not None:
        out.update({
            "dr2_id": match["ID"] + match.get("m_ID", ""), "dr2_type": match["Type"],
            "dr2_z": match["z"], "dr2_dist_mpc": match["Dist"], "dr2_logmbh": match["logMBH"],
            "dr2_mbh_method": match["Method"], "dr2_sep_arcsec": f"{match['sep_arcsec']:.2f}",
        })
    return out


def step_targets() -> None:
    bat, dr2 = load_catalogues()
    funnel: list[tuple[str, int]] = [("BAT 105-month sources", len(bat))]
    survivors: list[tuple[dict, dict | None]] = []
    rejected: list[dict] = []
    first_fail: dict[str, int] = {}
    if TAG == "v2":
        with open(config.TARGETS_DIR / "sample_v1.csv", newline="") as fh:
            v1_indices = {int(r["bat_index"]) for r in csv.DictReader(fh)}
        cut = lambda row, match: sample.catalogue_cut_v2(row, match, v1_indices)  # noqa: E731
    else:
        cut = sample.catalogue_cut
    for row in bat:
        calib = TAG == "v1" and row["counterpart"] in sample.CALIBRATION_COUNTERPARTS
        match = sample.match_dr2(row, dr2, CRITERIA.dr2_match_arcsec) if np.isfinite(
            row["ctpt_ra"]) else None
        fail = cut(row, match)
        if fail:
            first_fail[fail] = first_fail.get(fail, 0) + 1
            if fail != "S1 class" and not calib:
                rejected.append({**base_row(row, match), "failed_cut": fail})
            if not calib:
                continue
        survivors.append((row, match))
    remaining = len(bat)
    labels = ["S1 class", "S2 redshift", "S3 Galactic latitude", "S4 ecliptic latitude",
              "S5 BASS DR2 beamed/lensed"]
    if TAG == "v2":
        labels.insert(2, "E0 already in sample_v1")
    for label in labels:
        remaining -= first_fail.get(label, 0)
        suffix = ""
        if label.startswith("S2"):
            suffix = (f" ({sample.V2.z_min_exclusive:.3f} < z <= {sample.V2.z_max:.3f})"
                      if TAG == "v2" else f" (z <= {max(CRITERIA.z_grid):.3f})")
        funnel.append((label + suffix, remaining))

    # S6-S10 need the AllWISE Source Catalog around every surviving candidate.
    # The cones are fetched TARGET_BATCH positions per query (TAP table upload).
    fields: dict[int, dict[str, np.ndarray]] = {}
    for b in range(0, len(survivors), TARGET_BATCH):
        chunk = survivors[b:b + TARGET_BATCH]
        start = time.time()
        path = irsa.fetch_upload_csv(
            config.TABLE_ALLWISE_CAT, st.FIELD_COLUMNS,
            [(row["bat_index"], row["ctpt_ra"], row["ctpt_dec"]) for row, _ in chunk],
            TARGET_CONE_ARCSEC, config.CACHE_DIR / "allwise_cat" / f"{CAT}_targets_{b:04d}.csv",
            sort_keys=("pid", "designation"), manifest_path=config.MANIFEST_PATH)
        fields.update(irsa.split_by_pid(irsa.read_csv_table(path)))
        print(f"AllWISE cones {b + 1}-{b + len(chunk)} of {len(survivors)} "
              f"({time.time() - start:.0f}s)", flush=True)
    empty = {c: np.array([]) for c in st.FIELD_COLUMNS}
    passed: list[dict] = []
    allwise_fail: dict[str, int] = {}
    for n, (row, match) in enumerate(survivors, start=1):
        out = base_row(row, match)
        field = fields.get(row["bat_index"], empty)
        fail, i = sample.allwise_cut(field, row["ctpt_ra"], row["ctpt_dec"])
        if i >= 0:
            sep = sample.sep_arcsec(field["ra"][i], field["dec"][i], row["ctpt_ra"],
                                    row["ctpt_dec"])
            out.update({
                "allwise_designation": str(field["designation"][i]),
                "ra": f"{field['ra'][i]:.7f}", "dec": f"{field['dec'][i]:.7f}",
                "sep_ctpt_allwise_arcsec": f"{float(sep):.2f}",
                "w1mpro": f"{field['w1mpro'][i]:.3f}", "w2mpro": f"{field['w2mpro'][i]:.3f}",
                "ext_flg": int(field["ext_flg"][i]), "cc_flags": str(field["cc_flags"][i]),
                "nb": int(field["nb"][i]), "na": int(field["na"][i]),
                "w1rchi2": f"{field['w1rchi2'][i]:.2f}", "w2rchi2": f"{field['w2rchi2'][i]:.2f}",
            })
        catalogue_fail = cut(row, match)
        out["failed_cut"] = catalogue_fail or fail
        if fail and not catalogue_fail:  # an object leaves the funnel only once
            allwise_fail[fail] = allwise_fail.get(fail, 0) + 1
        if out["role"] == "calibration":
            if i < 0:
                raise RuntimeError(
                    f"calibration object {row['counterpart']} has no AllWISE match")
            passed.append(out)  # kept whatever the cuts say; failed_cut records them
        elif out["failed_cut"]:
            rejected.append(out)
        else:
            passed.append(out)
    for label in ("S6 no AllWISE counterpart", "S7 saturation margin", "S8 faint limit",
                  "S9 cc_flags / ph_qual", "S10 bright star nearby", "S10 blended neighbour"):
        remaining -= allwise_fail.get(label, 0)
        funnel.append((label, remaining))

    # S2 rule on the objects that passed everything else (calibration objects
    # that failed a cut are kept in the file but do not count).
    clean = [r for r in passed if not r["failed_cut"]]
    if TAG == "v2":  # no redshift rule: every object that passed is kept
        z_max, counts, selected = sample.V2.z_max, {}, list(passed)
        n_clean = len(clean)
    else:
        z_max, counts = sample.choose_z_max([float(r["z"]) for r in clean])
        selected = [r for r in passed if float(r["z"]) <= z_max or r["role"] == "calibration"]
        for r in passed:
            if r not in selected:
                rejected.append({**r, "failed_cut": f"S2 redshift rule (z > {z_max:.3f})"})
        n_clean = sum(1 for r in selected if not r["failed_cut"])
        funnel.append((f"S2 redshift rule: z <= {z_max:.3f}", n_clean))
    write_csv(config.TARGETS_DIR / f"{PREFIX}.csv", selected)
    write_csv(config.TARGETS_DIR / f"{PREFIX}_rejected.csv", rejected)
    summary = {
        "criteria": {k: getattr(CRITERIA, k) for k in CRITERIA.__dataclass_fields__},
        "funnel_targets": funnel,
        "count_by_z_grid_after_S1_S10": counts,
        "z_max": z_max,
        "n_selected_passing_all_cuts": n_clean,
        "n_calibration": sum(1 for r in selected if r["role"] == "calibration"),
        "calibration_failed_cuts": {r["counterpart"]: r["failed_cut"] for r in selected
                                    if r["role"] == "calibration"},
        "below_min_sample_size": n_clean < CRITERIA.min_sample_size,
    }
    (config.TARGETS_DIR / f"{PREFIX}_funnel.json").write_text(json.dumps(summary, indent=1) + "\n")
    for label, count in funnel:
        print(f"  {label}: {count}")
    print(f"  selected rows written: {len(selected)}")


def write_csv(path: Path, rows: list[dict]) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=SAMPLE_FIELDS, extrasaction="ignore",
                                lineterminator="\n")
        writer.writeheader()
        for r in sorted(rows, key=lambda r: int(r["bat_index"])):
            writer.writerow({k: r.get(k, "") for k in SAMPLE_FIELDS})


def read_sample() -> list[dict]:
    with open(config.TARGETS_DIR / f"{PREFIX}.csv", newline="") as fh:
        return list(csv.DictReader(fh))


def step_stars() -> None:
    rows = read_sample()
    stars: list[dict] = []
    funnels = {}
    faint_limit = STAR_W1_RANGE[1] + st.ISOLATION_DMAG
    empty = {c: np.array([]) for c in st.FIELD_COLUMNS}
    fields: dict[int, dict[str, np.ndarray]] = {}
    for n, r in enumerate(rows, start=1):
        ra, dec = float(r["ra"]), float(r["dec"])
        start = time.time()
        if int(r["bat_index"]) not in fields:
            # One query for STAR_FIELD_BATCH fields; only the current batch is
            # kept in memory.
            b = (n - 1) // STAR_FIELD_BATCH * STAR_FIELD_BATCH
            chunk = rows[b:b + STAR_FIELD_BATCH]
            path = irsa.fetch_upload_csv(
                config.TABLE_ALLWISE_CAT, st.FIELD_COLUMNS,
                [(int(c["bat_index"]), float(c["ra"]), float(c["dec"])) for c in chunk],
                st.FIELD_RADIUS_DEG * 3600.0 + 2 * st.ISOLATION_ARCSEC,
                config.CACHE_DIR / "allwise_cat" / f"{CAT}_fields_{b:04d}.csv",
                extra_where=f"t.w1mpro < {faint_limit:.1f}",
                sort_keys=("pid", "designation"), manifest_path=config.MANIFEST_PATH)
            fields = irsa.split_by_pid(irsa.read_csv_table(path))
        field = fields.get(int(r["bat_index"]), empty)
        # Known limitation (docs/SAMPLE.md section 6): star_candidates measures
        # distances with a small-angle formula that does not wrap at RA 0/360,
        # so in fields that straddle RA 0 only one side supplies candidates.
        # Left as it is so that this script reproduces the committed star list.
        ok, funnel, dist, pm = st.star_candidates(field, ra, dec, STAR_W1_RANGE)
        w1, w2 = field["w1mpro"], field["w2mpro"]
        rng = np.random.default_rng([config.SEED, SEED_TAG, int(r["bat_index"])])
        chosen: list[int] = []
        for lo, hi in STAR_MAG_BINS:
            idx = np.flatnonzero(ok & (w1 >= lo) & (w1 < hi))
            idx = idx[np.argsort(field["designation"][idx])]
            take = min(STARS_PER_BIN, idx.size)
            picked = rng.choice(idx, size=take, replace=False) if take else np.array([], int)
            funnel.append((f"bin {lo}-{hi}: candidates / drawn", f"{idx.size} / {take}"))
            chosen.extend(sorted(int(i) for i in picked))
        funnels[r["name"]] = funnel
        r["n_comparison_stars"] = len(chosen)
        for rank, i in enumerate(chosen, start=1):
            stars.append({
                "name": f"{r['name']}_s{rank:02d}", "field": r["name"],
                "allwise_designation": str(field["designation"][i]),
                "ra": f"{field['ra'][i]:.7f}", "dec": f"{field['dec'][i]:.7f}",
                "w1mpro": f"{w1[i]:.3f}", "w2mpro": f"{w2[i]:.3f}",
                "w1rchi2": f"{field['w1rchi2'][i]:.2f}", "w2rchi2": f"{field['w2rchi2'][i]:.2f}",
                "var_flg": str(field["var_flg"][i]),
                "j_ks": f"{field['j_m_2mass'][i] - field['k_m_2mass'][i]:.3f}",
                "pm_mas_yr": f"{pm[i]:.0f}",
                "dist_to_target_deg": f"{dist[i] / 3600.0:.3f}",
            })
        print(f"[{n}/{len(rows)}] {r['name']}: {len(chosen)} stars "
              f"({time.time() - start:.0f}s)", flush=True)

    # Objects whose light curves were already inspected in stages 1-4.
    seen = {analysis.load_ngc3786()["allwise_designation"]: "stage 1-4 (NGC 3786)"}
    seen.update({c["allwise_designation"]: "stage 2-4 control" for c in analysis.load_controls()})
    kept, dropped = [], []
    for r in rows:
        r["seen_before"] = seen.get(r["allwise_designation"], "")
        few = int(r["n_comparison_stars"]) < CRITERIA.min_comparison_stars
        if few and not r["failed_cut"]:
            r["failed_cut"] = "S11 comparison stars"
        if few and r["role"] != "calibration":
            dropped.append(r)
        else:
            kept.append(r)
    write_csv(config.TARGETS_DIR / f"{PREFIX}.csv", kept)
    if dropped:
        with open(config.TARGETS_DIR / f"{PREFIX}_rejected.csv", newline="") as fh:
            rejected = list(csv.DictReader(fh))
        write_csv(config.TARGETS_DIR / f"{PREFIX}_rejected.csv", rejected + dropped)
    kept_names = {r["name"] for r in kept}
    stars = [s for s in stars if s["field"] in kept_names]
    with open(config.TARGETS_DIR / f"{PREFIX}_stars.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(stars[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(stars)
    path = config.TARGETS_DIR / f"{PREFIX}_funnel.json"
    summary = json.loads(path.read_text())
    n_clean = sum(1 for r in kept if not r["failed_cut"])
    summary["funnel_targets"] = [f for f in summary["funnel_targets"]
                                 if not str(f[0]).startswith("S11")]
    summary["funnel_targets"].append([f"S11 at least {CRITERIA.min_comparison_stars} "
                                      "comparison stars", n_clean])
    summary["n_selected_passing_all_cuts"] = n_clean
    summary["below_min_sample_size"] = n_clean < CRITERIA.min_sample_size
    summary["star_seed"] = [config.SEED, SEED_TAG, "bat_index"]
    summary["star_w1_range"] = STAR_W1_RANGE
    summary["n_stars"] = len(stars)
    summary["stars_per_field"] = {str(k): int(sum(1 for r in kept
                                                  if int(r["n_comparison_stars"]) == k))
                                  for k in sorted({int(r["n_comparison_stars"]) for r in kept})}
    summary["seen_before"] = {r["name"]: r["seen_before"] for r in kept if r["seen_before"]}
    summary["star_funnel_per_field"] = funnels
    path.write_text(json.dumps(summary, indent=1) + "\n")
    print(f"kept {len(kept)} objects ({n_clean} passing all cuts), {len(stars)} stars")


def main() -> None:
    if len(sys.argv) != 2 or sys.argv[1] not in ("targets", "stars"):
        raise SystemExit(__doc__)
    config.TARGETS_DIR.mkdir(parents=True, exist_ok=True)
    (step_targets if sys.argv[1] == "targets" else step_stars)()


if __name__ == "__main__":
    main()
