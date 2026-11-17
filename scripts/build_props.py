"""Stage 5b: physical properties of the sample_v1 objects from published tables.

    python scripts/build_props.py

Eddington ratio and bolometric luminosity: BASS DR2 (Koss et al. 2022, VizieR
J/ApJS/261/2 table9, columns logEdd and logLbol; the file cached in stage 5a).
Column density: BASS DR1 X-ray catalogue (Ricci et al. 2017, ApJS 233, 17;
VizieR J/ApJS/233/17, tables 'agns' and 'Xraydata'), matched by position.
Nothing is recomputed. Output: targets/sample_v1_props.csv.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from build_sample import fetch_file  # noqa: E402
from irecho import config, sample, sample_io  # noqa: E402

LAMBDA_EDD_THRESHOLD = 0.02  # PREDICTIONS_5b.md section 6.1
DR1_MATCH_ARCSEC = 5.0
DR1_SOURCE = "J/ApJS/233/17"
OUT_PATH = config.TARGETS_DIR / "sample_v1_props.csv"
FIELDS = ["name", "role", "seen_before", "type_source", "type", "type_group", "z",
          "dr2_logmbh", "dr2_loglbol", "dr2_logedd", "above_threshold",
          "dr1_swift", "dr1_sep_arcsec", "dr1_lognh", "dr1_lognh_flag",
          "lambda_eff_nh", "above_threshold_nh"]
TYPE_1 = ("Sy1", "Sy1.2", "Sy1.5", "Sy1.8")
TYPE_2 = ("Sy1.9", "Sy2", "Sy2 candidate")


def type_group(label: str) -> str:
    if label in TYPE_1:
        return "1"
    if label in TYPE_2:
        return "2"
    return "other"


def _num(text: str) -> float:
    try:
        return float(text)
    except (TypeError, ValueError):
        return float("nan")


def main() -> None:
    cat_dir = config.CACHE_DIR / "catalogs"
    dr2 = sample.parse_vizier_tsv((cat_dir / "bass_dr2" / "J_ApJS_261_2_table9.tsv").read_text(
        errors="replace"))
    dr2_by_id = {}
    for r in dr2:
        dr2_by_id.setdefault(r["ID"].strip(), []).append(r)
    tables = {}
    for tab in ("agns", "Xraydata"):
        text = fetch_file(sample.BASS_DR2_URL,
                          {"-source": f"{DR1_SOURCE}/{tab}", "-out.max": "unlimited",
                           "-out.all": "1"},
                          cat_dir / "bass_dr1" / f"J_ApJS_233_17_{tab}.tsv",
                          f"catalogs/bass_dr1/J_ApJS_233_17_{tab}.tsv")
        tables[tab] = sample.parse_vizier_tsv(text)
    xray = {r["SWIFT"].strip(): r for r in tables["Xraydata"]}
    agn_ra = np.array([_num(r["RAJ2000"]) for r in tables["agns"]])
    agn_dec = np.array([_num(r["DEJ2000"]) for r in tables["agns"]])
    rows = []
    for t in sample_io.load_sample():
        row = {k: "" for k in FIELDS}
        row.update(name=t["name"], role=t["role"], seen_before=t["seen_before"], z=t["z"])
        label = t["dr2_type"] or t["bat_type"]
        row.update(type_source="BASS DR2" if t["dr2_type"] else "BAT 105-month", type=label,
                   type_group=type_group(label))
        if t["dr2_id"]:
            cands = dr2_by_id.get(t["dr2_id"].strip(), [])
            if cands:
                seps = [float(sample.sep_arcsec(t["ra"], t["dec"], _num(c["RAJ2000"]),
                                                _num(c["DEJ2000"]))) for c in cands]
                c = cands[int(np.argmin(seps))]
                row.update(dr2_logmbh=c["logMBH"].strip(), dr2_loglbol=c["logLbol"].strip(),
                           dr2_logedd=c["logEdd"].strip())
        log_edd = _num(row["dr2_logedd"])
        if np.isfinite(log_edd):
            row["above_threshold"] = int(10.0**log_edd >= LAMBDA_EDD_THRESHOLD)
        sep = sample.sep_arcsec(t["ra"], t["dec"], agn_ra, agn_dec)
        k = int(np.nanargmin(sep))
        if sep[k] <= DR1_MATCH_ARCSEC:
            swift = tables["agns"][k]["SWIFT"].strip()
            x = xray.get(swift)
            row.update(dr1_swift=swift, dr1_sep_arcsec=f"{float(sep[k]):.2f}")
            if x and np.isfinite(_num(x["logNH"])):
                row.update(dr1_lognh=x["logNH"].strip(), dr1_lognh_flag=x["f_logNH"].strip())
                lam_eff = max(LAMBDA_EDD_THRESHOLD, 10.0 ** _num(x["logNH"]) / 5e23)
                row["lambda_eff_nh"] = f"{lam_eff:.4g}"
                if np.isfinite(log_edd):
                    row["above_threshold_nh"] = int(10.0**log_edd >= lam_eff)
        rows.append(row)
    with open(OUT_PATH, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    main_rows = [r for r in rows if r["role"] == "science" and not r["seen_before"]]
    n = lambda cond: sum(1 for r in main_rows if cond(r))  # noqa: E731
    print(f"{len(rows)} objects, main sample {len(main_rows)}")
    print("  type 1 / 2 / other:", n(lambda r: r["type_group"] == "1"),
          n(lambda r: r["type_group"] == "2"), n(lambda r: r["type_group"] == "other"))
    print("  with logEdd:", n(lambda r: r["dr2_logedd"] != ""), " above / below threshold:",
          n(lambda r: r["above_threshold"] == 1), n(lambda r: r["above_threshold"] == 0))
    print("  with logNH:", n(lambda r: r["dr1_lognh"] != ""))


if __name__ == "__main__":
    main()
