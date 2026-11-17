"""Select the known-answer-test targets and freeze them in targets/.

Run once, before any light curve of these objects is inspected:

    python scripts/select_targets.py

Outputs
    targets/ngc3786.json          flare target (Son et al. 2022)
    targets/reference_stars.csv   non-variable reference stars (criteria R1-R13)
    targets/control_agn.csv       quiet-AGN control sample (criteria C1-C4)

The criteria are documented in docs/PIPELINE.md section 6 and are encoded
below as constants. The reference-star draw uses config.SEED.
"""

from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from irecho import config, irsa  # noqa: E402

# --- Flare target -----------------------------------------------------------
NGC3786_NED_NAME = "NGC 3786"

# --- Reference stars (R1-R13) -----------------------------------------------
FIELD_RADIUS_DEG = 1.5  # R1: within 1.5 deg of NGC 3786 (same visits as the target)
EXCLUDE_NEAR_TARGET_ARCMIN = 3.0  # R1: but not on the NGC 3786 / NGC 3788 pair
W1_RANGE = (9.5, 12.5)  # R2: unsaturated, single-exposure SNR of tens
STAR_COLOR_RANGE = (-0.15, 0.15)  # R3: photospheric W1-W2
MAX_RCHI2 = 3.0  # R10: point-like profile fit
MAX_J_KS = 1.0  # R11: stellar 2MASS colour
MAX_PM_MAS_YR = 150.0  # R12: stays inside the 3 arcsec cone for 14 yr
ISOLATION_ARCSEC = 15.0  # R13: no neighbour within 15 arcsec ...
ISOLATION_DMAG = 4.0  # ... brighter than W1_star + 4 mag
MAG_BINS = ((9.5, 10.5), (10.5, 11.5), (11.5, 12.5))
STARS_PER_BIN = 10

FIELD_COLUMNS = (
    "designation", "ra", "dec", "w1mpro", "w1sigmpro", "w2mpro", "w2sigmpro",
    "w1rchi2", "w2rchi2", "w1sat", "w2sat", "nb", "na", "cc_flags", "ext_flg",
    "var_flg", "ph_qual", "pmra", "pmdec", "n_2mass", "j_m_2mass", "k_m_2mass",
)

# --- Control AGN (C1-C4) ----------------------------------------------------
# C1: optically classified Seyfert 2 nuclei in nearby galaxies (obscured, so
#     the mid-IR is expected to respond slowly), chosen from the literature
#     before looking at any WISE light curve.
# C2: no changing-look event or mid-IR flare known to us at selection time
#     (not an exhaustive literature search; see docs/PIPELINE.md).
# C3: AllWISE W1 > 8.5 and W2 > 8.0 (below the saturation onset).
# C4: |b| > 20 deg (uncrowded field).
CONTROL_CANDIDATES = (
    "Mrk 573", "NGC 1386", "NGC 2273", "NGC 3081",
    "NGC 3393", "NGC 4941", "NGC 5347", "NGC 5728",
)
CONTROL_MATCH_ARCSEC = 3.0


def ned_lookup(name: str) -> dict:
    """Resolve a name with the NED Object Lookup service."""
    import requests

    for attempt in range(3):
        try:
            resp = requests.get(
                "https://ned.ipac.caltech.edu/srs/ObjectLookup",
                params={"name": name, "aliases": "false"}, timeout=60,
                headers={"User-Agent": "irecho/0.1 (python-requests)"},
            )
            resp.raise_for_status()
            data = resp.json()
            pref = data["Preferred"]
            return {
                "name": pref["Name"],
                "ra": float(pref["Position"]["RA"]),
                "dec": float(pref["Position"]["Dec"]),
                "z": pref.get("Redshift", {}).get("Value"),
            }
        except Exception:
            if attempt == 2:
                raise
            time.sleep(5)
    raise RuntimeError("unreachable")


def galactic_latitude_deg(ra_deg: float, dec_deg: float) -> float:
    """Galactic latitude from ICRS (J2000 pole: RA 192.85948, Dec 27.12825)."""
    ra, dec = np.deg2rad(ra_deg), np.deg2rad(dec_deg)
    ra_p, dec_p = np.deg2rad(192.85948), np.deg2rad(27.12825)
    sinb = np.sin(dec) * np.sin(dec_p) + np.cos(dec) * np.cos(dec_p) * np.cos(ra - ra_p)
    return float(np.rad2deg(np.arcsin(sinb)))


def sep_arcsec(ra1, dec1, ra2, dec2):
    """Small-angle separation in arcsec (adequate below a few degrees)."""
    cosd = np.cos(np.deg2rad(0.5 * (np.asarray(dec1) + np.asarray(dec2))))
    return 3600.0 * np.hypot((np.asarray(ra1) - np.asarray(ra2)) * cosd,
                             np.asarray(dec1) - np.asarray(dec2))


def slug(name: str) -> str:
    return name.lower().replace(" ", "")


def allwise_cone(name: str, ra: float, dec: float, radius_arcsec: float,
                 extra_where: str = "") -> dict[str, np.ndarray]:
    adql = irsa.cone_adql(config.TABLE_ALLWISE_CAT, FIELD_COLUMNS, ra, dec,
                          radius_arcsec, order_by="designation")
    if extra_where:
        adql = adql.replace(" ORDER BY", f" AND {extra_where} ORDER BY")
    path = config.CACHE_DIR / "allwise_cat" / f"{name}.csv"
    irsa.fetch_csv(adql, path, manifest_path=config.MANIFEST_PATH, timeout=900)
    return irsa.read_csv_table(path)


def select_ngc3786() -> dict:
    ned = ned_lookup(NGC3786_NED_NAME)
    cat = allwise_cone("ngc3786_cat", ned["ra"], ned["dec"], 10.0)
    sep = sep_arcsec(cat["ra"], cat["dec"], ned["ra"], ned["dec"])
    i = int(np.argmin(sep))
    target = {
        "name": "ngc3786",
        "ned_name": ned["name"], "ned_ra": ned["ra"], "ned_dec": ned["dec"],
        "ned_z": ned["z"],
        "allwise_designation": str(cat["designation"][i]),
        "ra": float(cat["ra"][i]), "dec": float(cat["dec"][i]),
        "sep_ned_allwise_arcsec": round(float(sep[i]), 3),
        "allwise_w1mpro": float(cat["w1mpro"][i]), "allwise_w2mpro": float(cat["w2mpro"][i]),
        "allwise_ext_flg": int(cat["ext_flg"][i]), "allwise_cc_flags": str(cat["cc_flags"][i]),
        "allwise_var_flg": str(cat["var_flg"][i]),
        "allwise_w1rchi2": float(cat["w1rchi2"][i]), "allwise_w2rchi2": float(cat["w2rchi2"][i]),
        "reference": "Son et al. 2022, ApJ 937, 3 (arXiv:2208.08062)",
    }
    (config.TARGETS_DIR / "ngc3786.json").write_text(json.dumps(target, indent=1) + "\n")
    return target


def star_candidates(field: dict[str, np.ndarray], center_ra: float, center_dec: float,
                    w1_range: tuple[float, float] = W1_RANGE):
    """Apply R1-R13 to an AllWISE field table. Returns (mask, funnel, dist_arcsec, pm)."""
    n = len(field["ra"])
    w1, w2 = field["w1mpro"], field["w2mpro"]
    var = field["var_flg"]
    var_ok = np.array(
        [len(v) >= 2 and v[0].isdigit() and v[1].isdigit() and int(v[0]) <= 5 and int(v[1]) <= 5
         for v in var]
    )
    dist = sep_arcsec(field["ra"], field["dec"], center_ra, center_dec)
    pm = np.hypot(field["pmra"], field["pmdec"])
    criteria = {
        "R1 field": (dist <= FIELD_RADIUS_DEG * 3600.0)
        & (dist > EXCLUDE_NEAR_TARGET_ARCMIN * 60.0),
        "R2 W1 range": (w1 >= w1_range[0]) & (w1 < w1_range[1]),
        "R3 W1-W2 colour": ((w1 - w2) > STAR_COLOR_RANGE[0]) & ((w1 - w2) < STAR_COLOR_RANGE[1]),
        "R4 cc_flags": np.array([c[:2] == "00" for c in field["cc_flags"]]),
        "R5 ext_flg": field["ext_flg"] == 0,
        "R6 nb, na": (field["nb"] == 1) & (field["na"] == 0),
        "R7 saturation": (field["w1sat"] == 0) & (field["w2sat"] == 0),
        "R8 var_flg": var_ok,
        "R9 ph_qual": np.array([p[:2] == "AA" for p in field["ph_qual"]]),
        "R10 rchi2": (field["w1rchi2"] < MAX_RCHI2) & (field["w2rchi2"] < MAX_RCHI2),
        "R11 2MASS": (field["n_2mass"] == 1)
        & ((field["j_m_2mass"] - field["k_m_2mass"]) <= MAX_J_KS),
        "R12 proper motion": pm < MAX_PM_MAS_YR,
    }
    ok = np.ones(n, dtype=bool)
    funnel = [("AllWISE sources in field query", n)]
    for label, mask in criteria.items():
        ok &= mask
        funnel.append((label, int(ok.sum())))

    # R13 isolation, evaluated only for the survivors.
    isolated = np.zeros(n, dtype=bool)
    for i in np.flatnonzero(ok):
        d = sep_arcsec(field["ra"], field["dec"], field["ra"][i], field["dec"][i])
        near = (d < ISOLATION_ARCSEC) & (d > 0) & (w1 < w1[i] + ISOLATION_DMAG)
        isolated[i] = not near.any()
    ok &= isolated
    funnel.append(("R13 isolation", int(ok.sum())))
    return ok, funnel, dist, pm


def select_reference_stars(center_ra: float, center_dec: float) -> list[dict]:
    faint_limit = W1_RANGE[1] + ISOLATION_DMAG
    field = allwise_cone(
        "field_ngc3786", center_ra, center_dec,
        FIELD_RADIUS_DEG * 3600.0 + 2 * ISOLATION_ARCSEC,
        extra_where=f"w1mpro < {faint_limit:.1f}",
    )
    ok, funnel, dist, pm = star_candidates(field, center_ra, center_dec)
    w1, w2, var = field["w1mpro"], field["w2mpro"], field["var_flg"]

    rng = np.random.default_rng(config.SEED)
    chosen: list[int] = []
    for lo, hi in MAG_BINS:
        idx = np.flatnonzero(ok & (w1 >= lo) & (w1 < hi))
        idx = idx[np.argsort(field["designation"][idx])]  # deterministic order
        take = min(STARS_PER_BIN, idx.size)
        picked = rng.choice(idx, size=take, replace=False) if take else np.array([], int)
        funnel.append((f"bin {lo}-{hi}: candidates / drawn", f"{idx.size} / {take}"))
        chosen.extend(sorted(int(i) for i in picked))

    stars = []
    for rank, i in enumerate(chosen, start=1):
        stars.append({
            "name": f"ref{rank:02d}",
            "allwise_designation": str(field["designation"][i]),
            "ra": f"{field['ra'][i]:.7f}", "dec": f"{field['dec'][i]:.7f}",
            "w1mpro": f"{w1[i]:.3f}", "w2mpro": f"{w2[i]:.3f}",
            "w1rchi2": f"{field['w1rchi2'][i]:.2f}", "w2rchi2": f"{field['w2rchi2'][i]:.2f}",
            "var_flg": str(var[i]),
            "j_ks": f"{field['j_m_2mass'][i] - field['k_m_2mass'][i]:.3f}",
            "pm_mas_yr": f"{pm[i]:.0f}",
            "dist_to_ngc3786_deg": f"{dist[i] / 3600.0:.3f}",
        })
    with open(config.TARGETS_DIR / "reference_stars.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(stars[0].keys()))
        writer.writeheader()
        writer.writerows(stars)
    (config.TARGETS_DIR / "reference_stars_funnel.json").write_text(
        json.dumps({"seed": config.SEED, "funnel": funnel}, indent=1) + "\n"
    )
    for label, count in funnel:
        print(f"  {label}: {count}")
    return stars


def select_controls() -> list[dict]:
    rows = []
    for cand in CONTROL_CANDIDATES:
        ned = ned_lookup(cand)
        cat = allwise_cone(f"{slug(cand)}_cat", ned["ra"], ned["dec"], 10.0)
        row = {"name": slug(cand), "ned_name": ned["name"],
               "ned_ra": f"{ned['ra']:.6f}", "ned_dec": f"{ned['dec']:.6f}",
               "ned_z": ned["z"],
               "gal_lat_deg": f"{galactic_latitude_deg(ned['ra'], ned['dec']):.1f}"}
        sep = sep_arcsec(cat["ra"], cat["dec"], ned["ra"], ned["dec"]) if len(cat["ra"]) else np.array([])
        near = np.flatnonzero(sep <= CONTROL_MATCH_ARCSEC)
        if near.size == 0:
            row.update({"selected": "no", "reason": "no AllWISE source within 3 arcsec"})
            rows.append(row)
            continue
        i = int(near[np.argmin(cat["w1mpro"][near])])
        w1, w2 = float(cat["w1mpro"][i]), float(cat["w2mpro"][i])
        reasons = []
        if not w1 > config.SATURATION_WARN_MAG["w1"]:
            reasons.append("C3 W1 too bright")
        if not w2 > config.SATURATION_WARN_MAG["w2"]:
            reasons.append("C3 W2 too bright")
        if abs(float(row["gal_lat_deg"])) <= 20.0:
            reasons.append("C4 low Galactic latitude")
        row.update({
            "allwise_designation": str(cat["designation"][i]),
            "ra": f"{cat['ra'][i]:.7f}", "dec": f"{cat['dec'][i]:.7f}",
            "sep_ned_allwise_arcsec": f"{sep[i]:.2f}",
            "w1mpro": f"{w1:.3f}", "w2mpro": f"{w2:.3f}",
            "ext_flg": int(cat["ext_flg"][i]), "cc_flags": str(cat["cc_flags"][i]),
            "var_flg": str(cat["var_flg"][i]),
            "selected": "no" if reasons else "yes", "reason": "; ".join(reasons),
        })
        rows.append(row)
        print(f"  {cand}: W1={w1:.2f} W2={w2:.2f} b={row['gal_lat_deg']} -> {row['selected']} {row['reason']}")
    fields = ["name", "ned_name", "ned_ra", "ned_dec", "ned_z", "gal_lat_deg",
              "allwise_designation", "ra", "dec", "sep_ned_allwise_arcsec", "w1mpro",
              "w2mpro", "ext_flg", "cc_flags", "var_flg", "selected", "reason"]
    with open(config.TARGETS_DIR / "control_agn.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return rows


def main() -> None:
    config.TARGETS_DIR.mkdir(parents=True, exist_ok=True)
    print("NGC 3786")
    target = select_ngc3786()
    print(f"  {target['allwise_designation']} ra={target['ra']} dec={target['dec']}")
    print("Reference stars")
    stars = select_reference_stars(target["ra"], target["dec"])
    print(f"  selected {len(stars)}")
    print("Control AGN")
    controls = select_controls()
    print(f"  selected {sum(r['selected'] == 'yes' for r in controls)} of {len(controls)}")


if __name__ == "__main__":
    main()
