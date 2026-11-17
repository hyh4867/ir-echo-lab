"""Independent cross-check of the NGC 3786 flare amplitude.

This script deliberately shares no code with src/irecho. It reads the cached
NEOWISE CSV with the standard library only and uses a cruder recipe:

  * no quality cuts (only rows with a measured magnitude and uncertainty),
  * visits split wherever two exposures are more than 100 days apart,
  * plain mean of magnitudes per visit (the pipeline uses a flux median
    after quality cuts and a 5-sigma clip),
  * baseline = mean of the visit means before 2020-06-01.

If the pipeline's amplitude came from its cuts, clipping or flux conversion,
this recipe would give a different number.

    python scripts/independent_check.py
"""

from __future__ import annotations

import csv
import json
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "data" / "cache" / "neowise" / "ngc3786.csv"
OUT_PATH = ROOT / "results" / "independent_check.json"
BASELINE_END_MJD = 59001.0   # 2020-06-01
PAPER_END_MJD = 59215.0      # 2021-01-01
GAP_DAYS = 100.0


def read_rows(path: Path) -> list[dict]:
    rows = []
    with path.open(newline="") as fh:
        for r in csv.DictReader(fh):
            try:
                row = {k: float(r[k]) for k in ("mjd", "w1mpro", "w1sigmpro", "w2mpro", "w2sigmpro")}
            except ValueError:
                continue  # a band without a measured magnitude or uncertainty
            rows.append(row)
    return sorted(rows, key=lambda r: r["mjd"])


def split_visits(rows: list[dict]) -> list[list[dict]]:
    visits: list[list[dict]] = [[rows[0]]]
    for prev, cur in zip(rows, rows[1:]):
        if cur["mjd"] - prev["mjd"] > GAP_DAYS:
            visits.append([])
        visits[-1].append(cur)
    return visits


def blackbody_check(base_w1: float, base_w2: float, flare_w1: float, flare_w2: float) -> dict:
    """Stage-3 quantities for the first flare visit, again without src/irecho:
    own constants, bisection instead of Brent, Hubble-law distance instead of
    the LCDM integral, no redshift terms, no differential correction."""
    import math

    h, c, k, sigma = 6.62607015e-34, 2.99792458e8, 1.380649e-23, 5.670374419e-8
    f0 = {"w1": 309.540, "w2": 171.787}  # Jy
    lam = {"w1": 3.3526e-6, "w2": 4.6028e-6}  # m

    def jy(mag: float, band: str) -> float:
        return f0[band] * 10 ** (-0.4 * mag)

    def planck(band: str, temp: float) -> float:
        nu = c / lam[band]
        return 2 * h * nu**3 / c**2 / math.expm1(h * nu / (k * temp))

    d1 = jy(flare_w1, "w1") - jy(base_w1, "w1")
    d2 = jy(flare_w2, "w2") - jy(base_w2, "w2")
    lo, hi = 300.0, 5000.0
    for _ in range(60):  # the flux ratio rises with temperature
        mid = 0.5 * (lo + hi)
        if planck("w1", mid) / planck("w2", mid) < d1 / d2:
            lo = mid
        else:
            hi = mid
    temp = 0.5 * (lo + hi)
    dist_m = (0.009083 * c / 1e3 / 67.4) * 3.0856775814913673e22  # cz / H0
    solid = d2 * 1e-26 / planck("w2", temp)
    lum_w = 4 * math.pi * dist_m**2 * solid * sigma * temp**4 / math.pi
    radius_m = math.sqrt(lum_w / (4 * math.pi * sigma * temp**4))
    return {
        "recipe": "mean magnitudes, no cuts, no differential correction, cz/H0 distance, "
                  "no redshift terms, bisection",
        "excess_w1_mjy": round(1e3 * d1, 3), "excess_w2_mjy": round(1e3 * d2, 3),
        "temperature_k": round(temp, 1),
        "l_dust_erg_s": float(f"{lum_w * 1e7:.4g}"),
        "r_bb_light_days": round(radius_m / (c * 86400.0), 2),
        "distance_mpc": round(dist_m / 3.0856775814913673e22, 2),
    }


def main() -> None:
    rows = read_rows(CSV_PATH)
    visits = split_visits(rows)
    table = [{
        "mjd": statistics.fmean(r["mjd"] for r in v),
        "n": len(v),
        "w1": statistics.fmean(r["w1mpro"] for r in v),
        "w2": statistics.fmean(r["w2mpro"] for r in v),
    } for v in visits]
    base = [t for t in table if t["mjd"] < BASELINE_END_MJD]
    flare = [t for t in table if BASELINE_END_MJD <= t["mjd"] < PAPER_END_MJD]
    assert len(flare) == 1, "expected exactly one visit between 2020-06-01 and 2021-01-01"
    f = flare[0]
    b1 = statistics.fmean(t["w1"] for t in base)
    b2 = statistics.fmean(t["w2"] for t in base)
    out = {
        "source_file": str(CSV_PATH.relative_to(ROOT)),
        "recipe": "no quality cuts, visits split at gaps > 100 d, mean magnitude per visit, "
                  "baseline = mean of visit means before MJD 59001 (2020-06-01)",
        "n_rows_used": len(rows),
        "n_visits": len(table),
        "n_baseline_visits": len(base),
        "baseline_w1_mag": round(b1, 4),
        "baseline_w2_mag": round(b2, 4),
        "baseline_w1_std_mag": round(statistics.stdev(t["w1"] for t in base), 4),
        "baseline_w2_std_mag": round(statistics.stdev(t["w2"] for t in base), 4),
        "flare_visit_mjd": round(f["mjd"], 3),
        "flare_visit_n_rows": f["n"],
        "flare_w1_mag": round(f["w1"], 4),
        "flare_w2_mag": round(f["w2"], 4),
        "w1_dmag": round(b1 - f["w1"], 4),
        "w2_dmag": round(b2 - f["w2"], 4),
        "color_baseline_mag": round(b1 - b2, 4),
        "color_flare_mag": round(f["w1"] - f["w2"], 4),
    }
    out["delta_color_mag"] = round(out["color_flare_mag"] - out["color_baseline_mag"], 4)
    out["physics"] = blackbody_check(b1, b2, f["w1"], f["w2"])
    OUT_PATH.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
