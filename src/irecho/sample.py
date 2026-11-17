"""Stage 5a: definition of the nearby hard-X-ray AGN sample (docs/SAMPLE.md).

Everything here is a pure function of catalogue rows, so the selection can be
unit-tested without network access. scripts/build_sample.py does the
downloads and writes targets/sample_v1*.

The criteria S1-S11 below were committed before the script was run on the
catalogues (docs/SAMPLE.md section 1 records what had been looked at by then).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .config import SATURATION_WARN_MAG

# --- parent catalogues --------------------------------------------------------
BAT105_URL = ("https://swift.gsfc.nasa.gov/bat_survey/bs105mon/data/"
              "BAT_105m_catalog_07jul2019.txt")
BASS_DR2_URL = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv"
BASS_DR2_PARAMS = {"-source": "J/ApJS/261/2/table9", "-out.max": "unlimited", "-out.all": "1"}

# Objects that were analysed in stages 1-4. They stay in the list whatever the
# cuts say, flagged role = 'calibration', and are never counted as discoveries.
CALIBRATION_COUNTERPARTS = ("NGC 3786", "NGC 2273")


@dataclass(frozen=True)
class SampleCriteria:
    """S1-S11. See docs/SAMPLE.md section 2 for the reason behind each value."""

    # S1 BAT 105-month class codes: 40 Seyfert I, 50 Seyfert II, 60 LINER.
    bat_classes: tuple[int, ...] = (40, 50, 60)
    # S2 redshift: 0 < z <= z_max, z_max chosen by the rule in choose_z_max().
    z_grid: tuple[float, ...] = (0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040, 0.045, 0.050)
    max_sample_size: int = 250
    min_sample_size: int = 150  # reported if not reached; no further action
    # S3 Galactic latitude.
    min_abs_gal_lat_deg: float = 20.0
    # S4 ecliptic latitude (near the ecliptic poles WISE coverage is nearly
    # continuous and a 'visit' is no longer a 1-2 day group of exposures).
    max_abs_ecl_lat_deg: float = 80.0
    # S5 BASS DR2 types that are removed when a DR2 match exists.
    dr2_match_arcsec: float = 5.0
    dr2_excluded_types: tuple[str, ...] = ("BZQ", "BZB", "BZG", "BZQ/Lense", "Sy1/Lense")
    # S6 AllWISE counterpart.
    allwise_match_arcsec: float = 3.0
    # S7 saturation margin (config.SATURATION_WARN_MAG, same as control cut C3).
    w1_bright_limit: float = SATURATION_WARN_MAG["w1"]
    w2_bright_limit: float = SATURATION_WARN_MAG["w2"]
    # S8 faint limit: the error model was checked on stars with W1 9.5-12.5.
    w1_faint_limit: float = 12.5
    # S9 AllWISE contamination flags in W1 and W2 must be '00'; ph_qual 'AA'.
    # S10 neighbours.
    bright_star_w1: float = 5.0
    bright_star_radius_arcsec: float = 300.0
    blend_radius_arcsec: float = 12.0
    blend_dmag: float = 2.5
    # S11 comparison stars (criteria R1-R13 of scripts/select_targets.py).
    min_comparison_stars: int = 6


CRITERIA = SampleCriteria()


# --- coordinates --------------------------------------------------------------
def galactic_latitude_deg(ra_deg: float, dec_deg: float) -> float:
    """Galactic latitude from ICRS (J2000 pole: RA 192.85948, Dec 27.12825)."""
    ra, dec = math.radians(ra_deg), math.radians(dec_deg)
    ra_p, dec_p = math.radians(192.85948), math.radians(27.12825)
    sinb = (math.sin(dec) * math.sin(dec_p)
            + math.cos(dec) * math.cos(dec_p) * math.cos(ra - ra_p))
    return math.degrees(math.asin(max(-1.0, min(1.0, sinb))))


def ecliptic_latitude_deg(ra_deg: float, dec_deg: float) -> float:
    """Ecliptic latitude for the J2000 mean obliquity (23.4392911 deg)."""
    ra, dec, eps = math.radians(ra_deg), math.radians(dec_deg), math.radians(23.4392911)
    sinb = math.sin(dec) * math.cos(eps) - math.cos(dec) * math.sin(eps) * math.sin(ra)
    return math.degrees(math.asin(max(-1.0, min(1.0, sinb))))


def sep_arcsec(ra1, dec1, ra2, dec2):
    """Angular separation in arcsec (haversine; valid at any separation)."""
    ra1, dec1, ra2, dec2 = (np.radians(np.asarray(x, dtype=float)) for x in (ra1, dec1, ra2, dec2))
    a = (np.sin(0.5 * (dec2 - dec1)) ** 2
         + np.cos(dec1) * np.cos(dec2) * np.sin(0.5 * (ra2 - ra1)) ** 2)
    return np.degrees(2.0 * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0)))) * 3600.0


# --- catalogue parsers --------------------------------------------------------
def _float_or_nan(text: str) -> float:
    text = text.strip()
    try:
        return float(text) if text else float("nan")
    except ValueError:
        return float("nan")


def parse_bat105(text: str) -> list[dict]:
    """Rows of the pipe-delimited BAT 105-month table (Oh et al. 2018)."""
    rows = []
    for line in text.splitlines():
        parts = line.split("|")
        if len(parts) < 22 or not parts[0].strip().isdigit():
            continue
        rows.append({
            "bat_index": int(parts[0]),
            "bat_name": parts[1].strip(),
            "counterpart": parts[5].strip(),
            "ctpt_ra": _float_or_nan(parts[7]),
            "ctpt_dec": _float_or_nan(parts[8]),
            "z": _float_or_nan(parts[17]),
            "bat_class": int(parts[20]) if parts[20].strip().isdigit() else -1,
            "bat_type": parts[21].strip(),
        })
    return rows


def parse_vizier_tsv(text: str) -> list[dict]:
    """Rows of one VizieR ``asu-tsv`` table (header, units, dashes, data)."""
    lines = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
    if len(lines) < 3:
        return []
    header = lines[0].split("\t")
    start = next(i for i, ln in enumerate(lines) if ln.startswith("---")) + 1
    return [dict(zip(header, (c.strip() for c in ln.split("\t")))) for ln in lines[start:]]


def match_dr2(row: dict, dr2: list[dict], radius_arcsec: float) -> dict | None:
    """Nearest BASS DR2 entry within ``radius_arcsec`` of the BAT counterpart."""
    if not dr2:
        return None
    ra = np.array([float(d["RAJ2000"]) for d in dr2])
    dec = np.array([float(d["DEJ2000"]) for d in dr2])
    sep = sep_arcsec(ra, dec, row["ctpt_ra"], row["ctpt_dec"])
    i = int(np.argmin(sep))
    if sep[i] > radius_arcsec:
        return None
    return {**dr2[i], "sep_arcsec": float(sep[i])}


# --- cuts ---------------------------------------------------------------------
def catalogue_cut(row: dict, dr2_match: dict | None, c: SampleCriteria = CRITERIA) -> str:
    """First failed criterion among S1-S5 ('' when all pass). S2 is evaluated
    against the top of the redshift grid; the final z_max is applied later."""
    if row["bat_class"] not in c.bat_classes:
        return "S1 class"
    if not (np.isfinite(row["z"]) and 0.0 < row["z"] <= max(c.z_grid)):
        return "S2 redshift"
    if not abs(galactic_latitude_deg(row["ctpt_ra"], row["ctpt_dec"])) > c.min_abs_gal_lat_deg:
        return "S3 Galactic latitude"
    if not abs(ecliptic_latitude_deg(row["ctpt_ra"], row["ctpt_dec"])) < c.max_abs_ecl_lat_deg:
        return "S4 ecliptic latitude"
    if dr2_match is not None and dr2_match.get("Type", "") in c.dr2_excluded_types:
        return "S5 BASS DR2 beamed/lensed"
    return ""


@dataclass(frozen=True)
class SampleV2Range:
    """Stage 6B: the independent sample keeps every criterion of sample_v1
    except the redshift range, and excludes the sample_v1 objects."""

    z_min_exclusive: float = 0.04
    z_max: float = 0.30


V2 = SampleV2Range()


def catalogue_cut_v2(row: dict, dr2_match: dict | None, v1_indices: set[int],
                     c: SampleCriteria = CRITERIA, r: SampleV2Range = V2) -> str:
    """S1-S5 for sample_v2: as catalogue_cut, with 0.04 < z <= 0.30 and the
    sample_v1 objects removed (E0)."""
    if row["bat_class"] not in c.bat_classes:
        return "S1 class"
    if not (np.isfinite(row["z"]) and r.z_min_exclusive < row["z"] <= r.z_max):
        return "S2 redshift"
    if row["bat_index"] in v1_indices:
        return "E0 already in sample_v1"
    if not abs(galactic_latitude_deg(row["ctpt_ra"], row["ctpt_dec"])) > c.min_abs_gal_lat_deg:
        return "S3 Galactic latitude"
    if not abs(ecliptic_latitude_deg(row["ctpt_ra"], row["ctpt_dec"])) < c.max_abs_ecl_lat_deg:
        return "S4 ecliptic latitude"
    if dr2_match is not None and dr2_match.get("Type", "") in c.dr2_excluded_types:
        return "S5 BASS DR2 beamed/lensed"
    return ""


def allwise_cut(field: dict[str, np.ndarray], ra: float, dec: float,
                c: SampleCriteria = CRITERIA) -> tuple[str, int]:
    """First failed criterion among S6-S10 for an AllWISE cone table around a
    candidate, and the row index of the adopted counterpart (-1 if none)."""
    if len(field["ra"]) == 0:
        return "S6 no AllWISE counterpart", -1
    sep = sep_arcsec(field["ra"], field["dec"], ra, dec)
    w1, w2 = field["w1mpro"], field["w2mpro"]
    near = np.flatnonzero((sep <= c.allwise_match_arcsec) & np.isfinite(w1))
    if near.size == 0:
        return "S6 no AllWISE counterpart", -1
    i = int(near[np.argmin(w1[near])])  # brightest W1 inside the match radius
    if not (w1[i] > c.w1_bright_limit and np.isfinite(w2[i]) and w2[i] > c.w2_bright_limit):
        return "S7 saturation margin", i
    if not w1[i] <= c.w1_faint_limit:
        return "S8 faint limit", i
    if str(field["cc_flags"][i])[:2] != "00" or str(field["ph_qual"][i])[:2] != "AA":
        return "S9 cc_flags / ph_qual", i
    d = sep_arcsec(field["ra"], field["dec"], field["ra"][i], field["dec"][i])
    others = np.arange(len(d)) != i
    if (others & (d <= c.bright_star_radius_arcsec) & (w1 < c.bright_star_w1)).any():
        return "S10 bright star nearby", i
    if (others & (d <= c.blend_radius_arcsec) & (w1 < w1[i] + c.blend_dmag)).any():
        return "S10 blended neighbour", i
    return "", i


def choose_z_max(redshifts, c: SampleCriteria = CRITERIA) -> tuple[float, dict[str, int]]:
    """S2 rule: the largest grid value whose count does not exceed
    ``max_sample_size``. ``redshifts`` are the objects that passed every other
    target-level cut. Returns (z_max, {grid value: count}). If even the
    smallest grid value gives more than the maximum, the smallest is used."""
    z = np.asarray(list(redshifts), dtype=float)
    counts = {f"{g:.3f}": int((z <= g).sum()) for g in c.z_grid}
    allowed = [g for g in c.z_grid if counts[f"{g:.3f}"] <= c.max_sample_size]
    return (max(allowed) if allowed else min(c.z_grid)), counts
