"""Stage 6A: other-wavelength light curves of the flare candidates.

ASAS-SN Sky Patrol v2 (V and g, large-aperture photometry that includes most of
the host) and the Swift-BAT 157-month survey light curves (14-195 keV,
2004-12 to 2017-12). Both are public; see docs/DATA.md.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from .config import CACHE_DIR

S6_DIR = CACHE_DIR / "s6"
JD_MINUS_MJD = 2400000.5
_FITS_TYPES = {"D": ">f8", "E": ">f4", "J": ">i4", "I": ">i2", "K": ">i8", "B": "u1", "L": "S1"}


def read_asassn(path: Path, max_mag_err: float = 0.1) -> dict[str, dict[str, np.ndarray]]:
    """Good-quality detections per filter ('V', 'g'): quality 'G', a measured
    magnitude brighter than the limit, and magnitude error below ``max_mag_err``."""
    out: dict[str, dict[str, list]] = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            try:
                mag, err, lim = float(r["mag"]), float(r["mag_err"]), float(r["limit"])
            except ValueError:
                continue
            if r["quality"] != "G" or not (0 < err <= max_mag_err) or not mag < lim:
                continue
            d = out.setdefault(r["phot_filter"], {"mjd": [], "mag": [], "err": [], "camera": []})
            d["mjd"].append(float(r["jd"]) - JD_MINUS_MJD)
            d["mag"].append(mag)
            d["err"].append(err)
            d["camera"].append(r["camera"])
    return {f: {k: np.array(v) for k, v in d.items()} for f, d in out.items()}


def read_fits_table(path: Path, hdu: int = 1) -> dict[str, np.ndarray]:
    """Minimal reader for a FITS binary-table extension (scalar and fixed-length
    vector columns of types D, E, J, I, K, B, A). Enough for the BAT survey
    light-curve files; not a general FITS reader."""
    raw = Path(path).read_bytes()
    pos, found = 0, -1
    while pos < len(raw):
        cards = {}
        while True:
            block = raw[pos:pos + 2880].decode("ascii", errors="replace")
            pos += 2880
            done = False
            for i in range(0, 2880, 80):
                card = block[i:i + 80]
                key = card[:8].strip()
                if key == "END":
                    done = True
                    break
                if card[8:10] == "= ":
                    cards[key] = card[10:].split("/")[0].strip().strip("'").strip()
            if done:
                break
        found += 1
        naxis = int(cards.get("NAXIS", 0))
        size = 0
        if naxis:
            size = abs(int(cards["BITPIX"])) // 8
            for k in range(1, naxis + 1):
                size *= int(cards[f"NAXIS{k}"])
            size += int(cards.get("PCOUNT", 0))
        padded = (size + 2879) // 2880 * 2880
        if found == hdu:
            dtype = []
            for k in range(1, int(cards["TFIELDS"]) + 1):
                form = cards[f"TFORM{k}"]
                code = form.lstrip("0123456789")[0]
                repeat = int(form[: len(form) - len(form.lstrip("0123456789"))] or 1)
                name = cards.get(f"TTYPE{k}", f"col{k}")
                if code == "A":
                    dtype.append((name, f"S{repeat}"))
                else:
                    dtype.append((name, _FITS_TYPES[code], (repeat,)) if repeat > 1
                                 else (name, _FITS_TYPES[code]))
            table = np.frombuffer(raw[pos:pos + int(cards["NAXIS1"]) * int(cards["NAXIS2"])],
                                  dtype=np.dtype(dtype))
            out = {name: np.array(table[name]) for name in table.dtype.names}
            out["_header"] = cards
            return out
        pos += padded
    raise ValueError(f"HDU {hdu} not found in {path}")


def weighted_mean(values, errors) -> tuple[float, float]:
    v, e = np.asarray(values, dtype=float), np.asarray(errors, dtype=float)
    ok = np.isfinite(v) & np.isfinite(e) & (e > 0)
    if not ok.any():
        return float("nan"), float("nan")
    w = 1.0 / e[ok] ** 2
    return float(np.sum(w * v[ok]) / w.sum()), float(w.sum() ** -0.5)


SWIFT_MET_MJD0 = 51910.0  # Swift mission elapsed time zero: 2001-01-01 00:00 UTC
OPT_MIN_BINS = 6
OPT_BIN_ERR_FLOOR = 0.02  # mag, as in stage 5b (an assumed value)


def optical_change(v: dict, mjd, mag, half_width: float = 91.3, min_points: int = 5) -> dict:
    """Optical change between the IR-faint and IR-bright thirds of the WISE
    visits that have an optical bin. Same definition as stage 5b Q3
    (stage5b.ztf_band), but without the requirement on the size of the IR change."""
    from .epochs import MEDIAN_SE_FACTOR
    from .stats import robust_std
    from .ztf import bin_to_visits

    bins = bin_to_visits(mjd, mag, v["mjd"], half_width, min_points)
    ok = np.flatnonzero(np.isfinite(bins["mag"]))
    out = {"n_points": int(np.asarray(mjd).size), "n_bins": int(ok.size),
           "first_bin_mjd": float(v["mjd"][ok[0]]) if ok.size else None,
           "last_bin_mjd": float(v["mjd"][ok[-1]]) if ok.size else None}
    if ok.size < OPT_MIN_BINS:
        out["status"] = "fewer than 6 visits with an optical bin"
        return out
    mean_ir = 0.5 * (v["w1_mag"][ok] + v["w2_mag"][ok])
    order = np.argsort(mean_ir, kind="stable")
    n3 = max(2, ok.size // 3)
    bright, faint = ok[order[:n3]], ok[order[-n3:]]

    def se(idx):
        m = bins["mag"][idx]
        err = np.hypot(bins["err"][idx], OPT_BIN_ERR_FLOOR)
        scatter = MEDIAN_SE_FACTOR * robust_std(m) / np.sqrt(m.size) if m.size > 1 else 0.0
        return float(max(scatter, np.sqrt(np.sum(err**2)) / m.size))

    dm = float(np.median(bins["mag"][faint]) - np.median(bins["mag"][bright]))
    err = float(np.hypot(se(faint), se(bright)))
    out.update(status="usable",
               ir_change_w2=float(np.median(v["w2_mag"][faint]) - np.median(v["w2_mag"][bright])),
               dm_opt=dm, dm_opt_err=err, upper_limit_3sigma=float(max(dm, 0.0) + 3.0 * err),
               bins_mjd=[float(x) for x in v["mjd"][ok]],
               bins_mag=[float(x) for x in bins["mag"][ok]])
    return out


def bat_yearly(table: dict) -> list[dict]:
    """Exposure-independent weighted mean of the monthly 14-195 keV rates per
    calendar year."""
    mjd = SWIFT_MET_MJD0 + table["TIME"] / 86400.0
    year = np.floor(2000.0 + (mjd - 51544.5) / 365.25).astype(int)
    rows = []
    for y in sorted(set(year.tolist())):
        sel = (year == y) & (table["RATE_ERR"] > 0)
        mean, err = weighted_mean(table["RATE"][sel], table["RATE_ERR"][sel])
        rows.append({"year": int(y), "n_months": int(sel.sum()), "rate": mean, "rate_err": err})
    return rows
