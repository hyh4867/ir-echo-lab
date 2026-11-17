"""Quality cuts on single-exposure measurements.

The cuts follow the NEOWISE Explanatory Supplement II.3.a/II.3.b and the
AllWISE Explanatory Supplement III.2.b. docs/PIPELINE.md quotes the passages.
"""

from __future__ import annotations

import numpy as np

from .config import QUALITY, QualityCuts

BAND_INDEX = {"w1": 0, "w2": 1}


def _flag_char(values: np.ndarray, index: int) -> np.ndarray:
    """Character ``index`` of each flag string ('?' when the string is short)."""
    return np.array([v[index] if len(v) > index else "?" for v in values], dtype=str)


def normalise(table: dict[str, np.ndarray], survey: str) -> dict[str, np.ndarray]:
    """Map the two IRSA schemas onto common column names.

    ``survey`` is 'neowise' (neowiser_p1bs_psd) or 'allwise' (allwise_p3as_mep).
    Columns that only exist in NEOWISE are filled with neutral values for
    AllWISE so that the same mask code runs on both.
    """
    n = len(table["mjd"])
    if survey == "neowise":
        out = {
            "mjd": table["mjd"],
            "ra": table["ra"],
            "dec": table["dec"],
            "frame": np.char.add(table["scan_id"], np.char.mod("%03d", table["frame_num"].astype(int))),
            "w1mag": table["w1mpro"], "w1err": table["w1sigmpro"], "w1rchi2": table["w1rchi2"],
            "w2mag": table["w2mpro"], "w2err": table["w2sigmpro"], "w2rchi2": table["w2rchi2"],
            "w1sat": table["w1sat"], "w2sat": table["w2sat"],
            "qual_frame": table["qual_frame"],
            "sso_flg": table["sso_flg"],
        }
    elif survey == "allwise":
        out = {
            "mjd": table["mjd"],
            "ra": table["ra"],
            "dec": table["dec"],
            "frame": table["frame_id"],
            "w1mag": table["w1mpro_ep"], "w1err": table["w1sigmpro_ep"],
            "w1rchi2": table["w1rchi2_ep"],
            "w2mag": table["w2mpro_ep"], "w2err": table["w2sigmpro_ep"],
            "w2rchi2": table["w2rchi2_ep"],
            "w1sat": np.zeros(n), "w2sat": np.zeros(n),
            "qual_frame": np.full(n, np.nan),
            "sso_flg": np.zeros(n),
        }
    else:
        raise ValueError(f"unknown survey {survey!r}")
    out["qi_fact"] = table["qi_fact"]
    out["saa_sep"] = table["saa_sep"]
    out["cc_flags"] = table["cc_flags"]
    out["moon_masked"] = table["moon_masked"]
    out["survey"] = np.full(n, survey)
    return out


def band_mask(obs: dict[str, np.ndarray], band: str,
              cuts: QualityCuts = QUALITY) -> tuple[np.ndarray, dict[str, int]]:
    """Boolean mask of usable measurements in one band.

    Returns the mask and, for each criterion, how many rows fail it when it is
    applied on its own (rows can fail several criteria).
    """
    i = BAND_INDEX[band]
    survey = obs["survey"]
    is_neo = survey == "neowise"
    mag, err = obs[f"{band}mag"], obs[f"{band}err"]

    tests: dict[str, np.ndarray] = {}
    # A null uncertainty marks an upper limit (SNR < 2), not a measurement.
    tests["measured"] = np.isfinite(mag) & np.isfinite(err) & (err > 0)
    tests["qual_frame"] = ~is_neo | (obs["qual_frame"] > cuts.min_qual_frame)
    tests["qi_fact"] = np.where(
        is_neo,
        obs["qi_fact"] > cuts.min_qi_fact_neowise,
        obs["qi_fact"] >= cuts.require_qi_fact_allwise,
    )
    tests["saa_sep"] = obs["saa_sep"] > cuts.min_saa_sep_deg
    if cuts.require_moon_unmasked:
        tests["moon_masked"] = _flag_char(obs["moon_masked"], i) == "0"
    if cuts.require_cc_clean:
        tests["cc_flags"] = _flag_char(obs["cc_flags"], i) == "0"
    if cuts.require_unsaturated:
        sat = obs[f"{band}sat"]
        tests["saturation"] = ~(np.isfinite(sat) & (sat > 0))
    if cuts.reject_sso:
        tests["sso_flg"] = ~(obs["sso_flg"] > 0)

    mask = np.ones(len(mag), dtype=bool)
    failed: dict[str, int] = {}
    for name, ok in tests.items():
        failed[name] = int((~ok).sum())
        mask &= ok
    failed["total_rows"] = int(len(mag))
    failed["kept"] = int(mask.sum())
    return mask, failed


def nearest_per_frame(obs: dict[str, np.ndarray], ra_deg: float,
                      dec_deg: float) -> np.ndarray:
    """Keep one detection per frame: the one closest to the target position.

    A frame can hold two detections inside the cone when the source is split by
    active deblending.
    """
    cosd = np.cos(np.deg2rad(dec_deg))
    sep = np.hypot((obs["ra"] - ra_deg) * cosd, obs["dec"] - dec_deg)
    keep = np.zeros(len(sep), dtype=bool)
    best: dict[str, int] = {}
    for idx, frame in enumerate(obs["frame"]):
        key = f"{obs['survey'][idx]}:{frame}"
        if key not in best or sep[idx] < sep[best[key]]:
            best[key] = idx
    keep[list(best.values())] = True
    return keep
