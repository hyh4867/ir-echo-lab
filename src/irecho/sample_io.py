"""Stage 5a: file locations and loaders for the sample_v1 fields.

One field = one sample object (position id 0) plus its comparison stars
(position ids 1..n, in the order of targets/sample_v1_stars.csv). The
single-exposure tables of a field are stored in one CSV per survey, written by
scripts/fetch_sample.py with a ``pid`` column.
"""

from __future__ import annotations

import os
import csv
from pathlib import Path

from . import config, irsa
from .lightcurve import LightCurve, build_lightcurve

# Which sample the stage-5a pipeline works on. "v1" is the 231-object sample of
# stage 5a; "v2" is the independent sample of stage 6B (IRECHO_SAMPLE=v2). The
# tag only changes file locations; every rule is the same.
SAMPLE_TAG = os.environ.get("IRECHO_SAMPLE", "v1")
if SAMPLE_TAG not in ("v1", "v2"):
    raise ValueError(f"unknown IRECHO_SAMPLE {SAMPLE_TAG!r}")
CACHE_KEY = {"v1": "s5", "v2": "s6b"}[SAMPLE_TAG]  # sub-directory of data/cache/
QA_PREFIX = {"v1": "qa_", "v2": "qa2_"}[SAMPLE_TAG]  # prefix of the results/ files
S5_DIR = config.CACHE_DIR / CACHE_KEY
SAMPLE_PATH = config.TARGETS_DIR / f"sample_{SAMPLE_TAG}.csv"
STARS_PATH = config.TARGETS_DIR / f"sample_{SAMPLE_TAG}_stars.csv"
FROZEN_DIR = config.PROJECT_ROOT / "data" / "frozen" / f"sample_{SAMPLE_TAG}"
FREEZE_PATH = config.PROJECT_ROOT / "data" / f"freeze_{SAMPLE_TAG}.json"
SURVEY_TABLES = {
    "neowise": (config.TABLE_NEOWISE, irsa.NEOWISE_COLUMNS, ("pid", "mjd", "source_id")),
    "allwise": (config.TABLE_ALLWISE_MEP, irsa.ALLWISE_MEP_COLUMNS, ("pid", "mjd", "frame_id")),
}


def load_sample() -> list[dict]:
    with open(SAMPLE_PATH, newline="") as fh:
        return [{**r, "ra": float(r["ra"]), "dec": float(r["dec"])} for r in csv.DictReader(fh)]


def load_stars() -> dict[str, list[dict]]:
    """Comparison stars grouped by field name, in file order."""
    out: dict[str, list[dict]] = {}
    with open(STARS_PATH, newline="") as fh:
        for r in csv.DictReader(fh):
            out.setdefault(r["field"], []).append(
                {**r, "ra": float(r["ra"]), "dec": float(r["dec"])})
    return out


def field_positions(target: dict, stars: list[dict]) -> list[tuple[int, float, float]]:
    """(pid, ra, dec): the target is 0, the stars are 1..n."""
    return [(0, target["ra"], target["dec"])] + [
        (k, s["ra"], s["dec"]) for k, s in enumerate(stars, start=1)]


def field_paths(name: str) -> dict[str, Path]:
    return {"neowise": S5_DIR / "neowise" / f"{name}.csv",
            "allwise": S5_DIR / "allwise_mep" / f"{name}.csv",
            "untimely": S5_DIR / "untimely" / f"{name}.csv"}


def fetch_field(target: dict, stars: list[dict], surveys=("neowise", "allwise"),
                refresh: bool = False) -> dict[str, Path]:
    """Download (or reuse) the single-exposure tables of one field: one TAP
    query per survey for the target and all of its comparison stars."""
    paths = field_paths(target["name"])
    for survey in surveys:
        table, columns, sort_keys = SURVEY_TABLES[survey]
        irsa.fetch_upload_csv(table, columns, field_positions(target, stars),
                              config.CONE_RADIUS_ARCSEC, paths[survey], sort_keys=sort_keys,
                              refresh=refresh, manifest_path=config.MANIFEST_PATH)
    return {s: paths[s] for s in surveys}


def field_is_cached(name: str, surveys=("neowise", "allwise")) -> bool:
    paths = field_paths(name)
    return all(paths[s].exists() for s in surveys)


def load_field(target: dict, stars: list[dict]) -> dict:
    """Light curves of one field from the cache: {'target': LightCurve or None,
    'stars': [LightCurve, ...], 'star_names': [...]}. An object without any
    usable single exposure is returned as None (target) or left out (stars)."""
    paths = field_paths(target["name"])
    parts = {}
    for survey in ("neowise", "allwise"):
        if not paths[survey].exists():
            raise FileNotFoundError(f"IRSA cache missing: {paths[survey]} "
                                    "(run scripts/fetch_sample.py)")
        parts[survey] = irsa.split_by_pid(irsa.read_csv_table(paths[survey]))

    def build(pid: int, name: str, ra: float, dec: float) -> LightCurve | None:
        tables = {s: parts[s][pid] for s in ("neowise", "allwise") if pid in parts[s]}
        if "neowise" not in tables:
            return None
        return build_lightcurve(name, ra, dec, tables)

    out = {"target": build(0, target["name"], target["ra"], target["dec"]),
           "stars": [], "star_names": [], "star_pids": [], "star_w1": []}
    for k, s in enumerate(stars, start=1):
        lc = build(k, s["name"], s["ra"], s["dec"])
        if lc is not None:
            out["stars"].append(lc)
            out["star_names"].append(s["name"])
            out["star_pids"].append(k)
            out["star_w1"].append(float(s["w1mpro"]))
    return out
