"""Regression tests for stage 5a on the cached sample_v1 data.

Skipped unless the stage-5a cache is present (scripts/fetch_sample.py and
scripts/fetch_untimely.py --sample). The tests never print or compare a target
magnitude: for targets they check counts, hashes and the blind contract of the
QA files. Expected numbers are the committed results/qa_*.json written by
scripts/run_qa.py and scripts/run_injection_stars.py from the same cache.
"""

from __future__ import annotations

import csv
import json
import subprocess

import numpy as np
import pytest

from irecho import (analysis, config, differential, errorfloor, irsa, qa, sample, sample_io,
                    untimely)
from irecho.analysis import visit_arrays
from irecho.config import BANDS

pytestmark = pytest.mark.cached

FREEZE_PATH = config.PROJECT_ROOT / "data" / "freeze_v1.json"
FROZEN_DIR = config.PROJECT_ROOT / "data" / "frozen" / "sample_v1"


def _cache_complete() -> bool:
    try:
        targets = sample_io.load_sample()
    except FileNotFoundError:
        return False
    return (FREEZE_PATH.exists() and len(targets) > 0
            and all(sample_io.field_is_cached(t["name"]) for t in targets))


if not _cache_complete():
    pytestmark = [pytest.mark.cached, pytest.mark.skip(
        reason="stage 5a cache absent: run scripts/fetch_sample.py and "
               "scripts/fetch_untimely.py --sample")]


def load_json(name: str) -> dict:
    return json.loads((config.RESULTS_DIR / name).read_text())


@pytest.fixture(scope="module")
def targets():
    return sample_io.load_sample()


@pytest.fixture(scope="module")
def stars():
    return sample_io.load_stars()


def test_sample_file_obeys_the_committed_criteria(targets, stars):
    funnel = json.loads((config.TARGETS_DIR / "sample_v1_funnel.json").read_text())
    c = sample.CRITERIA
    science = [t for t in targets if not t["failed_cut"]]
    assert len(science) == funnel["n_selected_passing_all_cuts"]
    assert c.min_sample_size <= len(science) <= c.max_sample_size
    assert sorted(t["counterpart"] for t in targets if t["role"] == "calibration") == [
        "NGC 2273", "NGC 3786"]
    for t in science:
        assert int(t["bat_class"]) in c.bat_classes
        assert 0.0 < float(t["z"]) <= funnel["z_max"]
        assert abs(float(t["gal_lat_deg"])) > c.min_abs_gal_lat_deg
        assert abs(float(t["ecl_lat_deg"])) < c.max_abs_ecl_lat_deg
        assert c.w1_bright_limit < float(t["w1mpro"]) <= c.w1_faint_limit
        assert float(t["w2mpro"]) > c.w2_bright_limit
        assert t["cc_flags"][:2] == "00"
        assert float(t["sep_ctpt_allwise_arcsec"]) <= c.allwise_match_arcsec
        assert t["dr2_type"] not in c.dr2_excluded_types
        assert int(t["n_comparison_stars"]) >= c.min_comparison_stars
        assert len(stars[t["name"]]) == int(t["n_comparison_stars"])
    names = [s["name"] for group in stars.values() for s in group]
    assert len(names) == len(set(names)) == funnel["n_stars"]
    for group in stars.values():
        for s in group:
            assert 9.5 <= float(s["w1mpro"]) < 11.5
            assert 3.0 / 60.0 < float(s["dist_to_target_deg"]) <= 1.5


def test_redshift_rule_reproduces_from_the_rejected_list(targets):
    """z_max is the largest grid value with at most 250 objects passing S1-S10."""
    funnel = json.loads((config.TARGETS_DIR / "sample_v1_funnel.json").read_text())
    with open(config.TARGETS_DIR / "sample_v1_rejected.csv", newline="") as fh:
        rejected = list(csv.DictReader(fh))
    z = [float(t["z"]) for t in targets if not t["failed_cut"]
         or t["failed_cut"] == "S11 comparison stars"]
    z += [float(r["z"]) for r in rejected
          if r["failed_cut"].startswith("S2 redshift rule") or r["failed_cut"].startswith("S11")]
    z_max, counts = sample.choose_z_max(z)
    assert z_max == funnel["z_max"]
    assert counts == funnel["count_by_z_grid_after_S1_S10"]


def test_cached_files_match_manifest_hashes(targets):
    manifest = json.loads(config.MANIFEST_PATH.read_text())
    keys = [k for k in manifest if k.startswith(("s5/", "allwise_cat/s5_", "catalogs/"))]
    assert len([k for k in keys if k.startswith("s5/neowise/")]) == len(targets)
    assert len([k for k in keys if k.startswith("s5/allwise_mep/")]) == len(targets)
    for key in keys:
        path = config.CACHE_DIR / key
        assert path.exists(), key
        assert irsa.sha256_file(path) == manifest[key]["sha256"], key


def test_upload_query_returns_the_same_rows_as_the_single_cone(targets):
    """NGC 3786 was downloaded with a single 3 arcsec cone in stage 2 and again
    as position 0 of its field in stage 5a."""
    t = next(t for t in targets if t["counterpart"] == "NGC 3786")
    old_target = analysis.load_ngc3786()
    assert (t["ra"], t["dec"]) == pytest.approx((old_target["ra"], old_target["dec"]), abs=1e-7)
    for survey, key in (("neowise", "source_id"), ("allwise", "frame_id")):
        new = irsa.split_by_pid(irsa.read_csv_table(sample_io.field_paths(t["name"])[survey]))[0]
        old = irsa.read_csv_table(irsa.cache_paths("ngc3786")[survey])
        assert sorted(new[key]) == sorted(old[key])
        mag = "w1mpro" if survey == "neowise" else "w1mpro_ep"
        a = new[mag][np.argsort(new[key], kind="stable")]
        b = old[mag][np.argsort(old[key], kind="stable")]
        assert np.array_equal(a, b, equal_nan=True)


def _leaves(obj, path=""):
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from _leaves(v, f"{path}/{k}")
    else:
        yield path, obj


def test_qa_targets_file_is_blind(targets):
    """No per-visit series and no target magnitude in results/qa_targets.json."""
    doc = load_json("qa_targets.json")
    assert set(doc) == {t["name"] for t in targets}
    allowed_mag_keys = ("median_pipeline_minus_untimely_mag",)
    for name, row in doc.items():
        for path, value in _leaves(row):
            if path == "/flags":
                assert all(isinstance(f, str) for f in value)
                continue
            assert not isinstance(value, (list, tuple)), f"{name}{path} is a series"
            leaf = path.rsplit("/", 1)[-1]
            if "mag" in leaf and "mmag" not in leaf:
                assert leaf in allowed_mag_keys, f"{name}{path}"
    # the other committed stage-5a result files do not name any target in a series either
    for row in csv.DictReader(open(config.RESULTS_DIR / "qa_stars.csv")):
        assert "_s" in row["star"]


def test_frozen_light_curves_are_not_in_git_and_match_the_freeze(targets):
    tracked = subprocess.run(["git", "ls-files", "data/frozen", "data/cache"], capture_output=True,
                             text=True, cwd=config.PROJECT_ROOT).stdout.strip()
    assert tracked == ""
    freeze = json.loads(FREEZE_PATH.read_text())
    qa_targets = load_json("qa_targets.json")
    ok = [n for n, r in qa_targets.items() if r["status"] == "ok"]
    assert sorted(freeze["light_curves"]) == sorted(ok)
    for name, rec in freeze["light_curves"].items():
        path = FROZEN_DIR / f"{name}.csv"
        assert path.exists(), name
        assert irsa.sha256_file(path) == rec["sha256"], name
        assert rec["n_visits"] == qa_targets[name]["n_visits_corrected"]


@pytest.mark.parametrize("which", ["first", "ngc3786", "last"])
def test_field_qa_and_frozen_file_rebuild_from_cache(targets, stars, which, tmp_path):
    """Recompute one field from the cache: the star chi-square, the counts and
    the frozen light-curve hash must equal the committed values."""
    t = {"first": targets[0], "last": targets[-1],
         "ngc3786": next(x for x in targets if x["counterpart"] == "NGC 3786")}[which]
    name = t["name"]
    expected = load_json("qa_fields.json")[name]
    trow = load_json("qa_targets.json")[name]
    floors = load_json("qa_error_floor.json")["floors_mmag"][name]
    field = sample_io.load_field(t, stars[name])
    sv = qa.star_tables(field["stars"])
    loo = differential.leave_one_out(sv, err_key=qa.ERR_KEY)
    offsets = differential.field_offsets(sv)
    rows = qa.star_qa(field["star_names"], field["star_w1"], sv, loo)
    got = qa.field_summary(rows, offsets)
    assert got["n_stars"] == expected["n_stars"]
    for b in BANDS:
        for key in ("chi2_red_median_before", "chi2_red_median_after", "offset_std_mmag"):
            assert got[b][key] == pytest.approx(expected[b][key], rel=1e-9)
    v = differential.apply_offsets(visit_arrays(field["target"], "neowise"), offsets,
                                   err_key=qa.ERR_KEY)
    counts = qa.target_counts(field["target"], v)
    for key in ("n_frames", "n_visits_corrected", "n_visits_allwise_both", "first_visit_date"):
        assert counts[key] == trow[key]
    pairs = qa.untimely_pairs(v, untimely.read_field(
        sample_io.field_paths(name)["untimely"]).get(0))
    assert pairs == trow["untimely"]
    final = errorfloor.apply_floor(v, {b: 1e-3 * floors[b] for b in BANDS})
    path = tmp_path / f"{name}.csv"
    qa.write_frozen(path, final)
    freeze = json.loads(FREEZE_PATH.read_text())
    assert irsa.sha256_file(path) == freeze["light_curves"][name]["sha256"]


def test_qa_summary_is_consistent_with_the_per_object_files(targets):
    summary = load_json("qa_summary.json")
    qa_targets = load_json("qa_targets.json")
    fields = load_json("qa_fields.json")
    with open(config.RESULTS_DIR / "qa_stars.csv", newline="") as fh:
        star_rows = list(csv.DictReader(fh))
    assert summary["n_sample"] == len(targets)
    assert summary["n_with_light_curve"] == sum(r["status"] == "ok" for r in qa_targets.values())
    assert summary["n_comparison_stars_with_data"] == len(star_rows)
    for b in BANDS:
        after = np.array([float(r[f"{b}_chi2_red_after"]) for r in star_rows])
        assert summary["stars"][b]["chi2_red_after"]["p50"] == pytest.approx(
            np.median(after), abs=1e-3)
        med = [f[b]["chi2_red_median_after"] for f in fields.values()]
        assert summary["fields"][b]["chi2_red_median_after"]["p50"] == pytest.approx(
            np.median(med), rel=1e-9)
        floors = load_json("qa_error_floor.json")
        own = [v[b] for n, v in floors["floors_mmag"].items() if floors["source"][n][b] == "own"]
        assert summary["error_floor"][b]["n_own"] == len(own)
        assert all(x >= 0 for x in own)


def test_error_floor_follows_the_documented_rule():
    """floor^2 = max(0, s^2 - s_ref^2), own value only with >= MIN_PAIRS pairs."""
    floors = load_json("qa_error_floor.json")
    qa_targets = load_json("qa_targets.json")
    with open(sample_io.SAMPLE_PATH, newline="") as fh:
        w1 = {r["name"]: float(r["w1mpro"]) for r in csv.DictReader(fh)}
    for name, row in floors["floors_mmag"].items():
        u = qa_targets[name]["untimely"]
        for b in BANDS:
            ref = errorfloor.reference_for(w1[name], floors["reference_scatter_of_stars"][b])
            if floors["source"][name][b] == "own":
                assert u[b]["n_pairs"] >= errorfloor.MIN_PAIRS
                expect = 1e3 * errorfloor.floor_from_scatter(1e-3 * u[b]["pair_scatter_mmag"], ref)
                assert row[b] == pytest.approx(expect, abs=0.011)
            else:
                assert u["status"] != "ok" or u[b]["n_pairs"] < errorfloor.MIN_PAIRS
                assert row[b] == pytest.approx(floors["fallback_floor_mmag"][b], abs=0.011)


def test_star_injection_table_is_complete_and_monotonic():
    doc = load_json("qa_injection_stars.json")
    with open(config.RESULTS_DIR / "qa_injection_completeness.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    grid = config.INJECTION
    per_group = (len(grid.amplitudes_mag) * len(grid.rise_years) * len(grid.baseline_visits))
    assert len(rows) == 3 * per_group  # all stars and the two W1 bins
    assert doc["n_trials_total"] == (sum(int(r["n_trials"]) for r in rows
                                         if r["star_w1_bin"] == "all")
                                     + doc["null"]["all"]["n_null_segments"])
    # larger injected steps are never recovered less often (instant rise, all stars)
    for base in grid.baseline_visits:
        comp = [float(r["completeness_a_or_b"]) for r in rows
                if r["star_w1_bin"] == "all" and int(r["baseline_visits"]) == base
                and float(r["rise_yr"]) == 0.0]
        assert len(comp) == len(grid.amplitudes_mag)
        assert all(b >= a - 0.01 for a, b in zip(comp, comp[1:]))
    null = doc["null"]["all"]
    assert null["flagged_a_or_b"] <= null["flagged_a"] + null["flagged_b"]
    assert null["n_stars"] == doc["n_stars"]
