"""Regression and known-answer tests on the cached IRSA data.

These need data/cache/ (scripts/fetch_data.py) and are skipped when it is
absent, because raw IRSA tables are not committed. They never touch the
network. Expected numbers come from results/*.json, which scripts/run_analysis.py
wrote from the same cache, and from Son et al. (2022) for the known answer.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

from irecho import analysis, config, irsa
from irecho.lightcurve import load_target, write_epochs_csv

pytestmark = pytest.mark.cached


def _cache_complete() -> bool:
    try:
        targets = analysis.all_targets()
    except FileNotFoundError:
        return False
    return all(p.exists() for t in targets for p in irsa.cache_paths(t["name"]).values())


if not _cache_complete():
    pytestmark = [pytest.mark.cached,
                  pytest.mark.skip(reason="IRSA cache absent: run scripts/fetch_data.py")]


def load_json(name: str) -> dict:
    return json.loads((config.RESULTS_DIR / name).read_text())


@pytest.fixture(scope="module")
def ngc3786_lc():
    t = analysis.load_ngc3786()
    return load_target(t["name"], t["ra"], t["dec"])


@pytest.fixture(scope="module")
def ngc3786(ngc3786_lc):
    return analysis.ngc3786_analysis(ngc3786_lc)


@pytest.fixture(scope="module")
def star_lcs():
    return [load_target(t["name"], t["ra"], t["dec"]) for t in analysis.load_reference_stars()]


@pytest.fixture(scope="module")
def control_lcs():
    return [load_target(t["name"], t["ra"], t["dec"]) for t in analysis.load_controls()]


# --- data integrity -----------------------------------------------------------
def test_cached_files_match_manifest_hashes():
    manifest = json.loads(config.MANIFEST_PATH.read_text())
    checked = 0
    for t in analysis.all_targets():
        for path in irsa.cache_paths(t["name"]).values():
            key = str(path.relative_to(config.CACHE_DIR))
            assert irsa.sha256_file(path) == manifest[key]["sha256"], key
            checked += 1
    assert checked == 2 * len(analysis.all_targets())


def test_target_counts():
    assert len(analysis.load_reference_stars()) >= 20
    assert len(analysis.load_controls()) >= 5


# --- (b) NGC 3786: known answer from Son et al. (2022) -------------------------
def test_ngc3786_flare_onset_is_bracketed_by_2020_visits(ngc3786):
    """Paper: the flare 'appeared to have occurred around mid-2020'."""
    p = ngc3786["paper_range"]
    assert p["n_flagged"] == 1
    assert p["onset_date"].startswith("2020-12")
    table = ngc3786["visit_table"]
    onset = next(i for i, row in enumerate(table) if row["flagged"])
    assert table[onset - 1]["date"].startswith("2020-05")
    assert not any(row["flagged"] for row in table[:onset])


def test_ngc3786_amplitude_matches_paper(ngc3786):
    """Paper: amplitude '~0.5 and 0.8 mag in W1 and W2, respectively'."""
    v = ngc3786["paper_range"]["onset_visit"]
    assert v["w1_dmag"] == pytest.approx(0.5, abs=0.15)
    assert v["w2_dmag"] == pytest.approx(0.8, abs=0.15)
    assert v["w2_dmag"] > v["w1_dmag"]


def test_ngc3786_becomes_redder_when_brighter(ngc3786):
    """Paper: 'W1-W2 became redder as it brightened in the W1-band'."""
    v = ngc3786["paper_range"]["onset_visit"]
    assert v["delta_color_mag"] > 5 * v["color_err_mag"]
    assert v["delta_color_mag"] > 3 * ngc3786["paper_range"]["baseline"]["color_std_mag"]


def test_ngc3786_matches_committed_summary(ngc3786):
    saved = load_json("ngc3786_summary.json")
    for rng in ("paper_range", "full_range"):
        assert ngc3786[rng]["flagged_dates"] == saved[rng]["flagged_dates"]
        for visit in ("onset_visit", "peak_w2_visit"):
            for key in ("w1_mag", "w2_mag", "w1_dmag", "w2_dmag", "color_mag",
                        "delta_color_mag", "w1_significance", "w2_significance"):
                assert ngc3786[rng][visit][key] == pytest.approx(saved[rng][visit][key], abs=1e-5)
        t_new = ngc3786[rng]["onset_visit"]["excess_blackbody"]["temperature_k"]
        t_old = saved[rng]["onset_visit"]["excess_blackbody"]["temperature_k"]
        assert t_new == pytest.approx(t_old, abs=1e-3)
    assert ngc3786["n_visits"] == saved["n_visits"]
    assert ngc3786["cut_report"] == saved["cut_report"]


def test_ngc3786_epoch_table_matches_committed_csv(ngc3786_lc, tmp_path):
    out = tmp_path / "ngc3786.csv"
    write_epochs_csv(ngc3786_lc, out)
    assert out.read_text() == (config.RESULTS_DIR / "epochs" / "ngc3786.csv").read_text()


def test_ngc3786_stays_below_saturation(ngc3786_lc):
    e = ngc3786_lc.epochs
    assert np.nanmin(e["w1_mag"]) > 8.0  # NEOWISE Expl. Supp. II.2: W1 ~ 8.0
    assert np.nanmin(e["w2_mag"]) > 7.5  # W2 ~ 7.5


# --- (a) reference stars --------------------------------------------------------
def test_reference_stars_match_committed_summary(star_lcs):
    result = analysis.reference_star_analysis(star_lcs)
    saved = load_json("reference_stars_summary.json")["summary"]
    assert result["summary"]["n_stars"] == saved["n_stars"] == 30
    for b in ("w1", "w2"):
        for model in analysis.ERROR_MODELS:
            for key in ("chi2_red_median", "pooled_floor_mag_neowise"):
                assert result["summary"][b][model][key] == pytest.approx(
                    saved[b][model][key], abs=1e-5)


def test_no_reference_star_is_flagged_as_a_flare(star_lcs):
    result = analysis.reference_star_analysis(star_lcs)
    saved = load_json("reference_stars_summary.json")["summary"]
    assert result["summary"]["n_stars_flagged_by_detector"] == saved["n_stars_flagged_by_detector"]
    assert result["summary"]["n_stars_flagged_by_detector"] == 0


def test_flare_epoch_is_not_a_field_wide_zero_point_shift(star_lcs, ngc3786):
    """The stars share NGC 3786's visits. At the flare visit their common
    offset must be tiny compared with the flare amplitude."""
    common = analysis.common_mode_analysis(star_lcs)
    flare_mjd = ngc3786["paper_range"]["onset_mjd"]
    k = int(np.argmin(np.abs(np.array(common["visit_mjd"]) - flare_mjd)))
    assert abs(common["visit_mjd"][k] - flare_mjd) < 15.0
    for b, amp in (("w1", ngc3786["paper_range"]["onset_visit"]["w1_dmag"]),
                   ("w2", ngc3786["paper_range"]["onset_visit"]["w2_dmag"])):
        assert abs(common[b]["offset_mag"][k]) < 0.05 * amp


# --- (c) controls ---------------------------------------------------------------
def test_controls_match_committed_summary(control_lcs):
    result = analysis.control_analysis(control_lcs)
    saved = load_json("controls_summary.json")
    assert result["n_objects"] == saved["n_objects"] == 8
    assert result["n_objects_flagged"] == saved["n_objects_flagged"]
    for new, old in zip(result["per_object"], saved["per_object"]):
        assert new["name"] == old["name"]
        assert new["flagged_dates"] == old["flagged_dates"]
        for b in ("w1", "w2"):
            assert new[f"{b}_max_significance"] == pytest.approx(
                old[f"{b}_max_significance"], abs=1e-5)
    assert result["injection_recovery_pooled"] == saved["injection_recovery_pooled"]


# --- independent cross-check (added after the first full run) -----------------
def test_independent_recipe_agrees_with_pipeline(ngc3786):
    """scripts/independent_check.py shares no code with src/irecho (no cuts,
    mean magnitudes, standard library only). Its amplitudes differ from the
    pipeline's mainly through the baseline definition (mean vs. median), so
    the tolerance is the baseline scatter, not the measurement error."""
    ind = load_json("independent_check.json")
    v = ngc3786["paper_range"]["onset_visit"]
    assert ind["n_visits"] == ngc3786["n_visits"]["neowise"]
    assert ind["w1_dmag"] == pytest.approx(v["w1_dmag"], abs=0.06)
    assert ind["w2_dmag"] == pytest.approx(v["w2_dmag"], abs=0.06)
    assert ind["flare_w1_mag"] == pytest.approx(v["w1_mag"], abs=0.02)
    assert ind["flare_w2_mag"] == pytest.approx(v["w2_mag"], abs=0.02)
    assert ind["delta_color_mag"] == pytest.approx(v["delta_color_mag"], abs=0.03)
