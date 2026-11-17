"""Regression tests for stages 0, 3 and 4 on the cached data.

Skipped unless the full cache is present (stage-2 targets, the 80 comparison
stars of the control fields and the unTimely extracts). Expected numbers are
results/stage34_*.json, written by scripts/run_stage34.py from the same cache.
Tolerances on physical agreement (unTimely, temperature range) were chosen
after the first full run and are documented as such in the report.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

from irecho import analysis, config, injection, irsa, stage34, untimely, ztf

pytestmark = pytest.mark.cached


def _cache_complete() -> bool:
    try:
        stars = analysis.load_field_stars()
        galaxies = [analysis.load_ngc3786()] + analysis.load_controls()
    except FileNotFoundError:
        return False
    return (len(stars) > 0
            and all(irsa.cache_paths(s["name"])["neowise"].exists() for s in stars)
            and all(untimely.cache_path(g["name"]).exists() for g in galaxies)
            and all(ztf.cache_path(n, b).exists() for n in ("ngc3786", "ngc2273") for b in "gr")
            and all(p.exists() for t in analysis.all_targets()
                    for p in irsa.cache_paths(t["name"]).values()))


if not _cache_complete():
    pytestmark = [pytest.mark.cached, pytest.mark.skip(
        reason="stage 3-4 cache absent: run scripts/fetch_data.py and scripts/fetch_untimely.py")]


def load_json(name: str) -> dict:
    return json.loads((config.RESULTS_DIR / name).read_text())


@pytest.fixture(scope="module")
def fields():
    return stage34.load_fields()


@pytest.fixture(scope="module")
def ngc_corrected(fields):
    return stage34.corrected_target(fields["ngc3786"])


def test_new_cache_files_match_manifest_hashes():
    manifest = json.loads(config.MANIFEST_PATH.read_text())
    paths = [irsa.cache_paths(s["name"])["neowise"] for s in analysis.load_field_stars()]
    paths += [untimely.cache_path(g["name"])
              for g in [analysis.load_ngc3786()] + analysis.load_controls()]
    paths += [ztf.cache_path(n, b) for n in ("ngc3786", "ngc2273") for b in "gr"]
    for path in paths:
        key = str(path.relative_to(config.CACHE_DIR))
        assert irsa.sha256_file(path) == manifest[key]["sha256"], key
    assert len(paths) == 80 + 9 + 4


def test_every_control_field_has_ten_comparison_stars(fields):
    assert len(fields["ngc3786"]["stars"]) == 30
    for name, field in fields.items():
        if name != "ngc3786":
            assert len(field["stars"]) == 10, name


# --- stage 0a ---------------------------------------------------------------------
def test_differential_correction_removes_the_w2_excess(fields):
    result = stage34.differential_validation(fields["ngc3786"]["stars"])
    saved = load_json("stage34_differential.json")
    for model in analysis.ERROR_MODELS:
        for b in ("w1", "w2"):
            for key in ("chi2_red_median_before", "chi2_red_median_after"):
                assert result["models"][model][b][key] == pytest.approx(
                    saved["models"][model][b][key], rel=1e-4)
        w2 = result["models"][model]["w2"]
        assert w2["chi2_red_median_after"] < w2["chi2_red_median_before"]
        assert w2["rms_median_mag_after"] < w2["rms_median_mag_before"]
    # The adopted model ends within 30% of 1 in both bands (it is not exactly 1:
    # the errors are overestimated, see README section 5.1).
    adopted = result["models"][stage34.ERR_KEY]
    for b in ("w1", "w2"):
        assert 0.7 < adopted[b]["chi2_red_median_after"] < 1.3


# --- stage 0b ---------------------------------------------------------------------
def test_untimely_excess_flux_agrees_with_the_pipeline(fields):
    from irecho.analysis import run_detector, visit_arrays
    res = run_detector(visit_arrays(fields["ngc3786"]["target"], "neowise"))
    flare = stage34.untimely_flare(fields["ngc3786"], res)
    saved = load_json("stage34_untimely.json")["ngc3786_flare"]
    assert flare["onset_matched"] and flare["n_baseline_visits_matched"] == 13
    for b, tol in (("w1", 0.05), ("w2", 0.10)):
        assert flare[b]["untimely_excess_flux_mjy"] == pytest.approx(
            saved[b]["untimely_excess_flux_mjy"], rel=1e-4)
        assert flare[b]["untimely_excess_flux_mjy"] == pytest.approx(
            flare[b]["pipeline_excess_flux_mjy"], rel=tol)
        # unTimely assigns less total flux, so its magnitude amplitude is larger.
        assert flare[b]["untimely_baseline_mag"] > flare[b]["pipeline_baseline_mag"]
        assert flare[b]["untimely_dmag"] > flare[b]["pipeline_dmag"]


# --- stage 3 ----------------------------------------------------------------------
def test_physics_matches_committed_results(ngc_corrected):
    z = json.loads((config.TARGETS_DIR / "ngc3786.json").read_text())["ned_z"]
    phys = stage34.physics_analysis(ngc_corrected, z)
    saved = load_json("stage34_physics_ngc3786.json")
    assert [r["date"] for r in phys["flagged_visits"]] == [
        r["date"] for r in saved["flagged_visits"]]
    assert len(phys["flagged_visits"]) == 8
    for new, old in zip(phys["flagged_visits"], saved["flagged_visits"]):
        for key in ("t_rest_k", "l_dust_erg_s", "r_bb_light_days", "t_rest_k_p16",
                    "r_bb_light_days_p84"):
            assert new[key] == pytest.approx(old[key], rel=1e-4)
        assert 800.0 < new["t_rest_k"] < 1400.0  # hot dust, below sublimation
        assert new["t_rest_k_p16"] < new["t_rest_k"] < new["t_rest_k_p84"]
    assert phys["luminosity_distance_mpc"] == pytest.approx(40.68, abs=0.01)


# --- stage 4 ----------------------------------------------------------------------
def test_detectors_on_real_light_curves(fields):
    det = stage34.detector_analysis(fields)
    saved = load_json("stage34_detectors.json")
    assert det["galaxies"]["ngc3786"]["a_n_flagged"] == 8
    assert det["galaxies"]["ngc3786"]["b"]["flagged"]
    assert det["stars"]["n"] == 110
    assert det["stars"]["n_flagged_a"] == saved["stars"]["n_flagged_a"] == 0
    assert det["stars"]["n_flagged_b"] == saved["stars"]["n_flagged_b"] == 0
    for name, g in det["galaxies"].items():
        assert g["a_flagged_dates"] == saved["galaxies"][name]["a_flagged_dates"], name
        assert g["b"]["flagged"] == saved["galaxies"][name]["b"]["flagged"], name
        if g["b"].get("w1_p") is not None:
            assert g["b"]["w1_p"] == pytest.approx(saved["galaxies"][name]["b"]["w1_p"], rel=1e-3)


def test_injection_counts_for_the_control_hosts(fields):
    """Detectors (a) and (b) are deterministic, so the control-galaxy rows of
    the committed table can be recomputed quickly (the ML column needs all
    hosts and is not recomputed here)."""
    hosts = [h for h in stage34.injection_hosts(fields) if h["kind"] == "control"]
    rows = injection.run_grid(hosts)
    saved = {(r["baseline_visits"], r["rise_years"], r["amp_mag"]): r
             for r in load_json("stage34_injection.json")["recovery"] if r["kind"] == "control"}
    grid = config.INJECTION
    for rise in grid.rise_years:
        for amp in grid.amplitudes_mag:
            sel = [r for r in rows if r["amp"] == amp and r["rise"] == rise]
            old = saved[("all", rise, amp)]
            assert len(sel) == old["n"]
            assert sum(r["a"] for r in sel) == old["a"]
            assert sum(r["b"] for r in sel) == old["b"]
    null = [r for r in rows if r["amp"] == 0.0]
    old_null = next(r for r in load_json("stage34_injection.json")["false_alarm"]
                    if r["kind"] == "control" and r["baseline_visits"] == "all")
    assert sum(r["a"] or r["b"] for r in null) == old_null["a_or_b"]
    assert np.isfinite([r["features"] for r in rows]).all()


def test_independent_physics_recipe_agrees(ngc_corrected):
    """scripts/independent_check.py recomputes the first flare visit with the
    standard library only (own constants, bisection, cz/H0 distance, raw
    photometry). Tolerances cover the known differences in recipe."""
    ind = load_json("independent_check.json")["physics"]
    first = load_json("stage34_physics_ngc3786.json")["flagged_visits"][0]
    assert ind["temperature_k"] == pytest.approx(first["t_obs_k"], rel=0.03)
    assert ind["l_dust_erg_s"] == pytest.approx(first["l_dust_erg_s"], rel=0.08)
    assert ind["r_bb_light_days"] == pytest.approx(first["r_bb_light_days"], rel=0.06)


def test_mean_brightening_matches_son2026(fields):
    """Second known answer: Son et al. (2026, Sec. 2) report that NGC 3786
    'has brightened by an average of 0.68 and 1.02 mag in the W1 and W2 bands
    ... compared to the pre-outburst state (MJD <~ 59000)'. Tolerance 0.03 mag
    (the paper does not say how the reference level was averaged)."""
    out = stage34.followup_comparisons(fields["ngc3786"])
    assert out["uncorrected"]["n_pre"] == 13 and out["uncorrected"]["n_post"] == 8
    assert out["uncorrected"]["w1_mean_minus_mean"] == pytest.approx(0.68, abs=0.03)
    assert out["uncorrected"]["w2_mean_minus_mean"] == pytest.approx(1.02, abs=0.03)
    # The 2010 colour is redder than the 2014-2020 baseline in both photometries.
    c = out["color_2010_vs_2014_2020"]
    assert min(c["untimely_2010_colors"]) > c["untimely_2014_2020_pre_outburst_median_color"] + 0.1
    assert min(c["pipeline_allwise_2010_colors"]) > c["pipeline_neowise_pre_outburst_median_color"] + 0.1


# --- docs/PREDICTIONS.md ------------------------------------------------------------
def test_hypothesis_verdicts(ngc_corrected):
    """The verdicts follow from the rules committed in docs/PREDICTIONS.md
    (0aefbe4). They are asserted here so that a change in the photometry that
    flips one of them cannot go unnoticed."""
    saved = load_json("stage34_hypotheses.json")
    phys = load_json("stage34_physics_ngc3786.json")
    hyp = stage34.hypothesis_tests(phys)
    opt = stage34.optical_comparison(ngc_corrected, phys["baseline"]["w1_flux_mjy"])
    for key in ("H1_temperature_mixing", "H2_silicate", "H3_dust_destruction"):
        assert hyp[key]["verdict"] == saved[key]["verdict"]
    assert hyp["H1_temperature_mixing"]["verdict"] == "undetermined"
    assert hyp["H2_silicate"]["verdict"] == "undetermined"
    assert hyp["H3_dust_destruction"]["verdict"] == "supported"
    assert opt["verdict"] == saved["alternative_central_source"]["verdict"] == "supported"
    for band in ("g", "r"):
        new, old = opt["bands"][band], saved["alternative_central_source"]["bands"][band]
        assert new["brightening_mag"] == pytest.approx(old["brightening_mag"], abs=1e-5)
        assert new["spearman_rho"] == pytest.approx(old["spearman_rho"], abs=1e-5)
        assert new["n_pairs_for_correlation"] == old["n_pairs_for_correlation"]
