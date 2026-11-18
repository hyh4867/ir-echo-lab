"""Stage 7b regression tests (post hoc supplement to the X-ray analysis).

They read committed results only and check that the corrected comparisons are
what the report says they are.
"""

import json
import math

import pytest

from irecho import config, xray

RES = config.RESULTS_DIR
pytestmark = pytest.mark.skipif(not (RES / "s7b_summary.json").exists(), reason="stage 7b results not present")


def load(name):
    return json.loads((RES / name).read_text())


def test_bat_common_band_uses_channels_1_to_4_selected_by_number():
    b = load("s7b_light.json")["bat_common_band"]
    lo, hi = zip(*b["ebounds_kev"])
    assert b["channels_used"] == xray.bat_channels_in_band(lo, hi, 14.0, 50.0) == [1, 2, 3, 4]
    assert b["channel_bounds_used_kev"] == [[14.0, 20.0], [20.0, 24.0], [24.0, 35.0], [35.0, 50.0]]
    free, fixed, full = (b[k] for k in ("bat_ch1_4_free_gamma", "bat_ch1_4_gamma_fixed_to_nustar", "bat_all_channels_integrated_14_50"))
    assert free["channels_noticed"] == fixed["channels_noticed"] == [1, 2, 3, 4]
    assert (free["n_free"], free["dof"]) == (2, 2) and (fixed["n_free"], fixed["dof"]) == (1, 3)
    assert full["channels_noticed"] == list(range(1, 9)) and full["dof"] == 6
    assert fixed["gamma"] == pytest.approx(b["nustar_14_50"]["gamma"])
    for rec in (free, fixed, full):
        r = rec["nustar_over_bat"]
        assert r["log_ratio"] == pytest.approx(b["nustar_14_50"]["log_flux_14_50"] - rec["log_flux_14_50"])
        assert r["lo"] < r["log_ratio"] < r["hi"]


def test_stage7_value_was_the_full_band_fit_and_is_kept_under_its_corrected_name():
    b = load("s7b_light.json")["bat_common_band"]["bat_all_channels_integrated_14_50"]
    old = load("s7_fits.json")["simple_checks"]
    assert b["gamma"] == pytest.approx(old["bat_powerlaw"]["gamma"], abs=1e-3)       # 2.04: the 14-195 keV index
    assert b["chi2"] == pytest.approx(old["bat_powerlaw"]["chi2"], abs=1e-3) and old["bat_powerlaw"]["dof"] == 6
    assert b["nustar_over_bat"]["log_ratio"] == pytest.approx(old["delta_14_50_powerlaw"]["value"], abs=1e-3)
    assert "not a fit of the 14-50 keV channels only" in load("s7_summary.json")["independent_checks"]["flux_14_50_powerlaw"][
        "description_corrected_in_7b"]
    s = load("s7b_summary.json")
    assert "all 8 channels" in s["flux_ratio_14_50"]["stage7_value_renamed"]["name"]
    for m in ("M0", "M2"):
        assert "not a model-independent" in s["c1_shape_tied"][m]["assumption"]


def test_xrt_constant_readings_differ_only_in_delta2_and_by_the_constant():
    f = load("s7b_fits.json")
    a, b = f["c2_m0_const_both_xrt"], f["c2_m0_const_2022_only"]
    assert a["statistic"] == pytest.approx(b["statistic"], abs=0.01)       # same fit, different reading
    assert a["const_xrt22"] == pytest.approx(b["const_xrt22"], abs=1e-3)
    assert b["delta2"]["value"] - a["delta2"]["value"] == pytest.approx(-math.log10(a["const_xrt22"]), abs=0.005)
    assert a["delta_nh"]["value"] == pytest.approx(b["delta_nh"]["value"], abs=1e-3)
    s = load("s7b_summary.json")["m0_delta2_delta_nh"]
    for key, rec in (("a_constant_on_both_xrt", a), ("b_constant_on_2022_xrt_only", b)):
        assert s[key]["delta2_90"] == pytest.approx([rec["delta2"]["lo"], rec["delta2"]["hi"]])
        assert s[key]["label_delta2_stat_only"] == xray.luminosity_label(s[key]["delta2"], *s[key]["delta2_90"])
        assert s[key]["label_delta_nh"] == xray.absorption_label(s[key]["delta_nh"], *s[key]["delta_nh_90"])
    env = load("s7b_summary.json")["m0_delta2_envelope"]
    assert env["lo"] == pytest.approx(a["delta2"]["lo"]) and env["hi"] == pytest.approx(max(b["delta2"]["hi"], s["stage7_xrt_constant_1"]["delta2_90"][1]))


def test_xrt_offset_is_also_there_with_the_reflection_model():
    m = load("s7b_fits.json")["m2_xrt_free"]
    assert m["const_xrt22_90"][1] < 1.0
    assert m["reference_const_fixed"]["statistic"] - m["statistic"] > 6.63
    assert any("lg10Lum = 36" in p for p in m["pegged"])        # the scattered component sits on its lower limit


def test_reload_in_a_new_session_and_energy_grid_reproduce_the_fits():
    f = load("s7b_fits.json")
    r = f["reload_m0"]
    assert r["returncode"] == 0 and r["statistic"] == pytest.approx(r["stage7_statistic"], abs=0.01)
    assert r["free_parameters"]["zTBabs.nH#3"] == pytest.approx(r["stage7_values"]["nh_1e22"], rel=1e-4)
    assert r["free_parameters"]["zpowerlw.PhoIndex#9"] == pytest.approx(r["stage7_values"]["gamma"], rel=1e-4)
    assert r["free_parameters"]["clumin.lg10Lum#8"] == pytest.approx(r["stage7_values"]["log_lint_2_10"], abs=1e-3)
    g = f["grid_3000_m2"]["difference"]
    assert abs(g["statistic"]) < 0.2 and abs(g["gamma"]) < 0.005 and abs(g["nh_fraction"]) < 0.005
    assert abs(g["refl_norm_fraction"]) < 0.01
    for key, rec in f["xcm"].items():
        if not isinstance(rec, dict) or "xcm" not in rec:
            continue
        assert rec["refit_statistic"] == pytest.approx(rec["stage7_statistic"], abs=0.05), key
        text = (config.PROJECT_ROOT / rec["xcm"]).read_text()
        assert "energies 0.1 300. 1500 log" in text and "/Users/" not in text and "\ncd " not in text
        assert (config.PROJECT_ROOT / rec["log"]).exists()
    d = f["saa_default_m0"]["difference"]
    assert abs(d["gamma"]) < 0.05 and abs(d["nh_fraction"]) < 0.1 and abs(d["flux_3_79_fraction"]) < 0.05


def test_xrt_corrections_and_the_2sxps_comparison():
    x = load("s7b_light.json")["xrt_corrections"]
    for key in ("xrt10", "xrt22"):
        r = x[key]
        assert r["reproduces_the_arf_used_in_stage_7"] is True
        for a, b, c in zip(r["factor_psf_and_expomap"], r["factor_psf_only"], r["factor_expomap_only"]):
            assert a == pytest.approx(b * c, rel=1e-6) and 0.5 < a < 1.0
        assert r["nearest_bad_column_pix_min"] <= 20          # a bad column crosses the region in at least one interval
    c = x["comparison_with_2sxps"]
    if "withheld" not in c["catalogue"]:                      # the catalogue row is not reproduced in the public copy
        assert c["catalogue"]["n_obs"] == 1 and c["catalogue"]["first_obs"].startswith("2010-12-29")
    for band in ("0.3-10", "2-10"):
        assert abs(c[band]["difference_sigma"]) < 1.0
        assert c[band]["own_rate_corrected"] == pytest.approx(c[band]["own_net_rate_raw"] / c[band]["correction_factor"])


def test_gti_overlap_chain_diagnostics_and_the_validation_record():
    d = load("s7b_light.json")
    g = d["gti_overlap"]
    assert g["overlap_with_fpma_s"] < 0.02 * g["xrt22_gti_total_s"]          # the two exposures are interleaved
    m = d["mcmc_diagnostics"]
    for key in ("delta2", "delta_nh"):
        assert m[key]["gelman_rubin_split"] < 1.05 and m[key]["effective_sample_size"] > 500
    assert m["delta2_autocorr_time_if_rows_were_walker_blocks"] < 2.0      # so rows do cycle through the walkers
    s = d["synthetic_validation_record"]
    assert s["binomial_se_at_90"] == pytest.approx(xray.binomial_se(0.9, 200))
    assert s["coverage_minus_nominal_in_sigma"]["gamma"] == pytest.approx(-2.12, abs=0.02)
    assert s["n_nh10_at_limit"] == s["n_gamma_at_limit"] == s["n_non_finite"] == 0
    assert 140 < s["xrt22_net_counts_median"] < 170


def test_hardness_ratio_is_a_poisson_lower_limit_and_the_truncated_value_is_gone():
    h = load("s7b_light.json")["hardness_ratio"]
    assert h["xrt22"]["soft_counts_source_region"] < h["xrt22"]["soft_background_expected"]
    assert 0.8 < h["xrt22"]["hr_lower_limit_90"] < 1.0
    assert h["xrt10"]["hr_p05_p50_p95"][0] < h["xrt10"]["hr_p05_p50_p95"][1] < h["xrt22"]["hr_lower_limit_90"]
    for key in ("xrt10", "xrt22"):
        assert "hardness_ratio" not in load("s7_summary.json")["xrt_count_rates"][key]
    frac = h["background_fraction_in_source_region"]
    assert frac["fpma"] < 0.1 and frac["fpmb"] < 0.1


def test_registered_verdict_of_stage_7_is_unchanged():
    assert load("s7b_summary.json")["registered_verdict_stage7"] == load("s7_summary.json")["registered_verdict"]
    assert load("s7_summary.json")["registered_verdict"]["verdict"] == "H1 지지, 가림 변화는 배제 못 함"
    plan = (config.PROJECT_ROOT / "docs" / "XRAY_PLAN.md").read_text()
    assert "## 9. 7b 사후 보완" in plan and "7단계 결과를 모두 본 뒤에 쓴 절" in plan
