"""Stage 7 regression tests on the committed X-ray results.

They need no HEASoft and no raw data: they check that the committed verdict
table follows from the committed fit results and the rules registered in
docs/XRAY_PLAN.md, and that the documents quote the committed numbers.
"""

import importlib.util
import json
import math

import pytest

from irecho import config, physics, xray

RES = config.RESULTS_DIR
pytestmark = pytest.mark.skipif(not (RES / "s7_summary.json").exists(), reason="stage 7 results not present")


def load(name):
    return json.loads((RES / name).read_text())


def _script(name):
    spec = importlib.util.spec_from_file_location(name, config.PROJECT_ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_observations_are_the_registered_ones_and_were_reprocessed():
    reg = load("s7_regions.json")
    assert reg["fpma"]["header"]["OBS_ID"] == reg["fpmb"]["header"]["OBS_ID"] == xray.OBS["nustar22"]
    assert reg["xrt10"]["header"]["OBS_ID"] == xray.OBS["xrt10"] and reg["xrt22"]["header"]["OBS_ID"] == xray.OBS["xrt22"]
    assert reg["xrt10"]["header"]["DATE-OBS"].startswith("2010-12-29")
    assert reg["xrt22"]["header"]["DATE-OBS"].startswith("2022-05-01")
    assert reg["fpma"]["header"]["DATE-OBS"].startswith("2022-05-01")
    for key in ("xrt10", "xrt22"):
        assert reg[key]["header"]["DATAMODE"] == "PHOTON"
    for key in reg:
        assert "V6.37.1" in reg[key]["header"]["SOFTVER"]
        assert reg[key]["centroid_offset_arcsec"] < 5.0          # the source is where the catalogue puts it
        assert not reg[key]["event_file"].startswith("/")
    assert reg["xrt10"]["exposure_s"] == pytest.approx(10309, abs=2)
    assert reg["xrt22"]["exposure_s"] == pytest.approx(6935, abs=2)
    assert reg["fpma"]["exposure_s"] == pytest.approx(22641, abs=2)
    for key in ("fpma", "fpmb"):
        b = reg[key]["background"]
        assert b["radius_arcsec"] == 100.0 and b["same_chip_fraction"] >= 0.95
        assert b["separation_arcsec"] >= 50 + 100 + 40
    assert reg["xrt10"]["contaminants_masked"] == [] and reg["xrt22"]["contaminants_masked"] == []


def test_data_document_lists_the_files_with_their_hashes():
    manifest = json.loads(config.MANIFEST_PATH.read_text())
    s7 = {k: v for k, v in manifest.items() if k.startswith("s7/")}
    assert len(s7) == 123
    doc = (config.PROJECT_ROOT / "docs" / "DATA.md").read_text()
    for obsid in xray.OBS.values():
        assert obsid in doc
    for key in ("s7/nustar/60160365002/event_uf/nu60160365002A_uf.evt.gz",
                "s7/swift/00040946001/xrt/event/sw00040946001xpcw3po_uf.evt.gz",
                "s7/swift/00081051001/xrt/event/sw00081051001xpcw3po_uf.evt.gz", "s7/bat157/SWIFT_J0929.7+6232.pha"):
        assert s7[key]["sha256"] in doc
    caldb = load("s7_caldb.json")
    assert caldb["nustar_subset"]["n_files"] == len(caldb["nustar_subset"]["files"]) == 94
    assert all(len(f["sha256"]) == 64 for f in caldb["nustar_subset"]["files"])
    assert {t["file"] for t in caldb["tarballs"]} >= {"goodfiles_swift_xrt.tar.Z", "goodfiles_swift_mis.tar.Z"}


def test_model_choice_follows_the_registered_aic_rule():
    fits, summary = load("s7_fits.json"), load("s7_summary.json")
    for model, rec in fits["fits_2022"].items():
        assert rec["aic"] == pytest.approx(rec["statistic"] + 2 * rec["n_free"])
        assert rec["spectra"] == ["xrt22", "fpma", "fpmb"]
    assert [fits["fits_2022"][m]["n_free"] for m in ("M0", "M1", "M2")] == [4, 5, 7]
    assert fits["chosen_model"] == summary["chosen_model"] == xray.choose_model(fits["aic_2022"])
    assert summary["aic_lowest_model"] == min(fits["aic_2022"], key=fits["aic_2022"].get)


def test_labels_and_verdict_follow_from_the_committed_fit_results():
    fits, summary = load("s7_fits.json"), load("s7_summary.json")
    for model, c in summary["by_model"].items():
        d1 = fits["c1"][model]["delta1"]
        assert c["delta1"]["value"] == pytest.approx(-math.log10(fits["c1"][model]["values"]["const_bat"]), abs=1e-6)
        lo, hi = xray.add_systematic(d1["value"], d1["lo"], d1["hi"], xray.SYS_NUSTAR_BAT_DEX)
        assert c["delta1"]["with_systematic_90"] == pytest.approx([lo, hi])
        assert c["delta1"]["label"] == xray.luminosity_label(d1["value"], lo, hi)
        v = fits["c2"][model]["values"]
        assert c["delta2"]["value"] == pytest.approx(v["log_lint_2_10"] - v["log_lint10_2_10"], abs=1e-6)
        assert c["delta_nh"]["value"] == pytest.approx(math.log10(v["nh_1e22"] / v["nh10_1e22"]), abs=1e-6)
        assert c["delta2"]["label"] == xray.luminosity_label(c["delta2"]["value"], *c["delta2"]["with_systematic_90"])
        assert c["delta_nh"]["label"] == xray.absorption_label(c["delta_nh"]["value"], *c["delta_nh"]["stat_90"],
                                                               measurable=c["delta_nh"]["measurable"])
        assert c["verdict"] == xray.overall_verdict(c["delta1"]["label"], c["delta2"]["label"], c["delta_nh"]["label"])
        if xray.LOG_LBOL_BASS is not None and c["energy_budget"]["f_need"] is not None:   # blank in the public copy
            assert c["energy_budget"]["f_need"] == pytest.approx(10 ** (xray.LOG_L_DUST_PEAK - xray.LOG_LBOL_BASS - c["delta1"]["value"]))
        assert c["energy_budget"]["label"] == "맞음"
    assert summary["registered_verdict"] == summary["by_model"][summary["chosen_model"]]["verdict"]


def test_profile_intervals_bracket_the_best_value_and_cross_the_90_percent_level():
    fits = load("s7_fits.json")
    for model in ("M0", "M1"):
        for key in ("delta2", "delta_nh"):
            p = fits["c2"][model][key]
            assert p["lo"] < p["value"] < p["hi"]
            assert max(p["delta_stat"]) > 2.706
            assert p["min_delta_stat"] > -0.5          # no much better minimum was found while profiling
        # an interval that excludes zero and a tied fit that is worse by more than 2.706 say the same thing
        tied = fits["c2"][model]["delta_stat_lum_tied"]
        assert (fits["c2"][model]["delta2"]["lo"] > 0) == (tied > 2.706)
        tied_nh = fits["c2"][model]["delta_stat_nh_tied"]
        d = fits["c2"][model]["delta_nh"]
        assert (d["lo"] > 0 or d["hi"] < 0) == (tied_nh > 2.706)


def test_reflection_model_leaves_the_2010_primary_unconstrained_and_is_reported_so():
    fits, summary = load("s7_fits.json"), load("s7_summary.json")
    c2 = fits["c2"]["M2"]
    assert c2["delta2"]["lo"] is None and c2["delta2"]["hi"] is None      # the profile found no bound
    assert c2["delta_stat_lum_tied"] < 2.706 and c2["delta_stat_nh_tied"] < 2.706
    m2 = summary["by_model"]["M2"]
    assert m2["delta2"]["label"] == xray.UNDECIDED and m2["delta_nh"]["label"] == xray.UNDECIDED
    assert m2["delta_nh"]["measurable"] is False
    assert m2["delta1"]["label"] == xray.HIGH          # the comparison with BAT does not need the 2010 spectrum


def test_luminosity_definition_checked_with_the_plain_normalisation():
    fits, chk = load("s7_fits.json"), load("s7_lum_check.json")
    m0 = fits["fits_2022"]["M0"]
    assert chk["statistic"] == pytest.approx(m0["statistic"], abs=0.05)
    lum = xray.zpowerlw_rest_luminosity(chk["norm"], chk["gamma"], xray.TARGET["z"],
                                        physics.luminosity_distance_mpc(xray.TARGET["z"]))
    assert math.log10(lum) == pytest.approx(chk["log_lint_2_10_from_norm"], abs=1e-9)
    # XSPEC's clumin (its own cosmology integral) against the project's distance and formula
    assert math.log10(lum) == pytest.approx(m0["values"]["log_lint_2_10"], abs=0.01)


def test_count_rates_recomputed_from_the_region_counts():
    reg, arf, summary = load("s7_regions.json"), load("s7_arf.json"), load("s7_summary.json")
    rates = _script("run_stage7").count_rates(reg, arf)
    for key in ("xrt10", "xrt22"):
        assert rates[key]["net_rate_2_10"] == pytest.approx(summary["xrt_count_rates"][key]["net_rate_2_10"])
        ratio = arf[key]["backscal_source"] / arf[key]["backscal_background"]
        assert ratio == pytest.approx(20.0 ** 2 / (110.0 ** 2 - 60.0 ** 2), rel=0.01)     # circle over annulus
        net = reg[key]["source_counts_2_10"] - reg[key]["background_counts_2_10"] * ratio
        assert rates[key]["net_counts_2_10"] == pytest.approx(net)
        assert rates[key]["net_counts_0p3_10"] >= xray.MIN_XRT_COUNTS
    r = summary["xrt_count_rates"]["ratio_2_10_2022_over_2010"]
    assert r["value"] == pytest.approx(rates["xrt22"]["net_rate_2_10"] / rates["xrt10"]["net_rate_2_10"])
    assert 0.95 < r["effective_area_ratio"] < 1.05


def test_independent_checks_point_the_same_way_as_the_registered_comparison():
    summary = load("s7_summary.json")
    chosen = summary["registered"]
    chk = summary["independent_checks"]
    # 14-50 keV power laws (no absorption model) against Delta1; XRT count rate ratio against Delta2
    assert (chk["flux_14_50_powerlaw"]["value"] > 0) == (chosen["delta1"]["value"] > 0)
    assert (chk["xrt_count_rate_ratio_2_10"]["log10"] > 0) == (chosen["delta2"]["value"] > 0)


def test_synthetic_validation_flags_follow_the_registered_pass_conditions():
    if not (RES / "s7_sim_validation.json").exists():
        pytest.skip("synthetic validation not present")
    sims = load("s7_sim_validation.json")
    assert sims["n_sims"] >= 200 and sims["model"] == load("s7_fits.json")["chosen_model"]
    b, cov = sims["bias"], sims["coverage_90"]
    assert sims["pass"]["nh"] == (abs(b["nh_fraction"]) <= 0.10 and 0.80 <= cov["nh"] <= 0.97)
    assert sims["pass"]["gamma"] == (abs(b["gamma"]) <= 0.1 and 0.80 <= cov["gamma"] <= 0.97)
    assert sims["pass"]["delta2"] == (abs(b["delta2_dex"]) <= 0.05 and 0.80 <= cov["delta2"] <= 0.97)


def test_no_local_paths_in_the_committed_stage_7_results():
    for path in RES.glob("s7_*.json"):
        text = path.read_text()
        assert "/Users/" not in text and "/var/folders" not in text, path.name
