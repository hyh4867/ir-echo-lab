"""Stage 7 unit tests: decision rules of docs/XRAY_PLAN.md and numerical helpers.

No HEASoft and no data are needed.
"""

import math

import numpy as np
import pytest

from irecho import physics, xray


def test_registered_numbers_are_the_ones_in_the_plan():
    plan = (xray.config.PROJECT_ROOT / "docs" / "XRAY_PLAN.md").read_text()
    assert xray.DELTA_L_MIN == 0.18 and "Δ ≥ +0.18" in plan
    assert xray.DELTA_L_SAME_BOUND == 0.30 and "hi < +0.30" in plan
    assert xray.DELTA_NH_MIN == 0.10 and "\\|ΔN\\| ≥ 0.10" in plan
    assert xray.DELTA_NH_SAME_BOUND == 0.20 and "hi < +0.20" in plan
    assert xray.SYS_NUSTAR_BAT_DEX == 0.06 and "0.06 dex(15%)" in plan
    assert xray.SYS_XRT_NUSTAR_DEX == 0.04 and "0.04 dex(10%)" in plan
    assert xray.AIC_STEP == 4.0 and "4 이상 줄어들 때만" in plan
    assert xray.MIN_XRT_COUNTS == 25 and "25 미만" in plan
    assert xray.LOG_L_DUST_PEAK == 43.11
    if xray.LOG_LBOL_BASS is not None:                      # blank in the public copy (catalogue value)
        assert 43.0 < xray.LOG_LBOL_BASS < 45.0             # log erg/s, not a luminosity in other units
    assert (xray.F_NEED_OK, xray.F_NEED_TIGHT) == (0.6, 1.0)
    for obsid in xray.OBS.values():
        assert obsid in plan


def test_luminosity_label_follows_section_6_1():
    f = xray.luminosity_label
    assert f(0.30, 0.10, 0.50) == xray.HIGH
    assert f(0.18, 0.01, 0.35) == xray.HIGH
    assert f(0.30, -0.05, 0.60) == xray.UNDECIDED      # interval reaches zero
    assert f(0.17, 0.05, 0.29) == xray.SAME            # significant but below a factor 1.5
    assert f(0.17, 0.05, 0.31) == xray.UNDECIDED       # upper bound outside +-0.30
    assert f(0.00, -0.29, 0.29) == xray.SAME
    assert f(0.00, -0.40, 0.40) == xray.UNDECIDED
    assert f(-0.25, -0.50, -0.02) == xray.LOW
    assert f(-0.25, -0.50, 0.02) == xray.UNDECIDED
    assert f(float("nan"), 0.0, 1.0) == xray.UNDECIDED
    inf = float("inf")
    assert f(0.50, 0.20, inf) == xray.HIGH             # open upper side: 'higher' needs only the lower bound
    assert f(0.10, -inf, 0.25) == xray.UNDECIDED       # open lower side: not 'same'
    assert f(-0.50, -inf, -0.20) == xray.LOW


def test_absorption_label_follows_section_6_2():
    f = xray.absorption_label
    assert f(0.15, 0.02, 0.30) == xray.CHANGED
    assert f(-0.40, -0.70, -0.10) == xray.CHANGED
    assert f(0.05, 0.01, 0.09) == xray.SAME            # excludes zero but smaller than 0.10 dex
    assert f(0.05, -0.10, 0.19) == xray.SAME
    assert f(0.05, -0.10, 0.25) == xray.UNDECIDED
    assert f(0.15, -0.02, 0.30) == xray.UNDECIDED
    assert f(0.15, 0.02, 0.30, measurable=False) == xray.UNDECIDED
    assert f(0.50, 0.20, float("inf")) == xray.CHANGED
    assert f(0.05, -float("inf"), 0.15) == xray.UNDECIDED


def test_overall_verdict_follows_section_6_3():
    v = lambda *a: xray.overall_verdict(*a)["verdict"]  # noqa: E731
    H, L, S, U, C = xray.HIGH, xray.LOW, xray.SAME, xray.UNDECIDED, xray.CHANGED
    assert v(H, S, S) == "H1 지지"
    assert v(U, H, S) == "H1 지지"
    assert v(H, U, U) == "H1 지지, 가림 변화는 배제 못 함"
    assert v(H, H, C) == "증광과 가림 변화 둘 다"
    assert v(S, S, C) == "H2 지지"
    assert v(U, S, C) == "H2 지지"
    assert v(H, L, S) == U                              # one comparison is 'lower'
    assert v(U, U, C) == U                              # H2 needs the 2010 luminosity to be 'same'
    assert v(S, S, S) == U
    assert "변화를 찾지 못함" in xray.overall_verdict(S, S, S)["note"]
    assert xray.overall_verdict(S, L, U)["note"] == "광도가 낮다는 비교가 있음"
    assert xray.overall_verdict(S, U, U)["note"] == "정밀도 부족"


def test_energy_budget_labels_and_the_stage_6a_ratio():
    if xray.LOG_LBOL_BASS is None:
        b = xray.energy_budget(0.0, log_lint_2022=42.5)
        assert b["f_need"] is None and b["label"] is None
        assert b["log_lbol_2022_from_kappa20"] == pytest.approx(42.5 + math.log10(20))
        pytest.skip("the BASS luminosity is blank in this copy (scripts/restore_catalog_columns.py fills it)")
    b = xray.energy_budget(0.0)
    ratio = 10 ** (xray.LOG_L_DUST_PEAK - xray.LOG_LBOL_BASS)
    assert b["f_need"] == pytest.approx(ratio) and b["label"] == "맞음"
    assert b["f_need"] == pytest.approx(0.15, abs=0.01)             # stage 6A: 0.15
    assert b["covering_factor_from_excess_only"] is None
    b = xray.energy_budget(math.log10(1.5))
    assert b["covering_factor_from_excess_only"] == pytest.approx(ratio / 0.5, rel=0.01)
    assert xray.energy_budget(-0.7)["label"] == "빠듯함"
    assert xray.energy_budget(-1.0)["label"] == "안 맞음"
    b = xray.energy_budget(0.0, log_lint_2022=42.5)
    assert b["log_lbol_2022_from_kappa20"] == pytest.approx(42.5 + math.log10(20))


def test_choose_model_steps_only_on_a_large_enough_aic_drop():
    assert xray.choose_model({"M0": 100.0, "M1": 97.0, "M2": 80.0}) == "M0"   # M1 not adopted, so M2 is not tried
    assert xray.choose_model({"M0": 100.0, "M1": 96.0, "M2": 93.0}) == "M1"
    assert xray.choose_model({"M0": 100.0, "M1": 90.0, "M2": 86.0}) == "M2"


def test_systematic_is_added_in_quadrature_at_90_percent():
    lo, hi = xray.add_systematic(0.2, 0.1, 0.35, 0.06)
    s = 1.645 * 0.06
    assert lo == pytest.approx(0.2 - math.hypot(0.1, s)) and hi == pytest.approx(0.2 + math.hypot(0.15, s))
    assert xray.add_systematic(0.2, 0.1, 0.35, 0.0) == pytest.approx((0.1, 0.35))


def test_powerlaw_flux_and_rest_frame_luminosity():
    # numerical integral of E N(E)
    e = np.linspace(2.0, 10.0, 200001)
    for gamma in (1.4, 1.8, 2.0, 2.3):
        num = np.trapezoid(3e-3 * e ** (1 - gamma), e) * 1.602176634e-9
        assert xray.powerlaw_band_flux(3e-3, gamma, 2.0, 10.0) == pytest.approx(num, rel=1e-6)
    # zpowerlw at z = 0 is a plain power law
    d = 100.0
    lum0 = xray.zpowerlw_rest_luminosity(1e-3, 1.8, 0.0, d)
    assert lum0 == pytest.approx(4 * math.pi * (d * 3.0856775814913673e24) ** 2
                                 * xray.powerlaw_band_flux(1e-3, 1.8, 2, 10))
    # with redshift: the observed band shrinks by (1+z) and so does the photon energy
    z = 0.0256
    lum = xray.zpowerlw_rest_luminosity(1e-3, 1.8, z, d)
    assert lum == pytest.approx(lum0 / (1 + z) ** 2, rel=1e-9)
    assert physics.luminosity_distance_mpc(z) == pytest.approx(116.07, abs=0.01)


def test_net_rate_hardness_ratio_and_coverage():
    rate, err = xray.net_rate(100, 50, 0.1, 1000.0)
    assert rate == pytest.approx(0.095) and err == pytest.approx(math.sqrt(100 + 0.5) / 1000)
    hr, hr_err = xray.hardness_ratio(0.03, 0.01, 0.003, 0.001)
    assert hr == pytest.approx(0.5)
    rng = np.random.default_rng(3)
    h, s = 0.03 + 0.003 * rng.standard_normal(200000), 0.01 + 0.001 * rng.standard_normal(200000)
    assert hr_err == pytest.approx(np.std((h - s) / (h + s)), rel=0.03)
    assert xray.interval_coverage(1.0, [0.5, 1.1, np.nan, 0.9], [1.5, 1.2, 2.0, 1.0]) == pytest.approx(2 / 3)


# --- stage 7b helpers -------------------------------------------------------------------------

BAT_LO = [14.0, 20.0, 24.0, 35.0, 50.0, 75.0, 100.0, 150.0]
BAT_HI = [20.0, 24.0, 35.0, 50.0, 75.0, 100.0, 150.0, 195.0]


def test_bat_channels_are_chosen_by_number_and_keep_the_edge_channels():
    assert xray.bat_channels_in_band(BAT_LO, BAT_HI, 14.0, 50.0) == [1, 2, 3, 4]      # channel 1 starts at 14.0
    assert xray.bat_channels_in_band(BAT_LO, BAT_HI, 14.0, 195.0) == list(range(1, 9))
    assert xray.bat_channels_in_band(BAT_LO, BAT_HI, 15.0, 50.0) == [2, 3, 4]         # only whole channels
    assert xray.bat_channels_in_band(BAT_LO, BAT_HI, 14.0, 60.0) == [1, 2, 3, 4]


def test_binomial_standard_error_of_a_coverage_fraction():
    assert xray.binomial_se(0.9, 200) == pytest.approx(0.0212, abs=1e-4)
    assert (0.90 - 0.855) / xray.binomial_se(0.9, 200) == pytest.approx(2.12, abs=0.01)


def test_autocorrelation_time_and_gelman_rubin_on_known_chains():
    rng = np.random.default_rng(11)
    white = rng.standard_normal((4000, 12))
    assert xray.autocorr_time(white) == pytest.approx(1.0, abs=0.15)
    assert xray.gelman_rubin(white) == pytest.approx(1.0, abs=0.01)
    phi = 0.9                                             # AR(1): tau = (1 + phi) / (1 - phi) = 19
    ar = np.zeros((20000, 12))
    for i in range(1, ar.shape[0]):
        ar[i] = phi * ar[i - 1] + rng.standard_normal(12)
    assert xray.autocorr_time(ar) == pytest.approx(19.0, rel=0.15)
    shifted = white + np.arange(12) * 1.0                 # walkers that sit at different values
    assert xray.gelman_rubin(shifted) > 2.0


def test_poisson_net_samples_handle_a_non_detection_and_a_detection():
    rng = np.random.default_rng(5)
    s = xray.poisson_net_samples(8, 196, 0.0471, 200000, rng)        # fewer counts than the background predicts
    assert s.min() >= 0.0
    assert 3.0 < np.percentile(s, 90) < 8.0
    h = xray.poisson_net_samples(111, 51, 0.0471, 200000, rng)       # clear detection
    assert np.median(h) == pytest.approx(111 - 51 * 0.0471, rel=0.03)
    assert np.std(h) == pytest.approx(math.sqrt(111), rel=0.08)


def test_gti_overlap():
    assert xray.gti_overlap_seconds([(0, 10), (20, 30)], [(5, 25)]) == 10.0
    assert xray.gti_overlap_seconds([(0, 10)], [(10, 20)]) == 0.0
