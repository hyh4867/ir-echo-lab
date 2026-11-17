import numpy as np
import pytest

from irecho.config import SlowRiseCriteria
from irecho.slowrise import (detect_slow_rise, mann_kendall_z, model_curve, theil_sen_slope)

MJD = 56700.0 + 182.625 * np.arange(21)
YEARS = (MJD - MJD[0]) / 365.25
SIGMA = 0.01


def noise(seed):
    rng = np.random.default_rng(seed)
    return 10.0 + rng.normal(0, SIGMA, MJD.size), 9.5 + rng.normal(0, SIGMA, MJD.size)


def run(w1, w2, **kw):
    err = np.full(MJD.size, SIGMA)
    return detect_slow_rise(MJD, w1, err, w2, err, **kw)


def ramp(start_index, amp, rise_years):
    return amp * np.clip((YEARS - YEARS[start_index]) / rise_years, 0.0, 1.0)


def test_mann_kendall_hand_computed():
    # n = 10 strictly increasing: S = 45, var = 125, Z = 44 / sqrt(125)
    assert mann_kendall_z(np.arange(10.0)) == pytest.approx(44 / np.sqrt(125))
    assert mann_kendall_z(-np.arange(10.0)) == pytest.approx(-44 / np.sqrt(125))
    assert mann_kendall_z(np.ones(10)) == 0.0


def test_theil_sen_exact_on_a_line_and_robust_to_an_outlier():
    y = 3.0 - 0.02 * YEARS
    assert theil_sen_slope(YEARS, y) == pytest.approx(-0.02)
    y[5] += 1.0
    assert theil_sen_slope(YEARS, y) == pytest.approx(-0.02, abs=2e-3)


def test_no_false_alarms_on_white_noise():
    assert sum(run(*noise(s)).flagged for s in range(200)) == 0


def test_f_test_alone_fires_rarely_on_white_noise():
    """Without the amplitude floor, the scan over onsets makes p <= 1e-3 in
    both bands rare but not impossible; it must stay below 1%."""
    hits = 0
    for s in range(300):
        r = run(*noise(s))
        hits += r.p["w1"] <= 1e-3 and r.p["w2"] <= 1e-3
    assert hits <= 3


def test_linear_ramp_is_flagged_with_its_onset():
    w1, w2 = noise(1)
    d = ramp(10, 0.3, 3.0)
    r = run(w1 - d, w2 - d)
    assert r.flagged and r.model == "hinge"
    assert abs(r.onset_index - 11) <= 1  # first rising visit is index 11
    assert r.dm_fit["w1"] == pytest.approx(0.3, abs=0.08)
    assert r.dm_cum["w2"] == pytest.approx(0.3, abs=0.03)
    assert r.mk_z["w1"] > 3


def test_whole_curve_linear_trend_is_flagged():
    w1, w2 = noise(2)
    d = 0.02 * YEARS  # 0.2 mag in 10 yr, from the first visit
    r = run(w1 - d, w2 - d)
    assert r.flagged
    assert r.theil_sen_mag_per_yr["w1"] == pytest.approx(0.02, abs=0.004)


def test_exponential_rise_is_flagged():
    w1, w2 = noise(3)
    f = 1.0 + 0.5 * np.exp((YEARS - YEARS[-1]) / 1.0)
    d = 2.5 * np.log10(f)
    r = run(w1 - d, w2 - d)
    assert r.flagged
    assert r.dm_fit["w1"] == pytest.approx(2.5 * np.log10(1.5), abs=0.08)


def test_single_band_rise_is_not_flagged():
    w1, w2 = noise(4)
    assert not run(w1 - ramp(8, 0.4, 3.0), w2).flagged


def test_dimming_is_not_flagged():
    w1, w2 = noise(5)
    d = ramp(8, 0.4, 3.0)
    r = run(w1 + d, w2 + d)
    assert not r.flagged
    assert r.mk_z["w1"] < -3


def test_significant_but_small_rise_fails_the_amplitude_floor():
    w1, w2 = noise(6)
    d = ramp(8, 0.06, 3.0)
    r = run(w1 - d, w2 - d)
    assert not r.flagged
    assert r.p["w1"] < 1e-3 and r.p["w2"] < 1e-3


def test_jump_in_the_last_visit_is_not_a_slow_rise():
    """One elevated visit at the end is detector (a)'s case: the cumulative
    change uses the median of the last three visits."""
    w1, w2 = noise(7)
    w1[-1] -= 0.6
    w2[-1] -= 0.6
    assert not run(w1, w2).flagged


def test_too_few_visits():
    err = np.full(7, SIGMA)
    r = detect_slow_rise(MJD[:7], np.full(7, 10.0), err, np.full(7, 9.5), err)
    assert not r.flagged and r.model is None


def test_nan_visits_are_ignored_and_order_is_checked():
    w1, w2 = noise(8)
    d = ramp(10, 0.3, 3.0)
    w1, w2 = w1 - d, w2 - d
    w1[3] = np.nan
    assert run(w1, w2).n_visits == 20
    err = np.full(MJD.size, SIGMA)
    with pytest.raises(ValueError):
        detect_slow_rise(MJD[::-1], w1[::-1], err, w2[::-1], err)


def test_thresholds_are_parameters():
    w1, w2 = noise(9)
    d = ramp(8, 0.06, 3.0)
    assert run(w1 - d, w2 - d, criteria=SlowRiseCriteria(min_cum_mag=0.03)).flagged


def test_model_curve_follows_the_ramp():
    w1, w2 = noise(10)
    d = ramp(10, 0.3, 3.0)
    r = run(w1 - d, w2 - d)
    fit = model_curve(MJD, w1 - d, np.full(MJD.size, SIGMA), r)
    assert np.std(fit - (10.0 - d)) < 0.03
