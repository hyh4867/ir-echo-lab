"""Stage 8 unit tests: null light curves, direct branch comparison, driver prescriptions."""

import math

import numpy as np
import pytest

from irecho import breath, hyst8


def _loop_curve(rng, n=22, loop=0.0, noise=0.01):
    """A smooth brightening and fading; `loop` shifts W1 on the rising branch."""
    t = np.linspace(0, 1, n)
    f2 = 1.0 + 0.8 * np.exp(-0.5 * ((t - 0.5) / 0.18) ** 2)
    rising = np.gradient(f2) > 0
    f1 = 0.7 + 0.6 * (f2 - 1.0) + loop * rising * (f2 - 1.0)
    s1, s2 = np.full(n, noise), np.full(n, noise)
    return f1 + noise * rng.standard_normal(n), s1, f2 + noise * rng.standard_normal(n), s2, rising


def test_branch_definition_is_the_one_of_stage_5b():
    rng = np.random.default_rng(1)
    f1, s1, f2, s2, _ = _loop_curve(rng)
    prep = breath._prepare(f1, s1, f2, s2, breath.SETTINGS)
    assert np.array_equal(hyst8.branch_rising(f1, f2), prep["rising"])


def test_eiv_fit_recovers_slope_and_offset_with_errors_in_both_variables():
    rng = np.random.default_rng(2)
    slopes, offsets, pulls = [], [], []
    for _ in range(300):
        x_true = rng.uniform(1.0, 3.0, 40)
        lab = (rng.random(40) < 0.5).astype(float)
        y_true = 0.5 + 0.8 * x_true + 0.05 * lab
        x = x_true + 0.05 * rng.standard_normal(40)
        y = y_true + 0.03 * rng.standard_normal(40)
        fit = hyst8.eiv_fit(x, y, np.full(40, 0.05), np.full(40, 0.03), 0.0, [lab])
        slopes.append(fit["theta"][1])
        offsets.append(fit["theta"][2])
        pulls.append((fit["theta"][2] - 0.05) / math.sqrt(fit["cov"][2, 2]))
    assert np.mean(slopes) == pytest.approx(0.8, abs=0.01)        # ordinary least squares would give about 0.79
    assert np.mean(offsets) == pytest.approx(0.05, abs=0.004)
    assert np.std(pulls) == pytest.approx(1.0, abs=0.15)
    assert hyst8.eiv_fit([1, 2], [1, 2], [0.1, 0.1], [0.1, 0.1]) is None


def test_null_truths_are_single_valued_curves_without_a_loop():
    rng = np.random.default_rng(3)
    f1, s1, f2, s2, rising = _loop_curve(rng, loop=0.3)
    lin = hyst8.truth_linear(f1, s1, f2, s2)
    assert np.allclose(np.diff(lin) / np.diff(f2), (lin[-1] - lin[0]) / (f2[-1] - f2[0]), rtol=1e-6)
    one_t = hyst8.truth_single_temperature(f1, s1, f2, s2)
    order = np.argsort(f2)
    assert np.all(np.diff(one_t[order]) >= -1e-12)                # monotonic in W2: one value of W1 per W2
    prep = breath._prepare(f1, s1, f2, s2, breath.SETTINGS)
    pos = f2 - prep["a"][1] > 0
    assert np.sum(one_t[pos] - prep["a"][0]) == pytest.approx(np.sum(f1[pos] - prep["a"][0]), rel=1e-3)
    assert hyst8.smooth3(np.array([1.0, 2.0, 6.0, 2.0])) == pytest.approx([1.5, 3.0, 10 / 3, 4.0])


def test_noise_has_the_requested_correlation_and_band_switch():
    rng = np.random.default_rng(4)
    t = np.zeros(200000)
    a, b = hyst8.add_noise(rng, t, np.ones(t.size), t, 2 * np.ones(t.size), r=0.6)
    assert np.corrcoef(a, b)[0, 1] == pytest.approx(0.6, abs=0.01)
    assert np.std(a) == pytest.approx(1.0, abs=0.01) and np.std(b) == pytest.approx(2.0, abs=0.02)
    a, b = hyst8.add_noise(rng, t, np.ones(t.size), t, np.ones(t.size), bands="w1")
    assert np.std(a) > 0.9 and np.all(b == 0)


def test_d_statistic_with_default_switches_equals_breath_measure():
    rng = np.random.default_rng(5)
    f1, s1, f2, s2, _ = _loop_curve(rng, loop=0.2, noise=0.004)
    m = breath.measure(f1, s1, f2, s2, with_in_place=False)
    d = hyst8.d_statistic(f1, s1, f2, s2)
    assert m["n_rising"] >= 4 and m["n_declining"] >= 4
    assert d["d"] == pytest.approx(m["d_hyst"], rel=1e-9) and d["sigma_d"] == pytest.approx(m["sigma_d"], rel=1e-9)
    assert d["measurable"] == m["hyst_measurable"]
    alt = hyst8.d_statistic(f1, s1, f2, s2, pivot="mean")
    assert np.isfinite(alt["d"]) and alt["d"] != pytest.approx(d["d"], rel=1e-6)


def test_direct_offset_is_zero_without_a_loop_and_recovers_an_injected_one():
    rng = np.random.default_rng(6)
    none, loop = [], []
    for _ in range(200):
        f1, s1, f2, s2, rising = _loop_curve(rng, loop=0.0, noise=0.006)
        r0 = hyst8.direct_offset(f1, s1, f2, s2, rising)
        f1, s1, f2, s2, rising = _loop_curve(rng, loop=0.15, noise=0.006)
        r1 = hyst8.direct_offset(f1, s1, f2, s2, rising)
        if r0["usable"]:
            none.append(r0["delta"] / r0["sigma_delta"])
        if r1["usable"]:
            loop.append(r1["delta"])
    assert len(none) > 150 and abs(np.mean(none)) < 0.25 and np.std(none) == pytest.approx(1.0, abs=0.25)
    assert np.median(loop) > 0.02                                 # rising branch has more W1: delta > 0
    r1 = hyst8.direct_offset(f1, s1, f2, s2, rising)
    assert r1["d_colour"] < 0 and r1["lo"] < r1["hi"] and r1["n_rise_in"] >= 3 and r1["n_fade_in"] >= 3


def test_direct_offset_refuses_branches_that_do_not_overlap():
    f2 = np.array([1.0, 1.1, 1.2, 1.3, 2.0, 2.1, 2.2, 2.3])
    rising = np.array([True] * 4 + [False] * 4)
    r = hyst8.direct_offset(0.5 * f2, np.full(8, 0.01), f2, np.full(8, 0.01), rising)
    assert not r["usable"] and math.isnan(r["delta"])


def test_driver_prescriptions_at_their_reference_points():
    log_lbol_ref = (90.596 + 23.0) / 2.5            # M_i = -23
    assert hyst8.absolute_i_mag(log_lbol_ref) == pytest.approx(-23.0)
    sigma, tau = hyst8.macleod_drw(log_lbol_ref, 9.0)
    assert sigma == pytest.approx(10 ** -0.51 / math.sqrt(2)) and tau == pytest.approx(10 ** 2.4)
    s2, t2 = hyst8.macleod_drw(log_lbol_ref - 0.4, 8.0)   # one magnitude fainter, ten times lighter
    assert s2 / sigma == pytest.approx(10 ** (0.131 - 0.18)) and t2 / tau == pytest.approx(10 ** (0.03 - 0.21))
    assert hyst8.burke_tau(8.0) == pytest.approx(107.0) and hyst8.burke_tau(9.0) == pytest.approx(107.0 * 10 ** 0.38)
