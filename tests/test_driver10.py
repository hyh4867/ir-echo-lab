"""Stage 10: the optical-driver fit (src/irecho/driver10.py) on made-up light curves."""

from __future__ import annotations

import numpy as np
import pytest

from irecho import driver10 as d10
from irecho import hyst8
from irecho import transfer9 as t9


@pytest.fixture(scope="module")
def setup():
    rng = np.random.default_rng(3)
    nights = np.sort(np.concatenate([58200 + 365.25 * k + rng.uniform(0, 230, 60) for k in range(8)]))
    visits = 56700 + 182.0 * np.arange(22) + rng.normal(0, 2, 22)
    s_opt = np.full(nights.size, 0.015)
    hyper = {"var": 0.02, "tau": 300.0, "jitter_var": 1e-8}
    prep = d10.prepare(nights, s_opt, hyper, visits, joint=True)
    return {"nights": nights, "visits": visits, "s_opt": s_opt, "hyper": hyper, "prep": prep,
            "grid": d10.truth_grid(nights, visits)}


def make(setup, lag1, lag2, shape="f1", alpha1=0.0, seed=0, noise=True):
    rng = np.random.default_rng(seed)
    grid, nights, visits = setup["grid"], setup["nights"], setup["visits"]
    mags = t9.drw_magnitudes(rng, grid.size, d10.TRUTH_STEP_DAYS, 300.0)
    ell = 10.0 ** (-0.4 * d10.lognormal_amp(0.15) * mags)
    d1, d2 = d10.linear_truth(ell, lag1, lag2, shape, alpha1)
    f1 = d10.scale_band(np.interp(visits, grid, d1), 0.08, 1.0)
    f2 = d10.scale_band(np.interp(visits, grid, d2), 0.08, 1.0)
    y = np.interp(nights, grid, ell / ell.mean())
    s = np.full(visits.size, 0.012)
    if noise:
        y = y + rng.normal(0, setup["s_opt"])
        f1, f2 = f1 + rng.normal(0, s), f2 + rng.normal(0, s)
    return y, f1, f2, s


def test_grid_points_stay_inside_the_window():
    pts = d10.grid_points()
    assert len(pts) == 564
    for shape, tau in pts:
        lo, hi = d10.support(shape, tau)
        assert lo >= -d10.SUPPORT_FUTURE_DAYS - 1e-9 and hi <= d10.SUPPORT_PAST_DAYS + 1e-9
    base = [tau for shape, tau in pts if shape == d10.BASE]
    assert min(base) == -65.0 and max(base) == 465.0


def test_kernel_rows_sum_to_one_and_have_the_mean_lag():
    g0, h, n = 0.0, 2.0, 1001
    grid = g0 + h * np.arange(n)
    for shape in d10.SHAPES:
        for tau in (-40.0, 0.0, 35.0, 250.0):
            c = d10.kernel_rows(g0, h, n, [1200.0], tau, shape)
            assert c.sum() == pytest.approx(1.0, abs=1e-9)
            assert float((c @ grid)[0]) == pytest.approx(1200.0 - tau, abs=1e-6)       # exact for a linear driver
    with pytest.raises(ValueError):
        d10.kernel_rows(g0, h, n, [100.0], 300.0, "f1")


def test_top_hat_average_is_exact_for_a_piecewise_linear_driver():
    g0, h, n = 0.0, 2.0, 501
    grid = g0 + h * np.arange(n)
    rng = np.random.default_rng(0)
    d = rng.normal(size=n)
    c = d10.kernel_rows(g0, h, n, [700.3], 100.0, "f1")        # window of lags [50, 150]
    fine = np.linspace(700.3 - 150.0, 700.3 - 50.0, 200001)
    assert float((c @ d)[0]) == pytest.approx(float(np.mean(np.interp(fine, grid, d))), abs=2e-5)


def test_gp_gain_reproduces_a_constant_and_shrinks_the_variance_at_the_nights(setup):
    nights, s_opt, hyper = setup["nights"], setup["s_opt"], setup["hyper"]
    grid = d10.driver_grid(nights[0], nights[-1])
    gain, p = d10.gp_driver(nights, s_opt, hyper, grid)
    assert np.allclose(gain @ np.full(nights.size, 3.7), 3.7)                       # flat prior on the mean
    var = np.diag(p)
    assert var.min() > 0
    near = np.abs(grid[:, None] - nights[None, :]).min(axis=1)
    assert np.median(var[near < 1.0]) < 0.2 * np.median(var[near > 40.0])        # seasonal gaps are less certain
    assert var.max() < 1.5 * hyper["var"]


def test_usable_visits_need_the_whole_window_and_no_long_gap(setup):
    nights, visits = setup["nights"], setup["visits"]
    ok = d10.usable_visits(visits, nights)
    assert ok.sum() == setup["prep"]["idx"].size >= d10.MIN_VISITS
    assert np.all(visits[ok] - d10.SUPPORT_PAST_DAYS >= nights[0]) and np.all(visits[ok] + d10.SUPPORT_FUTURE_DAYS <= nights[-1])
    holed = nights[(nights < 59300) | (nights > 59300 + 300)]
    ok2 = d10.usable_visits(visits, holed)
    assert ok2.sum() < ok.sum()


def test_noise_free_lags_are_recovered(setup):
    y, f1, f2, s = make(setup, 100.0, 140.0, noise=False)
    a = d10.fit_band(setup["prep"], y, f1, s)
    b = d10.fit_band(setup["prep"], y, f2, s)
    assert abs(a["tau"] - 100.0) <= 10.0 and abs(b["tau"] - 140.0) <= 10.0
    assert a["r_w"] > 0.95 and a["ds_c1"] > 50 and a["ds_01"] > 0
    assert a["tau_q16"] <= a["tau_q50"] <= a["tau_q84"]


def test_sample_of_noisy_fits_has_small_bias(setup):
    tau1, tau2 = [], []
    for seed in range(12):
        y, f1, f2, s = make(setup, 100.0, 140.0, seed=100 + seed)
        tau1.append(d10.fit_band(setup["prep"], y, f1, s, full=False)["tau"])
        tau2.append(d10.fit_band(setup["prep"], y, f2, s, full=False)["tau"])
    assert abs(np.median(tau1) - 100.0) <= 15.0 and abs(np.median(tau2) - 140.0) <= 15.0
    assert 15.0 <= np.median(np.array(tau2) - np.array(tau1)) <= 65.0


def test_undelayed_share_is_seen_by_the_model_with_alpha(setup):
    alphas, tau_plain, tau_alpha = [], [], []
    for seed in range(12):
        y, f1, f2, s = make(setup, 150.0, 190.0, alpha1=0.3, seed=200 + seed)
        a = d10.fit_band(setup["prep"], y, f1, s)
        alphas.append(a["alpha"])
        tau_plain.append(a["tau"])
        tau_alpha.append(a["tau_ma"])
        assert a["ds_1a"] >= -1e-6 or a["tau_ma"] != a["tau"]
    assert 0.15 <= np.median(alphas) <= 0.45
    assert np.median(tau_plain) < np.median(tau_alpha)            # the undelayed share pulls the plain lag down


def test_branch_offset_is_zero_for_a_fixed_kernel_and_never_makes_the_fit_worse(setup):
    y, f1, f2, s = make(setup, 100.0, 140.0, seed=5)
    rising = hyst8.branch_rising(f1, f2)
    a = d10.fit_band(setup["prep"], y, f1, s, rising)
    assert a["ds_1e"] >= -1e-6
    assert abs(a["eta"]) < 5.0 * a["sigma_eta"] + 0.02


def test_joint_fit_gives_an_interval_for_the_lag_difference(setup):
    y, f1, f2, s = make(setup, 100.0, 140.0, seed=7)
    a = d10.fit_band(setup["prep"], y, f1, s)
    b = d10.fit_band(setup["prep"], y, f2, s)
    j = d10.fit_joint(setup["prep"], y, f1, s, f2, s, 0.3, a, b)
    assert j["dtau_q025"] <= j["dtau_q16"] <= j["dtau_q50"] <= j["dtau_q84"] <= j["dtau_q975"]
    assert j["dtau_q025"] - 10.0 <= 40.0 <= j["dtau_q975"] + 10.0
    assert abs(j["tau1_joint"] - a["tau"]) <= 30.0


def test_colour_offset_of_the_model_has_the_sign_of_the_lag(setup):
    pred, obs = [], []
    for seed in range(30):
        y, f1, f2, s = make(setup, 100.0, 180.0, seed=300 + seed)
        a = d10.fit_band(setup["prep"], y, f1, s)
        b = d10.fit_band(setup["prep"], y, f2, s)
        off = d10.colour_offsets(f1, s, f2, s, 0.0, a, b)
        if off["usable_obs"] and off["usable_pred"]:
            obs.append(off["d_colour_obs"])
            pred.append(off["d_colour_pred"])
    assert len(obs) >= 5 and np.median(obs) < 0 and np.median(pred) < 0     # W2 lags W1: bluer when rising


def test_truth_kernels_have_unit_sum_and_the_stated_mean_lag():
    for shape in ("f0", "f05", "f1", "f2", "gwide", "exp"):
        k = d10.kernel_array(shape, 120.0)
        assert k.sum() == pytest.approx(1.0, abs=1e-9)
        tol = 4.0 if shape in ("gwide", "exp") else 0.6
        assert d10.kernel_mean_lag(k) == pytest.approx(120.0, abs=tol)
    assert d10.kernel_array("f1", 0.0).size == 1


def test_nightly_binning_and_clipping():
    mjd = np.array([58000.30, 58000.31, 58000.33, 58001.30, 58003.3, 58003.31])
    mag = np.array([15.0, 15.02, 14.98, 15.1, 15.0, 15.0])
    out = d10.nightly(mjd, mag, np.full(6, 0.01))
    assert out["t"].size == 3 and out["s"][0] == pytest.approx(out["s"][1] / np.sqrt(3), rel=0.15)
    t = np.arange(100.0)
    y = np.ones(100) + 0.01 * np.sin(t / 10.0)
    y[40] = 1.5
    keep = d10.clip_nights(t, y)
    assert not keep[40] and keep.sum() == 99


def test_literature_relations():
    w1, w2 = d10.lyu_lags_days(44.5)
    assert w1 == pytest.approx(115.0, abs=1.0) and w2 == pytest.approx(145.0, abs=1.0)
    m1, m2 = d10.mandal_lags_days(44.5)
    assert m1 == pytest.approx(156.0, abs=1.0) and m2 == pytest.approx(231.0, abs=1.5)
