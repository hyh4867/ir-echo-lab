"""Stage 9: band-to-band transfer (src/irecho/transfer9.py). Synthetic inputs only."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import quad

from irecho import config, hyst8
from irecho import transfer9 as t9

T = np.array([56643., 56831, 57008, 57193, 57367, 57557, 57726, 57924, 58090, 58288, 58450,
              58652, 58814, 59019, 59182, 59384, 59546, 59748, 59910, 60115, 60300, 60479])


def test_numbers_fixed_in_the_plan():
    plan = (config.PROJECT_ROOT / "docs" / "PREDICTIONS_9.md").read_text()
    assert t9.SUPPORT_PAST_DAYS == 350.0 and "350일" in plan
    assert t9.SUPPORT_FUTURE_DAYS == 150.0 and "150일" in plan
    assert t9.MAX_GAP_DAYS == 250.0 and "250일" in plan
    assert t9.TAU_STEP_DAYS == 5.0 and t9.WIDTHS_DAYS == (0.0, 100.0, 200.0, 300.0) and "0, 100, 200, 300일" in plan
    assert t9.MIN_VISITS_FIT == 10 and "10개 이상" in plan
    assert t9.WELL_DS == 9.0 and t9.WELL_HALF_WIDTH_DAYS == 45.0 and "ΔS(T0→T1) ≥ 9" in plan and "45일 이하" in plan
    assert set(t9.TRUTHS) == {"T0", "T1_20", "T1_40", "T1_80", "T1_160", "T1_lyu", "T2_40", "T2_80", "T0b", "T3_1yr", "T3_3yr"}
    for name in t9.TRUTHS:
        assert f"`{name}`" in plan


def test_interpolation_weights_are_exact_and_refuse_extrapolation():
    rng = np.random.default_rng(0)
    f = rng.normal(size=T.size)
    u = np.array([56643.0, 57000.0, 58333.3, 60479.0])
    assert np.allclose(t9.interp_weights(T, u) @ f, np.interp(u, T, f))
    for bad in (T[0] - 1.0, T[-1] + 1.0):
        with pytest.raises(ValueError):
            t9.interp_weights(T, [bad])
        with pytest.raises(ValueError):
            t9.antiderivative_weights(T, [bad])
    with pytest.raises(ValueError):
        t9.kernel_weights(T, T[:1], 10.0, 0.0)          # the first visit has no past


def test_window_average_is_the_exact_integral_of_the_interpolant():
    rng = np.random.default_rng(1)
    f = rng.normal(size=T.size)
    for tau, width in ((60.0, 200.0), (-40.0, 100.0), (175.0, 300.0)):
        c = t9.kernel_weights(T, T[5:8], tau, width)
        assert np.allclose(c.sum(axis=1), 1.0)
        for row, t_eval in zip(c, T[5:8]):
            lo, hi = t_eval - tau - width / 2, t_eval - tau + width / 2
            exact = quad(lambda x: np.interp(x, T, f), lo, hi, points=list(T[(T > lo) & (T < hi)]), limit=200)[0] / width
            assert row @ f == pytest.approx(exact, abs=1e-9)
    # a straight line in time is shifted by any symmetric window exactly
    line = 3.0 + 0.002 * (T - T[0])
    for width in (0.0, 100.0, 300.0):
        assert np.allclose(t9.kernel_weights(T, T[6:9], 80.0, width) @ line, 3.0 + 0.002 * (T[6:9] - 80.0 - T[0]))


def test_usable_visits_and_grid_respect_the_no_extrapolation_rule():
    use = t9.usable_visits(T)
    assert np.all(T[use] - t9.SUPPORT_PAST_DAYS >= T[0]) and np.all(T[use] + t9.SUPPORT_FUTURE_DAYS <= T[-1])
    assert not use[0] and not use[1] and not use[-1] and use[5]
    gappy = np.delete(T, 10)                               # one missing visit: a gap of about a year
    assert np.diff(gappy).max() > t9.MAX_GAP_DAYS
    use_g = t9.usable_visits(gappy)
    k = int(np.argmax(np.diff(gappy)))
    for i in np.nonzero(use_g)[0]:
        assert gappy[i] + t9.SUPPORT_FUTURE_DAYS <= gappy[k] or gappy[i] - t9.SUPPORT_PAST_DAYS >= gappy[k + 1]
    assert use_g.sum() < use.sum()
    pts = t9.grid()
    assert (0.0, 0.0) in pts and len(pts) == 284
    for tau, width in pts:
        assert tau - width / 2 >= -t9.SUPPORT_FUTURE_DAYS - 1e-9 and tau + width / 2 <= t9.SUPPORT_PAST_DAYS + 1e-9
    prep = t9.prepare(T)                                   # every weight matrix can be built: no extrapolation anywhere
    assert prep["c"].shape == (284, int(use.sum()), T.size) and np.allclose(prep["c"].sum(axis=2), 1.0)


def test_residual_covariance_against_monte_carlo():
    rng = np.random.default_rng(2)
    idx = np.nonzero(t9.usable_visits(T))[0]
    c = t9.kernel_weights(T, T[idx], 120.0, 100.0)
    s1 = rng.uniform(0.05, 0.15, T.size)
    s2 = rng.uniform(0.05, 0.15, T.size)
    r, b = 0.4, 1.3
    cov = t9.residual_covariance(c, idx, s1, s2, r, b)
    n = 200_000
    common = rng.normal(size=(n, T.size))
    e1 = s1 * (np.sqrt(r) * common + np.sqrt(1 - r) * rng.normal(size=(n, T.size)))
    e2 = s2 * (np.sqrt(r) * common + np.sqrt(1 - r) * rng.normal(size=(n, T.size)))
    resid = e2[:, idx] - b * e1 @ c.T
    emp = np.cov(resid.T)
    assert np.allclose(emp, cov, atol=4e-4)
    # with no lag the covariance is the effective variance of a straight-line fit with errors in both bands
    c0 = t9.kernel_weights(T, T[idx], 0.0, 0.0)
    cov0 = t9.residual_covariance(c0, idx, s1, s2, r, b)
    assert np.allclose(cov0, np.diag(s2[idx] ** 2 + b ** 2 * s1[idx] ** 2 - 2 * b * r * s1[idx] * s2[idx]))
    diag, q, xs = t9.covariance_parts(c, idx, s1, s2, r)
    assert np.allclose(np.diag(diag) + b ** 2 * q - b * xs, cov)


def _piecewise_object(rng, tau, noise=0.0):
    f1 = 10.0 + np.cumsum(rng.normal(size=T.size))
    s1 = np.full(T.size, 0.05)
    s2 = np.full(T.size, 0.05)
    f2 = np.full(T.size, np.nan)
    use = t9.usable_visits(T)
    f2[use] = 2.0 + 1.5 * (t9.kernel_weights(T, T[use], tau, 0.0) @ f1)
    f2[~use] = 2.0 + 1.5 * f1[~use]
    if noise:
        f1 = f1 + rng.normal(0, noise, T.size)
        f2 = f2 + rng.normal(0, noise, T.size)
    return f1, s1, f2, s2


def test_fit_recovers_a_shift_of_a_piecewise_linear_light_curve():
    rng = np.random.default_rng(3)
    f1, s1, f2, s2 = _piecewise_object(rng, 60.0)
    fit = t9.fit_object(T, f1, s1, f2, s2, 0.0)
    assert fit["usable"] and fit["tau_t1"] == 60.0 and fit["chi2_t1"] == pytest.approx(0.0, abs=1e-6)
    assert fit["a_t1"] == pytest.approx(2.0, abs=1e-6) and fit["b_t1"] == pytest.approx(1.5, abs=1e-6)
    assert fit["ds_01"] > 100 and fit["ds_12"] == pytest.approx(0.0, abs=1e-6)
    assert abs(fit["eta_t1"]) < 1e-6 and fit["ds_eta_t1"] == pytest.approx(0.0, abs=1e-6)
    assert fit["tau_lo"] <= 60.0 <= fit["tau_hi"] and not fit["tau_at_edge"] and fit["well_determined"]
    off = t9.offsets(f1, s1, f2, s2, 0.0, fit)
    if off["usable_obs"]:
        assert off["d_colour_t1"] == pytest.approx(off["d_colour_obs"], abs=1e-9)
        assert off["d_colour_t0"] == pytest.approx(0.0, abs=1e-9)     # no lag: W2 is a line in W1, no loop
    # no lag in the truth: T0 is as good as T1
    f1, s1, f2, s2 = _piecewise_object(rng, 0.0)
    fit0 = t9.fit_object(T, f1, s1, f2, s2, 0.0)
    assert fit0["tau_t1"] == 0.0 and fit0["ds_01"] == pytest.approx(0.0, abs=1e-6) and not fit0["well_determined"]
    # too few visits
    assert not t9.fit_object(T[:8], f1[:8], s1[:8], f2[:8], s2[:8])["usable"]


def test_nested_statistics_and_the_branch_offset():
    rng = np.random.default_rng(4)
    f1, s1, f2, s2 = _piecewise_object(rng, 40.0, noise=0.05)
    rising = hyst8.branch_rising(f1, f2)
    fit = t9.fit_object(T, f1, s1, f2, s2, 0.3, rising=rising)
    assert fit["s_t2"] <= fit["s_t1"] <= fit["s_t0"]
    assert fit["ds_eta_t0"] >= -1e-9 and fit["ds_eta_t1"] >= -1e-9 and fit["ds_eta_t2"] >= -1e-9
    # an offset added to W2 on the rising visits is found by eta with the right sign and size
    idx = fit["idx"]
    g2 = f2.copy()
    g2[rising] += 0.8
    fit_b = t9.fit_object(T, f1, s1, g2, s2, 0.3, rising=rising)
    assert fit_b["eta_t1"] == pytest.approx(0.8, abs=0.25) and fit_b["ds_eta_t1"] > 20
    assert fit_b["eta_colour_t1"] > 0                       # W2 higher on the rising branch: redder there
    assert idx.size == fit["n_used"] >= t9.MIN_VISITS_FIT


def test_direct_offset_with_given_visits_matches_the_default():
    rng = np.random.default_rng(5)
    f2 = 10 + np.sin(np.linspace(0, 6, T.size)) * 3
    f1 = 5 + 0.5 * f2 + rng.normal(0, 0.05, T.size)
    s = np.full(T.size, 0.05)
    d = hyst8.direct_offset(f1, s, f2, s, r=0.2)
    assert d["usable"]
    inside = (f2 >= d["lo"]) & (f2 <= d["hi"])
    e = hyst8.direct_offset(f1, s, f2, s, r=0.2, inside=inside)
    assert e["delta"] == pytest.approx(d["delta"]) and e["d_colour"] == pytest.approx(d["d_colour"])


def test_truth_models_have_the_stated_properties():
    rng = np.random.default_rng(6)
    dt = t9.GRID_STEP_DAYS
    n = t9.driver_length(T)
    mags = t9.drw_magnitudes(rng, n, dt, 400.0)
    assert abs(np.std(t9.drw_magnitudes(rng, 200_000, dt, 200.0)) - 1.0) < 0.05
    ell = 10.0 ** (-0.4 * 0.3 * mags)
    w1, w2 = t9.truth_fluxes("linear", ell, dt, 100.0, 0.0, None)
    assert np.allclose(w1, w2)                                              # T0: one response for both bands
    w1, w2 = t9.truth_fluxes("linear", ell, dt, 100.0, 40.0, None)
    k = int(40.0 / dt)
    d1 = (w1 - (1 - t9.DUST_FRACTION)) / t9.DUST_FRACTION
    d2 = (w2 - (1 - t9.DUST_FRACTION)) / t9.DUST_FRACTION
    ratio = d2[k + 50:] / d1[50:-k]                                        # W2 response = W1 response 40 d earlier
    assert np.allclose(ratio, ratio[0], rtol=1e-9)
    # T2: the mean extra delay of the tail model equals the stated delay
    step = np.zeros(4000)
    step[100:] = 1.0
    r1 = t9.boxcar_response(step, dt, 100.0)
    r2 = (1 - t9.TAIL_SHARE) * r1 + t9.TAIL_SHARE * t9.exponential_tail(r1, dt, 80.0 / t9.TAIL_SHARE)
    delay = np.sum(r1 - r2) * dt                                            # area between the two step responses
    assert delay == pytest.approx(80.0, rel=0.05)
    # amplitude: the noise-free W2 at the visits has the requested fractional rms
    f1o = 10 + rng.normal(size=T.size)
    f2o = 12 + 1.5 * rng.normal(size=T.size)
    s1 = np.full(T.size, 0.1)
    s2 = np.full(T.size, 0.12)
    target = np.sqrt(np.var(f2o) - np.mean(s2 ** 2)) / np.mean(f2o)
    for name, (kind, delay_days, tau_re) in t9.TRUTHS.items():
        v1, v2 = t9.mock_truth(mags, T, f1o, s1, f2o, s2, kind, 100.0, 40.0 if delay_days == "lyu" else delay_days, tau_re)
        assert np.std(v2) / np.mean(v2) == pytest.approx(target, rel=0.05), name
        assert np.median(v1) == pytest.approx(np.median(f1o)) and np.median(v2) == pytest.approx(np.median(f2o))
    v1, v2 = t9.mock_truth(mags, T, f1o, s1, f2o, s2, "linear", 100.0, 0.0)
    assert np.allclose(np.polyval(np.polyfit(v1, v2, 1), v1), v2)           # T0 at the visits: W2 is a line in W1
    w1_lag, w2_lag = t9.lyu_lags_days(np.log10(3.828e33 * 1e11))
    assert w1_lag == pytest.approx(10 ** 2.10) and w2_lag == pytest.approx(10 ** 2.20)


def test_sample_statistics_counts_and_medians():
    rows = []
    for tau, d_obs, d_pred in ((10.0, -0.02, -0.01), (30.0, -0.04, -0.03), (-5.0, 0.01, 0.0)):
        row = {"usable": True, "tau_t1": tau, "ds_01": 4.0, "ds_02": 5.0, "ds_12": 1.0, "well_determined": tau > 20,
               "width_t2": 100.0, "usable_obs": True, "comparable_obs": True, "d_colour_obs": d_obs}
        for name in ("t0", "t1", "t2"):
            row.update({f"ds_eta_{name}": 2.0, f"eta_colour_{name}": -0.01, f"usable_{name}": True, f"d_colour_{name}": d_pred})
        rows.append(row)
    rows.append({"usable": False})
    st = t9.sample_statistics(rows)
    assert st["n_fit"] == 3 and st["n_offset"] == 3 and st["n_well"] == 1
    assert st["median_tau"] == 10.0 and st["frac_tau_positive"] == pytest.approx(2 / 3)
    assert st["median_d_colour_obs"] == -0.02 and st["median_d_colour_res_t1"] == pytest.approx(-0.01)
    assert st["frac_res_bluer_t2"] == pytest.approx(2 / 3) and st["ratio_pred_obs_t1"] == pytest.approx(0.5)
