"""Stage 5b slope and hysteresis estimators on toy models (no data, no network)."""

import numpy as np
import pytest

from irecho import breath, photometry
from irecho.breath import ToyModel

MJD = 56650.0 + 182.6 * np.arange(21)
TINY = {"w1": 1e-4, "w2": 1e-4}  # 0.1 mmag: effectively noise-free


def _measure(model, sigma=0.3, ramp=1.5, err=TINY, seed=5, **kw):
    rng = np.random.default_rng(seed)
    v = breath.synthetic_visits(rng, model, MJD, err, sigma, ramp_mag=ramp, **kw)
    return breath.measure(v["w1_flux"], v["w1_flux_err"], v["w2_flux"], v["w2_flux_err"]), v


def test_blackbody_pieces_agree_with_photometry_module():
    for t in (500.0, 1000.0, 1500.0):
        assert float(breath.rho(t)) == pytest.approx(photometry.blackbody_flux_ratio(t), rel=1e-12)
        for b in ("w1", "w2"):
            lo, hi = breath.planck_rel(b, t * 0.999), breath.planck_rel(b, t * 1.001)
            numeric = (np.log(hi) - np.log(lo)) / (np.log(1.001) - np.log(0.999))
            assert float(breath.gamma(b, t)) == pytest.approx(numeric, rel=1e-5)


def test_beta_to_s_limits():
    assert breath.beta_to_s(0.0, 1200.0) == 0.0
    # beta = 0.306 at 1000 K corresponds to s close to 0.25 (gamma_2 = 3.27)
    assert breath.beta_to_s(0.306, 1000.0) == pytest.approx(0.25, abs=0.005)


def test_boundary_radius_follows_then_relaxes():
    t = np.arange(0.0, 4000.0, 10.0)
    ell = np.where(t < 1000.0, 1.0, np.where(t < 2000.0, 4.0, 1.0))
    instant = breath.boundary_radius(t, ell, 0.0)
    never = breath.boundary_radius(t, ell, float("inf"))
    slow = breath.boundary_radius(t, ell, 1.0)
    assert np.allclose(instant, np.sqrt(ell))
    assert np.allclose(never, np.maximum.accumulate(np.sqrt(ell)))
    assert slow[150] == pytest.approx(2.0)  # follows the rise at once
    k = int(np.argmin(np.abs(t - (2000.0 + 365.25))))  # one e-folding after the drop
    assert slow[k] == pytest.approx(1.0 + np.exp(-1.0), abs=0.02)


def test_boundary_model_has_constant_colour_and_flux_proportional_to_luminosity():
    t = np.arange(0.0, 3000.0, 10.0)
    ell = 1.0 + 0.5 * np.sin(t / 300.0)
    f = breath.toy_fluxes(ToyModel("boundary", 0.0, lag_days=0.0), t, ell)
    assert np.allclose(f["w1"] / f["w2"], breath.rho(1500.0))
    assert np.allclose(f["w2"], ell)


def test_fixed_area_model_follows_t_proportional_to_l_quarter():
    t = np.arange(0.0, 3000.0, 10.0)
    ell = np.linspace(0.5, 2.0, t.size)
    f = breath.toy_fluxes(ToyModel("fixed", t_fixed=1000.0, lag_days=0.0), t, ell)
    assert np.allclose(f["w1"] / f["w2"], breath.rho(1000.0 * ell**0.25))


def test_host_fraction_is_what_was_asked_for():
    t = np.arange(0.0, 100.0, 10.0)
    f = breath.toy_fluxes(ToyModel("boundary", 0.0, host_w1=0.6, lag_days=0.0), t, np.ones(t.size))
    dust_w1 = float(breath.rho(1500.0))
    assert (f["w1"][0] - dust_w1) / f["w1"][0] == pytest.approx(0.6)
    assert (f["w1"][0] - dust_w1) / (f["w2"][0] - 1.0) == pytest.approx(breath.HOST_F1_OVER_F2)


@pytest.mark.parametrize("host", [0.0, 0.5, 0.8])
def test_boundary_shift_gives_zero_slope_whatever_the_host(host):
    m, _ = _measure(ToyModel("boundary", 0.0, host_w1=host))
    assert m["status"] == "ok" and abs(m["beta"]) < 0.01
    assert m["t_pivot_k"] == pytest.approx(1500.0, rel=0.02)
    assert m["class"] == "boundary" and m["measurable"]


def test_heating_in_place_gives_r_of_one_without_host_and_more_with_host():
    m0, _ = _measure(ToyModel("fixed"))
    m1, _ = _measure(ToyModel("fixed", host_w1=0.5))
    assert m0["beta"] > 0.05 and m0["class"] == "in place"
    assert m0["r"] == pytest.approx(1.0, abs=0.15)
    assert m1["r"] > m0["r"] + 0.2


def test_heating_in_place_excess_matches_direct_planck():
    a2, t_ref, t_hot = 2.0, 900.0, 1100.0
    omega = a2 / breath.planck_rel("w2", t_ref)
    e2 = omega * (breath.planck_rel("w2", t_hot) - breath.planck_rel("w2", t_ref))
    e1 = omega * (breath.planck_rel("w1", t_hot) - breath.planck_rel("w1", t_ref))
    assert float(breath.heating_in_place_excess(e2, a2, t_ref)) == pytest.approx(float(e1))


def test_slope_does_not_depend_on_flux_units():
    _, v = _measure(ToyModel("fixed", host_w1=0.3), err={"w1": 0.02, "w2": 0.015})
    a = breath.measure(v["w1_flux"], v["w1_flux_err"], v["w2_flux"], v["w2_flux_err"])
    b = breath.measure(7.0 * v["w1_flux"], 7.0 * v["w1_flux_err"],
                       7.0 * v["w2_flux"], 7.0 * v["w2_flux_err"])
    assert b["beta"] == pytest.approx(a["beta"], abs=1e-4)
    assert b["sigma_beta"] == pytest.approx(a["sigma_beta"], rel=1e-2)


def test_noise_makes_small_changes_unmeasurable():
    m, _ = _measure(ToyModel("boundary", 0.0), sigma=0.05, ramp=0.0,
                    err={"w1": 0.02, "w2": 0.015})
    assert not m["measurable"]


def test_too_few_visits():
    m = breath.measure(np.ones(5), np.ones(5) * 0.01, np.ones(5), np.ones(5) * 0.01)
    assert m["status"] == "too few visits" and not m["measurable"]
    flat = breath.measure(np.ones(20), np.full(20, 0.01), np.ones(20), np.full(20, 0.01))
    assert flat["status"] == "too few used visits"


def test_classify_rules():
    assert breath.classify(0.00, 0.02, 0.10) == "boundary"
    assert breath.classify(0.09, 0.02, 0.10) == "in place"
    assert breath.classify(-0.06, 0.02, 0.10) == "negative"
    assert breath.classify(0.03, 0.02, 0.05) == "unclear"
    assert breath.classify(float("nan"), 0.02, 0.05) == ""


def test_slow_reformation_gives_positive_hysteresis_and_instant_gives_none():
    d_slow, d_fast = [], []
    for seed in range(12):
        m, _ = _measure(ToyModel("boundary", 2.0), sigma=0.6, ramp=0.0, seed=seed)
        if np.isfinite(m["d_hyst"]):
            d_slow.append(m["d_hyst"])
        m, _ = _measure(ToyModel("boundary", 0.0), sigma=0.6, ramp=0.0, seed=seed)
        if np.isfinite(m["d_hyst"]):
            d_fast.append(m["d_hyst"])
    assert len(d_slow) >= 8 and np.median(d_slow) > 0.05
    assert np.max(np.abs(d_fast)) < 0.02


def test_band_lag_mimics_hysteresis():
    d = []
    for seed in range(12):
        m, _ = _measure(ToyModel("boundary", 0.0, lag_w2_extra_days=90.0), sigma=0.6, ramp=0.0,
                        seed=seed)
        if np.isfinite(m["d_hyst"]):
            d.append(m["d_hyst"])
    assert np.median(d) > 0.1  # a confound: same sign as slow reformation


def test_synthetic_visits_are_reproducible():
    a = breath.synthetic_visits(np.random.default_rng(3), ToyModel(), MJD, TINY, 0.3)
    b = breath.synthetic_visits(np.random.default_rng(3), ToyModel(), MJD, TINY, 0.3)
    assert np.array_equal(a["w1_flux"], b["w1_flux"])
