"""Stage 5a error floor (synthetic data only; no target light curve is used)."""

import numpy as np
import pytest

from irecho import errorfloor as ef
from irecho.stats import constant_chi2


def _visits(n, rng, w1_err=0.01, w2_err=0.012):
    w1 = 10.0 + rng.normal(0, w1_err, n)
    w2 = 9.5 + rng.normal(0, w2_err, n)
    return {"mjd": 57000.0 + 180.0 * np.arange(n),
            "w1_mag": w1, "w1_err": np.full(n, w1_err),
            "w2_mag": w2, "w2_err": np.full(n, w2_err),
            "w1_flux_jy": 309.54 * 10 ** (-0.4 * w1), "w2_flux_jy": 171.787 * 10 ** (-0.4 * w2),
            "w1_flux_err_jy": np.zeros(n), "w2_flux_err_jy": np.zeros(n),
            "color": w1 - w2, "color_err": np.zeros(n)}


def test_pair_scatter_ignores_source_variability():
    # A strongly variable source measured twice with independent 10 and 15 mmag errors:
    # the difference only knows about the measurement errors.
    rng = np.random.default_rng(1)
    n = 4000
    truth = 10.0 + 0.5 * np.sin(np.arange(n) / 30.0) + np.linspace(0, -1.0, n)
    a = truth + rng.normal(0, 0.010, n)
    b = truth + 0.3 + rng.normal(0, 0.015, n)
    ps = ef.pair_scatter(a, b)
    assert ps.n_pairs == n and ps.usable
    assert ps.scatter_mag == pytest.approx(np.hypot(0.010, 0.015), rel=0.05)
    assert ps.median_diff_mag == pytest.approx(-0.3, abs=0.002)


def test_pair_scatter_is_robust_and_handles_missing():
    rng = np.random.default_rng(2)
    a = 10.0 + rng.normal(0, 0.01, 400)
    b = 10.0 + rng.normal(0, 0.01, 400)
    clean = ef.pair_scatter(a, b).scatter_mag
    b_bad = b.copy()
    b_bad[:8] += 1.0  # 2 % of the coadds are wrong by a magnitude
    assert ef.pair_scatter(a, b_bad).scatter_mag == pytest.approx(clean, rel=0.1)
    b_nan = b.copy()
    b_nan[::2] = np.nan
    assert ef.pair_scatter(a, b_nan).n_pairs == 200
    few = ef.pair_scatter(a[:5], b[:5])
    assert few.n_pairs == 5 and not few.usable and np.isfinite(few.scatter_mag)
    one = ef.pair_scatter(a[:1], b[:1])
    assert np.isnan(one.scatter_mag) and not one.usable
    none = ef.pair_scatter([], [])
    assert none.n_pairs == 0 and np.isnan(none.median_diff_mag)


def test_reference_scatter_bins_and_lookup():
    scat = [0.010, 0.012, 0.014, 0.020, 0.022, np.nan]
    mags = [9.6, 10.0, 10.4, 10.6, 11.4, 11.0]
    ref = ef.reference_scatter(scat, mags)
    assert ref["n_stars"] == [3, 2] and ref["scatter_mag"] == pytest.approx([0.012, 0.021])
    assert ef.reference_for(8.9, ref) == pytest.approx(0.012)  # brighter than every star
    assert ef.reference_for(10.5, ref) == pytest.approx(0.021)
    assert ef.reference_for(12.4, ref) == pytest.approx(0.021)  # fainter than every star
    # an empty bin falls back to the median of all stars
    ref2 = ef.reference_scatter([0.01, 0.02, 0.03], [9.6, 9.8, 10.0])
    assert ref2["n_stars"] == [3, 0] and ref2["scatter_mag"] == pytest.approx([0.02, 0.02])


def test_floor_from_scatter():
    assert ef.floor_from_scatter(0.05, 0.03) == pytest.approx(0.04)
    assert ef.floor_from_scatter(0.02, 0.03) == 0.0  # never negative
    assert np.isnan(ef.floor_from_scatter(np.nan, 0.03))
    assert np.isnan(ef.floor_from_scatter(0.05, np.nan))


def test_recovers_an_injected_floor():
    # Pipeline: 10 mmag statistical + 30 mmag hidden visit-level error.
    # unTimely: 12 mmag. Stars: the same without the hidden term.
    rng = np.random.default_rng(3)
    n = 3000
    truth = 10.0 + rng.normal(0, 0.2, n)
    star = ef.pair_scatter(truth + rng.normal(0, 0.010, n), truth + rng.normal(0, 0.012, n))
    gal = ef.pair_scatter(truth + rng.normal(0, np.hypot(0.010, 0.030), n),
                          truth + rng.normal(0, 0.012, n))
    floor = ef.floor_from_scatter(gal.scatter_mag, star.scatter_mag)
    assert floor == pytest.approx(0.030, rel=0.08)


def test_adopt_floors_uses_own_value_or_sample_median():
    per_object = {"a": {"w1": 0.02, "w2": 0.01}, "b": {"w1": 0.04, "w2": 0.03},
                  "c": {"w1": 0.06, "w2": np.nan}, "d": {"w1": 0.50, "w2": 0.50}}
    usable = {"a": {"w1": True, "w2": True}, "b": {"w1": True, "w2": True},
              "c": {"w1": True, "w2": True}, "d": {"w1": False, "w2": False}}
    floors, fallback = ef.adopt_floors(per_object, usable)
    assert fallback == {"w1": pytest.approx(0.04), "w2": pytest.approx(0.02)}
    assert floors["a"]["w1"] == 0.02 and floors["a"]["w1_source"] == "own"
    assert floors["c"]["w2"] == pytest.approx(0.02) and floors["c"]["w2_source"] == "sample median"
    # too few pairs: the object's own (unreliable) 0.5 mag is not used
    assert floors["d"]["w1"] == pytest.approx(0.04) and floors["d"]["w1_source"] == "sample median"


def test_apply_floor_adds_in_quadrature_and_keeps_magnitudes():
    rng = np.random.default_rng(4)
    v = _visits(20, rng)
    out = ef.apply_floor(v, {"w1": 0.03, "w2": 0.04})
    assert np.array_equal(out["w1_mag"], v["w1_mag"]) and np.array_equal(out["mjd"], v["mjd"])
    assert out["w1_err"] == pytest.approx(np.hypot(0.01, 0.03))
    assert out["w2_err"] == pytest.approx(np.hypot(0.012, 0.04))
    assert out["color_err"] == pytest.approx(np.hypot(out["w1_err"], out["w2_err"]))
    assert out["w1_flux_err_jy"] == pytest.approx(
        v["w1_flux_jy"] * np.log(10) / 2.5 * out["w1_err"])
    assert out["w1_floor"] == pytest.approx(0.03)
    assert v["w1_err"] == pytest.approx(0.01)  # input untouched
    assert ef.apply_floor(v, {"w1": 0.0, "w2": 0.0})["w1_err"] == pytest.approx(v["w1_err"])


def test_apply_floor_refuses_missing_floor():
    v = _visits(5, np.random.default_rng(5))
    with pytest.raises(ValueError, match="invalid error floor"):
        ef.apply_floor(v, {"w1": float("nan"), "w2": 0.01})
    with pytest.raises(ValueError, match="invalid error floor"):
        ef.apply_floor(v, {"w1": -0.01, "w2": 0.01})


def test_floor_restores_chi2_of_a_constant_source():
    # Constant source whose true visit error is 10 mmag (quoted) + 30 mmag (hidden).
    rng = np.random.default_rng(6)
    n = 2000
    v = _visits(n, rng, w1_err=0.010, w2_err=0.010)
    v["w1_mag"] = v["w1_mag"] + rng.normal(0, 0.030, n)
    before = constant_chi2(v["w1_mag"], v["w1_err"]).chi2_red
    after = constant_chi2(*[ef.apply_floor(v, {"w1": 0.030, "w2": 0.0})[k]
                            for k in ("w1_mag", "w1_err")]).chi2_red
    assert before == pytest.approx(10.0, rel=0.1)
    assert after == pytest.approx(1.0, rel=0.08)
