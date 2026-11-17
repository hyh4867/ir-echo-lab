import numpy as np
import pytest

from irecho import stats


def test_robust_std_gaussian():
    rng = np.random.default_rng(3)
    assert stats.robust_std(rng.normal(0, 2.0, 20000)) == pytest.approx(2.0, rel=0.03)


def test_robust_std_ignores_outlier_and_nan():
    assert stats.robust_std([1.0, 2.0, 3.0, 4.0, 5.0, 1000.0, np.nan]) == pytest.approx(
        1.4826 * 1.5)
    assert np.isnan(stats.robust_std([1.0]))


def test_constant_chi2_hand_computed():
    # values 1, 2, 3 with unit errors: mean 2, chi2 = 1 + 0 + 1 = 2, dof = 2
    r = stats.constant_chi2([1.0, 2.0, 3.0], [1.0, 1.0, 1.0])
    assert (r.n, r.dof) == (3, 2)
    assert r.weighted_mean == pytest.approx(2.0)
    assert r.chi2 == pytest.approx(2.0)
    assert r.chi2_red == pytest.approx(1.0)
    assert r.p_value == pytest.approx(np.exp(-1.0))  # chi2 sf with 2 dof = exp(-x/2)
    assert r.rms == pytest.approx(1.0)


def test_constant_chi2_uses_inverse_variance_weights():
    r = stats.constant_chi2([0.0, 10.0], [1.0, 3.0])
    assert r.weighted_mean == pytest.approx(1.0)  # (0*1 + 10/9) / (1 + 1/9)
    assert r.chi2 == pytest.approx(1.0 + 81.0 / 9.0)


def test_constant_chi2_drops_invalid_points():
    r = stats.constant_chi2([1.0, np.nan, 3.0, 5.0], [1.0, 1.0, np.nan, 0.0])
    assert r.n == 1 and np.isnan(r.chi2_red)


def test_reduced_chi2_is_one_for_correct_errors():
    rng = np.random.default_rng(11)
    values = [stats.constant_chi2(rng.normal(10.0, 0.02, 21), np.full(21, 0.02)).chi2_red
              for _ in range(400)]
    assert np.mean(values) == pytest.approx(1.0, abs=0.05)


def test_error_floor_recovers_injected_systematic():
    rng = np.random.default_rng(5)
    quoted, extra = 0.01, 0.03
    values = rng.normal(10.0, np.hypot(quoted, extra), 4000)
    floor = stats.error_floor(values, np.full(values.size, quoted))
    assert floor == pytest.approx(extra, rel=0.05)


def test_error_floor_zero_when_errors_are_generous():
    rng = np.random.default_rng(6)
    values = rng.normal(10.0, 0.01, 200)
    assert stats.error_floor(values, np.full(200, 0.05)) == 0.0


def test_pooled_error_floor():
    rng = np.random.default_rng(7)
    sets_v = [rng.normal(m, np.hypot(0.01, 0.02), 21) for m in np.linspace(9.5, 12.5, 200)]
    sets_e = [np.full(21, 0.01)] * 200
    assert stats.pooled_error_floor(sets_v, sets_e) == pytest.approx(0.02, rel=0.08)
