"""Stage 6B (B4): the derivative-regression lag estimator on exact inputs."""

import numpy as np
import pytest

from irecho import bandlag


def smooth_object(lag_days, seed, n=21, noise=0.0):
    rng = np.random.default_rng(seed)
    t = 56650.0 + 182.6 * np.arange(n)
    phase, period = rng.uniform(0, 2 * np.pi), rng.uniform(2500.0, 4000.0)
    f = lambda tt: 1.0 + 0.3 * np.sin(2 * np.pi * tt / period + phase)  # noqa: E731
    f1 = f(t)
    f2 = 0.5 + 1.7 * f(t - lag_days) + rng.normal(0.0, noise, n)
    return bandlag.object_terms(t, f1, f2, np.full(n, max(noise, 1e-3)))


def test_derivative_of_a_line():
    t = np.arange(0.0, 1000.0, 100.0)
    g = bandlag.derivative(t, 3.0 + 0.02 * t)
    assert np.isnan(g[0]) and np.isnan(g[-1]) and np.allclose(g[1:-1], 0.02)


@pytest.mark.parametrize("lag", [0.0, 40.0, 120.0])
def test_recovers_a_lag_from_slow_noise_free_light_curves(lag):
    objs = [smooth_object(lag, s) for s in range(12)]
    res = bandlag.stacked_lag(objs)
    assert res["lag_days"] == pytest.approx(lag, abs=12.0) and not res["at_grid_edge"]


def test_scaled_lag_and_noise():
    scale = np.linspace(0.5, 2.0, 30)
    objs = [smooth_object(50.0 * s, k, noise=0.003) for k, s in enumerate(scale)]
    res = bandlag.stacked_lag(objs, scale=scale)
    assert res["lag_days"] == pytest.approx(50.0, abs=3.0 * res["lag_err_days"] + 6.0)
    assert res["lag_err_days"] > 0 and res["n_objects"] == 30


def test_too_few_visits_gives_none():
    t = 56650.0 + 182.6 * np.arange(6)
    assert bandlag.object_terms(t, np.ones(6), np.ones(6), np.full(6, 0.01)) is None
