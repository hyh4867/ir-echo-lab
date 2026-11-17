import numpy as np
import pytest

from irecho.config import FlareCriteria
from irecho.flare import detect_flares


def quiet_curve(n=20, sigma=0.02, seed=0):
    rng = np.random.default_rng(seed)
    mjd = 56700.0 + 182.0 * np.arange(n)
    w1 = 10.3 + rng.normal(0, sigma, n)
    w2 = 10.0 + rng.normal(0, sigma, n)
    err = np.full(n, sigma)
    return mjd, w1, err.copy(), w2, err.copy()


def test_recovers_injected_flare_onset_and_amplitude():
    mjd, w1, e1, w2, e2 = quiet_curve(seed=1)
    w1[13:] -= np.array([0.5, 0.45, 0.4, 0.3, 0.25, 0.2, 0.15])
    w2[13:] -= np.array([0.8, 0.75, 0.6, 0.5, 0.4, 0.3, 0.25])
    res = detect_flares(mjd, w1, e1, w2, e2)
    assert res.detected
    assert res.onset_index == 13
    assert res.onset_mjd == pytest.approx(mjd[13])
    amp1, t1 = res.peak("w1")
    amp2, t2 = res.peak("w2")
    assert amp1 == pytest.approx(0.5, abs=0.05) and t1 == pytest.approx(mjd[13])
    assert amp2 == pytest.approx(0.8, abs=0.05) and t2 == pytest.approx(mjd[13])
    assert not res.flagged[:13].any()


def test_flagged_visits_do_not_contaminate_the_baseline():
    mjd, w1, e1, w2, e2 = quiet_curve(n=24, seed=2)
    w1[10:] -= 0.6
    w2[10:] -= 0.6  # long plateau: later visits must still be measured against quiescence
    res = detect_flares(mjd, w1, e1, w2, e2)
    assert res.flagged[10:].all()
    assert res.dmag["w1"][-1] == pytest.approx(0.6, abs=0.05)
    assert res.baseline_mag["w1"][-1] == pytest.approx(10.3, abs=0.02)


def test_no_false_positives_on_gaussian_noise():
    flagged = sum(detect_flares(*quiet_curve(seed=s)).flagged.sum() for s in range(300))
    assert flagged == 0


def test_single_band_brightening_is_not_a_flare():
    mjd, w1, e1, w2, e2 = quiet_curve(seed=3)
    w1[12] -= 0.8  # e.g. a W1-only artefact
    res = detect_flares(mjd, w1, e1, w2, e2)
    assert not res.detected
    assert res.significance["w1"][12] > 5 and res.significance["w2"][12] < 5
    assert np.isnan(res.onset_mjd) and np.isnan(res.peak("w1")[0])


def test_significant_but_small_brightening_fails_amplitude_floor():
    mjd, w1, e1, w2, e2 = quiet_curve(sigma=0.005, seed=4)
    w1[15] -= 0.1
    w2[15] -= 0.1  # 14 sigma but below the 0.2 mag floor
    res = detect_flares(mjd, w1, e1, w2, e2)
    assert res.significance["w1"][15] > 5 and not res.detected
    loose = detect_flares(mjd, w1, e1, w2, e2, FlareCriteria(min_amp_mag=0.05))
    assert loose.flagged[15]


def test_dimming_is_ignored():
    mjd, w1, e1, w2, e2 = quiet_curve(seed=5)
    w1[10:] += 0.7
    w2[10:] += 0.7
    assert not detect_flares(mjd, w1, e1, w2, e2).detected


def test_needs_minimum_baseline():
    mjd, w1, e1, w2, e2 = quiet_curve(seed=6)
    w1[2] -= 1.0
    w2[2] -= 1.0  # only two earlier visits: not evaluated
    res = detect_flares(mjd, w1, e1, w2, e2)
    assert np.isnan(res.significance["w1"][2]) and not res.flagged[2]


def test_nan_visits_are_skipped_everywhere():
    mjd, w1, e1, w2, e2 = quiet_curve(seed=7)
    w2[5] = np.nan
    w1[14] -= 0.5
    w2[14] -= 0.8
    res = detect_flares(mjd, w1, e1, w2, e2)
    assert res.onset_index == 14
    assert np.isnan(res.dmag["w1"][5]) and not res.flagged[5]


def test_baseline_scatter_floor_is_the_visit_error():
    # identical baseline magnitudes would give zero robust scatter
    mjd = 56700.0 + 182.0 * np.arange(8)
    w1 = np.full(8, 10.3)
    w2 = np.full(8, 10.0)
    err = np.full(8, 0.03)
    res = detect_flares(mjd, w1, err, w2, err)
    assert res.baseline_scatter["w1"][-1] == pytest.approx(0.03)
    assert np.all(res.significance["w1"][4:] == 0)


def test_requires_time_order():
    mjd, w1, e1, w2, e2 = quiet_curve(seed=8)
    with pytest.raises(ValueError):
        detect_flares(mjd[::-1], w1, e1, w2, e2)
