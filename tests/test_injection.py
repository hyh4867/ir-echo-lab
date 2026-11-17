import numpy as np
import pytest

from irecho import injection as inj
from irecho import mlcompare
from irecho.config import InjectionGrid

MJD = 56700.0 + 182.625 * np.arange(21)


def host(seed, name="h", kind="star", sigma=0.008):
    rng = np.random.default_rng(seed)
    return {"name": name, "kind": kind, "visits": {
        "mjd": MJD.copy(),
        "w1_mag": 10.0 + rng.normal(0, sigma, 21), "w1_err": np.full(21, sigma),
        "w2_mag": 9.6 + rng.normal(0, sigma, 21), "w2_err": np.full(21, sigma)}}


def test_segment_onsets():
    assert inj.segment_onsets(21, 4, 7) == list(range(4, 15))
    assert inj.segment_onsets(21, 12, 7) == [12, 13, 14]
    assert inj.segment_onsets(10, 4, 7) == []


def test_cut_segment_layout():
    seg = inj.cut_segment(host(1)["visits"], onset=9, baseline_visits=4, post_visits=7)
    assert seg["mjd"].size == 11
    assert seg["mjd"][0] == MJD[5] and seg["mjd"][4] == MJD[9]


def test_inject_step_and_ramp():
    seg = inj.cut_segment(host(1)["visits"], 9, 4, 7)
    step = inj.inject(seg, 4, 0.3, 0.0)
    assert step["w1_mag"][:4] == pytest.approx(seg["w1_mag"][:4])
    assert seg["w1_mag"][4:] - step["w1_mag"][4:] == pytest.approx(0.3)
    assert seg["w2_mag"][4:] - step["w2_mag"][4:] == pytest.approx(0.3)
    slow = inj.inject(seg, 4, 0.3, 3.0)
    added = seg["w1_mag"] - slow["w1_mag"]
    assert added[:4] == pytest.approx(0.0)
    assert added[4] == pytest.approx(0.3 * 0.5 / 3.0, rel=1e-3)  # half a year into a 3-yr rise
    assert added[-1] == pytest.approx(0.3)  # 3.5 yr after the last baseline visit
    assert np.all(np.diff(added) >= 0)
    assert seg["w1_mag"] is not slow["w1_mag"]


def test_detector_a_only_counts_flags_after_the_onset():
    seg = inj.cut_segment(host(2)["visits"], 12, 8, 7)
    assert inj.run_detectors(seg, None) == {"a": False, "b": False, "q1": False}
    assert inj.run_detectors(inj.inject(seg, 8, 0.5, 0.0), 8)["a"]
    early = {k: v.copy() for k, v in seg.items()}
    for b in ("w1", "w2"):
        early[f"{b}_mag"][5] -= 0.6  # a spike before the nominal onset
    assert inj.run_detectors(early, None)["a"]  # counted as a false alarm
    assert not inj.run_detectors(early, 8)["a"]  # not counted as a recovery


def test_slow_ramp_needs_detector_b():
    seg = inj.cut_segment(host(3)["visits"], 12, 8, 7)
    hit = inj.run_detectors(inj.inject(seg, 8, 0.2, 3.0), 8)
    assert hit["b"] and not hit["a"]


def test_q1_class_needs_a_fitted_rise_of_0p2_mag_in_both_bands():
    seg = inj.cut_segment(host(3)["visits"], 12, 8, 7)
    small = inj.run_detectors(inj.inject(seg, 8, 0.12, 3.0), 8)
    large = inj.run_detectors(inj.inject(seg, 8, 0.5, 3.0), 8)
    assert small["b"] and not small["q1"]
    assert large["b"] and large["q1"]


def test_features_shape_and_response():
    seg = inj.cut_segment(host(4)["visits"], 9, 4, 7)
    f0 = inj.features(seg)
    f1 = inj.features(inj.inject(seg, 4, 0.5, 0.0))
    assert f0.shape == (len(inj.FEATURE_NAMES),) == (18,)
    assert np.isfinite(f0).all() and np.isfinite(f1).all()
    names = inj.FEATURE_NAMES
    assert f1[names.index("w1_max_brightening_mag")] == pytest.approx(0.5, abs=0.05)
    assert f1[names.index("w2_cum_change_mag")] == pytest.approx(0.5, abs=0.05)
    assert f1[names.index("color_change_mag")] == pytest.approx(
        f0[names.index("color_change_mag")], abs=1e-9)  # same amplitude in both bands


def test_grid_row_count_and_ml_columns():
    grid = InjectionGrid(amplitudes_mag=(0.05, 0.5), rise_years=(0.0, 3.0),
                         baseline_visits=(4, 12), post_visits=7)
    hosts = [host(10 + i, name=f"h{i}") for i in range(10)]
    rows = inj.run_grid(hosts, grid)
    assert len(rows) == 10 * (11 + 3) * (1 + 2 * 2)
    assert sum(r["amp"] == 0.0 for r in rows) == 10 * 14
    big = [r for r in rows if r["amp"] == 0.5 and r["rise"] == 0.0]
    assert all(r["a"] for r in big)
    assert not any(r["a"] or r["b"] for r in rows if r["amp"] == 0.0)

    rows = mlcompare.evaluate(rows)
    again = mlcompare.evaluate([dict(r) for r in rows])
    assert [r["ml_matched"] for r in rows] == [r["ml_matched"] for r in again]  # seeded
    null = [r for r in rows if r["amp"] == 0.0]
    assert np.mean([r["ml_far05"] for r in null]) < 0.25
    assert np.mean([r["ml_far05"] for r in big]) > 0.9
