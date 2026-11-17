import numpy as np
import pytest

from irecho import epochs
from irecho.config import EpochBinning
from irecho.photometry import mag_to_flux_jy


def test_assign_epochs_splits_on_gaps():
    mjd = np.array([0.0, 0.5, 1.0, 180.0, 180.7, 365.0])
    assert epochs.assign_epochs(mjd, gap_days=30).tolist() == [0, 0, 0, 1, 1, 2]


def test_assign_epochs_keeps_input_order_when_unsorted():
    mjd = np.array([180.0, 0.0, 365.0, 0.5, 180.7])
    assert epochs.assign_epochs(mjd, gap_days=30).tolist() == [1, 0, 2, 0, 1]


def test_assign_epochs_gap_threshold_is_exclusive():
    assert epochs.assign_epochs([0.0, 30.0], gap_days=30).tolist() == [0, 0]
    assert epochs.assign_epochs([0.0, 30.1], gap_days=30).tolist() == [0, 1]
    assert epochs.assign_epochs([]).size == 0


def test_summarise_visit_median_and_formal_error():
    mag = np.array([10.00, 10.02, 9.98, 10.01, 9.99, 10.00, 10.03, 9.97, 10.00])
    err = np.full(mag.size, 0.02)
    s = epochs.summarise_visit(np.arange(mag.size, dtype=float), mag, err, "w1")
    assert s["n"] == 9 and s["n_clipped"] == 0
    assert s["mag"] == pytest.approx(10.00, abs=1e-9)
    assert s["flux"] == pytest.approx(mag_to_flux_jy(10.0, "w1"))
    assert s["mjd"] == pytest.approx(4.0)
    # standard error of the median: sqrt(pi/2) * 0.02 / sqrt(9)
    assert s["err_formal"] == pytest.approx(np.sqrt(np.pi / 2) * 0.02 / 3.0, rel=1e-3)
    assert s["err"] == pytest.approx(max(s["err_formal"], s["err_scatter"]))
    assert s["scatter"] == pytest.approx(np.std(mag, ddof=1))


def test_summarise_visit_adopts_scatter_when_larger_than_formal():
    rng = np.random.default_rng(4)
    mag = 10.0 + rng.normal(0, 0.08, 14)  # real scatter four times the quoted error
    s = epochs.summarise_visit(np.zeros(14), mag, np.full(14, 0.02), "w2")
    assert s["err_scatter"] > s["err_formal"]
    assert s["err"] == pytest.approx(s["err_scatter"])


def test_summarise_visit_clips_cosmic_ray():
    mag = np.array([10.0, 10.01, 9.99, 10.0, 10.02, 9.98, 10.0, 10.01, 9.99, 8.0])
    s = epochs.summarise_visit(np.zeros(mag.size), mag, np.full(mag.size, 0.02), "w1")
    assert s["n_clipped"] == 1 and s["n"] == 9
    assert s["scatter"] < 0.02


def test_summarise_visit_needs_min_frames():
    s = epochs.summarise_visit(np.zeros(4), np.full(4, 10.0), np.full(4, 0.02), "w1")
    assert s["n"] == 4 and np.isnan(s["mag"]) and np.isnan(s["err"])
    s = epochs.summarise_visit(np.zeros(4), np.full(4, 10.0), np.full(4, 0.02), "w1",
                               EpochBinning(min_frames=3))
    assert s["mag"] == pytest.approx(10.0)


def test_error_model_calibration_on_gaussian_noise():
    """For pure Gaussian noise with correct catalogue errors, the scatter of the
    visit medians should match the formal error (within sampling noise)."""
    rng = np.random.default_rng(8)
    sigma, n = 0.02, 13
    meds, formal = [], []
    for _ in range(3000):
        s = epochs.summarise_visit(np.zeros(n), 10.0 + rng.normal(0, sigma, n),
                                   np.full(n, sigma), "w1")
        meds.append(s["mag"])
        formal.append(s["err_formal"])
    ratio = np.std(meds) / np.mean(formal)
    # the asymptotic sqrt(pi/2) factor is slightly conservative at n = 13
    assert 0.85 < ratio < 1.02


def test_bin_band_shapes_and_missing_visits():
    mjd = np.concatenate([np.linspace(0, 1, 8), np.linspace(180, 181, 2), np.linspace(360, 361, 8)])
    ids = epochs.assign_epochs(mjd)
    out = epochs.bin_band(mjd, np.full(mjd.size, 11.0), np.full(mjd.size, 0.03), ids, "w1", 3)
    assert out["mag"].shape == (3,)
    assert out["n"].tolist() == [8, 2, 8]
    assert np.isnan(out["mag"][1]) and out["mag"][0] == pytest.approx(11.0)
