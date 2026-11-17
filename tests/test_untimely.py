import numpy as np
import pytest

from irecho import untimely as ut
from irecho.photometry import mag_to_flux_jy


def row(epoch, band, flux, sep=0.2, primary=1, mjd=None):
    return {"EPOCH": epoch, "band": band, "flux": flux, "dflux": flux * 0.001,
            "sep_arcsec": sep, "primary": primary, "MJDMEAN": mjd or 56800.0 + 182 * epoch,
            "fracflux": 1.0, "flags_info": 66}


def test_nanomaggy_zero_point_matches_the_pipeline():
    flux_nmgy = 10.0 ** ((22.5 - 10.0) / 2.5)  # a 10th-magnitude source
    for b in ("w1", "w2"):
        assert ut.nanomaggy_to_jy(flux_nmgy, b) == pytest.approx(mag_to_flux_jy(10.0, b))


def test_light_curve_picks_nearest_primary_detection():
    rows = [row(0, 1, 1e5), row(0, 2, 2e5),
            row(1, 1, 1.2e5, sep=3.0), row(1, 1, 5e3, sep=0.4),  # galaxy split in two
            row(1, 1, 9e9, sep=0.1, primary=0),  # duplicate from an overlapping tile
            row(1, 2, 2.4e5), row(2, 1, 1e5, sep=9.0)]  # too far
    lc = ut.light_curve(rows)
    assert list(lc["epoch"]) == [0, 1, 2]
    assert lc["w1_mag"][0] == pytest.approx(22.5 - 2.5 * 5)
    assert lc["w1_mag"][1] == pytest.approx(22.5 - 2.5 * np.log10(5e3))
    assert lc["w1_n_det"][1] == 2
    assert np.isnan(lc["w1_mag"][2]) and lc["w1_n_det"][2] == 0
    assert lc["w1_err"][0] == pytest.approx(2.5 / np.log(10) * 0.001)


def test_match_visits():
    a = np.array([100.0, 280.0, 460.0])
    b = np.array([102.0, 300.0, 461.0, 640.0])
    assert ut.match_visits(a, b) == [(0, 0), (2, 2)]
