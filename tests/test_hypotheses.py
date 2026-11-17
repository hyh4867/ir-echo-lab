"""Decision rules of docs/PREDICTIONS.md and the ZTF binning, on synthetic input."""

import numpy as np
import pytest

from irecho import stage34, ztf

YEARS = 2020.9 + 0.5 * np.arange(8)


def phys(temps, half_width=30.0):
    return {"flagged_visits": [
        {"date": f"visit{i}", "year": float(y), "t_rest_k": float(t),
         "t_rest_k_p16": float(t - half_width), "t_rest_k_p84": float(t + half_width)}
        for i, (y, t) in enumerate(zip(YEARS, temps))]}


def verdicts(temps, **kw):
    h = stage34.hypothesis_tests(phys(temps, **kw))
    return (h["H1_temperature_mixing"]["verdict"], h["H2_silicate"]["verdict"],
            h["H3_dust_destruction"]["verdict"])


def test_steady_cooling_supports_h1_and_rejects_h2():
    temps = 1400.0 - 120.0 * (YEARS - YEARS[0])  # 1400 K -> 980 K
    h1, h2, _ = verdicts(temps)
    assert h1 == "supported"
    assert h2 == "rejected"  # not constant even with the 70 K systematic


def test_constant_temperature_supports_h2_and_rejects_h1():
    rng = np.random.default_rng(1)
    temps = 1100.0 + rng.normal(0, 10.0, 8)
    h1, h2, h3 = verdicts(temps)
    assert h1 == "rejected"  # a decline of 100 K is excluded
    assert h2 == "supported"
    assert h3 == "supported"


def test_constant_but_outside_the_range_rejects_h2():
    assert verdicts(np.full(8, 1400.0) + np.linspace(-5, 5, 8))[1] == "rejected"


def test_scatter_between_statistical_and_systematic_is_undetermined():
    temps = 1100.0 + np.array([90, -90, 60, -60, 90, -90, 60, -60.0])
    h1, h2, _ = verdicts(temps)
    assert h2 == "undetermined"  # rejected with 30 K errors, not with 70 K added
    assert h1 == "undetermined"  # scatter too large to exclude a 100 K decline


def test_h3_limits():
    hot = np.array([1800.0, 1750.0, 1200, 1200, 1200, 1200, 1200, 1200])
    near = np.array([1450.0, 1400.0, 1200, 1200, 1200, 1200, 1200, 1200])
    cool = np.array([1000.0, 950.0, 1300, 1000, 1000, 1000, 1000, 1000])
    assert verdicts(hot)[2] == "rejected"
    assert verdicts(near)[2] == "undetermined"
    assert verdicts(cool)[2] == "supported"
    h3 = stage34.hypothesis_tests(phys(cool))["H3_dust_destruction"]
    assert h3["rerise_dates_2sigma_statistical"] == ["visit2"]
    assert h3["rerise_dates_2sigma_with_systematic"] == ["visit2"]


def test_h1_numbers_on_an_exact_line():
    temps = 1200.0 - 50.0 * (YEARS - YEARS[0])
    h1 = stage34.hypothesis_tests(phys(temps))["H1_temperature_mixing"]
    assert h1["slope_k_per_yr"] == pytest.approx(-50.0)
    assert h1["delta_t_k"] == pytest.approx(-175.0)
    assert h1["mann_kendall_z"] < -2.5


def test_ztf_binning_to_visits():
    visits = np.array([58000.0, 58180.0, 58360.0])
    mjd = np.concatenate([58000 + np.arange(-20, 20, 4.0), 58180 + np.arange(0, 12, 4.0),
                          [58275.0]])  # 10 points, 3 points, 1 point 85/95 d from visits
    mag = np.concatenate([np.full(10, 15.0), np.full(3, 15.5), [16.0]])
    mag[0] = 18.0  # an outlier does not move the median
    bins = ztf.bin_to_visits(mjd, mag, visits)
    assert list(bins["n"]) == [10, 3, 1]
    assert bins["mag"][0] == pytest.approx(15.0)
    assert np.isnan(bins["mag"][1]) and np.isnan(bins["mag"][2])  # fewer than 5 points
    assert bins["mjd"][0] == pytest.approx(58000.0 - 2.0)


def test_ztf_read_clean_keeps_one_object_id(tmp_path):
    path = tmp_path / "x_g.csv"
    rows = ["oid,mjd,mag,magerr,catflags"]
    rows += [f"1,{58000 + i},15.0,0.01,0" for i in range(6)]
    rows += [f"1,{58010 + i},15.0,0.01,32768" for i in range(3)]
    rows += [f"2,{58020 + i},15.3,0.01,0" for i in range(2)]
    path.write_text("\n".join(rows) + "\n")
    lc = ztf.read_clean(path)
    assert lc["oid"] == "1" and lc["n_rows_total"] == 11
    assert lc["n_rows_clean"] == 8 and lc["n_rows_used"] == 6
    assert np.all(np.diff(lc["mjd"]) > 0)
