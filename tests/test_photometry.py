import numpy as np
import pytest

from irecho import photometry as ph
from irecho.config import F0_JY, LAMBDA_ISO_UM


def test_zero_magnitude_gives_zero_point():
    assert ph.mag_to_flux_jy(0.0, "w1") == pytest.approx(309.540)
    assert ph.mag_to_flux_jy(0.0, "w2") == pytest.approx(171.787)


def test_five_magnitudes_is_factor_100():
    assert ph.mag_to_flux_jy(5.0, "w1") == pytest.approx(F0_JY["w1"] / 100.0)


def test_mag_flux_round_trip():
    mags = np.array([8.5, 10.279, 14.0])
    for band in ("w1", "w2"):
        assert ph.flux_jy_to_mag(ph.mag_to_flux_jy(mags, band), band) == pytest.approx(mags)


def test_error_propagation_small_error_limit():
    flux = ph.mag_to_flux_jy(10.0, "w1")
    ferr = ph.mag_err_to_flux_err(flux, 0.01)
    # 0.01 mag is 0.921 per cent in flux
    assert ferr / flux == pytest.approx(0.00921034, rel=1e-5)
    assert ph.flux_err_to_mag_err(flux, ferr) == pytest.approx(0.01)
    # agrees with a finite difference of the conversion
    numeric = (ph.mag_to_flux_jy(10.0 - 0.005, "w1") - ph.mag_to_flux_jy(10.0 + 0.005, "w1"))
    assert ferr == pytest.approx(numeric, rel=1e-4)


def test_color_and_error():
    color, err = ph.color_w1_w2(10.3, 0.03, 10.0, 0.04)
    assert color == pytest.approx(0.3)
    assert err == pytest.approx(0.05)


def test_planck_matches_wien_and_rayleigh_jeans_limits():
    lam = 3.3526
    nu = ph.C_LIGHT / (lam * 1e-6)
    # Rayleigh-Jeans: B_nu -> 2 nu^2 k T / c^2 for h nu << k T
    t_hot = 1e7
    assert ph.planck_nu(lam, t_hot) == pytest.approx(
        2 * nu**2 * ph.K_BOLTZ * t_hot / ph.C_LIGHT**2, rel=1e-3)
    # Wien: B_nu -> (2 h nu^3 / c^2) exp(-h nu / k T) for h nu >> k T
    t_cold = 200.0
    wien = 2 * ph.H_PLANCK * nu**3 / ph.C_LIGHT**2 * np.exp(-ph.H_PLANCK * nu / (ph.K_BOLTZ * t_cold))
    assert ph.planck_nu(lam, t_cold) == pytest.approx(wien, rel=1e-6)


def test_blackbody_ratio_is_monotonic_in_temperature():
    temps = np.linspace(200, 5000, 200)
    ratios = np.array([ph.blackbody_flux_ratio(t) for t in temps])
    assert np.all(np.diff(ratios) > 0)


@pytest.mark.parametrize("temperature", [300.0, 800.0, 1200.0, 1500.0, 1800.0, 3000.0])
def test_blackbody_temperature_round_trip(temperature):
    ratio = ph.blackbody_flux_ratio(temperature)
    assert ph.blackbody_temperature(ratio) == pytest.approx(temperature, rel=1e-6)


def test_blackbody_temperature_out_of_range_is_nan():
    assert np.isnan(ph.blackbody_temperature(-1.0))
    assert np.isnan(ph.blackbody_temperature(float("nan")))
    # bluer than the Rayleigh-Jeans limit (lambda2/lambda1)^2 cannot be a blackbody
    rj_limit = (LAMBDA_ISO_UM["w2"] / LAMBDA_ISO_UM["w1"]) ** 2
    assert np.isnan(ph.blackbody_temperature(rj_limit * 1.01))


def test_color_to_flux_ratio():
    assert ph.color_to_flux_ratio(0.0) == pytest.approx(309.540 / 171.787)
    # 2.5 mag redder -> W2 ten times brighter relative to W1
    assert ph.color_to_flux_ratio(2.5) == pytest.approx(309.540 / 171.787 / 10.0)
