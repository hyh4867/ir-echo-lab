import pytest

from irecho import theory


def test_parsec_in_light_days():
    assert theory.PC_IN_LIGHT_DAYS == pytest.approx(1191.3, rel=1e-4)


def test_sublimation_radius_reference_values():
    # Netzer (2015) Eq. 1 at L = 1e46 erg/s and the reference temperatures
    assert theory.sublimation_radius_pc(1e46, "graphite") == pytest.approx(0.5)
    assert theory.sublimation_radius_pc(1e46, "silicate") == pytest.approx(1.3)


def test_sublimation_radius_scales_as_sqrt_luminosity():
    r1 = theory.sublimation_radius_pc(1e43)
    r2 = theory.sublimation_radius_pc(1e45)
    assert r2 / r1 == pytest.approx(10.0)


def test_sublimation_radius_temperature_and_anisotropy():
    hot = theory.sublimation_radius_pc(1e46, "graphite", t_sub_k=1500.0)
    assert hot == pytest.approx(0.5 * (1800.0 / 1500.0) ** 2.6)
    assert theory.sublimation_radius_pc(1e46, f_theta=0.5) == pytest.approx(0.25)


def test_lyu2019_lag_normalisation_and_slope():
    l11 = 1e11 * theory.L_SUN_ERG_S
    assert theory.lyu2019_lag_days(l11, "w1") == pytest.approx(10**2.10)
    assert theory.lyu2019_lag_days(l11, "w2") == pytest.approx(10**2.20)
    assert theory.lyu2019_lag_days(100 * l11, "w1") / theory.lyu2019_lag_days(l11, "w1") == (
        pytest.approx(100**0.47))


def test_light_travel_days():
    assert theory.light_travel_days(1.0) == pytest.approx(theory.PC_IN_LIGHT_DAYS)
