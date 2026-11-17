import numpy as np
import pytest

from irecho import physics
from irecho.config import LAMBDA_ISO_UM
from irecho.photometry import planck_nu


def sphere_fluxes_jy(t_k, radius_m, dist_m):
    """F_nu = pi B_nu(T) (R/D)^2 for a blackbody sphere, in Jy."""
    return [np.pi * planck_nu(LAMBDA_ISO_UM[b], t_k) * (radius_m / dist_m) ** 2 / physics.JY
            for b in ("w1", "w2")]


@pytest.mark.parametrize("t_k,r_ld", [(800.0, 10.0), (1100.0, 25.0), (1500.0, 60.0)])
def test_blackbody_round_trip(t_k, r_ld):
    d_mpc = 40.0
    f1, f2 = sphere_fluxes_jy(t_k, r_ld * physics.LIGHT_DAY_M, d_mpc * physics.MPC_M)
    q = physics.blackbody_quantities(f1, f2, 0.0, d_mpc)
    assert q["t_rest_k"] == pytest.approx(t_k, rel=1e-6)
    assert q["r_bb_light_days"] == pytest.approx(r_ld, rel=1e-6)
    l_expected = 4 * np.pi * (r_ld * physics.LIGHT_DAY_M) ** 2 * physics.SIGMA_SB * t_k**4 * 1e7
    assert q["l_dust_erg_s"] == pytest.approx(l_expected, rel=1e-6)


def test_redshift_raises_rest_temperature_only():
    f1, f2 = sphere_fluxes_jy(1000.0, 20 * physics.LIGHT_DAY_M, 40 * physics.MPC_M)
    q0 = physics.blackbody_quantities(f1, f2, 0.0, 40.0)
    q1 = physics.blackbody_quantities(f1, f2, 0.01, 40.0)
    assert q1["t_obs_k"] == pytest.approx(q0["t_obs_k"])
    assert q1["t_rest_k"] == pytest.approx(1.01 * q0["t_rest_k"])
    assert q1["l_dust_erg_s"] == pytest.approx(q0["l_dust_erg_s"])
    assert q1["r_bb_light_days"] == pytest.approx(q0["r_bb_light_days"] / 1.01**2)


def test_scalings():
    """Doubling both fluxes doubles L and multiplies R by sqrt(2); L scales as D^2."""
    q1 = physics.blackbody_quantities(0.010, 0.012, 0.0, 40.0)
    q2 = physics.blackbody_quantities(0.020, 0.024, 0.0, 40.0)
    q3 = physics.blackbody_quantities(0.010, 0.012, 0.0, 80.0)
    assert q2["t_rest_k"] == pytest.approx(q1["t_rest_k"])
    assert q2["l_dust_erg_s"] == pytest.approx(2 * q1["l_dust_erg_s"])
    assert q2["r_bb_light_days"] == pytest.approx(np.sqrt(2) * q1["r_bb_light_days"])
    assert q3["l_dust_erg_s"] == pytest.approx(4 * q1["l_dust_erg_s"])


def test_non_positive_excess_is_nan():
    for f1, f2 in [(-0.01, 0.01), (0.01, 0.0), (np.nan, 0.01)]:
        assert np.isnan(physics.blackbody_quantities(f1, f2, 0.0, 40.0)["t_rest_k"])


def test_luminosity_distance_low_z_limit():
    z = 0.001
    hubble = physics.C_LIGHT / 1e3 * z / physics.H0_KM_S_MPC
    assert physics.luminosity_distance_mpc(z) == pytest.approx(hubble, rel=2e-3)
    assert physics.luminosity_distance_mpc(0.009083) == pytest.approx(40.68, abs=0.01)


def test_flux_correction_table_and_clamp():
    assert physics.flux_correction(800.0, "w1") == pytest.approx(0.9884)
    assert physics.flux_correction(1131.0, "w2") == pytest.approx(0.9877)
    assert physics.flux_correction(5000.0, "w1") == pytest.approx(0.9801)
    mid = physics.flux_correction(950.0, "w1")
    assert 0.9801 < mid < 0.9884


def test_color_correction_shifts_temperature_by_about_one_percent():
    q0 = physics.blackbody_quantities(0.0146, 0.0173, 0.009, 40.7)
    q1 = physics.blackbody_quantities(0.0146, 0.0173, 0.009, 40.7, color_correct=True)
    assert abs(q1["t_rest_k"] / q0["t_rest_k"] - 1) < 0.02


def test_monte_carlo_brackets_point_estimate():
    rng = np.random.default_rng(5)
    flux = {"w1": 0.0362, "w2": 0.0317}
    base = {"w1": 0.0216, "w2": 0.0144}
    err = {"w1": 0.0003, "w2": 0.0003}
    mc = physics.monte_carlo(flux, err, base, err, 0.009, 40.7, 2000, rng)
    q = physics.blackbody_quantities(flux["w1"] - base["w1"], flux["w2"] - base["w2"], 0.009, 40.7)
    assert mc["valid_draw_fraction"] == 1.0
    for key in ("t_rest_k", "l_dust_erg_s", "r_bb_light_days"):
        assert mc[f"{key}_p16"] < q[key] < mc[f"{key}_p84"]
