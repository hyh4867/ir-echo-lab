"""Stage 6B: sample_v2 cut, the lag-luminosity reference and the decision rule."""

import sys
from pathlib import Path

import pytest

from irecho import sample

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import run_d_lum_power as power  # noqa: E402
import run_stage6b as s6  # noqa: E402


def row(**kw):
    base = {"bat_class": 50, "z": 0.08, "ctpt_ra": 150.0, "ctpt_dec": 30.0, "bat_index": 999}
    return {**base, **kw}


def test_sample_v2_redshift_range_and_exclusion():
    cut = lambda r, m=None: sample.catalogue_cut_v2(r, m, {5})  # noqa: E731
    assert cut(row()) == ""
    assert cut(row(z=0.04)) == "S2 redshift"  # the lower edge belongs to sample_v1
    assert cut(row(z=0.0401)) == "" and cut(row(z=0.30)) == ""
    assert cut(row(z=0.301)) == "S2 redshift" and cut(row(z=float("nan"))) == "S2 redshift"
    assert cut(row(bat_class=80)) == "S1 class"
    assert cut(row(bat_index=5)) == "E0 already in sample_v1"
    assert cut(row(ctpt_ra=266.4, ctpt_dec=-29.0)) == "S3 Galactic latitude"
    assert cut(row(), {"Type": "BZQ"}) == "S5 BASS DR2 beamed/lensed"
    assert (sample.V2.z_min_exclusive, sample.V2.z_max) == (0.04, 0.30)


def test_lyu_lags_reference_values():
    log_l11 = 11.0 + 33.583  # 1e11 L_sun in erg/s (log10 3.828e33 = 33.583)
    w1, w2 = power.lyu_lags(log_l11)
    assert w1 == pytest.approx(10 ** 2.10, rel=1e-3) and w2 == pytest.approx(10 ** 2.20, rel=1e-3)
    w1, w2 = power.lyu_lags(log_l11 + 1.0)
    assert w2 - w1 == pytest.approx(75.1, abs=0.2)
    assert set(power.MODELS) == {"model1", "model1_half", "model1_double", "model2"}


@pytest.mark.parametrize("rho, q1, q2, n, expected", [
    (0.05, 0.30, 0.25, 40, "사용자 예측 지지"),
    (0.45, 0.30, 0.25, 40, "사용자 예측 기각"),
    (0.27, 0.30, 0.25, 40, "판정불가"),       # positive but weaker than model 1
    (0.25, 0.17, 0.38, 20, "판정불가"),       # overlap of the two models at small n
    (0.10, 0.17, 0.38, 20, "사용자 예측 지지"),
    (0.50, 0.17, 0.38, 20, "사용자 예측 기각"),
    (0.60, 0.17, 0.38, 19, "판정불가 (표본 부족)"),
])
def test_decision_rule(rho, q1, q2, n, expected):
    assert s6.verdict(rho, q1, q2, n) == expected


def test_threshold_interpolation():
    by_n = {"20": {"rho_p5": 0.1}, "40": {"rho_p5": 0.3}}
    assert s6.interp_threshold(by_n, 30, "rho_p5") == pytest.approx(0.2)
    assert s6.interp_threshold(by_n, 10, "rho_p5") == pytest.approx(0.1)  # clamped at the grid edge
