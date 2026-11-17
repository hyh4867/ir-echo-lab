import numpy as np
import pytest

from irecho import differential as d
from irecho.photometry import mag_to_flux_jy
from irecho.stats import constant_chi2

MJD = 56700.0 + 182.0 * np.arange(20)


def make_star(rng, common, mag=10.5, sigma=0.006, shift_days=0.0):
    out = {"mjd": MJD + shift_days}
    for b, scale in (("w1", 1.0), ("w2", 1.5)):
        m = mag + scale * common + rng.normal(0, sigma, MJD.size)
        out[f"{b}_mag"] = m
        for key in ("err", "err_formal", "err_scatter"):
            out[f"{b}_{key}"] = np.full(MJD.size, sigma)
        out[f"{b}_flux_jy"] = mag_to_flux_jy(m, b)
        out[f"{b}_flux_err_jy"] = out[f"{b}_flux_jy"] * 0.921 * sigma
    return out


@pytest.fixture
def field():
    rng = np.random.default_rng(11)
    common = rng.normal(0, 0.01, MJD.size)  # shared zero-point error per visit
    stars = [make_star(rng, common, shift_days=rng.uniform(-2, 2)) for _ in range(20)]
    target = make_star(rng, common, mag=10.2)
    return common, stars, target


def test_offsets_recover_the_common_mode(field):
    common, stars, _ = field
    offs = d.field_offsets(stars)
    assert offs.mjd.size == MJD.size
    truth = common - np.median(common)
    resid = offs.offset["w1"] - np.median(offs.offset["w1"]) - truth
    assert np.abs(resid).max() < 5 * offs.offset_err["w1"].max()
    assert offs.offset_err["w1"].mean() == pytest.approx(1.2533 * 0.006 / np.sqrt(20), rel=0.35)
    assert np.all(offs.n_stars["w2"] == 20)


def test_correction_removes_shared_scatter(field):
    _, stars, target = field
    corrected = d.apply_offsets(target, d.field_offsets(stars))
    assert np.std(target["w2_mag"]) > 2 * np.std(corrected["w2_mag"])
    assert np.all(corrected["w1_err"] > target["w1_err"])  # offset error is added
    flux_ratio = corrected["w1_flux_jy"] / mag_to_flux_jy(corrected["w1_mag"], "w1")
    assert flux_ratio == pytest.approx(1.0)
    assert corrected["color"] == pytest.approx(corrected["w1_mag"] - corrected["w2_mag"])


def test_leave_one_out_gives_chi2_near_one(field):
    _, stars, _ = field
    before = np.median([constant_chi2(s["w2_mag"], s["w2_err"]).chi2_red for s in stars])
    after = np.median([constant_chi2(s["w2_mag"], s["w2_err"]).chi2_red
                       for s in d.leave_one_out(stars)])
    assert before > 3.0
    assert after == pytest.approx(1.0, abs=0.35)


def test_floor_and_error_model_choice(field):
    _, stars, target = field
    target = dict(target, w1_err_scatter=np.full(MJD.size, 0.002))
    out = d.apply_offsets(target, d.field_offsets(stars), err_key="err_scatter",
                          floor_mag={"w1": 0.01})
    assert np.all(out["w1_err"] >= np.hypot(0.002, 0.01))
    assert np.all(out["w1_err"] < 0.0125)


def test_visits_without_a_field_visit_are_dropped(field):
    _, stars, target = field
    late = dict(target, mjd=target["mjd"].copy())
    late["mjd"][-1] += 60.0  # no star visit within 15 days
    out = d.apply_offsets(late, d.field_offsets(stars))
    assert out["mjd"].size == MJD.size - 1


def test_too_few_stars_gives_no_offset(field):
    _, stars, target = field
    offs = d.field_offsets(stars[:2])
    assert np.isnan(offs.offset["w1"]).all()
    assert d.apply_offsets(target, offs)["mjd"].size == 0


def test_a_variable_star_does_not_bias_the_median(field):
    common, stars, _ = field
    stars = [dict(s) for s in stars]
    stars[0]["w1_mag"] = stars[0]["w1_mag"] + np.linspace(0, 0.5, MJD.size)
    clean = d.field_offsets(stars[1:]).offset["w1"]
    dirty = d.field_offsets(stars).offset["w1"]
    assert np.abs(dirty - clean).max() < 0.004
