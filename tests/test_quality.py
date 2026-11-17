import numpy as np
import pytest

from conftest import DEC0, RA0, make_allwise_table, make_neowise_table
from irecho import quality
from irecho.config import QualityCuts


@pytest.fixture
def obs():
    return quality.normalise(make_neowise_table([56700.0], n_frames=10), "neowise")


def test_clean_table_passes_everything(obs):
    for band in ("w1", "w2"):
        mask, failed = quality.band_mask(obs, band)
        assert mask.all()
        assert failed["kept"] == failed["total_rows"] == 10


@pytest.mark.parametrize(
    "column, value, criterion, bands_hit",
    [
        ("qual_frame", 0.0, "qual_frame", ("w1", "w2")),
        ("qi_fact", 0.0, "qi_fact", ("w1", "w2")),
        ("saa_sep", 3.0, "saa_sep", ("w1", "w2")),
        ("saa_sep", -10.0, "saa_sep", ("w1", "w2")),
        ("sso_flg", 1.0, "sso_flg", ("w1", "w2")),
        ("moon_masked", "10", "moon_masked", ("w1",)),
        ("moon_masked", "01", "moon_masked", ("w2",)),
        ("cc_flags", "D000", "cc_flags", ("w1",)),
        ("cc_flags", "0h00", "cc_flags", ("w2",)),
        ("w1sat", 0.05, "saturation", ("w1",)),
        ("w2err", np.nan, "measured", ("w2",)),  # upper limit: null uncertainty
        ("w1mag", np.nan, "measured", ("w1",)),
    ],
)
def test_each_cut_rejects_only_the_affected_band(obs, column, value, criterion, bands_hit):
    obs[column] = obs[column].copy()
    obs[column][3] = value
    for band in ("w1", "w2"):
        mask, failed = quality.band_mask(obs, band)
        if band in bands_hit:
            assert not mask[3] and mask.sum() == 9
            assert failed[criterion] == 1
        else:
            assert mask.all()


def test_qi_fact_half_is_kept_for_neowise_but_not_allwise(obs):
    obs["qi_fact"] = np.full(10, 0.5)
    assert quality.band_mask(obs, "w1")[0].all()
    allwise = quality.normalise(make_allwise_table([55340.0], n_frames=10), "allwise")
    allwise["qi_fact"] = np.full(10, 0.5)
    assert not quality.band_mask(allwise, "w1")[0].any()


def test_allwise_schema_is_normalised_and_passes():
    allwise = quality.normalise(make_allwise_table([55340.0], n_frames=10), "allwise")
    assert set(allwise) == set(quality.normalise(make_neowise_table([56700.0]), "neowise"))
    for band in ("w1", "w2"):
        assert quality.band_mask(allwise, band)[0].all()
    # four-character AllWISE moon mask: third character (W3) must not matter
    allwise["moon_masked"] = np.full(10, "0010")
    assert quality.band_mask(allwise, "w2")[0].all()
    allwise["moon_masked"] = np.full(10, "0100")
    assert not quality.band_mask(allwise, "w2")[0].any()


def test_cuts_can_be_switched_off(obs):
    obs["cc_flags"] = np.full(10, "DD00")
    assert not quality.band_mask(obs, "w1")[0].any()
    assert quality.band_mask(obs, "w1", QualityCuts(require_cc_clean=False))[0].all()


def test_unknown_survey_raises():
    with pytest.raises(ValueError):
        quality.normalise(make_neowise_table([56700.0]), "spitzer")


def test_nearest_per_frame_keeps_closest_detection(obs):
    # duplicate frame 0: the copy is 2 arcsec away and must be dropped
    dup = {k: np.concatenate([v, v[:1]]) for k, v in obs.items()}
    dup["ra"][-1] = RA0 + 2.0 / 3600.0
    dup["dec"][0] = DEC0
    dup["ra"][0] = RA0
    keep = quality.nearest_per_frame(dup, RA0, DEC0)
    assert keep.sum() == 10 and keep[0] and not keep[-1]
