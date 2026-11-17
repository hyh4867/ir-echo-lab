"""Stage 5a sample definition: parsers, coordinates and cuts (no network)."""

import numpy as np
import pytest

from irecho import sample
from irecho.sample import CRITERIA, SampleCriteria

BAT_TEXT = """
    |             BAT_NAME|     RA|    DEC|    SNR|              COUNTERPART_NAME|                     OTHER_NAME| CTPT_RA|CTPT_DEC|    FLUX| FLUX_LO| FLUX_HI|CONTA|GAMM|GAMM|GAMM|  CHI_SQ_R|REDSHIF|  LUM|ASSOC_STREN|        CL2|                 TYPE
    |                     |    deg|    deg|       |                              |                               |     deg|     deg|10^-12er|10^-12er|10^-12er|     |    |    |    |          |       |log[e|           |           |
   1|   SWIFT J0001.0-0708|  0.228| -7.164|   7.73|       2MASX J00004876-0709117|                               |  0.2032| -7.1532|   13.49|   11.10|   16.01|  0.0|2.21|1.80|2.70|      1.30| 0.0375|43.64|           |         50|Sy1.9
1630|   SWIFT J2343.0-5920|355.738|-59.325|   5.37|APMUKS(BJ) B233953.42-593202.5|                               |355.6517|-59.2565|    6.97|    4.80|    8.75|  0.0|2.49|1.80|3.53|      0.50|       |     |           |         70|Unknown AGN
"""

VIZIER_TEXT = (
    "#comment\n"
    "recno\tID\tm_ID\tRAJ2000\tDEJ2000\tType\tz\tDist\tlogMBH\tMethod\n"
    " \t \t \tdeg\tdeg\t \t \tMpc\t[Msun]\t \n"
    "-----\t----\t-\t-----------\t-------------\t---------\t-----------\t-------\t-----\t--------\n"
    "    1\t   1\t \t000.2032337\t-07.153222700\tSy1.9    \t 0.03749562\t  165.2\t 7.61\tVdisp   \n"
    "    2\t   9\t \t010.0000000\t+20.000000000\tBZQ      \t 0.50000000\t 2800.0\t     \t        \n"
)


def test_parse_bat105_reads_fields_and_blank_redshift():
    rows = sample.parse_bat105(BAT_TEXT)
    assert [r["bat_index"] for r in rows] == [1, 1630]
    assert rows[0]["counterpart"] == "2MASX J00004876-0709117"
    assert rows[0]["ctpt_ra"] == pytest.approx(0.2032)
    assert rows[0]["z"] == pytest.approx(0.0375)
    assert rows[0]["bat_class"] == 50 and rows[0]["bat_type"] == "Sy1.9"
    assert np.isnan(rows[1]["z"]) and rows[1]["bat_class"] == 70


def test_parse_vizier_tsv_and_match():
    dr2 = sample.parse_vizier_tsv(VIZIER_TEXT)
    assert len(dr2) == 2 and dr2[0]["Type"] == "Sy1.9" and dr2[1]["logMBH"] == ""
    row = sample.parse_bat105(BAT_TEXT)[0]
    match = sample.match_dr2(row, dr2, 5.0)
    assert match["ID"] == "1" and match["sep_arcsec"] < 1.0
    far = dict(row, ctpt_dec=row["ctpt_dec"] + 0.01)  # 36 arcsec away
    assert sample.match_dr2(far, dr2, 5.0) is None


def test_galactic_and_ecliptic_latitude_known_points():
    # North Galactic pole and Galactic centre (Sgr A*).
    assert sample.galactic_latitude_deg(192.85948, 27.12825) == pytest.approx(90.0, abs=1e-4)
    assert sample.galactic_latitude_deg(266.41683, -29.00781) == pytest.approx(-0.05, abs=0.02)
    # North ecliptic pole (RA 18h, Dec 90 - obliquity) and the vernal equinox.
    assert sample.ecliptic_latitude_deg(270.0, 66.5607089) == pytest.approx(90.0, abs=1e-4)
    assert sample.ecliptic_latitude_deg(0.0, 0.0) == pytest.approx(0.0, abs=1e-9)
    assert sample.ecliptic_latitude_deg(90.0, 23.4392911) == pytest.approx(0.0, abs=1e-6)


def test_sep_arcsec_small_and_large():
    assert float(sample.sep_arcsec(10.0, 0.0, 10.0, 1.0 / 3600.0)) == pytest.approx(1.0, rel=1e-6)
    assert float(sample.sep_arcsec(10.0, 60.0, 10.0 + 2.0 / 3600.0, 60.0)) == pytest.approx(
        1.0, rel=1e-4)
    assert float(sample.sep_arcsec(0.0, 0.0, 180.0, 0.0)) == pytest.approx(180 * 3600.0)


def _row(**kw):
    base = {"bat_class": 50, "z": 0.02, "ctpt_ra": 150.0, "ctpt_dec": 30.0}
    return {**base, **kw}


def test_catalogue_cut_order():
    assert sample.catalogue_cut(_row(), None) == ""
    assert sample.catalogue_cut(_row(bat_class=80), None) == "S1 class"
    assert sample.catalogue_cut(_row(z=float("nan")), None) == "S2 redshift"
    assert sample.catalogue_cut(_row(z=0.0501), None) == "S2 redshift"
    assert sample.catalogue_cut(_row(z=0.05), None) == ""
    # Galactic plane (Sgr A*) and north ecliptic pole.
    assert sample.catalogue_cut(_row(ctpt_ra=266.4, ctpt_dec=-29.0), None) == "S3 Galactic latitude"
    assert sample.catalogue_cut(_row(ctpt_ra=270.0, ctpt_dec=66.56), None) == "S4 ecliptic latitude"
    assert sample.catalogue_cut(_row(), {"Type": "BZQ"}) == "S5 BASS DR2 beamed/lensed"
    assert sample.catalogue_cut(_row(), {"Type": "Sy2"}) == ""


def _field(rows):
    keys = ("ra", "dec", "w1mpro", "w2mpro")
    out = {k: np.array([r[k] for r in rows], dtype=float) for k in keys}
    out["cc_flags"] = np.array([r.get("cc_flags", "0000") for r in rows])
    out["ph_qual"] = np.array([r.get("ph_qual", "AAAA") for r in rows])
    return out


def _src(dra_arcsec=0.0, ddec_arcsec=0.0, w1=10.0, w2=9.5, **kw):
    return {"ra": 150.0 + dra_arcsec / 3600.0 / np.cos(np.deg2rad(30.0)),
            "dec": 30.0 + ddec_arcsec / 3600.0, "w1mpro": w1, "w2mpro": w2, **kw}


def test_allwise_cut_pass_and_counterpart_choice():
    # Two sources inside 3 arcsec: the brighter one is adopted.
    fail, i = sample.allwise_cut(_field([_src(1.0, 0.0, w1=13.0, w2=13.0), _src(0.5, 0.5)]),
                                 150.0, 30.0, SampleCriteria(blend_radius_arcsec=0.0))
    assert fail == "" and i == 1


@pytest.mark.parametrize("rows, expected", [
    ([_src(5.0, 0.0)], "S6 no AllWISE counterpart"),
    ([_src(w1=8.4)], "S7 saturation margin"),
    ([_src(w2=7.9)], "S7 saturation margin"),
    ([_src(w2=float("nan"))], "S7 saturation margin"),
    ([_src(w1=12.6, w2=12.0)], "S8 faint limit"),
    ([_src(cc_flags="0h00")], "S9 cc_flags / ph_qual"),
    ([_src(ph_qual="ABAA")], "S9 cc_flags / ph_qual"),
    ([_src(), _src(250.0, 0.0, w1=4.5, w2=4.5)], "S10 bright star nearby"),
    ([_src(), _src(10.0, 0.0, w1=12.0, w2=12.0)], "S10 blended neighbour"),
])
def test_allwise_cut_failures(rows, expected):
    assert sample.allwise_cut(_field(rows), 150.0, 30.0)[0] == expected


def test_allwise_cut_ignores_harmless_neighbours():
    rows = [_src(),
            _src(10.0, 0.0, w1=12.6, w2=12.6),  # inside 12 arcsec but > 2.5 mag fainter
            _src(20.0, 0.0, w1=10.5, w2=10.5),  # bright but outside 12 arcsec
            _src(0.0, 310.0, w1=4.0, w2=4.0)]  # very bright but outside 300 arcsec
    assert sample.allwise_cut(_field(rows), 150.0, 30.0) == ("", 0)
    assert sample.allwise_cut(_field([]), 150.0, 30.0) == ("S6 no AllWISE counterpart", -1)


def test_choose_z_max_rule():
    c = SampleCriteria(z_grid=(0.01, 0.02, 0.03), max_sample_size=5)
    z = [0.005, 0.008, 0.012, 0.015, 0.019, 0.021, 0.025, 0.029]
    z_max, counts = sample.choose_z_max(z, c)
    assert counts == {"0.010": 2, "0.020": 5, "0.030": 8}
    assert z_max == 0.02  # largest grid value with at most 5 objects
    assert sample.choose_z_max(z, SampleCriteria(z_grid=(0.01, 0.02, 0.03),
                                                 max_sample_size=100))[0] == 0.03
    # Even the smallest grid value is too many: fall back to the smallest.
    assert sample.choose_z_max(z, SampleCriteria(z_grid=(0.01, 0.02), max_sample_size=1))[0] == 0.01


def test_default_criteria_are_the_documented_ones():
    assert CRITERIA.bat_classes == (40, 50, 60)
    assert max(CRITERIA.z_grid) == 0.05 and CRITERIA.max_sample_size == 250
    assert (CRITERIA.w1_bright_limit, CRITERIA.w2_bright_limit) == (8.5, 8.0)
    assert CRITERIA.min_abs_gal_lat_deg == 20.0 and CRITERIA.max_abs_ecl_lat_deg == 80.0
