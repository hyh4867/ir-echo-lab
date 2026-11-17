import csv

import numpy as np
import pytest

from conftest import DEC0, RA0, make_allwise_table, make_neowise_table
from irecho.flare import detect_flares
from irecho.lightcurve import EPOCH_COLUMNS, build_lightcurve, mjd_to_year, write_epochs_csv


def build(neowise_visits, allwise_visits, **kwargs):
    tables = {
        "neowise": make_neowise_table(neowise_visits, **kwargs),
        "allwise": make_allwise_table(allwise_visits),
    }
    return build_lightcurve("synthetic", RA0, DEC0, tables)


def test_visits_from_both_surveys_are_binned(neowise_visits, allwise_visits):
    lc = build(neowise_visits, allwise_visits)
    assert lc.n_epochs == 22
    assert lc.epochs["survey"].tolist() == ["allwise"] * 2 + ["neowise"] * 20
    assert lc.select("neowise").sum() == 20
    assert np.all(np.diff(lc.epochs["mjd"]) > 100)
    assert np.all(lc.epochs["w1_n"] == 12)
    assert lc.epochs["w1_mag"] == pytest.approx(np.full(22, 10.3), abs=0.03)
    assert lc.epochs["w2_mag"] == pytest.approx(np.full(22, 10.0), abs=0.03)
    assert lc.epochs["color"] == pytest.approx(np.full(22, 0.3), abs=0.04)
    assert lc.epochs["color_err"] == pytest.approx(
        np.hypot(lc.epochs["w1_err"], lc.epochs["w2_err"]))


def test_cut_report_counts_rejections(neowise_visits, allwise_visits):
    table = make_neowise_table(neowise_visits)
    table["saa_sep"][:12] = 1.0  # first visit taken near the SAA: all frames rejected
    lc = build_lightcurve("synthetic", RA0, DEC0, {"neowise": table})
    assert lc.cut_report["neowise"]["w1"]["saa_sep"] == 12
    assert lc.epochs["w1_n"][0] == 0 and np.isnan(lc.epochs["w1_mag"][0])
    assert np.isfinite(lc.epochs["w1_mag"][1:]).all()
    # the visit still exists in the table and keeps a time stamp
    assert lc.n_epochs == 20 and np.isfinite(lc.epochs["mjd"][0])


def test_end_to_end_flare_recovery_from_single_exposures(neowise_visits, allwise_visits):
    dmag = np.zeros(20)
    dmag[13:] = -np.array([0.5, 0.45, 0.4, 0.3, 0.25, 0.2, 0.15])
    lc = build(neowise_visits, allwise_visits, dmag=dmag)
    neo = lc.select("neowise")
    res = detect_flares(lc.epochs["mjd"][neo], lc.epochs["w1_mag"][neo], lc.epochs["w1_err"][neo],
                        lc.epochs["w2_mag"][neo], lc.epochs["w2_err"][neo])
    assert res.onset_index == 13
    assert res.peak("w1")[0] == pytest.approx(0.5, abs=0.03)


def test_write_epochs_csv_round_trip(tmp_path, neowise_visits, allwise_visits):
    lc = build(neowise_visits, allwise_visits)
    path = tmp_path / "out" / "epochs.csv"
    write_epochs_csv(lc, path)
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 22 and tuple(rows[0].keys()) == EPOCH_COLUMNS
    assert float(rows[5]["w1_mag"]) == pytest.approx(lc.epochs["w1_mag"][5], abs=1e-4)
    assert rows[0]["survey"] == "allwise"


def test_empty_input_raises():
    with pytest.raises(ValueError):
        build_lightcurve("none", RA0, DEC0, {})


def test_mjd_to_year():
    assert mjd_to_year(51544.5) == pytest.approx(2000.0)
    assert mjd_to_year(59000.0) == pytest.approx(2020.41, abs=0.01)
