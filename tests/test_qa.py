"""Stage 5a blind QA: synthetic fields only."""

import csv
import json

import numpy as np
import pytest

from conftest import DEC0, RA0, make_allwise_table, make_neowise_table
from irecho import differential, errorfloor, irsa, qa, sample_io
from irecho.analysis import visit_arrays
from irecho.lightcurve import build_lightcurve

VISITS = 56700.0 + 182.0 * np.arange(20)
ZP = 0.03 * np.sin(np.arange(20))  # zero-point offset shared by the whole field


def star(seed, w1=10.3):
    return build_lightcurve(f"s{seed}", RA0, DEC0, {
        "neowise": make_neowise_table(VISITS, w1=w1, w2=w1, sigma=0.02, seed=seed, dmag=ZP)})


def target(dmag=None, **kw):
    return build_lightcurve("t", RA0, DEC0, {
        "neowise": make_neowise_table(VISITS, seed=99, dmag=ZP if dmag is None else dmag, **kw),
        "allwise": make_allwise_table(np.array([55340.0, 55530.0]))})


def untimely_rows(mjds, w1_mags, w2_mags, sep=0.5):
    rows = []
    for e, (t, m1, m2) in enumerate(zip(mjds, w1_mags, w2_mags)):
        for band, m in ((1, m1), (2, m2)):
            rows.append({"sep_arcsec": sep, "band": band, "EPOCH": e, "primary": 1,
                         "MJDMEAN": float(t), "flux": 10 ** (-0.4 * (m - 22.5)),
                         "dflux": 10 ** (-0.4 * (m - 22.5)) * 0.005, "fracflux": 0.98,
                         "flags_info": 0})
    return rows


def test_star_qa_shows_the_correction_working():
    stars = [star(s) for s in range(1, 9)]
    sv = qa.star_tables(stars)
    loo = differential.leave_one_out(sv, err_key=qa.ERR_KEY)
    rows = qa.star_qa([lc.name for lc in stars], [10.3] * 8, sv, loo)
    before = np.median([r["w1_chi2_red_before"] for r in rows])
    after = np.median([r["w1_chi2_red_after"] for r in rows])
    assert before > 5.0  # a 21 mmag rms zero-point term against ~6 mmag visit errors
    assert 0.5 < after < 1.8
    assert all(r["n_visits"] == 20 and r["n_visits_corrected"] == 20 for r in rows)
    offsets = differential.field_offsets(sv)
    summary = qa.field_summary(rows, offsets)
    assert summary["n_stars"] == 8 and summary["n_field_visits"] == 20
    assert summary["w1"]["offset_std_mmag"] == pytest.approx(1e3 * np.std(ZP, ddof=1), rel=0.15)
    assert summary["w1"]["min_stars_per_visit"] == 8
    orow = qa.offset_rows("f", offsets)
    assert len(orow) == 20 and orow[0]["field"] == "f"
    assert [r["w1_offset_mmag"] for r in orow] == pytest.approx(
        1e3 * (ZP - np.median(ZP)), abs=8.0)


def _numbers(obj):
    """All leaf values of a nested dict/list."""
    if isinstance(obj, dict):
        for v in obj.values():
            yield from _numbers(v)
    elif isinstance(obj, (list, tuple)):
        for v in obj:
            yield from _numbers(v)
    else:
        yield obj


def test_target_counts_are_blind():
    # A target with a 1 mag brightening: none of the returned numbers may carry it.
    dmag = ZP.copy()
    dmag[12:] -= 1.0
    stars = [star(s) for s in range(1, 9)]
    offsets = differential.field_offsets(qa.star_tables(stars))
    lc = target(dmag=dmag)
    v = differential.apply_offsets(visit_arrays(lc, "neowise"), offsets, err_key=qa.ERR_KEY)
    out = qa.target_counts(lc, v)
    flat = target(dmag=ZP)
    v0 = differential.apply_offsets(visit_arrays(flat, "neowise"), offsets, err_key=qa.ERR_KEY)
    ref = qa.target_counts(flat, v0)
    assert out["n_visits_neowise_both"] == 20 and out["n_visits_allwise_both"] == 2
    assert out["n_visits_corrected"] == 20 and out["n_frames"] == {"neowise": 240, "allwise": 24}
    assert out["first_visit_date"] == "2014-02-12" and out["max_gap_days"] == pytest.approx(182.0)
    assert out["max_visit_duration_days"] == pytest.approx(1.5)
    assert out["w1"]["saturated_frame_fraction"] == 0.0
    assert out["w1"]["n_visits_brighter_than_saturation_margin"] == 0
    # no per-visit series anywhere in the output (only the nested cut report is a dict tree)
    assert not any(isinstance(x, (list, np.ndarray)) for x in out.values())
    json.dumps(out)  # plain numbers only
    # the varying and the flat target give the same counts and (nearly) the same error sizes
    for key in ("n_visits_corrected", "n_frames", "first_visit_date", "max_gap_days"):
        assert out[key] == ref[key]
    assert out["w1"]["median_stat_err_mmag"] == pytest.approx(ref["w1"]["median_stat_err_mmag"],
                                                              rel=0.3)
    # and no number in the output is anywhere near the target's magnitudes (9.3-10.4)
    leaves = [x for x in _numbers({k: v for k, v in out.items() if k != "cut_report"})
              if isinstance(x, float)]
    assert not any(8.0 < x < 11.0 for x in leaves)


def test_target_counts_report_saturation():
    table = make_neowise_table(VISITS, w1=8.2, w2=7.8, seed=5)
    table["w1sat"][:60] = 0.02  # a quarter of the frames have saturated pixels
    lc = build_lightcurve("bright", RA0, DEC0, {"neowise": table})
    stars = [build_lightcurve(f"s{s}", RA0, DEC0, {"neowise": make_neowise_table(VISITS, seed=s)})
             for s in range(1, 6)]
    v = differential.apply_offsets(visit_arrays(lc, "neowise"),
                                   differential.field_offsets(qa.star_tables(stars)),
                                   err_key=qa.ERR_KEY)
    out = qa.target_counts(lc, v)
    assert out["w1"]["saturated_frame_fraction"] == pytest.approx(0.25)
    # the five visits whose W1 frames are all saturated have no W1 point, so they
    # are not two-band visits; the other fifteen are all above the margin
    assert out["n_visits_neowise_w1"] == 15 and out["n_visits_neowise_w2"] == 20
    assert v["mjd"].size == 15
    assert out["w1"]["n_visits_brighter_than_saturation_margin"] == 15
    assert out["w2"]["n_visits_brighter_than_saturation_margin"] == 15


def test_untimely_pairs_measure_only_the_difference():
    rng = np.random.default_rng(3)
    truth = 10.0 - np.linspace(0.0, 1.0, 20)  # a strongly brightening source
    v = {"mjd": VISITS.copy(), "w1_mag": truth + rng.normal(0, 0.01, 20),
         "w2_mag": truth - 0.5 + rng.normal(0, 0.01, 20),
         "w1_err": np.full(20, 0.01), "w2_err": np.full(20, 0.01)}
    # unTimely covers the first 14 visits, 0.2 mag brighter (more host light), 15 mmag noise
    rows = untimely_rows(VISITS[:14] + 0.7, truth[:14] - 0.2 + rng.normal(0, 0.015, 14),
                         truth[:14] - 0.7 + rng.normal(0, 0.015, 14))
    out = qa.untimely_pairs(v, rows)
    assert out["status"] == "ok" and out["n_untimely_epochs"] == 14
    assert out["w1"]["n_pairs"] == 14 and out["w1"]["usable"]
    assert out["w1"]["median_pipeline_minus_untimely_mag"] == pytest.approx(0.2, abs=0.02)
    assert 5.0 < out["w1"]["pair_scatter_mmag"] < 40.0  # ~18 mmag, not the 1 mag trend
    assert out["w1"]["max_detections_within_5arcsec"] == 1
    assert qa.untimely_pairs(v, None) == {"status": "no unTimely detection"}
    assert qa.untimely_pairs(v, []) == {"status": "no unTimely detection"}
    few = qa.untimely_pairs(v, rows[:8])  # four epochs
    assert few["w1"]["n_pairs"] == 4 and not few["w1"]["usable"]


def test_write_frozen_round_trip(tmp_path):
    stars = [star(s) for s in range(1, 9)]
    offsets = differential.field_offsets(qa.star_tables(stars))
    v = differential.apply_offsets(visit_arrays(target(), "neowise"), offsets, err_key=qa.ERR_KEY)
    final = errorfloor.apply_floor(v, {"w1": 0.02, "w2": 0.03})
    path = tmp_path / "frozen" / "t.csv"
    qa.write_frozen(path, final)
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    assert tuple(rows[0].keys()) == qa.FROZEN_COLUMNS and len(rows) == 20
    k = 7
    assert float(rows[k]["w1_mag"]) == pytest.approx(final["w1_mag"][k], abs=1e-5)
    total = np.sqrt(float(rows[k]["w1_err_stat"]) ** 2 + float(rows[k]["w1_offset_err"]) ** 2
                    + float(rows[k]["w1_floor"]) ** 2)
    assert float(rows[k]["w1_err"]) == pytest.approx(total, abs=2e-5)
    assert float(rows[k]["w2_floor"]) == pytest.approx(0.03)
    # writing the same table twice gives the same bytes (the freeze hash is reproducible)
    first = irsa.sha256_file(path)
    qa.write_frozen(path, final)
    assert irsa.sha256_file(path) == first


# --- field files -------------------------------------------------------------------
def _write_field_csv(path, tables: dict[int, dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    header = ["pid"] + list(next(iter(tables.values())).keys())
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(header)
        for pid, table in tables.items():
            for i in range(len(table["mjd"])):
                writer.writerow([pid] + [table[c][i] for c in header[1:]])


def test_load_field_splits_positions(tmp_path, monkeypatch):
    monkeypatch.setattr(sample_io, "S5_DIR", tmp_path)
    tgt = {"name": "bat9999", "ra": RA0, "dec": DEC0}
    stars = [{"name": f"bat9999_s{k:02d}", "ra": RA0, "dec": DEC0, "w1mpro": f"{10 + 0.1 * k}"}
             for k in (1, 2, 3)]
    assert sample_io.field_positions(tgt, stars) == [
        (0, RA0, DEC0), (1, RA0, DEC0), (2, RA0, DEC0), (3, RA0, DEC0)]
    assert not sample_io.field_is_cached("bat9999")
    paths = sample_io.field_paths("bat9999")
    # star 2 has no NEOWISE rows at all; star 3 has no AllWISE rows
    _write_field_csv(paths["neowise"], {
        0: make_neowise_table(VISITS, seed=1), 1: make_neowise_table(VISITS, seed=2),
        3: make_neowise_table(VISITS[:15], seed=3)})
    _write_field_csv(paths["allwise"], {
        0: make_allwise_table(np.array([55340.0, 55530.0])),
        1: make_allwise_table(np.array([55340.0]), seed=4)})
    assert sample_io.field_is_cached("bat9999")
    field = sample_io.load_field(tgt, stars)
    assert field["target"].n_epochs == 22
    assert field["star_names"] == ["bat9999_s01", "bat9999_s03"]
    assert field["star_pids"] == [1, 3] and field["star_w1"] == pytest.approx([10.1, 10.3])
    assert field["stars"][0].n_epochs == 21 and field["stars"][1].n_epochs == 15
    assert field["stars"][1].select("allwise").sum() == 0
