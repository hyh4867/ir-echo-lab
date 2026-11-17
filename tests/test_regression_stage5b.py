"""Stage 5b regression tests on the frozen light curves and committed results.

Skipped when the frozen files are not on disk. These tests were written after
the first run of scripts/run_stage5b.py; they check that the committed tables
follow from the frozen data and the committed rules, not that the science is
right.
"""

import csv
import json

import numpy as np
import pytest

from irecho import breath, config, irsa, sample_io
from irecho import stage5b as s5

RES = config.RESULTS_DIR
FREEZE_PATH = config.PROJECT_ROOT / "data" / "freeze_v1.json"
pytestmark = pytest.mark.skipif(
    not (s5.FROZEN_DIR.exists() and (RES / "s5b_objects.csv").exists()),
    reason="frozen light curves or stage 5b results not present")


def rows(name="s5b_objects.csv"):
    with open(RES / name, newline="") as fh:
        return list(csv.DictReader(fh))


def summary():
    return json.loads((RES / "s5b_summary.json").read_text())


def test_frozen_files_still_match_the_freeze():
    freeze = json.loads(FREEZE_PATH.read_text())
    assert len(freeze["light_curves"]) == 231
    for name, rec in freeze["light_curves"].items():
        assert irsa.sha256_file(s5.FROZEN_DIR / f"{name}.csv") == rec["sha256"], name
    assert summary()["freeze_check"]["n_mismatch"] == 0


def test_sample_split_follows_the_predictions_document():
    r = rows()
    assert len(r) == 231
    separate = sorted(x["counterpart"] for x in r if x["sample"] == "separate")
    assert separate == ["NGC 2273", "NGC 3081", "NGC 3393", "NGC 3786", "NGC 4941", "NGC 5728"]
    s = summary()["floor"]
    assert s["main_all"]["n"] == 225
    assert s["main_clean"]["n"] == sum(1 for x in r if x["sample"] == "main" and not x["qa_flags"])


def test_settings_are_the_registered_ones():
    st = breath.SETTINGS
    assert (st.n_ref, st.min_snr, st.min_used, st.sigma_scale, st.separation) == (3, 3.0, 6, 1.0, 3.0)
    assert (st.hyst_min_each, st.hyst_max_sigma) == (4, 0.10)
    val = json.loads((RES / "s5b_validation.json").read_text())
    assert val["settings"] == {k: getattr(st, k) for k in st.__dataclass_fields__}
    assert val["d_min"] == pytest.approx(0.0245) == summary()["d_min"]
    assert s5.LAMBDA_EDD_THRESHOLD == 0.02 and s5.MIN_PER_GROUP == 5 and s5.MIN_HYSTERESIS == 8


@pytest.mark.parametrize("with_floor, table", [(True, "s5b_objects.csv"),
                                               (False, "s5b_objects_nofloor.csv")])
def test_object_rows_rebuild_from_the_frozen_files(with_floor, table):
    r = {x["name"]: x for x in rows(table)}
    picks = [x for x in r.values() if x["brightened"] == "1"][:4] + list(r.values())[:3]
    for x in picks:
        v = s5.read_frozen(x["name"], with_floor)
        d = s5.run_detectors(v)
        assert int(d["q1"]) == int(x["q1"]) and int(d["b_flagged"]) == int(x["b_flagged"])
        assert int(d["a_flagged"]) == int(x["a_flagged"])
        if d["brightened"]:
            c = s5.colour_change(v, d["start_index"])
            if x["delta_color"] == "":  # fewer than three visits after the baseline
                assert np.isnan(c["delta_color"])
            else:
                assert c["delta_color"] == pytest.approx(float(x["delta_color"]), rel=1e-4,
                                                         abs=1e-6)
        m = breath.measure(v["w1_flux_jy"], v["w1_flux_err_jy"], v["w2_flux_jy"], v["w2_flux_err_jy"])
        assert m["status"] == x["sl_status"]
        if m["status"] == "ok":
            assert m["beta"] == pytest.approx(float(x["sl_beta"]), rel=1e-3, abs=1e-5)
            assert m["sigma_beta"] == pytest.approx(float(x["sl_sigma_beta"]), rel=1e-3)


def test_floor_only_changes_the_errors():
    a = s5.read_frozen("bat0001", True)
    b = s5.read_frozen("bat0001", False)
    assert np.array_equal(a["w1_mag"], b["w1_mag"]) and np.array_equal(a["w2_flux_jy"], b["w2_flux_jy"])
    assert np.all(b["w1_err"] <= a["w1_err"] + 1e-9)


def test_summary_counts_follow_from_the_object_table():
    r = [x for x in rows() if x["sample"] == "main"]
    s = summary()["floor"]["main_all"]
    q1 = [x["name"] for x in r if x["q1"] == "1"]
    assert s["q1"]["n_slow_rise"] == len(q1) and s["q1"]["objects"] == q1
    lo, hi = s5.wilson(len(q1), len(r))
    assert s["q1"]["ci95"] == pytest.approx([lo, hi]) and s["q1"]["verdict"] == s5.verdict_q1(lo, hi)
    br = [x for x in r if x["brightened"] == "1" and x["delta_color"]]
    red = sum(1 for x in br if float(x["delta_color"]) > 0)
    assert (s["q2"]["n_brightened"], s["q2"]["n_redder"]) == (len(br), red)
    assert s["q2"]["verdict"] == s5.verdict_majority(red, len(br))[0]
    use = [x for x in r if x["brightened"] == "1" and x["ztf_usable"] == "1"]
    acc = sum(1 for x in use if x["ztf_verdict"] == "accompanies")
    assert (s["q3"]["n_usable"], s["q3"]["n_accompanies"]) == (len(use), acc)
    assert s["q3"]["verdict"] == s5.verdict_majority(acc, len(use))[0]
    meas = [x for x in r if x["sl_measurable"] == "1"]
    assert s["h_breath"]["n_measurable"] == len(meas)
    hm = [float(x["sl_d_hyst"]) for x in r if x["sl_hyst_measurable"] == "1"]
    assert s["hysteresis"]["n"] == len(hm)
    assert s["hysteresis"]["median_d"] == pytest.approx(float(np.median(hm)), rel=1e-4)


def test_final_verdict_is_split_when_the_two_samples_disagree():
    for mode in ("floor", "nofloor"):
        s = summary()[mode]
        for key, v in s["verdicts"].items():
            assert v["main_all"] == s["main_all"][key]["verdict"]
            assert v["main_clean"] == s["main_clean"][key]["verdict"]
            same = v["main_all"].split(" (")[0] == v["main_clean"].split(" (")[0]
            assert (v["final"] == "갈림") == (not same)
            if same:
                assert v["final"].split(" (")[0] == v["main_all"].split(" (")[0]


def test_measurable_objects_satisfy_the_criterion():
    for x in rows():
        if x["sl_measurable"] == "1":
            assert int(float(x["sl_n_used"])) >= 6
            assert 3.0 * float(x["sl_sigma_beta"]) <= float(x["sl_beta_in_place_min"]) * (1 + 1e-4)
        if x["sl_hyst_measurable"] == "1":
            assert min(int(float(x["sl_n_rising"])), int(float(x["sl_n_declining"]))) >= 4
            assert float(x["sl_sigma_d"]) <= 0.10 * (1 + 1e-4)


def test_props_table_threshold_column():
    with open(config.TARGETS_DIR / "sample_v1_props.csv", newline="") as fh:
        props = list(csv.DictReader(fh))
    assert len(props) == 231 == len(sample_io.load_sample())
    for p in props:
        if p["dr2_logedd"]:
            assert int(p["above_threshold"]) == int(10 ** float(p["dr2_logedd"]) >= 0.02)
        else:
            assert p["above_threshold"] == ""


def test_q1_completeness_table_is_complete():
    path = RES / "s5b_q1_completeness.csv"
    if not path.exists():
        pytest.skip("Q1 completeness not yet computed")
    with open(path, newline="") as fh:
        t = list(csv.DictReader(fh))
    assert len(t) == 4 * 7 * 3
    assert {float(x["floor_mmag"]) for x in t} == {0.0, 10.0, 20.0, 50.0}
    for x in t:
        assert 0.0 <= float(x["ci95_low"]) <= float(x["completeness_q1"]) <= float(x["ci95_high"]) <= 1.0
        assert int(x["n_q1"]) <= int(x["n_b"]) <= int(x["n_trials"])
