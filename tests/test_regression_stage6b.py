"""Stage 6B regression tests on the frozen sample_v2 light curves and the
committed results.

Skipped when the frozen files are not on disk. Written after the first run of
scripts/run_stage6b.py; they check that the committed tables follow from the
frozen data and the registered rules, not that the science is right.
"""

import csv
import importlib.util
import json

import numpy as np
import pytest

from irecho import breath, config, irsa
from irecho import stage5b as s5

RES = config.RESULTS_DIR
FROZEN_V2 = config.PROJECT_ROOT / "data" / "frozen" / "sample_v2"
FREEZE_V2 = config.PROJECT_ROOT / "data" / "freeze_v2.json"
pytestmark = pytest.mark.skipif(
    not (FROZEN_V2.exists() and (RES / "s6b_objects.csv").exists()),
    reason="frozen sample_v2 light curves or stage 6B results not present")


def _script(name):
    spec = importlib.util.spec_from_file_location(name, config.PROJECT_ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def rows():
    with open(RES / "s6b_objects.csv", newline="") as fh:
        return list(csv.DictReader(fh))


def summary():
    return json.loads((RES / "s6b_summary.json").read_text())


def targets(tag):
    with open(config.TARGETS_DIR / f"sample_{tag}.csv", newline="") as fh:
        return list(csv.DictReader(fh))


def test_frozen_v2_files_still_match_the_freeze():
    freeze = json.loads(FREEZE_V2.read_text())
    assert len(freeze["light_curves"]) == 210 and freeze["code_tree_clean"] is True
    for name, rec in freeze["light_curves"].items():
        assert irsa.sha256_file(FROZEN_V2 / f"{name}.csv") == rec["sha256"], name
    assert summary()["freeze_check"]["n_mismatch"] == 0


def test_sample_v2_is_disjoint_from_sample_v1_and_in_the_registered_redshift_range():
    v1, v2 = targets("v1"), targets("v2")
    assert len(v2) == 210
    assert not {t["name"] for t in v1} & {t["name"] for t in v2}
    assert not {t["counterpart"] for t in v1} & {t["counterpart"] for t in v2}
    z = np.array([float(t["z"]) for t in v2])
    assert z.min() > 0.040 and z.max() <= 0.30
    assert max(float(t["z"]) for t in v1) <= 0.040


def test_object_table_covers_the_sample_and_the_counts_match_the_summary():
    r, s = rows(), summary()
    assert len(r) == s["n_with_light_curve"] == s["n_sample"] == 210
    meas = [x for x in r if x["hyst_measurable"] == "1"]
    assert len(meas) == s["n_measurable"] == s["d_replication"]["n"]
    assert sum(1 for x in meas if x["x"]) == s["n_measurable_with_lbol"] == s["d_vs_luminosity"]["n"]
    st = breath.SETTINGS
    for x in meas:
        assert int(x["n_rising"]) >= st.hyst_min_each and int(x["n_declining"]) >= st.hyst_min_each
        assert float(x["sigma_d"]) <= st.hyst_max_sigma


def test_d_of_three_objects_recomputed_from_the_frozen_files():
    meas = [x for x in rows() if x["hyst_measurable"] == "1"]
    for x in (meas[0], meas[len(meas) // 2], meas[-1]):
        v = s5.read_frozen(x["name"], frozen_dir=FROZEN_V2)
        m = breath.measure(v["w1_flux_jy"], v["w1_flux_err_jy"], v["w2_flux_jy"], v["w2_flux_err_jy"],
                           with_in_place=False)
        assert m["hyst_measurable"]
        assert m["d_hyst"] == pytest.approx(float(x["d_hyst"]), rel=1e-4, abs=1e-6)
        assert m["sigma_d"] == pytest.approx(float(x["sigma_d"]), rel=1e-4, abs=1e-6)


def test_predictor_is_the_registered_one():
    for x in rows():
        if x["x"]:
            assert float(x["x"]) == pytest.approx(
                np.log10(1.0 + float(x["z"])) + 0.5 * float(x["log_lbol"]), abs=2e-4)


def test_summary_statistics_follow_from_the_object_table():
    from scipy.stats import spearmanr
    s = summary()
    meas = [x for x in rows() if x["hyst_measurable"] == "1"]
    d = np.array([float(x["d_hyst"]) for x in meas])
    assert s["d_replication"]["median_d"] == pytest.approx(float(np.median(d)), abs=1e-5)
    use = [x for x in meas if x["x"]]
    # the table stores x to six significant figures, which makes a few ties in the ranks
    rho = spearmanr([float(x["x"]) for x in use], [float(x["d_hyst"]) for x in use]).statistic
    assert s["d_vs_luminosity"]["rho"] == pytest.approx(float(rho), abs=5e-3)


def test_checks_file_repeats_the_registered_numbers_and_is_marked_as_post_hoc():
    c = json.loads((RES / "s6b_checks.json").read_text())
    rel = summary()["d_vs_luminosity"]
    assert "after the registered result was seen" in c["note"]
    assert c["registered"]["n"] == rel["n"]
    assert c["registered"]["rho"] == pytest.approx(rel["rho"], abs=5e-3)
    assert c["registered"]["verdict_by_the_registered_rule"] == rel["verdict"]
    assert c["v1_against_v2"]["n_v1"] == 71 and c["v1_against_v2"]["n_v2"] == summary()["n_measurable"]


def test_verdict_follows_the_registered_rule_and_the_committed_thresholds():
    run = _script("run_stage6b")
    rel = summary()["d_vs_luminosity"]
    final = json.loads((RES / "s6b_power_final.json").read_text())
    assert final["stage"] == "final"
    for key, field in (("model1", "verdict"), ("model1_half", "verdict_if_lag_half"),
                       ("model1_double", "verdict_if_lag_double")):
        q1 = run.interp_threshold(final["models"][key]["by_n"], rel["n"], "rho_p5")
        assert rel["q1"][key] == pytest.approx(q1, abs=1e-9)
        assert rel[field] == run.verdict(rel["rho"], q1, rel["q2"], rel["n"])
    assert run.MIN_OBJECTS == 20


def test_final_power_file_used_the_real_conditions_of_every_object():
    final = json.loads((RES / "s6b_power_final.json").read_text())
    power = _script("run_d_lum_power")
    objects = power.load_objects()
    cond = power.load_conditions("final", [o["name"] for o in objects])
    assert len(cond) == len(objects) == final["n_objects_with_lbol"]
    assert all(final["models"][m]["n_objects_per_realisation"] == len(objects) for m in final["models"])
    # 'final' used the visits of the field; the frozen light curve of a target can have
    # up to three fewer (found by this test; the 'exact' run checks that it does not matter)
    extra = [len(cond[o["name"]][0]) - len(s5.read_frozen(o["name"], frozen_dir=FROZEN_V2)["mjd"])
             for o in objects]
    assert min(extra) == 0 and max(extra) <= 3
    assert sum(1 for e in extra if e > 0) == 24
    exact = power.load_conditions("exact", [o["name"] for o in objects])
    assert all(len(exact[o["name"]][0]) == len(s5.read_frozen(o["name"], frozen_dir=FROZEN_V2)["mjd"])
               for o in objects[:20])


def test_band_lag_estimator_was_applied_only_if_its_validation_passed():
    val = json.loads((RES / "s6b_bandlag_validation_final.json").read_text())
    applied = (RES / "s6b_bandlag.json").exists()
    assert applied == bool(val["pass"]["all"])
