"""Stage 8 regression tests on the committed results (validation of the hysteresis statistic).

The labels must follow from the stored numbers and the criteria fixed in
docs/VALIDATION_8.md. Tests that recompute from the frozen light curves are
skipped when those are not on disk.
"""

import csv
import importlib.util
import json

import numpy as np
import pytest
from scipy.stats import binomtest, wilcoxon

from irecho import config, hyst8

RES = config.RESULTS_DIR
FROZEN_V2 = config.PROJECT_ROOT / "data" / "frozen" / "sample_v2"
pytestmark = pytest.mark.skipif(not (RES / "s8_direct.json").exists(), reason="stage 8 results not present")


def load(name):
    return json.loads((RES / name).read_text())


def rows():
    with open(RES / "s8_direct_objects.csv", newline="") as fh:
        return list(csv.DictReader(fh))


def _script(name):
    spec = importlib.util.spec_from_file_location(name, config.PROJECT_ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_constants_are_the_ones_fixed_in_the_plan():
    plan = (config.PROJECT_ROOT / "docs" / "VALIDATION_8.md").read_text()
    assert hyst8.DUST_FRACTION == 0.5 and "먼지의 몫 (null b) | 0.5" in plan
    assert hyst8.MIN_EACH_DIRECT == 3 and hyst8.MIN_WIDTH_SIGMA == 3.0 and hyst8.MAX_SIGMA_COLOUR == 0.03
    assert "구간마다 공통 구간 안 방문 ≥ 3, 폭 ≥ 3 × W2 오차 중앙값, σ(Δ색) ≤ 0.03 mag" in plan
    assert "A = −0.51, B = −0.479, C = 0.131, D = 0.18" in plan and hyst8.MACLEOD["sf_inf"] == (-0.51, -0.479, 0.131, 0.18)
    assert "A = 2.4, B = 0.17, C = 0.03, D = 0.21" in plan and hyst8.MACLEOD["tau"] == (2.4, 0.17, 0.03, 0.21)
    assert "이미 본 것" in plan and "94개 모두 D > 0" in plan


def test_error_correlation_was_estimated_for_every_object():
    for tag, n in (("v1", 231), ("v2", 210)):
        c = load(f"s8_error_correlation_{tag}.json")
        assert c["n_targets"] == c["n_targets_with_r"] == n and c["min_pairs"] == 6
        r = [o["r"] for o in c["objects"].values()]
        assert c["target_r_median"] == pytest.approx(float(np.median(r)))
        assert 0.2 < c["target_r_median"] < 0.6           # the errors of the two bands are positively correlated


def test_null_labels_follow_the_criteria_and_no_variant_is_biased():
    null = load("s8_null.json")
    assert null["settings"]["n_null"] == 200 and null["settings"]["n_variant"] == 100
    for tag, obs in (("v1", 0.228), ("v2", 0.268)):
        for kind in ("linear", "single_t"):
            for variant, s in null["samples"][tag]["nulls"][kind].items():
                if not isinstance(s, dict) or "label" not in s:
                    continue
                q = s["median_d_p2p5_p16_p50_p84_p97p5_p99"]
                rej = s["wilcoxon_rejection_rate"] or 0.0
                unbiased = abs(q[2]) <= 0.03 and 0.45 <= s["frac_positive_mean"] <= 0.55 and rej <= 0.10
                assert (s["label"] == "치우침 없음") == unbiased, (tag, kind, variant)
                assert s["label"] == "치우침 없음", (tag, kind, variant)
                assert s["n_realisations"] == (200 if variant == "primary" else 100)
            p = null["samples"][tag]["nulls"][kind]["primary"]
            assert p["observed_median_d"] == obs and obs > p["median_d_p2p5_p16_p50_p84_p97p5_p99"][5]
            assert p["observed_vs_null"] == "null로 설명되지 않음"
            assert p["direct"]["usable_under_null"] is True
            # the null reproduces the number of measurable objects of the real analysis within 15%
            n_obs = {"v1": 71, "v2": 94}[tag]
            assert abs(p["n_measurable_p5_p50_p95"][1] - n_obs) <= 0.15 * n_obs


def test_direct_statistics_follow_from_the_object_table():
    direct = load("s8_direct.json")
    for tag in ("v1", "v2"):
        sub = [r for r in rows() if r["sample"] == tag and r["comparable"] == "1"]
        st = direct["samples"][tag]["direct"]
        delta = np.array([float(r["delta"]) for r in sub])
        col = np.array([float(r["d_colour"]) for r in sub])
        assert st["n"] == len(sub) and st["n_delta_positive"] == int(np.sum(delta > 0))
        assert np.all((delta > 0) == (col < 0))                    # more W1 on the rise = bluer on the rise
        assert st["median_d_colour_mag"] == pytest.approx(float(np.median(col)), abs=2e-6)
        assert st["binomial_p_one_sided"] == pytest.approx(binomtest(st["n_delta_positive"], st["n"], 0.5, alternative="greater").pvalue, rel=1e-6)
        for r in sub:
            assert int(r["n_rise_in"]) >= 3 and int(r["n_fade_in"]) >= 3 and float(r["sigma_colour"]) <= 0.03
            assert float(r["lo"]) < float(r["hi"])
        both = [r for r in sub if r["hyst_measurable"] == "1"]
        wd = direct["samples"][tag]["with_d"]
        assert wd["n_both"] == len(both) == {"v1": 71, "v2": 94}[tag]
        assert wd["n_same_sign"] == sum((float(r["delta"]) > 0) == (float(r["d_hyst"]) > 0) for r in both)


def test_direct_labels_follow_the_registered_rule():
    direct, null = load("s8_direct.json"), load("s8_null.json")
    labels = []
    for tag in ("v1", "v2"):
        st = direct["samples"][tag]["direct"]
        nd = [null["samples"][tag]["nulls"][k]["primary"]["direct"]["frac_delta_positive_p2p5_p50_p97p5_p99"] for k in ("linear", "single_t")]
        p99 = max(x[3] for x in nd)
        assert st["null_frac_positive_p99"] == pytest.approx(p99)
        ci = st["median_d_colour_ci95"]
        supports = st["frac_delta_positive"] >= 0.65 and st["binomial_p_one_sided"] < 0.01 and ci[1] < 0 and st["frac_delta_positive"] > p99
        assert (st["label"] == "비대칭 지지") == supports
        labels.append(st["label"])
        ut = direct["samples"][tag]["untimely"]
        assert (ut["label"] == "같은 방향") == (ut["frac_delta_positive"] > 0.5 and ut["binomial_p_one_sided"] < 0.05)
        sp = direct["samples"][tag]["scan_parity"]
        assert (sp["label"] == "scan 방향과 무관") == (sp["relative_change"] <= 0.30)
        assert "added after the plan" in direct["samples"][tag]["not_in_plan_with_linear_time_term"]["note"]
    assert direct["overall"] == ("지지" if all(x == "비대칭 지지" for x in labels) else "사라짐" if all(x == "사라짐" for x in labels) else "약해짐")


@pytest.mark.skipif(not FROZEN_V2.exists(), reason="frozen sample_v2 light curves not present")
def test_direct_offset_of_three_objects_recomputed_from_the_frozen_files():
    r8 = _script("run_stage8")
    objs = {o["name"]: o for o in r8.load_sample("v2")}
    sub = [r for r in rows() if r["sample"] == "v2" and r["comparable"] == "1"]
    for r in (sub[0], sub[len(sub) // 2], sub[-1]):
        o = objs[r["name"]]
        d = hyst8.direct_offset(o["f1"], o["s1"], o["f2"], o["s2"], r=o["r"])
        assert d["comparable"]
        assert d["delta"] == pytest.approx(float(r["delta"]), rel=1e-4)
        assert d["d_colour"] == pytest.approx(float(r["d_colour"]), rel=1e-4, abs=1e-7)
        assert (d["n_rise_in"], d["n_fade_in"]) == (int(r["n_rise_in"]), int(r["n_fade_in"]))


def test_lag_model_verdict_follows_from_the_configurations():
    lag = load("s8_lag_model.json")
    assert set(lag["configs"]) == {f"{d}_scatter{s:.2f}" for d in ("macleod", "macleod_sf_burke_tau") for s in (0.0, 0.15, 0.30)}
    flags = []
    for c in lag["configs"].values():
        assert c["n_realisations"] == 200
        if "rho_p5_p50_p95" in c:
            assert c["outside"] == (lag["rho_observed"] < c["rho_p5_p50_p95"][0])
        flags.append(c.get("outside"))
    expect = ("판정불가 (mock의 측정 가능 수 부족)" if any(f is None for f in flags) else
              "여전히 바깥" if all(flags) else "바깥이 아님" if not any(flags) else "가정에 달림")
    assert lag["verdict"] == expect and lag["rho_observed"] == 0.129


def test_selection_numbers_and_quoted_p_values():
    s = load("s8_selection.json")
    assert s["v1_types"]["by_type_group"] == {"2": 135, "1": 88, "other": 2}
    outside = s["v1_types"]["objects_outside_type_1_and_2"]
    if isinstance(outside, list):                             # the names and types are withheld in the public copy
        assert len(outside) == 2 and all(o["type_group"] == "other" for o in outside)
    q = s["v1_q1_against_q2_q3"]
    assert (q["n_q1"], q["n_brightened_for_q2_q3"], q["n_in_both"]) == (16, 20, 16) and len(q["q2_not_in_q1"]) == 4
    g = s["v2_measurable_with_and_without_lbol"]
    assert g["with_lbol"]["n"] == 70 and g["without_lbol"]["n"] == 24
    with open(RES / "s6b_objects.csv", newline="") as fh:
        d2 = np.array([float(r["d_hyst"]) for r in csv.DictReader(fh) if r["hyst_measurable"] == "1"])
    assert s["p_values"]["v2"]["p"] == pytest.approx(wilcoxon(d2, alternative="greater").pvalue, rel=1e-9)
    assert s["p_values"]["v2"]["p"] == pytest.approx(1.9e-17, rel=0.05) and s["p_values"]["v1"]["p"] == pytest.approx(2.2e-13, rel=0.05)
    assert s["p_values"]["v2"]["n_positive"] == 94 and s["p_values"]["v1"]["n_positive"] == 68


def test_no_local_paths_in_the_committed_stage_8_results():
    for path in list(RES.glob("s8_*.json")) + [RES / "s8_direct_objects.csv"]:
        text = path.read_text()
        assert "/Users/" not in text and "/var/folders" not in text, path.name
