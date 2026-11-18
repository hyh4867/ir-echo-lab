"""Stage 9: the committed results are consistent with the rules of docs/PREDICTIONS_9.md.

Most tests read only committed files. The test that refits objects needs the frozen light curves.
"""

from __future__ import annotations

import csv
import json
import sys

import numpy as np
import pytest

from irecho import config
from irecho import transfer9 as t9

RES = config.RESULTS_DIR
TAGS = ("v1", "v2")
DRIVERS = ("fast", "slow")
FLAG_COLUMNS = ("usable", "tau_at_edge", "well_determined")


def load(name: str):
    path = RES / name
    if not path.exists():
        pytest.skip(f"{name} not present")
    return json.loads(path.read_text())


def object_rows() -> dict[str, list[dict]]:
    path = RES / "s9_objects.csv"
    if not path.exists():
        pytest.skip("s9_objects.csv not present")
    out: dict[str, list[dict]] = {t: [] for t in TAGS}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            row = {"name": r["name"]}
            for k, v in r.items():
                if k in ("sample", "name"):
                    continue
                if k in FLAG_COLUMNS or k.startswith(("usable_", "comparable_")):
                    row[k] = v == "1"
                else:
                    row[k] = float(v) if v != "" else float("nan")
            out[r["sample"]].append(row)
    return out


def inside(v, rng):
    return rng[0] <= v <= rng[1]


def test_identifiability_flags_follow_the_criteria_of_section_3_2():
    ident = load("s9_identifiability.json")
    assert set(ident["configs"]) == {f"{t}_{d}" for t in TAGS for d in DRIVERS}
    flags = {k: [] for k in ("lag_recovered", "t0_vs_t1", "t1_vs_t2", "fixed_vs_t3")}
    for rec in ident["configs"].values():
        assert rec["n_samples"] == 50 and set(rec["truths"]) == set(t9.TRUTHS)
        c = rec["criteria"]
        ok = []
        for name, delay in (("T1_40", 40.0), ("T1_80", 80.0)):
            b = c[f"lag_bias_{name}"]
            good = abs(b["median_error_days"]) <= max(10.0, 0.25 * delay) and b["half_width_16_84_days"] <= 15.0
            assert b["pass"] == good
            ok.append(good)
        assert c["lag_recovered"] == all(ok)
        p = c["t0_vs_t1_power"]
        assert c["t0_vs_t1"] == (max(p["median_tau"]["T1_lyu"], p["mean_ds_01"]["T1_lyu"]) >= 0.8)
        assert c["t1_vs_t2"] == (min(c["t1_vs_t2_power"].values()) >= 0.8)
        p3 = c["fixed_vs_t3_power"]
        assert c["fixed_vs_t3"] == all(max(p3[k][name] for k in p3) >= 0.8 for name in ("T3_1yr", "T3_3yr"))
        for k in flags:
            flags[k].append(c[k])
        # with no lag in the truth the fitted lags scatter around zero
        assert abs(rec["truths"]["T0"]["median_tau"]["median"]) <= 5.0
        assert abs(rec["truths"]["T0"]["median_d_colour_obs"]["median"]) < 0.003
    assert ident["overall"] == {k: all(v) for k, v in flags.items()}


def test_object_table_respects_the_grid_and_carries_no_catalogue_values():
    rows = object_rows()
    assert len(rows["v1"]) == 225 and len(rows["v2"]) == 210
    with open(RES / "s9_objects.csv", newline="") as fh:
        header = next(csv.reader(fh))
    assert not {"type", "z", "log_lbol", "loglbol", "type_group"} & set(header)
    for tag in TAGS:
        for r in rows[tag]:
            if not r["usable"]:
                continue
            assert r["n_used"] >= t9.MIN_VISITS_FIT
            assert -t9.SUPPORT_FUTURE_DAYS <= r["tau_t1"] <= t9.SUPPORT_PAST_DAYS
            assert r["width_t2"] in t9.WIDTHS_DAYS
            assert r["tau_t2"] - r["width_t2"] / 2 >= -t9.SUPPORT_FUTURE_DAYS - 1e-6
            assert r["tau_t2"] + r["width_t2"] / 2 <= t9.SUPPORT_PAST_DAYS + 1e-6
            assert r["tau_lo"] <= r["tau_t1"] <= r["tau_hi"]
            assert r["ds_01"] >= -1e-6 and r["ds_12"] >= -1e-6 and r["ds_02"] == pytest.approx(r["ds_01"] + r["ds_12"], abs=2e-4, rel=1e-4)
            for k in ("t0", "t1", "t2"):
                assert not np.isfinite(r[f"ds_eta_{k}"]) or r[f"ds_eta_{k}"] >= -1e-6
            well = r["ds_01"] >= t9.WELL_DS and r["tau_half_width"] <= t9.WELL_HALF_WIDTH_DAYS and not r["tau_at_edge"]
            assert r["well_determined"] == well


def test_sample_statistics_are_those_of_the_object_table():
    rows = object_rows()
    tr = load("s9_transfer.json")
    for tag in TAGS:
        st = t9.sample_statistics(rows[tag])
        ref = tr["samples"][tag]["statistics"]
        assert st["n_fit"] == ref["n_fit"] and st["n_offset"] == ref["n_offset"] and st["n_well"] == ref["n_well"]
        for k in ("median_tau", "frac_tau_positive", "mean_ds_01", "mean_ds_12", "mean_ds_eta_t2", "median_d_colour_obs",
                  "median_d_colour_pred_t1", "median_d_colour_pred_t2", "median_d_colour_res_t2", "frac_res_bluer_t2"):
            assert st[k] == pytest.approx(ref[k], rel=2e-4, abs=2e-6), (tag, k)


def test_refit_of_frozen_light_curves_reproduces_the_table():
    frozen = config.PROJECT_ROOT / "data" / "frozen" / "sample_v2"
    if not frozen.exists() or not (RES / "s9_objects.csv").exists() or not (RES / "s5b_objects.csv").exists():
        pytest.skip("frozen sample_v2 light curves not present")
    sys.path.insert(0, str(config.PROJECT_ROOT / "scripts"))
    import run_stage9 as r9

    rows = {r["name"]: r for r in object_rows()["v2"]}
    for o in r9.load("v2")[:6]:
        fit = t9.fit_object(o["mjd"], o["f1"], o["s1"], o["f2"], o["s2"], o["r"], prep=o["prep"])
        ref = rows[o["name"]]
        assert fit["usable"] == ref["usable"]
        if fit["usable"]:
            assert fit["tau_t1"] == ref["tau_t1"] and fit["width_t2"] == ref["width_t2"]
            assert fit["ds_01"] == pytest.approx(ref["ds_01"], rel=1e-4, abs=1e-4)
            off = t9.offsets(o["f1"], o["s1"], o["f2"], o["s2"], o["r"], fit)
            if off["usable_obs"]:
                assert off["d_colour_t2"] == pytest.approx(ref["d_colour_t2"], rel=1e-4, abs=1e-7)
            # no visit outside the window that can be predicted without extrapolation
            t = o["mjd"]
            assert np.all(t[fit["idx"]] - t9.SUPPORT_PAST_DAYS >= t[0]) and np.all(t[fit["idx"]] + t9.SUPPORT_FUTURE_DAYS <= t[-1])


def test_verdicts_follow_the_rules_of_sections_4_and_5():
    v = load("s9_verdict.json")
    ident = load("s9_identifiability.json")
    assert v["identifiability"] == ident["overall"]
    labels = {}
    for tag in TAGS:
        rec = v["samples"][tag]
        assert rec["n_null_samples"]["fitted"] == {d: 50 for d in DRIVERS}
        for cmp_ in rec["model_comparison"].values():
            assert cmp_["simpler_model_not_enough"] == all(cmp_["observed"] > cmp_["null"][d]["null_max"] for d in DRIVERS)
        for k in ("t1", "t2"):
            e = rec["explained"][k]
            bands = e["null_fitted"]
            in_m = [inside(e["m_res_mag"], bands[d]["m_res"]) for d in DRIVERS]
            in_f = [inside(e["f_res_bluer"], bands[d]["f_res"]) for d in DRIVERS]
            expect = "설명됨" if all(in_m) and all(in_f) else ("남음" if not any(in_m) else "판정불가")
            assert e["label"] == expect
            assert e["m_res_mag"] == pytest.approx(rec["statistics"][f"median_d_colour_res_{k}"])
            for d in DRIVERS:
                assert e["remaining_mag"][d] == pytest.approx(e["m_res_mag"] - bands[d]["m_res_median"])
        labels[tag] = {k: rec["explained"][k]["label"] for k in ("t1", "t2")}
        assert rec["explained"]["label_sample"] == (labels[tag]["t1"] if labels[tag]["t1"] == labels[tag]["t2"] else "모형에 달림")
    overall = v["explained_overall"]
    s1, s2 = (v["samples"][t]["explained"]["label_sample"] for t in TAGS)
    expect = s1 if s1 == s2 else "표본에 따라 다름"
    assert overall["label_by_rule"] == expect
    if expect == "설명됨" and not ident["overall"]["fixed_vs_t3"]:
        assert "재형성과는 구별하지 못함" in overall["wording_after_b1"]
    else:
        assert overall["wording_after_b1"] == expect
    b3 = v["lag_luminosity"]
    k, rho = b3["observed"]["combined"]["k"], b3["observed"]["combined"]["rho"]
    per = []
    for d in DRIVERS:
        lyu, zero = b3["null"]["T1_lyu"][d], b3["null"]["T0"][d]
        if inside(k, zero["k_null_range"]):
            per.append("지연 없음과 구별 안 됨")
        elif inside(k, lyu["k_null_range"]) and inside(rho, lyu["rho_null_range"]):
            per.append("문헌 관계와 맞음")
        elif not inside(k, lyu["k_null_range"]):
            per.append("문헌 관계와 다름")
        else:
            per.append("판정불가")
    assert b3["label_by_rule"] == (per[0] if len(set(per)) == 1 else "판정불가")
    testable = ident["overall"]["lag_recovered"] and ident["overall"]["t0_vs_t1"]
    assert b3["label_after_b1"] == (b3["label_by_rule"] if testable else "이 자료로 시험할 수 없음")


def test_mock5_flags_follow_section_7():
    m = load("s9_mock5.json")
    assert len(m["configs"]) == 6
    for rec in m["configs"].values():
        assert rec["n_samples"] == 100
        n, d = rec["n_measurable_p5_p50_p95"], rec["median_d_p5_p50_p95"]
        rep = n[0] <= m["observed"]["n_measurable_with_lbol"] <= n[2] and d[0] <= m["observed"]["median_d"] <= d[2]
        assert rec["reproduces_the_sample"] == rep
        if "outside" in rec:
            assert rec["outside"] == (rec["fraction_at_or_below_observed"] < 0.05)


def test_stage9_result_files_hold_no_local_paths():
    paths = sorted(RES.glob("s9_*"))
    if not paths:
        pytest.skip("no stage 9 results")
    for path in paths:
        text = path.read_text()
        assert "/Users/" not in text and "/var/folders" not in text, path.name
