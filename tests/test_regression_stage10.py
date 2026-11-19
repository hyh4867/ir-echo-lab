"""Stage 10: the committed results are consistent with the rules of docs/PREDICTIONS_10.md.

Most tests read only committed files. The tests that open light curves need the frozen ZTF nights
(data/frozen/ztf_v3/) and the frozen WISE light curves; they are skipped without them.
"""

from __future__ import annotations

import csv
import hashlib
import json
import sys

import numpy as np
import pytest

from irecho import config
from irecho import driver10 as d10

RES = config.RESULTS_DIR
DRIVERS = ("fast", "slow")
FROZEN = config.PROJECT_ROOT / "data" / "frozen" / "ztf_v3"
CATALOGUE_KEYS = {"type", "type_group", "z", "log_lbol", "loglbol", "lambda_edd", "log_nh", "log_mbh"}


def load(name: str):
    path = RES / name
    if not path.exists():
        pytest.skip(f"{name} not present")
    return json.loads(path.read_text())


def object_rows() -> list[dict]:
    path = RES / "s10_objects.csv"
    if not path.exists():
        pytest.skip("s10_objects.csv not present")
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def num(v) -> float:
    return float(v) if v != "" else float("nan")


def inside(v, b) -> bool:
    return b["lo"] <= v <= b["hi"]


def test_sample_table_follows_the_criteria_of_section_2():
    s = load("s10_sample.json")
    summary, objs = s["summary"], s["objects"]
    c = summary["criteria"]
    steps = summary["steps"]
    order = ("type1", "ztf_object", "light_curve", "nights", "visits", "optical_variable", "ir_variable")
    assert all(steps[a] >= steps[b] for a, b in zip(order, order[1:]))
    assert steps["type1"] == len(objs) == 234 and steps["ir_variable"] == summary["n_usable"]
    assert summary["stop"] == (summary["n_usable"] < 30)
    n_usable = 0
    for rec in objs.values():
        passing = [b for b in ("g", "r") if rec["bands"][b].get("pass_optical")]
        assert rec["driver_band"] == (passing[0] if passing else None)
        for b in rec["bands"].values():
            if b.get("pass_optical"):
                assert b["n_nights"] >= c["min_nights"] and b["n_visits_usable"] >= c["min_visits"]
                assert b["var_n_bins"] >= c["min_bins"] and b["var_f30"] >= c["min_f30"] and b["var_f30_over_sigma"] >= c["min_f30_sigma"]
                assert b["pass"] == (b["p_const_w1"] < c["p_ir_variable"] and b["p_const_w2"] < c["p_ir_variable"])
                assert c["flag_rho"] == 0.3 and (b["flag_airmass"] is None or b["flag_airmass"] == (
                    abs(b["rho_airmass"]) >= c["flag_rho"] and b["p_airmass"] < c["flag_p"]))
        if rec["usable"]:
            n_usable += 1
            assert rec["bands"][rec["driver_band"]]["pass"]
    assert n_usable == summary["n_usable"]
    assert (c["support_past_days"], c["support_future_days"], c["max_gap_days"]) == (
        d10.SUPPORT_PAST_DAYS, d10.SUPPORT_FUTURE_DAYS, d10.MAX_GAP_DAYS)


def test_result_files_carry_no_catalogue_values_and_no_local_paths():
    paths = sorted(RES.glob("s10_*"))
    if not paths:
        pytest.skip("no stage 10 results")
    for path in paths:
        text = path.read_text()
        assert "/Users/" not in text and "/var/folders" not in text and "/private/tmp" not in text, path.name
    with open(RES / "s10_objects.csv", newline="") as fh:
        header = next(csv.reader(fh))
    assert not CATALOGUE_KEYS & set(header)

    def walk(node):
        if isinstance(node, dict):
            assert not ("name" in node and CATALOGUE_KEYS & set(node)), sorted(node)
            for v in node.values():
                walk(v)
        elif isinstance(node, list):
            for v in node:
                walk(v)

    for path in paths:
        if path.suffix == ".json":
            walk(json.loads(path.read_text()))
    sample = load("s10_sample.json")
    for rec in sample["objects"].values():
        assert not CATALOGUE_KEYS & set(rec)


def test_freeze_lists_every_light_curve_and_matches_the_files_when_present():
    path = config.PROJECT_ROOT / "data" / "freeze_v3.json"
    if not path.exists():
        pytest.skip("freeze_v3.json not present")
    freeze = json.loads(path.read_text())
    sample = load("s10_sample.json")["objects"]
    assert freeze["n_objects"] == len(sample) and freeze["n_usable"] == sum(r["usable"] for r in sample.values())
    checked = 0
    for name, rec in freeze["light_curves"].items():
        assert rec["usable"] == sample[name]["usable"] and rec["driver_band"] == sample[name]["driver_band"]
        for band, b in rec["bands"].items():
            assert len(b["raw_sha256"]) == 64
            f = FROZEN / f"{name}_{band}.csv"
            if b["sha256"] and f.exists():
                assert hashlib.sha256(f.read_bytes()).hexdigest() == b["sha256"]
                checked += 1
    if not FROZEN.exists():
        pytest.skip("frozen ZTF nights not present (hashes not compared)")
    assert checked > 0


def test_identifiability_flags_follow_the_criteria_of_section_5_1():
    ident = load("s10_identifiability.json")
    assert ident["n_realisations"] == 50 and set(ident["drivers"]) == set(DRIVERS)
    names = ("lag_recovered", "dtau_recovered", "detected", "alpha_recovered", "alpha_power_ok_b20", "width_recovered",
             "fixed_vs_history", "relations_separated")
    per = {k: [] for k in names}
    for rec in ident["drivers"].values():
        c, tab = rec["criteria"], rec["truths"]
        assert len(tab) == 21
        ok = []
        for name in ("A_100", "A_200", "A_lyu"):
            for b in (1, 2):
                x = c[f"lag_bias_{name}_w{b}"]
                good = (x["median_error_days"] is not None and abs(x["median_error_days"]) <= max(10.0, 0.10 * x["injected_median_days"])
                        and x["half_width_16_84_days"] <= 15.0)
                assert x["pass"] == good
                ok.append(good)
        assert c["lag_recovered"] == all(ok)
        ok = []
        for name in ("A_50", "A_100", "A_200"):
            x = c[f"dtau_bias_{name}"]
            good = x["median_error_days"] is not None and abs(x["median_error_days"]) <= 8.0 and x["half_width_16_84_days"] <= 10.0
            assert x["pass"] == good
            ok.append(good)
        assert c["dtau_recovered"] == (all(ok) and c["dtau_power_lyu_vs_a0"] >= 0.8)
        assert c["detected"] == (c["detection_power"] >= 0.8)
        bias = [c[f"alpha_bias_{n}"]["pass"] for n in ("A_lyu", "B_20", "B_40")]
        for n in ("A_lyu", "B_20", "B_40"):
            x = c[f"alpha_bias_{n}"]
            assert x["pass"] == (x["median_alpha"] is not None and abs(x["median_alpha"] - x["injected"]) <= 0.1)
        assert c["alpha_power_ok_b40"] == (max(c["alpha_power_B_40"].values()) >= 0.8)
        assert c["alpha_power_ok_b20"] == (max(c["alpha_power_B_20"].values()) >= 0.8)
        assert c["alpha_recovered"] == (all(bias) and c["alpha_power_ok_b40"])
        assert c["width_recovered"] == all(max(v.values()) >= 0.8 for v in c["width_power"].values())
        assert c["fixed_vs_history"] == all(max(c["fixed_vs_history_power"][n].values()) >= 0.8 for n in ("C_1yr", "C_3yr"))
        assert c["relations_separated"] == (max(c["relations_power"].values()) >= 0.8)
        for k in names:
            per[k].append(c[k])
        # with the infrared made from an independent driver few objects pass the detection in both bands
        assert tab["X0"]["all"]["frac_detected"]["median"] <= 0.05
    assert ident["overall"] == {k: all(v) for k, v in per.items()}


def test_object_table_respects_the_grid_and_the_detection_rule():
    rows = object_rows()
    sample = load("s10_sample.json")
    assert len(rows) == sample["summary"]["n_usable"]
    n_d = 0
    for r in rows:
        assert int(r["n_used"]) >= d10.MIN_VISITS and int(r["n_nights"]) >= d10.MIN_NIGHTS
        assert sample["objects"][r["name"]]["driver_band"] == r["driver_band"]
        edge = False
        for b in ("w1", "w2"):
            tau = num(r[f"tau_{b}"])
            assert -65.0 <= tau <= 465.0 and tau % 5.0 == 0.0
            assert num(r[f"tau_q16_{b}"]) <= num(r[f"tau_q50_{b}"]) <= num(r[f"tau_q84_{b}"])
            assert num(r[f"tau_ma_{b}"]) >= d10.MA_MIN_TAU_DAYS
            assert num(r[f"ds_1a_{b}"]) > -50 and -1e-6 <= num(r[f"ds_12_{b}"])
            assert not np.isfinite(num(r[f"ds_1e_{b}"])) or num(r[f"ds_1e_{b}"]) >= -1e-6
            assert (r[f"detected_{b}"] == "1") == (int(r["n_swaps"]) >= 10 and num(r[f"r_w_{b}"]) > num(r[f"null_rw_max_{b}"]))
            edge |= r[f"tau_at_edge_{b}"] == "1"
            assert (r[f"tau_at_edge_{b}"] == "1") == (tau in (-65.0, 465.0))
        assert num(r["dtau"]) == pytest.approx(num(r["tau_w2"]) - num(r["tau_w1"]))
        assert num(r["dtau_q025"]) <= num(r["dtau_q16"]) <= num(r["dtau_q50"]) <= num(r["dtau_q84"]) <= num(r["dtau_q975"])
        in_d = r["detected_w1"] == "1" and r["detected_w2"] == "1" and not edge
        assert (r["in_d"] == "1") == in_d
        n_d += in_d
    fit = load("s10_fit.json")
    assert fit["n_objects"] == len(rows) and fit["n_d"] == n_d == fit["statistics"]["d"]["n"]
    d = [r for r in rows if r["in_d"] == "1"]
    if len(d) >= 3:
        st = fit["statistics"]["d"]
        assert st["med_tau1"] == pytest.approx(float(np.median([num(r["tau_w1"]) for r in d])))
        assert st["med_tau2"] == pytest.approx(float(np.median([num(r["tau_w2"]) for r in d])))
        assert st["med_dtau"] == pytest.approx(float(np.median([num(r["dtau"]) for r in d])))


def test_verdicts_follow_the_rules_of_sections_5_2_and_6():
    v = load("s10_verdict.json")
    ident = load("s10_identifiability.json")
    flags = ident["overall"]
    assert v["identifiability"] == flags
    small = v["n_d"] < 20
    not_testable, too_small = "이 자료로 시험할 수 없음", "판정불가(표본 부족)"

    def final(rule, needs):
        if not all(flags[k] for k in needs):
            return not_testable
        return too_small if small else rule

    q1 = v["q1"]
    ins = [q1["q"] is not None and "lo" in q1["null_A_s9"][d] and inside(q1["q"], q1["null_A_s9"][d]) for d in DRIVERS]
    known = q1["q"] is not None and all("lo" in q1["null_A_s9"][d] for d in DRIVERS)
    rule = "맞음" if all(ins) else ("다름" if known and not any(ins) else "판정불가")
    assert q1["label_by_rule"] == rule and q1["label"] == final(rule, ("dtau_recovered", "detected"))
    for rel in ("lyu", "mandal"):
        rec = v["q2"][rel]
        keys = (f"k_{rel}_w1", f"k_{rel}_w2", "rho_w1", "rho_w2")
        all_in = all(x is True for k in keys for x in rec["inside"][k])
        k_out = any(all(x is False for x in rec["inside"][k]) for k in keys[:2])
        rule = "맞음" if all_in else ("다름" if k_out else "판정불가")
        assert rec["label_by_rule"] == rule and rec["label"] == final(rule, ("lag_recovered", "detected"))
    q3 = v["q3"]
    assert q3["label"] == final(q3["label_by_rule"], ("alpha_recovered", "detected"))
    assert q3["explains_1p4"] == final(q3["explains_1p4_by_rule"], ("alpha_recovered", "dtau_recovered", "detected"))
    obs, nul = q3["observed"], q3["null_fitted"]
    keys = ("med_alpha_w1", "mean_ds_1a_w1")
    if all("hi" in nul[d][k] for d in DRIVERS for k in keys) and all(obs[k] is not None for k in keys):
        ab = all(obs[k] > nul[d][k]["hi"] for d in DRIVERS for k in keys)
        ins3 = all(inside(obs[k], nul[d][k]) for d in DRIVERS for k in keys)
        assert q3["label_by_rule"] == ("지연 없는 성분 있음" if ab else ("근거 없음" if ins3 else "판정불가"))
    q4 = v["q4"]
    assert q4["attempted"] == flags["fixed_vs_history"]
    if not flags["fixed_vs_history"]:
        assert q4["label"] == "이 자료로 가를 수 없음"
        assert v["offsets"]["wording"] == "inter-band lag 모형과 모순되지 않음 (재형성과는 구별하지 못함)"
    for k, rec in v["predictions"].items():
        if "predicted" in rec and rec["observed"] is not None:
            assert rec["inside"] == (rec["predicted"][0] <= rec["observed"] <= rec["predicted"][1]), k


def test_null_with_the_fitted_lags_has_fifty_samples_per_driver():
    nf = load("s10_null_fitted.json")
    assert nf["n_realisations"] == 50 and set(nf["drivers"]) == set(DRIVERS)
    for rec in nf["drivers"].values():
        assert rec["all"]["med_tau1"]["n"] == 50


def test_refit_of_frozen_light_curves_reproduces_the_table():
    rows = object_rows()
    if not FROZEN.exists() or not (config.PROJECT_ROOT / "data" / "frozen" / "sample_v2").exists():
        pytest.skip("frozen light curves not present")
    sys.path.insert(0, str(config.PROJECT_ROOT / "scripts"))
    import run_stage10 as r10

    objs = {o["name"]: o for o in r10.load()}
    for r in rows[:3]:
        o = objs[r["name"]]
        prep = d10.prepare(o["opt_t"], o["opt_s"], o["hyper"], o["mjd"])
        assert prep["idx"].size == int(r["n_used"])
        # no visit whose window leaves the optical light curve
        t = o["mjd"][prep["idx"]]
        assert np.all(t - d10.SUPPORT_PAST_DAYS >= o["opt_t"][0]) and np.all(t + d10.SUPPORT_FUTURE_DAYS <= o["opt_t"][-1])
        for b, f, s in (("w1", o["f1"], o["s1"]), ("w2", o["f2"], o["s2"])):
            fit = d10.fit_band(prep, o["opt_y"], f, s, full=False)
            assert fit["tau"] == num(r[f"tau_{b}"])
            assert fit["r_w"] == pytest.approx(num(r[f"r_w_{b}"]), rel=1e-4, abs=1e-5)
