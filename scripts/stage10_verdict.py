"""Stage 10: the verdicts of docs/PREDICTIONS_10.md (sections 5.2 and 6) from the result files.

    python scripts/stage10_verdict.py

Reads results/s10_identifiability.json, s10_fit.json, s10_null_fitted.json and s10_objects.csv
and writes results/s10_verdict.json. Post hoc comparisons are marked as such and do not enter a label.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import run_stage10 as r10  # noqa: E402
from irecho import driver10 as d10  # noqa: E402

RES = r10.RES
DRIVERS = tuple(r10.DRIVERS)
STAGE9_BIAS = 0.87            # stage 9, B3: fitted lag / lag put in (k of the T1_lyu null: 0.79-0.96)
MIN_D = 20
NOT_TESTABLE = "이 자료로 시험할 수 없음"
TOO_SMALL = "판정불가(표본 부족)"
# Mandal et al. 2024, Table 5 (BCES, ICCF), without the accretion-disc correction ('NON ADC'): post hoc
MANDAL_NON_ADC = {"w1": (0.38, 2.61), "w2": (0.35, 2.72)}
PREDICTED = {"med_dtau": (30.0, 80.0), "q": (0.7, 1.6), "k_lyu_w1": (0.6, 1.3), "k_lyu_w2": (0.7, 1.5),
             "k_mandal_w1": (0.4, 1.0), "k_mandal_w2": (0.4, 1.0), "slope_w1": (0.2, 0.6), "slope_w2": (0.2, 0.6),
             "med_alpha_w1": (0.05, 0.35), "frac_z_alpha_w1": (0.10, 0.40), "k_alpha": (0.8, 1.3)}


def inside(v, b: dict) -> bool | None:
    if v is None or not np.isfinite(v) or "lo" not in b:
        return None
    return bool(b["lo"] <= v <= b["hi"])


def above(v, b: dict) -> bool | None:
    if v is None or not np.isfinite(v) or "hi" not in b:
        return None
    return bool(v > b["hi"])


def main() -> None:
    ident = json.loads((RES / "s10_identifiability.json").read_text())
    fit = json.loads((RES / "s10_fit.json").read_text())
    fitted = json.loads((RES / "s10_null_fitted.json").read_text())["drivers"]
    flags = ident["overall"]
    real = fit["statistics"]["d"]
    n_d = fit["n_d"]
    null = lambda driver, truth, key, s="d": ident["drivers"][driver]["truths"][truth][s].get(key, {})  # noqa: E731
    out = {"n_objects": fit["n_objects"], "n_d": n_d, "identifiability": flags, "d_too_small": n_d < MIN_D}

    def finish(label: str, needs: tuple[str, ...]) -> str:
        if not all(flags[k] for k in needs):
            return NOT_TESTABLE
        return TOO_SMALL if n_d < MIN_D else label

    # ---- Q1
    q = real.get("q")
    bands = {d: null(d, "A_s9", "q") for d in DRIVERS}
    ins = [inside(q, bands[d]) for d in DRIVERS]
    rule = "맞음" if all(x is True for x in ins) else ("다름" if all(x is False for x in ins) else "판정불가")
    out["q1"] = {"n": real.get("n_q1"), "q": q, "d_days": real.get("d"), "rho": real.get("rho_q1"),
                 "med_dtau_days": real.get("med_dtau"), "med_tau9_days": real.get("med_tau9"),
                 "q_bootstrap_ci95": fit["bootstrap_q"]["d"]["q_ci95"], "med_dtau_bootstrap_ci95": fit["bootstrap_q"]["d"]["med_dtau_ci95"],
                 "q_with_stage9_bias_corrected": q * STAGE9_BIAS if q is not None else None,
                 "null_A_s9": bands, "null_A0": {d: null(d, "A0", "med_dtau") for d in DRIVERS},
                 "q_over_null_median": {d: q / bands[d]["median"] if q is not None and bands[d].get("median") else None for d in DRIVERS},
                 "all_objects": {k: fit["statistics"]["all"].get(k) for k in ("n_q1", "q", "d", "rho_q1", "med_dtau", "med_tau9")},
                 "label_by_rule": rule, "label": finish(rule, ("dtau_recovered", "detected"))}
    # ---- Q2
    out["q2"] = {"relations_separated": flags["relations_separated"]}
    for rel, truth in (("lyu", "A_lyu"), ("mandal", "A_mandal")):
        keys = (f"k_{rel}_w1", f"k_{rel}_w2", "rho_w1", "rho_w2")
        ins = {k: [inside(real.get(k), null(d, truth, k)) for d in DRIVERS] for k in keys}
        k_out = [all(x is False for x in ins[k]) for k in keys[:2]]
        if all(all(x is True for x in v) for v in ins.values()):
            rule = "맞음"
        elif any(k_out):
            rule = "다름"
        else:
            rule = "판정불가"
        out["q2"][rel] = {"observed": {k: real.get(k) for k in keys + ("slope_w1", "slope_w2", "n_lbol")},
                          "null": {d: {k: null(d, truth, k) for k in keys + ("slope_w1", "slope_w2")} for d in DRIVERS},
                          "corrected_multiple": {k: {d: real[k] / null(d, truth, k)["median"] if real.get(k) is not None and null(d, truth, k).get("median") else None
                                                     for d in DRIVERS} for k in keys[:2]},
                          "inside": ins, "label_by_rule": rule, "label": finish(rule, ("lag_recovered", "detected"))}
    # ---- Q3
    keys = ("med_alpha_w1", "mean_ds_1a_w1")
    ab = {k: [above(real.get(k), fitted[d]["d"].get(k, {})) for d in DRIVERS] for k in keys}
    ins = {k: [inside(real.get(k), fitted[d]["d"].get(k, {})) for d in DRIVERS] for k in keys}
    if all(all(x is True for x in v) for v in ab.values()):
        rule = "지연 없는 성분 있음"
    elif all(all(x is True for x in v) for v in ins.values()):
        rule = "근거 없음"
    else:
        rule = "판정불가"
    k10, ka = real.get("k10"), real.get("k_alpha")
    n10 = {d: null(d, "A_lyu", "k10") for d in DRIVERS}
    if rule == "지연 없는 성분 있음" and all(above(k10, n10[d]) for d in DRIVERS) and all(inside(ka, n10[d]) for d in DRIVERS):
        rule14 = "설명함"
    elif all(above(ka, n10[d]) for d in DRIVERS):
        rule14 = "설명하지 못함"
    else:
        rule14 = "판정불가"
    note14 = None
    if all(inside(k10, n10[d]) for d in DRIVERS):
        note14 = "광학 driver로 잰 Δτ는 Lyu의 차와 구별되지 않음"
    out["q3"] = {"observed": {k: real.get(k) for k in ("med_alpha_w1", "mean_ds_1a_w1", "frac_z_alpha_w1", "med_alpha_w2",
                                                        "mean_ds_1a_w2", "frac_z_alpha_w2", "med_tau_ma_w1", "med_dtau_alpha",
                                                        "med_dtau_alpha_both", "q_alpha", "k9", "k10", "k_alpha")},
                 "null_fitted": {d: {k: fitted[d]["d"].get(k, {}) for k in ("med_alpha_w1", "mean_ds_1a_w1", "frac_z_alpha_w1",
                                                                          "med_alpha_w2", "mean_ds_1a_w2", "frac_z_alpha_w2")}
                                 for d in DRIVERS},
                 "null_A_lyu_k10": n10, "null_A_lyu_k_alpha": {d: null(d, "A_lyu", "k_alpha") for d in DRIVERS},
                 "label_by_rule": rule, "label": finish(rule, ("alpha_recovered", "detected")),
                 "explains_1p4_by_rule": rule14, "explains_1p4_note": note14,
                 "explains_1p4": finish(rule14, ("alpha_recovered", "dtau_recovered", "detected")),
                 "alpha_detectable": "α ≥ 0.4만 검출 가능" if flags["alpha_recovered"] and not flags["alpha_power_ok_b20"] else None}
    # ---- Q4
    q4 = {"attempted": bool(flags["fixed_vs_history"])}
    obs4 = {k: real.get(k) for k in ("m_res", "mean_ds_1e_w1", "mean_ds_1e_w2", "f_res_bluer")}
    n4 = {d: {k: fitted[d]["d"].get(k, {}) for k in obs4} for d in DRIVERS}
    if flags["fixed_vs_history"]:
        ok_all = all(inside(obs4[k], n4[d][k]) for d in DRIVERS for k in ("m_res", "mean_ds_1e_w1", "mean_ds_1e_w2"))
        out_m = all(inside(obs4["m_res"], n4[d]["m_res"]) is False for d in DRIVERS)
        q4["label"] = "고정 kernel로 설명됨" if ok_all else ("남음" if out_m else "판정불가")
    else:
        q4["label"] = "이 자료로 가를 수 없음"
    q4.update(observed_reported_only=obs4, null_fitted=n4)
    out["q4"] = q4
    # ---- the direct comparison of stage 8 predicted by the two kernels
    out["offsets"] = {"observed": {k: real.get(k) for k in ("n_offset", "m_obs", "m_pred", "m_res", "ratio_pred_obs", "f_res_bluer", "frac_obs_bluer")},
                      "null_fitted": {d: {k: fitted[d]["d"].get(k, {}) for k in ("m_obs", "m_pred", "m_res", "ratio_pred_obs", "f_res_bluer")}
                                      for d in DRIVERS},
                      "all_objects": {k: fit["statistics"]["all"].get(k) for k in ("n_offset", "m_obs", "m_pred", "m_res", "ratio_pred_obs")},
                      "wording": "inter-band lag 모형과 모순되지 않음 (재형성과는 구별하지 못함)" if not flags["fixed_vs_history"] else None}
    # ---- the registered numbers against what came out
    pred = {}
    for k, (lo, hi) in PREDICTED.items():
        v = real.get(k)
        pred[k] = {"predicted": [lo, hi], "observed": v, "inside": bool(v is not None and np.isfinite(v) and lo <= v <= hi)}
    pred["rho_positive"] = {"observed": [real.get("rho_w1"), real.get("rho_w2")],
                            "inside": bool((real.get("rho_w1") or -1) > 0 and (real.get("rho_w2") or -1) > 0)}
    pred["criterion_6_not_passed"] = {"observed": not flags["fixed_vs_history"], "inside": not flags["fixed_vs_history"]}
    out["predictions"] = pred
    # ---- post hoc: the Mandal relation without the accretion-disc correction; a shift of the luminosity scale
    objs = {o["name"]: o for o in r10.load()}
    rows = []
    with open(RES / "s10_objects.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            o = objs[r["name"]]
            if r["in_d"] == "1" and o["log_lbol"] is not None:
                rows.append((o["log_lbol"], o["z"], float(r["tau_w1"]), float(r["tau_w2"])))
    post = {"n": len(rows)}
    if rows:
        logl, z, t1, t2 = (np.array(v) for v in zip(*rows))
        x = logl - d10.MANDAL["log_l0"]
        for b, tau in (("w1", t1), ("w2", t2)):
            a, c = MANDAL_NON_ADC[b]
            post[f"k_mandal_non_adc_{b}"] = float(np.median(tau / (1.0 + z) / 10.0 ** (c + a * x)))
            # Mandal et al. divide by (1 + z)^0.38 (section 4.3) instead of (1 + z)
            post[f"k_mandal_{b}_with_their_redshift_factor"] = float(np.median(
                tau / (1.0 + z) ** 0.38 / 10.0 ** (d10.MANDAL[b][1] + d10.MANDAL[b][0] * x)))
        for shift in (-0.3, 0.3):
            lyu = np.array([d10.lyu_lags_days(v + shift) for v in logl])
            man = np.array([d10.mandal_lags_days(v + shift) for v in logl])
            post[f"shift_{shift:+.1f}_dex"] = {"k_lyu_w1": float(np.median(t1 / (1 + z) / lyu[:, 0])), "k_lyu_w2": float(np.median(t2 / (1 + z) / lyu[:, 1])),
                                               "k_mandal_w1": float(np.median(t1 / (1 + z) / man[:, 0])), "k_mandal_w2": float(np.median(t2 / (1 + z) / man[:, 1]))}
    out["post_hoc"] = post
    (RES / "s10_verdict.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    brief = {"n_d": n_d, "q1": out["q1"]["label"], "q2_lyu": out["q2"]["lyu"]["label"], "q2_mandal": out["q2"]["mandal"]["label"],
             "q3": out["q3"]["label"], "q3_1p4": out["q3"]["explains_1p4"], "q4": out["q4"]["label"]}
    print(json.dumps(brief, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
