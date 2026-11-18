"""Stage 9: model comparison, verdicts and lag against luminosity (docs/PREDICTIONS_9.md, sections 4-6).

    python scripts/run_stage9.py verdict

Reads results/s9_objects.csv, results/s9_transfer.json, results/s9_identifiability.json and the
synthetic samples in data/cache/s9/ ('ident' and 'fitted'). Writes results/s9_verdict.json.
Only aggregate numbers are written; no catalogue value of a single object.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import mannwhitneyu, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import run_stage9 as r9  # noqa: E402
from irecho import config  # noqa: E402
from irecho import transfer9 as t9  # noqa: E402

RES = config.RESULTS_DIR
TYPE1 = ("Sy1", "Sy1.2", "Sy1.5", "Sy1.8")
TYPE2 = ("Sy1.9", "Sy2", "Sy2 candidate")
# Mandal et al. 2024 (arXiv:2403.01885), Table 5, BCES, ICCF lags, with the accretion-disc correction:
# log(R / 1 light-day) = beta + alpha log(L_bol / 10^45.7 erg/s). Not in the plan: a post hoc comparison.
MANDAL = {"w1": (0.39, 2.66), "w2": (0.33, 2.76), "log_l0": 45.7}


def mandal_delay_rest(log_lbol: float) -> float:
    x = log_lbol - MANDAL["log_l0"]
    return 10.0 ** (MANDAL["w2"][1] + MANDAL["w2"][0] * x) - 10.0 ** (MANDAL["w1"][1] + MANDAL["w1"][0] * x)


def real_rows() -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {t: [] for t in r9.TAGS}
    with open(RES / "s9_objects.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            row = {"name": r["name"]}
            for k, v in r.items():
                if k in ("sample", "name"):
                    continue
                row[k] = (v == "1") if k.startswith(("usable", "comparable")) or k in ("tau_at_edge", "well_determined") else (
                    float(v) if v != "" else float("nan"))
            out[r["sample"]].append(row)
    return out


def type_group(o: dict) -> str:
    if "type_group" in o:
        return o["type_group"] if o["type_group"] in ("1", "2") else "other"
    return "1" if o["type"] in TYPE1 else ("2" if o["type"] in TYPE2 else "other")


def label_explained(m_res: float, f_res: float, bands: dict) -> str:
    inside_m = all(b["m_res"][0] <= m_res <= b["m_res"][1] for b in bands.values())
    inside_f = all(b["f_res"][0] <= f_res <= b["f_res"][1] for b in bands.values())
    outside_m = all(not (b["m_res"][0] <= m_res <= b["m_res"][1]) for b in bands.values())
    if inside_m and inside_f:
        return "설명됨"
    if outside_m:
        return "남음"
    return "판정불가"


def lagl_combined(objs_by_tag: dict, tau_by_tag: dict, delay_fn=None) -> dict:
    """k and rho for the union of the two samples (objects with a luminosity and a usable fit)."""
    rest, lum, ratio = [], [], []
    for tag in r9.TAGS:
        for o, tau in zip(objs_by_tag[tag], tau_by_tag[tag]):
            if tau is None or o["log_lbol"] is None or not np.isfinite(tau):
                continue
            t = tau / (1.0 + o["z"])
            rest.append(t)
            lum.append(o["log_lbol"])
            ratio.append(t / (o["lyu_delay_rest"] if delay_fn is None else delay_fn(o["log_lbol"])))
    return {"n": len(rest), "k": float(np.median(ratio)), "rho": float(spearmanr(lum, rest).statistic),
            "median_tau_rest_days": float(np.median(rest))}


def main() -> None:
    ident = json.loads((RES / "s9_identifiability.json").read_text())
    transfer = json.loads((RES / "s9_transfer.json").read_text())
    rows = real_rows()
    objs = {tag: r9.load(tag) for tag in r9.TAGS}
    out = {"note": "stage 9: verdicts by the rules of docs/PREDICTIONS_9.md; earlier verdicts are not changed",
           "identifiability": ident["overall"], "samples": {}}
    reduce = {"lag_not_recovered": not ident["overall"]["lag_recovered"], "t0_t1_not_separated": not ident["overall"]["t0_vs_t1"],
              "t1_t2_not_separated": not ident["overall"]["t1_vs_t2"], "fixed_t3_not_separated": not ident["overall"]["fixed_vs_t3"]}
    out["interpretation_reduced_by_b1"] = reduce
    labels = {}
    for tag in r9.TAGS:
        st = transfer["samples"][tag]["statistics"]
        rec = {"statistics": st, "model_comparison": {}, "explained": {}}
        t0 = {d: r9.read_samples("ident", tag, d) for d in r9.DRIVERS}
        fitted = {d: r9.read_samples("fitted", tag, d) for d in r9.DRIVERS}
        rec["n_null_samples"] = {"T0": {d: len(v) for d, v in t0.items()}, "fitted": {d: len(v) for d, v in fitted.items()}}
        # --- 4.1 every model against T0, and against each other
        for name, key, null, truth in (("T1 against T0", "mean_ds_01", t0, "T0"), ("T2 against T0", "mean_ds_02", t0, "T0"),
                                       ("T3 against T0", "mean_ds_eta_t0", t0, "T0"), ("T2 against T1", "mean_ds_12", fitted, "fitted"),
                                       ("T3 against T1", "mean_ds_eta_t1", fitted, "fitted"),
                                       ("T3 against T2", "mean_ds_eta_t2", fitted, "fitted")):
            per = {}
            for d in r9.DRIVERS:
                v = r9.stat(null[d], truth, key)
                per[d] = {"null_median": float(np.nanmedian(v)), "null_max": float(np.nanmax(v)),
                          "null_range": r9.band(v)}
            rec["model_comparison"][name] = {"statistic": key, "observed": st[key], "null_truth": truth, "null": per,
                                             "simpler_model_not_enough": bool(all(st[key] > per[d]["null_max"] for d in per))}
        # --- 4.2 explained by a fixed kernel?
        for k in ("t1", "t2"):
            bands = {}
            for d in r9.DRIVERS:
                bands[d] = {"m_res": r9.band(r9.stat(fitted[d], "fitted", f"median_d_colour_res_{k}")),
                            "f_res": r9.band(r9.stat(fitted[d], "fitted", f"frac_res_bluer_{k}")),
                            "m_res_median": float(np.nanmedian(r9.stat(fitted[d], "fitted", f"median_d_colour_res_{k}"))),
                            "ratio_median": float(np.nanmedian(r9.stat(fitted[d], "fitted", f"ratio_pred_obs_{k}"))),
                            "ratio": r9.band(r9.stat(fitted[d], "fitted", f"ratio_pred_obs_{k}")),
                            "m_obs": r9.band(r9.stat(fitted[d], "fitted", "median_d_colour_obs")),
                            "m_obs_median": float(np.nanmedian(r9.stat(fitted[d], "fitted", "median_d_colour_obs"))),
                            "m_pred_median": float(np.nanmedian(r9.stat(fitted[d], "fitted", f"median_d_colour_pred_{k}")))}
            m_res, f_res = st[f"median_d_colour_res_{k}"], st[f"frac_res_bluer_{k}"]
            lab = label_explained(m_res, f_res, bands)
            rec["explained"][k] = {
                "label": lab, "m_obs_mag": st["median_d_colour_obs"], "m_pred_mag": st[f"median_d_colour_pred_{k}"],
                "m_res_mag": m_res, "f_res_bluer": f_res, "ratio_pred_obs": st[f"ratio_pred_obs_{k}"], "null_fitted": bands,
                "remaining_mag": {d: m_res - bands[d]["m_res_median"] for d in bands},
                "offsets": transfer["samples"][tag][f"offsets_{k}"],
                "n_res_beyond_3sigma": {"observed": transfer["samples"][tag][f"offsets_{k}"]["n_res_beyond_3sigma"],
                                        "null_fitted": {d: {"median": float(np.nanmedian(r9.stat(fitted[d], "fitted", f"n_res_beyond_3sigma_{k}"))),
                                                            "null_range": r9.band(r9.stat(fitted[d], "fitted", f"n_res_beyond_3sigma_{k}"))}
                                                        for d in r9.DRIVERS}}}
        labels[tag] = {k: rec["explained"][k]["label"] for k in ("t1", "t2")}
        rec["explained"]["label_sample"] = labels[tag]["t1"] if labels[tag]["t1"] == labels[tag]["t2"] else "모형에 달림"
        # --- size comparison: which input delay gives the observed m_obs (B1 samples)
        size = {}
        for d in r9.DRIVERS:
            truths = ident["configs"][f"{tag}_{d}"]["truths"]
            delays = np.array([0.0, 20.0, 40.0, 80.0, 160.0])
            m = np.array([truths[n]["median_d_colour_obs"]["median"] for n in ("T0", "T1_20", "T1_40", "T1_80", "T1_160")])
            obs = st["median_d_colour_obs"]
            order = np.argsort(m)
            size[d] = {"delay_days": delays.tolist(), "m_obs_mag": m.tolist(),
                       "delay_giving_observed_m_obs_days": float(np.interp(obs, m[order], delays[order])),
                       "m_obs_of_T1_lyu_mag": truths["T1_lyu"]["median_d_colour_obs"]["median"],
                       "m_obs_of_T3_1yr_mag": truths["T3_1yr"]["median_d_colour_obs"]["median"],
                       "m_obs_of_T3_3yr_mag": truths["T3_3yr"]["median_d_colour_obs"]["median"]}
        rec["size_comparison"] = {"observed_m_obs_mag": st["median_d_colour_obs"], "fitted_median_tau_days": st["median_tau"],
                                  "median_lyu_delay_observed_frame_days": float(np.median([o["lyu_delay"] for o in objs[tag]])),
                                  "by_driver": size}
        # --- B4 (exploratory): split by type
        split = {}
        for grp in ("1", "2"):
            sel = [r for o, r in zip(objs[tag], rows[tag]) if type_group(o) == grp]
            s = t9.sample_statistics(sel)
            split[grp] = {k: s.get(k) for k in ("n_fit", "n_offset", "median_tau", "frac_tau_positive", "median_d_colour_obs",
                                                "median_d_colour_pred_t2", "median_d_colour_res_t2", "frac_res_bluer_t2",
                                                "ratio_pred_obs_t2", "n_well")}
        tau1 = [r["tau_t1"] for o, r in zip(objs[tag], rows[tag]) if r["usable"] and type_group(o) == "1"]
        tau2 = [r["tau_t1"] for o, r in zip(objs[tag], rows[tag]) if r["usable"] and type_group(o) == "2"]
        split["mann_whitney_p_tau_two_sided"] = float(mannwhitneyu(tau1, tau2, alternative="two-sided").pvalue)
        split["note"] = "exploratory, post hoc comparison; no verdict"
        rec["by_type_exploratory"] = split
        out["samples"][tag] = rec
    # section 4.2: the label of a sample is that of T1 and T2 when they agree, else '모형에 달림';
    # the two samples are then compared by their sample labels
    s1, s2 = (out["samples"][t]["explained"]["label_sample"] for t in r9.TAGS)
    lab = s1 if s1 == s2 else "표본에 따라 다름"
    wording = lab
    if lab == "설명됨" and reduce["fixed_t3_not_separated"]:
        wording = "inter-band lag 모형과 모순되지 않음 (재형성과는 구별하지 못함)"
    out["explained_overall"] = {"label_by_rule": lab, "wording_after_b1": wording, "by_sample": labels}

    # --- B3: lag against luminosity (primary: the two samples together)
    tau_real = {tag: [r["tau_t1"] if r["usable"] else None for r in rows[tag]] for tag in r9.TAGS}
    b3 = {"observed": {"combined": lagl_combined(objs, tau_real)}, "null": {}}
    for tag in r9.TAGS:
        b3["observed"][tag] = r9.lag_luminosity(objs[tag], [r if r["usable"] else None for r in rows[tag]])
    for truth in ("T1_lyu", "T0"):
        b3["null"][truth] = {}
        for d in r9.DRIVERS:
            s1, s2 = r9.read_samples("ident", "v1", d), r9.read_samples("ident", "v2", d)
            ks, rhos = [], []
            for a, b in zip(s1, s2):
                c = lagl_combined(objs, {"v1": a[truth]["tau"], "v2": b[truth]["tau"]})
                ks.append(c["k"])
                rhos.append(c["rho"])
            b3["null"][truth][d] = {"k_median": float(np.median(ks)), "k_null_range": r9.band(ks), "rho_median": float(np.median(rhos)),
                                    "rho_null_range": r9.band(rhos),
                                    "by_sample": {tag: {"k_null_range": ident["configs"][f"{tag}_{d}"]["truths"][truth]["lag_luminosity"]["k_null_range"],
                                                        "rho_null_range": ident["configs"][f"{tag}_{d}"]["truths"][truth]["lag_luminosity"]["rho_null_range"]}
                                                  for tag in r9.TAGS}}
    k, rho = b3["observed"]["combined"]["k"], b3["observed"]["combined"]["rho"]

    def inside(v, rng):
        return rng[0] <= v <= rng[1]

    per_driver = []
    for d in r9.DRIVERS:
        lyu, zero = b3["null"]["T1_lyu"][d], b3["null"]["T0"][d]
        if inside(k, zero["k_null_range"]):
            per_driver.append("지연 없음과 구별 안 됨")
        elif inside(k, lyu["k_null_range"]) and inside(rho, lyu["rho_null_range"]):
            per_driver.append("문헌 관계와 맞음")
        elif not inside(k, lyu["k_null_range"]):
            per_driver.append("문헌 관계와 다름")
        else:
            per_driver.append("판정불가")
    label = per_driver[0] if len(set(per_driver)) == 1 else "판정불가"
    b3["label_by_rule"] = label
    b3["label_by_driver"] = dict(zip(r9.DRIVERS, per_driver))
    b3["label_after_b1"] = "이 자료로 시험할 수 없음" if (reduce["lag_not_recovered"] or reduce["t0_t1_not_separated"]) else label
    b3["calibrated_factor"] = {d: k / b3["null"]["T1_lyu"][d]["k_median"] for d in r9.DRIVERS}
    b3["post_hoc_mandal_2024"] = {"observed": lagl_combined(objs, tau_real, mandal_delay_rest),
                                  "relation": "log(R/1 lt-day) = beta + alpha log(L_bol/10^45.7): W1 alpha 0.39 beta 2.66, W2 alpha 0.33 beta 2.76",
                                  "note": "not in the plan; no null distribution was made for this relation"}
    # selection: well-determined lags against the rest
    sel = {}
    for tag in r9.TAGS:
        well, rest = [], []
        for o, r in zip(objs[tag], rows[tag]):
            if not r["usable"]:
                continue
            snr = float(np.std(o["f2"]) / np.median(o["s2"]))
            (well if r["well_determined"] else rest).append({"log_lbol": o["log_lbol"], "z": o["z"], "snr": snr, "type": type_group(o),
                                                              "tau": r["tau_t1"]})

        def summary(g):
            lum = [x["log_lbol"] for x in g if x["log_lbol"] is not None]
            return {"n": len(g), "n_with_lbol": len(lum), "median_log_lbol": float(np.median(lum)) if lum else None,
                    "median_z": float(np.median([x["z"] for x in g])) if g else None,
                    "median_w2_rms_over_error": float(np.median([x["snr"] for x in g])) if g else None,
                    "n_type1": sum(x["type"] == "1" for x in g), "n_type2": sum(x["type"] == "2" for x in g),
                    "median_tau_days": float(np.median([x["tau"] for x in g])) if g else None,
                    "frac_tau_positive": float(np.mean([x["tau"] > 0 for x in g])) if g else None}

        rec = {"well_determined": summary(well), "rest": summary(rest)}
        if len(well) >= 5 and len(rest) >= 5:
            rec["mann_whitney_p_w2_snr"] = float(mannwhitneyu([x["snr"] for x in well], [x["snr"] for x in rest]).pvalue)
            lw = [x["log_lbol"] for x in well if x["log_lbol"] is not None]
            lr = [x["log_lbol"] for x in rest if x["log_lbol"] is not None]
            if len(lw) >= 5 and len(lr) >= 5:
                rec["mann_whitney_p_log_lbol"] = float(mannwhitneyu(lw, lr).pvalue)
        sel[tag] = rec
    b3["selection"] = sel
    out["lag_luminosity"] = b3
    (RES / "s9_verdict.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    print("explained:", json.dumps(out["explained_overall"], ensure_ascii=False))
    for tag in r9.TAGS:
        rec = out["samples"][tag]
        print(tag, {k: v["simpler_model_not_enough"] for k, v in rec["model_comparison"].items()})
        for k in ("t1", "t2"):
            e = rec["explained"][k]
            print("  ", k, e["label"], "m_obs %.4f m_pred %.4f m_res %.4f f_res %.3f ratio %.2f" % (
                e["m_obs_mag"], e["m_pred_mag"], e["m_res_mag"], e["f_res_bluer"], e["ratio_pred_obs"]),
                {d: (b["m_res"], b["f_res"]) for d, b in e["null_fitted"].items()})
    print("B3:", b3["label_by_rule"], "->", b3["label_after_b1"], json.dumps(b3["observed"]), json.dumps(b3["calibrated_factor"]))


if __name__ == "__main__":
    main()
