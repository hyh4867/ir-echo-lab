"""Stage 10: sample statistics, the identifiability summary, the real fit and the checks.

Called through scripts/run_stage10.py (ident-summary, fit, fitted, fitted-summary, checks).
"""

from __future__ import annotations

import csv
import json
import time

import numpy as np
from scipy.stats import spearmanr, theilslopes

import run_stage10 as r10
from irecho import config, hyst8
from irecho import driver10 as d10

RES = r10.RES
K = {k: i for i, k in enumerate(r10.FIT_KEYS)}
KO = {k: i for i, k in enumerate(r10.OFF_KEYS)}
FIXED_KERNEL_TRUTHS = ("Z0", "A0", "A_20", "A_50", "A_100", "A_200", "A_300", "A_lyu", "A_mandal", "A_s9", "G_lyu", "E_lyu",
                       "W_thin", "W_wide", "B_20", "B_40", "B_both")
OBJECT_COLUMNS = ("name", "sample", "driver_band", "n_nights", "n_used", "detected_w1", "detected_w2", "in_d",
                  "null_rw_max_w1", "null_rw_max_w2", "n_swaps")


def catalogue(objs: list[dict]) -> dict:
    nan = float("nan")
    return {"z": np.array([o["z"] for o in objs]),
            "logl": np.array([o["log_lbol"] if o["log_lbol"] is not None else nan for o in objs]),
            "lyu": np.array([o["lyu_rest"] if o["lyu_rest"] is not None else (nan, nan) for o in objs]),
            "mandal": np.array([o["mandal_rest"] if o["mandal_rest"] is not None else (nan, nan) for o in objs]),
            "tau9": np.array([o["tau9"] if o["tau9"] is not None else nan for o in objs]),
            "f30": np.array([o["f30"] for o in objs]),
            "dof": np.array([max(o["idx"].size - 3, 1) for o in objs], dtype=float),
            "flagged": np.array([bool(o["flag_seeing"]) or bool(o["flag_airmass"]) for o in objs])}


def _med(x) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return float(np.median(x)) if x.size else float("nan")


def _mean(x) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return float(np.mean(x)) if x.size else float("nan")


def _rho(x, y) -> float:
    ok = np.isfinite(x) & np.isfinite(y)
    return float(spearmanr(x[ok], y[ok]).statistic) if ok.sum() >= 8 and np.ptp(x[ok]) > 0 and np.ptp(y[ok]) > 0 else float("nan")


def set_stats(f, off, sel, cat, lags=None) -> dict:
    """Section 4 of the plan for the objects in ``sel``. f: (n_obj, 2, keys); off: (n_obj, keys)."""
    out = {"n": int(sel.sum())}
    if sel.sum() < 3:
        return out
    tau = f[sel][:, :, K["tau"]].astype(float)
    dtau = tau[:, 1] - tau[:, 0]
    z1 = 1.0 + cat["z"][sel]
    out.update(med_tau1=_med(tau[:, 0]), med_tau2=_med(tau[:, 1]), med_dtau=_med(dtau),
               frac_dtau_positive=float(np.mean(dtau > 0)))
    for key in ("ds_c1", "ds_01", "ds_12", "ds_1g", "ds_1a", "ds_1e", "var_share", "gap_share"):
        for b in (0, 1):
            out[f"mean_{key}_w{b + 1}"] = _mean(f[sel][:, b, K[key]])
    for b in (0, 1):
        # chi2 of M1 per degree of freedom (not in the plan; added after the first real fit to compare
        # the size of the residuals of the real and the synthetic light curves)
        out[f"med_chi2_red_w{b + 1}"] = _med(f[sel][:, b, K["chi2_m1"]] / cat["dof"][sel])
        out[f"med_alpha_w{b + 1}"] = _med(f[sel][:, b, K["alpha"]])
        out[f"frac_z_alpha_w{b + 1}"] = float(np.mean(f[sel][:, b, K["z_alpha"]] >= 2.0))
        out[f"med_tau_ma_w{b + 1}"] = _med(f[sel][:, b, K["tau_ma"]])
        for shape in ("f0", "f05", "f2", "mg"):
            out[f"med_tau_{shape}_w{b + 1}"] = _med(f[sel][:, b, K[f"tau_{shape}"]])
    for shape in ("f0", "f05", "f2", "mg"):
        out[f"med_dtau_{shape}"] = _med(f[sel][:, 1, K[f"tau_{shape}"]] - f[sel][:, 0, K[f"tau_{shape}"]])
    dtau_a = tau[:, 1] - f[sel][:, 0, K["tau_ma"]].astype(float)
    out["med_dtau_alpha"] = _med(dtau_a)
    out["med_dtau_alpha_both"] = _med(f[sel][:, 1, K["tau_ma"]].astype(float) - f[sel][:, 0, K["tau_ma"]].astype(float))
    # colour offsets (stage 8 direct comparison, measured and predicted by the two fixed kernels)
    o = off[sel]
    ok = (o[:, KO["usable_obs"]] > 0.5) & (o[:, KO["comparable_obs"]] > 0.5) & (o[:, KO["usable_pred"]] > 0.5)
    out["n_offset"] = int(ok.sum())
    if ok.sum() >= 3:
        obs, pred = o[ok, KO["d_colour_obs"]].astype(float), o[ok, KO["d_colour_pred"]].astype(float)
        out.update(m_obs=_med(obs), m_pred=_med(pred), m_res=_med(obs - pred), f_res_bluer=float(np.mean(obs - pred < 0)),
                   frac_obs_bluer=float(np.mean(obs < 0)),
                   ratio_pred_obs=float(np.median(pred) / np.median(obs)) if np.median(obs) != 0 else float("nan"))
    # lag against luminosity
    logl = cat["logl"][sel]
    rest = tau / z1[:, None]
    for b in (0, 1):
        for rel in ("lyu", "mandal"):
            out[f"k_{rel}_w{b + 1}"] = _med(rest[:, b] / cat[rel][sel][:, b])
        out[f"rho_w{b + 1}"] = _rho(logl, rest[:, b])
        pos = np.isfinite(logl) & (rest[:, b] > 0)
        out[f"slope_w{b + 1}"] = float(theilslopes(np.log10(rest[pos, b]), logl[pos]).slope) if pos.sum() >= 8 else float("nan")
    out["n_lbol"] = int(np.isfinite(logl).sum())
    d_lyu = cat["lyu"][sel][:, 1] - cat["lyu"][sel][:, 0]
    out["k10"] = _med(dtau / z1 / d_lyu)
    out["k_alpha"] = _med(dtau_a / z1 / d_lyu)
    # Q1: against the stage 9 lag of the same objects
    t9v = cat["tau9"][sel]
    has = np.isfinite(t9v)
    out["n_q1"] = int(has.sum())
    if has.sum() >= 3:
        m9 = float(np.median(t9v[has]))
        out.update(med_tau9=m9, q=float(np.median(dtau[has]) / m9) if m9 != 0 else float("nan"),
                   d=float(np.median(dtau[has] - t9v[has])), rho_q1=_rho(t9v[has], dtau[has]),
                   q_alpha=float(np.median(dtau_a[has]) / m9) if m9 != 0 else float("nan"),
                   k9=_med(t9v[has] / z1[has] / d_lyu[has]))
    if lags is not None:
        lg = lags[sel].astype(float)
        out.update(bias_tau1=out["med_tau1"] - _med(lg[:, 0]), bias_tau2=out["med_tau2"] - _med(lg[:, 1]),
                   bias_dtau=out["med_dtau"] - _med(lg[:, 1] - lg[:, 0]), injected_tau1=_med(lg[:, 0]),
                   injected_tau2=_med(lg[:, 1]),
                   rms_tau1=float(np.sqrt(np.mean((tau[:, 0] - lg[:, 0]) ** 2))),
                   rms_tau2=float(np.sqrt(np.mean((tau[:, 1] - lg[:, 1]) ** 2))),
                   rms_dtau=float(np.sqrt(np.mean((dtau - (lg[:, 1] - lg[:, 0])) ** 2))),
                   mad_dtau=_med(np.abs(dtau - (lg[:, 1] - lg[:, 0]))))
        for b in (0, 1):
            lo, hi = f[sel][:, b, K["tau_q16"]], f[sel][:, b, K["tau_q84"]]
            out[f"coverage_tau{b + 1}"] = float(np.mean((lo <= lg[:, b] + 2.5) & (hi >= lg[:, b] - 2.5)))
    return out


def sample_stats(f, off, det, cat, lags=None, extra_sel=None) -> dict:
    """Statistics of the detected set D and of all fitted objects."""
    usable = np.isfinite(f[:, 0, K["tau"]]) & np.isfinite(f[:, 1, K["tau"]])
    edge = (f[:, :, K["tau_at_edge"]] > 0.5).any(axis=1)
    in_d = usable & det.all(axis=1) & ~edge
    if extra_sel is not None:
        usable, in_d = usable & extra_sel, in_d & extra_sel
    out = {"d": set_stats(f, off, in_d, cat, lags), "all": set_stats(f, off, usable, cat, lags)}
    out["d"]["frac_detected"] = out["all"]["frac_detected"] = float(in_d.sum() / max(usable.sum(), 1))
    for b in (0, 1):
        out["all"][f"frac_detected_w{b + 1}"] = float((det[:, b] & usable).sum() / max(usable.sum(), 1))
    out["all"]["frac_edge"] = float((edge & usable).sum() / max(usable.sum(), 1))
    return out, in_d


def null_threshold(fits, truths) -> np.ndarray:
    """Per object and band: the third largest r_w of the X0 realisations (section 3 of the plan)."""
    rw = fits[:, truths.index("X0"), :, :, K["r_w"]].astype(float)          # obj, real, band
    rw = np.where(np.isfinite(rw), rw, -np.inf)
    return np.sort(rw, axis=1)[:, -3, :]


def band(values) -> dict:
    """Median, 16-84 per cent and the null range (smallest and largest value removed)."""
    v = np.sort(np.asarray([x for x in values if x is not None and np.isfinite(x)], dtype=float))
    if v.size < 5:
        return {"n": int(v.size)}
    return {"n": int(v.size), "median": float(np.median(v)), "p16": float(np.percentile(v, 16)), "p84": float(np.percentile(v, 84)),
            "lo": float(v[1]), "hi": float(v[-2])}


def frac_outside(values, ref: dict, side: str = "both") -> float:
    v = np.asarray([x for x in values if x is not None and np.isfinite(x)], dtype=float)
    if v.size == 0 or "lo" not in ref:
        return float("nan")
    if side == "above":
        return float(np.mean(v > ref["hi"]))
    return float(np.mean((v < ref["lo"]) | (v > ref["hi"])))


def realisation_stats(kind: str, driver: str, thr=None) -> tuple[dict, list[str], np.ndarray]:
    """Per truth model: the list of sample statistics over the realisations."""
    objs = r10.load()
    cat = catalogue(objs)
    pack = r10.read_samples(kind, driver)
    truths = pack["truths"]
    thr = null_threshold(pack["fits"], truths) if thr is None else thr
    out = {}
    for j, name in enumerate(truths):
        rows = []
        for i in range(pack["fits"].shape[2]):
            f = pack["fits"][:, j, i]
            det = f[:, :, K["r_w"]] > thr
            st, _ = sample_stats(f, pack["offs"][:, j, i], det, cat, pack["lags"][:, j])
            rows.append(st)
        out[name] = rows
    return out, truths, thr, pack


def ident_summary() -> None:
    summary = {"n_objects": len(r10.load()), "n_realisations": r10.N_IDENT, "drivers": {}, "criteria": {}}
    crit_all = {}
    thresholds = {}
    for driver in r10.DRIVERS:
        stats, truths, thr, pack = realisation_stats("ident", driver)
        thresholds[driver] = thr
        rec = {"truths": {}}
        keys = sorted({k for rows in stats.values() for r in rows for k in r["d"]})
        tab = {}
        for name in truths:
            tab[name] = {s: {k: band([r[s].get(k) for r in stats[name]]) for k in
                             sorted({k for r in stats[name] for k in r[s]})} for s in ("d", "all")}
        rec["truths"] = tab
        vals = lambda name, key, s="d": [r[s].get(key) for r in stats[name]]  # noqa: E731
        c = {}
        # 1. lags of each band at the sample level
        ok1 = []
        for name in ("A_100", "A_200", "A_lyu"):
            for b in (1, 2):
                bias = tab[name]["d"].get(f"bias_tau{b}", {})
                inj = tab[name]["d"].get(f"injected_tau{b}", {}).get("median", float("nan"))
                good = bool("median" in bias and abs(bias["median"]) <= max(10.0, 0.10 * inj)
                            and 0.5 * (bias["p84"] - bias["p16"]) <= 15.0)
                c[f"lag_bias_{name}_w{b}"] = {"median_error_days": bias.get("median"), "injected_median_days": inj,
                                             "half_width_16_84_days": 0.5 * (bias["p84"] - bias["p16"]) if "p84" in bias else None,
                                             "pass": good}
                ok1.append(good)
        c["lag_recovered"] = all(ok1)
        # 2. the lag difference
        ok2 = []
        for name in ("A_50", "A_100", "A_200"):
            bias = tab[name]["d"].get("bias_dtau", {})
            good = bool("median" in bias and abs(bias["median"]) <= 8.0 and 0.5 * (bias["p84"] - bias["p16"]) <= 10.0)
            c[f"dtau_bias_{name}"] = {"median_error_days": bias.get("median"),
                                      "half_width_16_84_days": 0.5 * (bias["p84"] - bias["p16"]) if "p84" in bias else None,
                                      "pass": good}
            ok2.append(good)
        c["dtau_power_lyu_vs_a0"] = frac_outside(vals("A_lyu", "med_dtau"), tab["A0"]["d"].get("med_dtau", {}), "above")
        c["dtau_recovered"] = bool(all(ok2) and c["dtau_power_lyu_vs_a0"] >= 0.8)
        # 3. detection
        c["detection_power"] = frac_outside(vals("A_lyu", "frac_detected", "all"), tab["X0"]["all"].get("frac_detected", {}), "above")
        c["detected"] = bool(c["detection_power"] >= 0.8)
        # 4. the undelayed share
        a_ok = []
        for name, inj in (("A_lyu", 0.0), ("B_20", 0.2), ("B_40", 0.4)):
            m = tab[name]["d"].get("med_alpha_w1", {}).get("median")
            good = bool(m is not None and abs(m - inj) <= 0.1)
            c[f"alpha_bias_{name}"] = {"median_alpha": m, "injected": inj, "pass": good}
            a_ok.append(good)
        for name in ("B_20", "B_40", "B_both"):
            c[f"alpha_power_{name}"] = {k: frac_outside(vals(name, k), tab["A_lyu"]["d"].get(k, {}), "above")
                                        for k in ("med_alpha_w1", "mean_ds_1a_w1")}
        c["alpha_power_ok_b40"] = bool(max(c["alpha_power_B_40"].values()) >= 0.8)
        c["alpha_power_ok_b20"] = bool(max(c["alpha_power_B_20"].values()) >= 0.8)
        c["alpha_recovered"] = bool(all(a_ok) and c["alpha_power_ok_b40"])
        # 5. the width
        c["width_power"] = {name: {b: frac_outside(vals(name, f"mean_ds_12_w{b}"), tab["A_lyu"]["d"].get(f"mean_ds_12_w{b}", {}), "above")
                                   for b in (1, 2)} for name in ("W_thin", "W_wide")}
        c["width_recovered"] = bool(all(max(v.values()) >= 0.8 for v in c["width_power"].values()))
        # 6. fixed kernel against the history-dependent response
        p6 = {}
        for name in ("C_0", "C_1yr", "C_3yr"):
            p6[name] = {}
            for key in ("m_res", "mean_ds_1e_w1", "mean_ds_1e_w2"):
                v = np.array([x if x is not None else np.nan for x in vals(name, key)], dtype=float)
                out_all = np.isfinite(v)
                for ref in FIXED_KERNEL_TRUTHS:
                    rb = tab[ref]["d"].get(key, {})
                    out_all &= (v < rb["lo"]) | (v > rb["hi"]) if "lo" in rb else False
                p6[name][key] = float(np.mean(out_all))
        c["fixed_vs_history_power"] = p6
        c["fixed_vs_history"] = bool(all(max(p6[name].values()) >= 0.8 for name in ("C_1yr", "C_3yr")))
        # 7. the two literature relations
        c["relations_power"] = {f"w{b}": frac_outside(vals("A_mandal", f"k_lyu_w{b}"), tab["A_lyu"]["d"].get(f"k_lyu_w{b}", {}))
                                for b in (1, 2)}
        c["relations_power_reverse"] = {f"w{b}": frac_outside(vals("A_lyu", f"k_mandal_w{b}"), tab["A_mandal"]["d"].get(f"k_mandal_w{b}", {}))
                                        for b in (1, 2)}
        c["relations_separated"] = bool(max(c["relations_power"].values()) >= 0.8)
        # per-object accuracy and coverage (no criterion)
        acc = {}
        for name in ("A_100", "A_lyu", "A_mandal", "G_lyu", "E_lyu"):
            acc[name] = {s: {k: tab[name][s].get(k, {}).get("median") for k in
                             ("rms_tau1", "rms_tau2", "rms_dtau", "mad_dtau", "coverage_tau1", "coverage_tau2", "n")}
                         for s in ("d", "all")}
        if driver == "fast":
            j = truths.index("A_lyu")
            jt = pack["joint"][:, :10].astype(float)                           # obj, the 10 realisations with a joint fit, 5
            true = (pack["lags"][:, j, 1] - pack["lags"][:, j, 0]).astype(float)[:, None]
            ok = np.isfinite(jt[:, :, 0])
            f = pack["fits"][:, j, :10]
            det = (f[:, :, :, K["r_w"]] > thr[:, None, :]).all(axis=2) & ~(f[:, :, :, K["tau_at_edge"]] > 0.5).any(axis=2)
            for label, sel in (("all", ok), ("d", ok & det)):
                acc[f"joint_dtau_{label}"] = {
                    "n": int(sel.sum()),
                    "coverage_16_84": float(np.mean(((jt[:, :, 0] <= true + 5.0) & (jt[:, :, 2] >= true - 5.0))[sel])),
                    "coverage_2p5_97p5": float(np.mean(((jt[:, :, 3] <= true + 5.0) & (jt[:, :, 4] >= true - 5.0))[sel])),
                    "median_half_width_16_84_days": float(np.median((0.5 * (jt[:, :, 2] - jt[:, :, 0]))[sel]))}
        c["per_object"] = acc
        rec["criteria"] = c
        summary["drivers"][driver] = rec
        crit_all[driver] = c
        del keys
    flags = ("lag_recovered", "dtau_recovered", "detected", "alpha_recovered", "alpha_power_ok_b20", "width_recovered",
             "fixed_vs_history", "relations_separated")
    summary["overall"] = {k: bool(all(crit_all[d][k] for d in r10.DRIVERS)) for k in flags}
    (RES / "s10_identifiability.json").write_text(json.dumps(summary, indent=1) + "\n")
    np.save(r10.CKPT / "ident_thresholds.npy", np.stack([thresholds[d] for d in r10.DRIVERS]))
    print(json.dumps(summary["overall"], indent=1))
    for d in r10.DRIVERS:
        c = crit_all[d]
        print(d, {k: v for k, v in c.items() if k not in ("per_object",)})


# ----------------------------------------------------------------------------- section 6: the real light curves
def choose_swaps(objs: list[dict]) -> list[list[int]]:
    """For every object 20 other objects (random order) whose optical light curve gives at least
    MIN_VISITS usable visits at this object's visit times."""
    rng = np.random.default_rng([config.SEED, r10.SEED_TAG, 77])
    out = []
    for i, o in enumerate(objs):
        picks = []
        for j in rng.permutation(len(objs)):
            if j == i:
                continue
            if d10.usable_visits(o["mjd"], objs[j]["opt_t"]).sum() >= d10.MIN_VISITS:
                picks.append(int(j))
            if len(picks) == r10.N_SWAPS:
                break
        out.append(picks)
    return out


def fit_real() -> None:
    objs = r10.load()
    swaps = choose_swaps(objs)
    users = {j: [i for i, picks in enumerate(swaps) if j in picks] for j in range(len(objs))}
    n = len(objs)
    fits = np.full((n, 2, len(r10.FIT_KEYS)), np.nan)
    offs = np.full((n, len(r10.OFF_KEYS)), np.nan)
    joint, extra = [None] * n, [None] * n
    null_rw = [[[], []] for _ in range(n)]
    null_tau = [[[], []] for _ in range(n)]
    rev = np.full((n, 2, len(r10.FIT_KEYS)), np.nan)
    examples = {}
    t0 = time.time()
    for j, o in enumerate(objs):
        state = d10.gp_state(o["opt_t"], o["opt_s"], o["hyper"])
        prep = d10.prepare(o["opt_t"], o["opt_s"], o["hyper"], o["mjd"], joint=True, state=state)
        a, b, off, fa, fb = r10.fit_pair(o, prep, o["opt_y"], o["f1"], o["f2"])
        fits[j, 0], fits[j, 1], offs[j] = a, b, off
        joint[j] = d10.fit_joint(prep, o["opt_y"], o["f1"], o["s1"], o["f2"], o["s2"], o["r"], fa, fb)
        extra[j] = {"n_used": int(prep["idx"].size), "a_w1": fa["a"], "b_w1": fa["b"], "a_w2": fb["a"], "b_w2": fb["b"],
                    "tau_lo_w1": fa["tau_lo"], "tau_hi_w1": fa["tau_hi"], "tau_lo_w2": fb["tau_lo"], "tau_hi_w2": fb["tau_hi"],
                    "shape_m2_w1": fa["shape_m2"], "shape_m2_w2": fb["shape_m2"], "tau_m2_w1": fa["tau_m2"], "tau_m2_w2": fb["tau_m2"]}
        examples[o["name"]] = {"idx": prep["idx"].tolist(), "pred_w1": fa["pred"].tolist(), "pred_w2": fb["pred"].tolist(),
                               "pred_ma_w1": fa["pred_ma"].tolist()}
        for i in users[j]:                                   # this object's driver against other objects' visits
            oi = objs[i]
            pi = d10.prepare(o["opt_t"], None, None, oi["mjd"], shapes=(d10.BASE,), state=state)
            for bnd, (f, s) in enumerate(((oi["f1"], oi["s1"]), (oi["f2"], oi["s2"]))):
                r = d10.fit_band(pi, o["opt_y"], f, s, full=False)
                null_rw[i][bnd].append(r["r_w"])
                null_tau[i][bnd].append(r["tau"])
        # time reversal of the same driver
        t_rev = (o["opt_t"][0] + o["opt_t"][-1]) - o["opt_t"][::-1]
        prep_r = d10.prepare(t_rev, o["opt_s"][::-1], o["hyper"], o["mjd"])
        if prep_r["W"] is not None:
            ar, br, _, _, _ = r10.fit_pair(o, prep_r, o["opt_y"][::-1], o["f1"], o["f2"])
            rev[j, 0], rev[j, 1] = ar, br
        if (j + 1) % 10 == 0 or j + 1 == n:
            print(f"{time.time() - t0:6.0f} s  {j + 1}/{n}", flush=True)
    det = np.array([[fits[i, b, K["r_w"]] > max(null_rw[i][b]) if len(null_rw[i][b]) >= 10 else False for b in (0, 1)]
                    for i in range(n)])
    cat = catalogue(objs)
    stats, in_d = sample_stats(fits, offs, det, cat)
    np.savez_compressed(r10.CKPT / "fit_real.npz", fits=fits, offs=offs, det=det, rev=rev,
                        null_rw=np.array([[np.pad(np.array(x[b], dtype=float), (0, r10.N_SWAPS - len(x[b])), constant_values=np.nan)
                                           for b in (0, 1)] for x in null_rw]),
                        null_tau=np.array([[np.pad(np.array(x[b], dtype=float), (0, r10.N_SWAPS - len(x[b])), constant_values=np.nan)
                                            for b in (0, 1)] for x in null_tau]))
    (r10.CKPT / "fit_examples.json").write_text(json.dumps(examples) + "\n")
    # object table: no catalogue values (type, z, L_bol), observed frame only
    with open(RES / "s10_objects.csv", "w", newline="") as fh:
        cols = list(OBJECT_COLUMNS) + [f"{k}_w{b + 1}" for b in (0, 1) for k in r10.FIT_KEYS] + list(r10.OFF_KEYS) + [
            "dtau", "dtau_q025", "dtau_q16", "dtau_q50", "dtau_q84", "dtau_q975", "tau1_joint", "tau2_joint", "p_dtau_positive",
            "tau_rev_w1", "tau_rev_w2", "rw_rev_w1", "rw_rev_w2"] + sorted(extra[0])
        wr = csv.writer(fh)
        wr.writerow(cols)

        def fmt(v):
            if isinstance(v, (bool, np.bool_)):
                return int(v)
            if isinstance(v, str):
                return v
            return "" if v is None or not np.isfinite(v) else f"{float(v):.6g}"

        for i, o in enumerate(objs):
            row = [o["name"], o["sample"], o["band"], o["opt_t"].size, extra[i]["n_used"], det[i, 0], det[i, 1], in_d[i],
                   max(null_rw[i][0]) if null_rw[i][0] else None, max(null_rw[i][1]) if null_rw[i][1] else None, len(null_rw[i][0])]
            row += [fits[i, b, K[k]] for b in (0, 1) for k in r10.FIT_KEYS] + [offs[i, KO[k]] for k in r10.OFF_KEYS]
            jt = joint[i]
            row += [fits[i, 1, K["tau"]] - fits[i, 0, K["tau"]], jt["dtau_q025"], jt["dtau_q16"], jt["dtau_q50"], jt["dtau_q84"],
                    jt["dtau_q975"], jt["tau1_joint"], jt["tau2_joint"], jt["p_dtau_positive"],
                    rev[i, 0, K["tau"]], rev[i, 1, K["tau"]], rev[i, 0, K["r_w"]], rev[i, 1, K["r_w"]]]
            row += [extra[i][k] for k in sorted(extra[0])]
            wr.writerow([fmt(v) for v in row])
    # controls
    pooled = np.array([x for i in range(n) for b in (0, 1) for x in null_rw[i][b]])
    pooled_tau = np.array([x for i in range(n) for b in (0, 1) for x in null_tau[i][b]])
    rev_ok = np.isfinite(rev[:, 0, K["tau"]])
    rev_det = np.array([[rev[i, b, K["r_w"]] > max(null_rw[i][b]) if len(null_rw[i][b]) >= 10 else False for b in (0, 1)]
                        for i in range(n)])
    controls = {
        "swapped": {"n_fits": int(pooled.size), "r_w_p50_p90_p95_p99": [float(np.percentile(pooled, q)) for q in (50, 90, 95, 99)],
                    "tau_median": float(np.median(pooled_tau)), "tau_p16_p84": [float(np.percentile(pooled_tau, q)) for q in (16, 84)],
                    "frac_tau_at_edge": float(np.mean((pooled_tau <= -65 + 1e-9) | (pooled_tau >= 465 - 1e-9))),
                    "n_swaps_min": int(min(len(x[0]) for x in null_rw))},
        "real_r_w": {f"w{b + 1}": [float(np.percentile(fits[:, b, K["r_w"]], q)) for q in (10, 50, 90)] for b in (0, 1)},
        "time_reversed": {"n": int(rev_ok.sum()),
                          "frac_detected_both": float((rev_det.all(axis=1) & rev_ok).sum() / max(rev_ok.sum(), 1)),
                          "med_tau_w1": _med(rev[:, 0, K["tau"]]), "med_tau_w2": _med(rev[:, 1, K["tau"]]),
                          "med_dtau": _med(rev[:, 1, K["tau"]] - rev[:, 0, K["tau"]]),
                          "med_alpha_w1": _med(rev[:, 0, K["alpha"]]),
                          "frac_z_alpha_w1": float(np.mean(rev[rev_ok, 0, K["z_alpha"]] >= 2.0)),
                          "r_w_median": [_med(rev[:, b, K["r_w"]]) for b in (0, 1)]},
    }
    # bootstrap interval of q (objects resampled)
    rng = np.random.default_rng([config.SEED, r10.SEED_TAG, 78])
    boot = {}
    for label, sel in (("d", in_d), ("all", np.isfinite(fits[:, 0, K["tau"]]))):
        ok = np.nonzero(sel & np.isfinite(cat["tau9"]))[0]
        dt = fits[ok, 1, K["tau"]] - fits[ok, 0, K["tau"]]
        t9v = cat["tau9"][ok]
        qs, ms = [], []
        for _ in range(2000):
            pick = rng.integers(0, ok.size, ok.size) if ok.size else np.zeros(0, int)
            if ok.size and np.median(t9v[pick]) != 0:
                qs.append(np.median(dt[pick]) / np.median(t9v[pick]))
                ms.append(np.median(dt[pick]))
        boot[label] = {"n": int(ok.size), "q_ci95": [float(np.percentile(qs, 2.5)), float(np.percentile(qs, 97.5))] if qs else None,
                       "med_dtau_ci95": [float(np.percentile(ms, 2.5)), float(np.percentile(ms, 97.5))] if ms else None}
    # post hoc (not in the plan): bootstrap intervals of the statistics of the detected set
    post = {}
    d_idx = np.nonzero(in_d)[0]
    if d_idx.size >= 5:
        keys = ("med_tau1", "med_tau2", "med_dtau", "k_lyu_w1", "k_lyu_w2", "k_mandal_w1", "k_mandal_w2", "slope_w1", "slope_w2",
                "rho_w1", "rho_w2", "med_alpha_w1", "med_alpha_w2", "med_dtau_alpha", "k10", "k9", "q")
        draws = {k: [] for k in keys}
        for _ in range(2000):
            pick = d_idx[rng.integers(0, d_idx.size, d_idx.size)]
            sub = {k: v[pick] for k, v in cat.items()}
            st = set_stats(fits[pick], offs[pick], np.ones(pick.size, dtype=bool), sub)
            for k in keys:
                v = st.get(k)
                if v is not None and np.isfinite(v):
                    draws[k].append(v)
        post = {k: [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))] for k, v in draws.items() if len(v) > 100}
    jq = np.array([[joint[i][k] for k in ("dtau_q16", "dtau_q50", "dtau_q84")] for i in range(n)])
    out = {"n_objects": n, "n_d": int(in_d.sum()), "statistics": stats, "controls": controls, "bootstrap_q": boot,
           "bootstrap_ci95_post_hoc": post,
           "joint": {"median_half_width_16_84_days_d": float(np.median(0.5 * (jq[in_d, 2] - jq[in_d, 0]))) if in_d.any() else None,
                     "median_half_width_16_84_days_all": float(np.median(0.5 * (jq[:, 2] - jq[:, 0]))),
                     "frac_p_dtau_positive_above_0p975_d": float(np.mean([joint[i]["p_dtau_positive"] > 0.975 for i in np.nonzero(in_d)[0]])) if in_d.any() else None},
           "subsets": {}}
    for label, sel in (("flag_free", ~cat["flagged"]), ("strong_optical", cat["f30"] >= np.median(cat["f30"])),
                       ("sample_v1", np.array([o["sample"] == "v1" for o in objs])),
                       ("sample_v2", np.array([o["sample"] == "v2" for o in objs])),
                       ("driver_g", np.array([o["band"] == "g" for o in objs])),
                       ("z_le_0p04", cat["z"] <= 0.04), ("z_gt_0p04", cat["z"] > 0.04)):
        st, d_sub = sample_stats(fits, offs, det, cat, extra_sel=sel)
        out["subsets"][label] = {"n": int(sel.sum()), "n_d": int(d_sub.sum()), "statistics": st}
    (RES / "s10_fit.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({"n": n, "n_d": int(in_d.sum()), "d": stats["d"]}, indent=1))


def load_fitted() -> None:
    with open(RES / "s10_objects.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            r10._FITTED[r["name"]] = (float(r["tau_w1"]), float(r["tau_w2"]))


def fitted_summary() -> None:
    """Null ranges of the sample statistics when the truth is two fixed kernels with the fitted lags."""
    thr = np.load(r10.CKPT / "ident_thresholds.npy")
    out = {}
    for di, driver in enumerate(r10.DRIVERS):
        stats, truths, _, _ = realisation_stats("fitted", driver, thr=thr[di])
        rows = stats["fitted"]
        out[driver] = {s: {k: band([r[s].get(k) for r in rows]) for k in sorted({k for r in rows for k in r[s]})}
                       for s in ("d", "all")}
    (RES / "s10_null_fitted.json").write_text(json.dumps({"n_realisations": r10.N_IDENT, "drivers": out}, indent=1) + "\n")
    print("written", RES / "s10_null_fitted.json")


def checks() -> None:
    """Section 7: other band, straight lines instead of the Gaussian process, hyper-parameters refitted."""
    objs = r10.load()
    sample = json.loads(r10.SAMPLE_PATH.read_text())["objects"]
    real = np.load(r10.CKPT / "fit_real.npz")
    fits, det = real["fits"], real["det"]
    cat = catalogue(objs)
    n = len(objs)
    lin = np.full((n, 2, len(r10.FIT_KEYS)), np.nan)
    lin_off = np.full((n, len(r10.OFF_KEYS)), np.nan)
    other = np.full((n, 2, len(r10.FIT_KEYS)), np.nan)
    t0 = time.time()
    for j, o in enumerate(objs):
        prep = d10.prepare(o["opt_t"], o["opt_s"], None, o["mjd"], method="linear")
        a, b, off, _, _ = r10.fit_pair(o, prep, o["opt_y"], o["f1"], o["f2"])
        lin[j, 0], lin[j, 1], lin_off[j] = a, b, off
        if o["other_band"]:
            ob = r10.read_nights(o["name"], o["other_band"])
            hyper = sample[o["name"]]["bands"][o["other_band"]]["hyper"]
            prep2 = d10.prepare(ob["t"], ob["s"], hyper, o["mjd"])
            if prep2["W"] is not None:
                a2, b2, _, _, _ = r10.fit_pair(o, prep2, ob["y"], o["f1"], o["f2"])
                other[j, 0], other[j, 1] = a2, b2
        if (j + 1) % 20 == 0 or j + 1 == n:
            print(f"{time.time() - t0:6.0f} s  {j + 1}/{n}", flush=True)
    st_lin, d_lin = sample_stats(lin, lin_off, det, cat)
    both = np.isfinite(other[:, 0, K["tau"]]) & np.isfinite(fits[:, 0, K["tau"]])
    d_real = np.isfinite(fits[:, 0, K["tau"]]) & det.all(axis=1) & ~(fits[:, :, K["tau_at_edge"]] > 0.5).any(axis=1)
    diff = {}
    for label, sel in (("all", both), ("d", both & d_real)):
        dt_main = fits[sel, 1, K["tau"]] - fits[sel, 0, K["tau"]]
        dt_other = other[sel, 1, K["tau"]] - other[sel, 0, K["tau"]]
        diff[label] = {"n": int(sel.sum()),
                       "med_tau1_main": _med(fits[sel, 0, K["tau"]]), "med_tau1_other": _med(other[sel, 0, K["tau"]]),
                       "med_tau2_main": _med(fits[sel, 1, K["tau"]]), "med_tau2_other": _med(other[sel, 1, K["tau"]]),
                       "med_dtau_main": _med(dt_main), "med_dtau_other": _med(dt_other),
                       "median_abs_diff_tau1": _med(np.abs(fits[sel, 0, K["tau"]] - other[sel, 0, K["tau"]])),
                       "median_abs_diff_tau2": _med(np.abs(fits[sel, 1, K["tau"]] - other[sel, 1, K["tau"]])),
                       "median_abs_diff_dtau": _med(np.abs(dt_main - dt_other)),
                       "rho_tau1": _rho(fits[sel, 0, K["tau"]], other[sel, 0, K["tau"]]),
                       "rho_tau2": _rho(fits[sel, 1, K["tau"]], other[sel, 1, K["tau"]])}
    # hyper-parameters refitted on every synthetic realisation (A_lyu, fast driver, 5 realisations)
    refit = {"fixed": [], "refit": []}
    for index in range(5):
        rows = {"fixed": [], "refit": []}
        for k, o in enumerate(objs):
            y, f1, f2, lag1, lag2 = r10.synth_cases(o, k, "fast", index, ["A_lyu"])["A_lyu"]
            if not (np.all(f1 > 0) and np.all(f2 > 0)):
                continue
            for label, hyper in (("fixed", o["hyper"]), ("refit", d10.fit_drw(o["opt_t"], y, o["opt_s"]))):
                prep = d10.prepare(o["opt_t"], o["opt_s"], hyper, o["mjd"], shapes=(d10.BASE,))
                ta = d10.fit_band(prep, y, f1, o["s1"], full=False)["tau"]
                tb = d10.fit_band(prep, y, f2, o["s2"], full=False)["tau"]
                rows[label].append((ta - lag1, tb - lag2, (tb - ta) - (lag2 - lag1)))
        for label in rows:
            arr = np.array(rows[label])
            refit[label].append([float(np.median(arr[:, c])) for c in range(3)])
        print(f"{time.time() - t0:6.0f} s  refit realisation {index + 1}/5", flush=True)
    out = {"linear_interpolation": {"statistics": st_lin, "n_d_with_the_main_detection_flags": int(d_lin.sum()),
                                    "median_abs_diff_tau1": _med(np.abs(lin[:, 0, K["tau"]] - fits[:, 0, K["tau"]])),
                                    "median_abs_diff_tau2": _med(np.abs(lin[:, 1, K["tau"]] - fits[:, 1, K["tau"]]))},
           "other_band": diff,
           "hyper_refit_A_lyu_fast": {"columns": ["median error tau1", "median error tau2", "median error dtau"],
                                      "all_objects_fixed": refit["fixed"], "all_objects_refit": refit["refit"]}}
    (RES / "s10_checks.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({"other_band": diff, "refit": out["hyper_refit_A_lyu_fast"]}, indent=1))


def dispatch(cmd: str, workers: int) -> None:
    if cmd == "ident-summary":
        ident_summary()
    elif cmd == "fit":
        if not (RES / "s10_identifiability.json").exists():
            raise SystemExit("the identifiability summary must exist before a real light curve is fitted")
        fit_real()
    elif cmd == "fitted":
        load_fitted()
        r10.run_samples("fitted", ["fitted"], workers)
    elif cmd == "fitted-summary":
        fitted_summary()
    elif cmd == "checks":
        checks()
    else:
        raise SystemExit(__doc__)
