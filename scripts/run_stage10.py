"""Stage 10: an optical (ZTF) driver for the W1 and W2 lags (docs/PREDICTIONS_10.md).

    python scripts/run_stage10.py build                 # nightly light curves, criteria, flags, freeze_v3
    python scripts/run_stage10.py ident --workers 3     # section 5: synthetic samples (checkpoint per object)
    python scripts/run_stage10.py ident-summary
    python scripts/run_stage10.py fit                   # section 6: the real light curves, swapped drivers
    python scripts/run_stage10.py fitted --workers 3    # null samples with the fitted lags
    python scripts/run_stage10.py checks                # section 7: sensitivity

The verdicts are made by scripts/stage10_verdict.py. Nothing is fitted to a real pair of optical
and infrared light curves before `ident-summary` has been written.
"""

from __future__ import annotations

import csv
import json
import multiprocessing
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.stats import chi2 as chi2_dist
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import fetch_ztf10 as fz  # noqa: E402
import run_stage9 as r9  # noqa: E402
from irecho import config, hyst8, irsa  # noqa: E402
from irecho import driver10 as d10  # noqa: E402
from irecho import transfer9 as t9  # noqa: E402
from irecho.stats import robust_std  # noqa: E402

RES = config.RESULTS_DIR
CKPT = config.CACHE_DIR / "s10"
FROZEN = config.PROJECT_ROOT / "data" / "frozen" / "ztf_v3"
FREEZE_PATH = config.PROJECT_ROOT / "data" / "freeze_v3.json"
SAMPLE_PATH = RES / "s10_sample.json"
SEED_TAG = 10
DRIVERS = {"fast": 200.0, "slow": 1000.0}
N_IDENT = 50
N_SWAPS = 20
BANDS = ("g", "r")
# section 2 of the plan
MAX_MAGERR = 0.2
BIN_DAYS = 30.0
MIN_BINS = 20
MIN_F30 = 0.03
MIN_F30_SIGMA = 3.0
BIN_ERROR_FLOOR = 0.005
P_IR_VARIABLE = 0.01
FLAG_RHO = 0.3
FLAG_P = 1e-3
RESID_HALF_WINDOW = 15.0
TRUTHS = {
    # name: (kind, shape, lag rule, alpha1, alpha2, tau_re)
    "X0": ("independent", "f1", "lyu", 0.0, 0.0, None),
    "Z0": ("linear", "f1", "zero", 0.0, 0.0, None),
    "A0": ("linear", "f1", "lyu_same", 0.0, 0.0, None),
    "A_20": ("linear", "f1", 20.0, 0.0, 0.0, None),
    "A_50": ("linear", "f1", 50.0, 0.0, 0.0, None),
    "A_100": ("linear", "f1", 100.0, 0.0, 0.0, None),
    "A_200": ("linear", "f1", 200.0, 0.0, 0.0, None),
    "A_300": ("linear", "f1", 300.0, 0.0, 0.0, None),
    "A_lyu": ("linear", "f1", "lyu", 0.0, 0.0, None),
    "A_mandal": ("linear", "f1", "mandal", 0.0, 0.0, None),
    "A_s9": ("linear", "f1", "s9", 0.0, 0.0, None),
    "G_lyu": ("linear", "gwide", "lyu", 0.0, 0.0, None),
    "E_lyu": ("linear", "exp", "lyu", 0.0, 0.0, None),
    "W_thin": ("linear", "f0", "lyu", 0.0, 0.0, None),
    "W_wide": ("linear", "f2", "lyu", 0.0, 0.0, None),
    "B_20": ("linear", "f1", "lyu", 0.2, 0.0, None),
    "B_40": ("linear", "f1", "lyu", 0.4, 0.0, None),
    "B_both": ("linear", "f1", "lyu", 0.3, 0.15, None),
    "C_0": ("boundary", "f1", "lyu", 0.0, 0.0, 0.0),
    "C_1yr": ("boundary", "f1", "lyu", 0.0, 0.0, 1.0),
    "C_3yr": ("boundary", "f1", "lyu", 0.0, 0.0, 3.0),
}
FIXED_DELTA_DAYS = 40.0
FIT_KEYS = ("tau", "tau_q16", "tau_q50", "tau_q84", "tau_at_edge", "r_w", "ds_c1", "ds_01", "ds_12", "ds_1g", "ds_1a",
            "ds_1e", "alpha", "z_alpha", "tau_ma", "eta", "tau_f0", "tau_f05", "tau_f2", "tau_mg", "var_share",
            "gap_share", "b", "chi2_m1")
OFF_KEYS = ("d_colour_obs", "d_colour_pred", "sigma_colour_obs", "usable_obs", "usable_pred", "comparable_obs")
_OBJ: list = []
_FITTED: dict = {}


# ----------------------------------------------------------------------------- section 2: the optical light curves
def wise_objects() -> dict[str, dict]:
    """Frozen WISE light curves of the type 1 objects, with the catalogue values of stage 9."""
    out = {}
    for tag in r9.TAGS:
        for o in r9.load(tag):
            if o["type"] in fz.TYPE1:
                out[o["name"]] = dict(o, sample=tag)
    return out


def read_points(path: Path) -> dict | None:
    """Rules 1-3 of section 2.3: catflags == 0, the object ID with the most such rows, 0 < magerr < 0.2."""
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    good = [r for r in rows if r["catflags"] == "0"]
    if not good:
        return None
    count: dict[str, int] = {}
    for r in good:
        count[r["oid"]] = count.get(r["oid"], 0) + 1
    oid = max(sorted(count), key=lambda k: count[k])
    keep = sorted((r for r in good if r["oid"] == oid and 0.0 < float(r["magerr"]) < MAX_MAGERR), key=lambda r: float(r["mjd"]))
    if not keep:
        return None
    r0 = keep[0]
    return {"n_rows_total": len(rows), "n_rows_clean": len(good), "n_rows_used": len(keep), "oid": oid,
            "key": (int(r0["field"]), int(r0["ccdid"], 0), int(r0["qid"], 0), r0["filtercode"]),
            "mjd": np.array([float(r["mjd"]) for r in keep]), "mag": np.array([float(r["mag"]) for r in keep]),
            "magerr": np.array([float(r["magerr"]) for r in keep]),
            "airmass": np.array([float(r["airmass"]) for r in keep]),
            "expid": np.array([int(r["expid"]) for r in keep])}


def seeing_of(points: dict) -> np.ndarray | None:
    path = fz.meta_path(*points["key"])
    if not path.exists():
        return None
    with open(path, newline="") as fh:
        table = {int(r["expid"]): float(r["seeing"]) for r in csv.DictReader(fh) if r["seeing"] not in ("", "null")}
    return np.array([table.get(int(e), np.nan) for e in points["expid"]])


def local_residuals(mjd, mag, half_window: float = RESID_HALF_WINDOW, min_other: int = 5) -> np.ndarray:
    out = np.full(mjd.size, np.nan)
    lo = np.searchsorted(mjd, mjd - half_window, side="left")
    hi = np.searchsorted(mjd, mjd + half_window, side="right")
    for i in range(mjd.size):
        if hi[i] - lo[i] - 1 >= min_other:
            others = np.concatenate([mag[lo[i]:i], mag[i + 1:hi[i]]])
            out[i] = mag[i] - np.median(others)
    return out


def spurious_flags(points: dict) -> dict:
    """Section 2.5: Spearman correlation of the local magnitude residual with seeing and airmass."""
    res = local_residuals(points["mjd"], points["mag"])
    out = {"n_resid": int(np.isfinite(res).sum())}
    for name, x in (("seeing", seeing_of(points)), ("airmass", points["airmass"])):
        if x is None:
            out[f"rho_{name}"], out[f"p_{name}"], out[f"flag_{name}"] = None, None, None
            continue
        ok = np.isfinite(res) & np.isfinite(x)
        if ok.sum() < 30:
            out[f"rho_{name}"], out[f"p_{name}"], out[f"flag_{name}"] = None, None, None
            continue
        sp = spearmanr(x[ok], res[ok])
        out[f"rho_{name}"], out[f"p_{name}"] = float(sp.statistic), float(sp.pvalue)
        out[f"flag_{name}"] = bool(abs(sp.statistic) >= FLAG_RHO and sp.pvalue < FLAG_P)
        out[f"n_{name}"] = int(ok.sum())
    return out


def binned_variability(t, y) -> dict:
    """Rule (c) of section 2.4: excess fractional rms of 30-day bins."""
    k = np.floor((t - t[0]) / BIN_DAYS).astype(int)
    vals, errs = [], []
    for b in np.unique(k):
        sel = k == b
        if sel.sum() >= 3:
            v = float(np.median(y[sel]))
            vals.append(v)
            errs.append(max(1.2533 * robust_std(y[sel]) / np.sqrt(sel.sum()), BIN_ERROR_FLOOR * v))
    out = {"n_bins": len(vals), "f30": None, "f30_over_sigma": None, "pass": False}
    if len(vals) < 3:
        return out
    vals, errs = np.array(vals), np.array(errs)
    excess = max(float(np.var(vals) - np.mean(errs ** 2)), 0.0)
    f30 = float(np.sqrt(excess) / np.mean(vals))
    ratio = float(np.sqrt(excess) / np.median(errs))
    out.update(f30=f30, f30_over_sigma=ratio,
               **{"pass": bool(len(vals) >= MIN_BINS and f30 >= MIN_F30 and ratio >= MIN_F30_SIGMA)})
    return out


def ir_variable(o: dict, idx) -> dict:
    out = {}
    for b, f, s in (("w1", o["f1"], o["s1"]), ("w2", o["f2"], o["s2"])):
        y, e = f[idx], s[idx]
        w = 1.0 / e ** 2
        chi2 = float(np.sum(w * (y - np.sum(w * y) / w.sum()) ** 2))
        out[f"p_const_{b}"] = float(chi2_dist.sf(chi2, idx.size - 1))
    out["pass"] = bool(out["p_const_w1"] < P_IR_VARIABLE and out["p_const_w2"] < P_IR_VARIABLE)
    return out


def frozen_path(name: str, band: str) -> Path:
    return FROZEN / f"{name}_{band}.csv"


def write_nights(path: Path, t, y, s) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["mjd", "flux_rel", "flux_err"])
        for a, b, c in zip(t, y, s):
            wr.writerow([f"{a:.5f}", f"{b:.6f}", f"{c:.6f}"])
    return irsa.sha256_file(path)


def read_nights(name: str, band: str) -> dict:
    with open(frozen_path(name, band), newline="") as fh:
        rows = list(csv.DictReader(fh))
    return {"t": np.array([float(r["mjd"]) for r in rows]), "y": np.array([float(r["flux_rel"]) for r in rows]),
            "s": np.array([float(r["flux_err"]) for r in rows])}


def build() -> None:
    """Sections 2.3-2.6. Looks at the optical and the infrared light curves separately; no fit of one to the other."""
    wise = wise_objects()
    targets = fz.targets()
    have = fz.have_objects(targets)
    records, freeze = {}, {}
    steps = {"type1": len(targets), "ztf_object": 0, "light_curve": 0, "nights": 0, "visits": 0, "optical_variable": 0,
             "ir_variable": 0}
    t0 = time.time()
    for n, tg in enumerate(targets, start=1):
        name = tg["name"]
        o = wise[name]
        rec = {"sample": tg["sample"], "has_ztf_object": any((name, b) in have for b in BANDS), "bands": {}}
        for band in BANDS:
            path = fz.curve_path(name, band)
            b = {"file": None}
            rec["bands"][band] = b
            if path is None:
                continue
            side = json.loads(path.with_suffix(".csv.meta.json").read_text())
            pts = read_points(path)
            b.update(file=str(path.relative_to(config.CACHE_DIR)), raw_sha256=irsa.sha256_file(path),
                     retrieved_utc=side.get("retrieved_utc"), n_rows_total=side.get("n_rows"))
            if pts is None:
                b["n_nights"] = 0
                continue
            night = d10.nightly(pts["mjd"], pts["mag"], pts["magerr"])
            keep = d10.clip_nights(night["t"], night["y"])
            t, y, s = night["t"][keep], night["y"][keep], night["s"][keep]
            b.update(oid=pts["oid"], n_rows_clean=pts["n_rows_clean"], n_rows_used=pts["n_rows_used"],
                     n_nights=int(t.size), n_nights_clipped=int((~keep).sum()), mjd_first=float(t[0]), mjd_last=float(t[-1]),
                     mag_ref=night["mag_ref"], sha256=write_nights(frozen_path(name, band), t, y, s))
            b.update(spurious_flags(pts))
            idx = np.nonzero(d10.usable_visits(o["mjd"], t))[0]
            var = binned_variability(t, y)
            b.update(n_visits_usable=int(idx.size), longest_gap_days=float(np.diff(t).max()) if t.size > 1 else None,
                     frac_rms_nightly=d10.excess_frac_rms(y, s), **{f"var_{k}": v for k, v in var.items()})
            b["pass_nights"] = bool(t.size >= d10.MIN_NIGHTS)
            b["pass_visits"] = bool(b["pass_nights"] and idx.size >= d10.MIN_VISITS)
            b["pass_optical"] = bool(b["pass_visits"] and var["pass"])
            if b["pass_optical"]:
                b.update(ir_variable(o, idx))
                b["hyper"] = d10.fit_drw(t, y, s)
        order = [band for band in BANDS if rec["bands"][band].get("pass_optical")]
        rec["driver_band"] = order[0] if order else None
        rec["usable"] = bool(rec["driver_band"] and rec["bands"][rec["driver_band"]]["pass"])
        bands = rec["bands"].values()
        steps["ztf_object"] += rec["has_ztf_object"]
        steps["light_curve"] += any(b["file"] for b in bands)
        steps["nights"] += any(b.get("pass_nights") for b in bands)
        steps["visits"] += any(b.get("pass_visits") for b in bands)
        steps["optical_variable"] += bool(order)
        steps["ir_variable"] += rec["usable"]
        records[name] = rec
        if n % 20 == 0 or n == len(targets):
            print(f"{time.time() - t0:6.0f} s  {n}/{len(targets)}  usable so far {steps['ir_variable']}", flush=True)
    usable = [k for k, r in records.items() if r["usable"]]
    flags = {"seeing": 0, "airmass": 0, "seeing_unknown": 0}
    for k in usable:
        b = records[k]["bands"][records[k]["driver_band"]]
        flags["seeing"] += bool(b.get("flag_seeing"))
        flags["airmass"] += bool(b.get("flag_airmass"))
        flags["seeing_unknown"] += b.get("flag_seeing") is None
    summary = {"steps": steps, "n_usable": len(usable), "driver_band": {b: sum(records[k]["driver_band"] == b for k in usable) for b in BANDS},
               "flags_among_usable": flags, "n_type2_excluded": None, "stop": len(usable) < 30,
               "criteria": {"max_magerr": MAX_MAGERR, "min_nights": d10.MIN_NIGHTS, "min_visits": d10.MIN_VISITS,
                            "support_past_days": d10.SUPPORT_PAST_DAYS, "support_future_days": d10.SUPPORT_FUTURE_DAYS,
                            "max_gap_days": d10.MAX_GAP_DAYS, "min_f30": MIN_F30, "min_f30_sigma": MIN_F30_SIGMA,
                            "min_bins": MIN_BINS, "p_ir_variable": P_IR_VARIABLE, "flag_rho": FLAG_RHO, "flag_p": FLAG_P}}
    n_all = {tag: len(r9.load(tag)) for tag in r9.TAGS}
    summary["n_type2_excluded"] = int(sum(n_all.values()) - len(targets))
    SAMPLE_PATH.write_text(json.dumps({"summary": summary, "objects": records}, indent=1, sort_keys=True) + "\n")
    for name, rec in records.items():
        freeze[name] = {"driver_band": rec["driver_band"], "usable": rec["usable"],
                        "bands": {band: {k: b.get(k) for k in ("file", "raw_sha256", "retrieved_utc", "oid", "n_nights", "sha256")}
                                  for band, b in rec["bands"].items() if b["file"]}}
    commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, cwd=config.PROJECT_ROOT).stdout.strip()
    clean = subprocess.run(["git", "status", "--porcelain", "--", "src", "scripts"], capture_output=True, text=True,
                           cwd=config.PROJECT_ROOT).stdout.strip() == ""
    FREEZE_PATH.write_text(json.dumps({"frozen_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                                       "code_commit": commit, "code_tree_clean": clean,
                                       "source": "IRSA ZTF light-curve service (public data releases), radius 3 arcsec",
                                       "n_objects": len(freeze), "n_usable": len(usable), "light_curves": freeze},
                                      indent=1, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=1))


def update_flags() -> None:
    """Section 2.5 again for every light curve, once more exposure metadata (seeing) has arrived.
    Only the flag fields of results/s10_sample.json change; the freeze and the criteria do not."""
    data = json.loads(SAMPLE_PATH.read_text())
    for name, rec in data["objects"].items():
        for band, b in rec["bands"].items():
            if not b.get("n_nights"):
                continue
            pts = read_points(config.CACHE_DIR / b["file"])
            b.update(spurious_flags(pts))
    usable = [k for k, r in data["objects"].items() if r["usable"]]
    flags = {"seeing": 0, "airmass": 0, "seeing_unknown": 0}
    for k in usable:
        b = data["objects"][k]["bands"][data["objects"][k]["driver_band"]]
        flags["seeing"] += bool(b.get("flag_seeing"))
        flags["airmass"] += bool(b.get("flag_airmass"))
        flags["seeing_unknown"] += b.get("flag_seeing") is None
    data["summary"]["flags_among_usable"] = flags
    SAMPLE_PATH.write_text(json.dumps(data, indent=1, sort_keys=True) + "\n")
    print(json.dumps(flags))


# ----------------------------------------------------------------------------- the objects that are fitted
def load() -> list[dict]:
    """Usable objects: WISE visits, the frozen optical nights of the driver band, hyper-parameters,
    and the lags that the truth models need (observed frame)."""
    if _OBJ:
        return _OBJ
    sample = json.loads(SAMPLE_PATH.read_text())["objects"]
    wise = wise_objects()
    s9 = {}
    with open(RES / "s9_objects.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            if r["usable"] == "1":
                s9[r["name"]] = float(r["tau_t1"])
    objs = []
    for name in sorted(sample):
        rec = sample[name]
        if not rec["usable"]:
            continue
        band = rec["driver_band"]
        b = rec["bands"][band]
        o = dict(wise[name], band=band, hyper=b["hyper"], frac_rms_opt=b["frac_rms_nightly"], f30=b["var_f30"],
                 flag_seeing=b.get("flag_seeing"), flag_airmass=b.get("flag_airmass"), tau9=s9.get(name))
        o.update({f"opt_{k}": v for k, v in read_nights(name, band).items()})
        o["other_band"] = next((x for x in BANDS if x != band and rec["bands"][x].get("pass_optical")), None)
        objs.append(o)
    lyu = np.array([d10.lyu_lags_days(o["log_lbol"]) for o in objs if o["log_lbol"] is not None])
    man = np.array([d10.mandal_lags_days(o["log_lbol"]) for o in objs if o["log_lbol"] is not None])
    zs = np.array([o["z"] for o in objs if o["log_lbol"] is not None])
    med = {"lyu": np.median(lyu * (1.0 + zs[:, None]), axis=0), "mandal": np.median(man * (1.0 + zs[:, None]), axis=0)}
    tau9_med = float(np.median([o["tau9"] for o in objs if o["tau9"] is not None]))
    for o in objs:
        has = o["log_lbol"] is not None
        for rel, fn in (("lyu", d10.lyu_lags_days), ("mandal", d10.mandal_lags_days)):
            o[f"{rel}_rest"] = fn(o["log_lbol"]) if has else None
            o[f"{rel}_obs"] = tuple((1.0 + o["z"]) * v for v in o[f"{rel}_rest"]) if has else tuple(float(v) for v in med[rel])
        o["tau9_truth"] = o["tau9"] if o["tau9"] is not None else tau9_med
        o["idx"] = np.nonzero(d10.usable_visits(o["mjd"], o["opt_t"]))[0]
        # amplitudes of the synthetic light curves: the infrared rms over all visits (as in stage 9),
        # the optical rms of the nights without the fitted white noise
        o["frac_rms_ir"] = (d10.excess_frac_rms(o["f1"], o["s1"]), d10.excess_frac_rms(o["f2"], o["s2"]))
        floor = 0.5 * float(np.median(o["opt_s"]))
        o["frac_rms_opt"] = float(np.sqrt(max(o["frac_rms_opt"] ** 2 - o["hyper"]["jitter_var"], floor ** 2)))
    _OBJ.extend(objs)
    return objs


def truth_lags(o: dict, name: str) -> tuple[float, float]:
    rule = TRUTHS[name][2] if name != "fitted" else "fitted"
    if rule == "fitted":
        return _FITTED[o["name"]]
    if rule == "zero":
        return 0.0, 0.0
    if rule == "lyu":
        return o["lyu_obs"]
    if rule == "mandal":
        return o["mandal_obs"]
    if rule == "lyu_same":
        return o["lyu_obs"][0], o["lyu_obs"][0]
    if rule == "s9":
        return o["lyu_obs"][0], max(o["lyu_obs"][0] + o["tau9_truth"], 0.0)
    return float(rule), float(rule) + FIXED_DELTA_DAYS


def slim(fit: dict) -> list[float]:
    if not fit.get("usable"):
        return [np.nan] * len(FIT_KEYS)
    vals = dict(fit)
    vals["ds_1g"] = fit["s_m1"] - fit.get("s_mg", np.nan)
    return [float(vals[k]) for k in FIT_KEYS]


def fit_pair(o: dict, prep: dict, y_opt, f1, f2, full: bool = True) -> tuple[list, list, list, dict, dict]:
    rising = hyst8.branch_rising(f1, f2)
    a = d10.fit_band(prep, y_opt, f1, o["s1"], rising, full=full)
    b = d10.fit_band(prep, y_opt, f2, o["s2"], rising, full=full)
    off = d10.colour_offsets(f1, o["s1"], f2, o["s2"], o["r"], a, b, rising)
    return slim(a) if full else None, slim(b) if full else None, [float(off[k]) for k in OFF_KEYS], a, b


# ----------------------------------------------------------------------------- section 5: synthetic samples
def synth_cases(o: dict, k: int, driver: str, index: int, truths: list[str]) -> dict:
    """Noisy optical nights and W1, W2 visits of every truth model for one object and one realisation.
    The driver and the noise are shared between the truth models."""
    seed = [config.SEED, SEED_TAG, list(DRIVERS).index(driver), index, k]
    grid = d10.truth_grid(o["opt_t"], o["mjd"])
    mags = t9.drw_magnitudes(np.random.default_rng(seed + [0]), grid.size, d10.TRUTH_STEP_DAYS, DRIVERS[driver])
    ell = 10.0 ** (-0.4 * d10.lognormal_amp(o["frac_rms_opt"]) * mags)
    rng = np.random.default_rng(seed + [1])
    z_opt = rng.normal(size=o["opt_t"].size) * np.sqrt(o["opt_s"] ** 2 + o["hyper"]["jitter_var"])
    med1, med2 = float(np.median(o["f1"])), float(np.median(o["f2"]))
    rms1, rms2 = o["frac_rms_ir"]
    out = {}

    def optical(e, host: float = 0.0):
        y = np.interp(o["opt_t"], grid, e)
        y = host + (1.0 - host) * y / float(np.mean(y))
        return y / float(np.median(y)) + z_opt

    for name in truths:
        kind, shape, _, a1, a2, tau_re = TRUTHS[name] if name != "fitted" else ("linear", "f1", None, 0.0, 0.0, None)
        lag1, lag2 = truth_lags(o, name)
        if kind == "boundary":
            e, w1, w2 = d10.boundary_truth(mags, grid, o["mjd"], rms2, lag1, tau_re)
            at = np.interp(o["opt_t"], grid, e)
            cv = float(np.std(at) / np.mean(at))
            host = float(np.clip(1.0 - o["frac_rms_opt"] / cv, 0.0, 0.95)) if cv > 0 else 0.0
            y = optical(e, host)
            v1, v2 = np.interp(o["mjd"], grid, w1), np.interp(o["mjd"], grid, w2)
            t1, t2 = v1 * med1 / float(np.median(v1)), v2 * med2 / float(np.median(v2))
        else:
            src = ell
            if kind == "independent":
                m2 = t9.drw_magnitudes(np.random.default_rng(seed + [3]), grid.size, d10.TRUTH_STEP_DAYS, DRIVERS[driver])
                src = 10.0 ** (-0.4 * d10.lognormal_amp(o["frac_rms_opt"]) * m2)
            d1, d2 = d10.linear_truth(src, lag1, lag2, shape, a1, a2)
            y = optical(ell)
            t1 = d10.scale_band(np.interp(o["mjd"], grid, d1), rms1, med1)
            t2 = d10.scale_band(np.interp(o["mjd"], grid, d2), rms2, med2)
        f1, f2 = hyst8.add_noise(np.random.default_rng(seed + [2]), t1, o["s1"], t2, o["s2"], o["r"])
        out[name] = (y, f1, f2, lag1, lag2)
    return out


def object_samples(args) -> str:
    """All realisations of all truth models for one object and one driver; one checkpoint file."""
    kind, k, driver, truths, n = args
    objs = load()
    o = objs[k]
    path = CKPT / kind / f"{driver}_{o['name']}.npz"
    prep = d10.prepare(o["opt_t"], o["opt_s"], o["hyper"], o["mjd"], joint=(kind == "ident"))
    fits = np.full((len(truths), n, 2, len(FIT_KEYS)), np.nan, dtype=np.float32)
    offs = np.full((len(truths), n, len(OFF_KEYS)), np.nan, dtype=np.float32)
    lags = np.zeros((len(truths), 2), dtype=np.float32)
    joint = np.full((n, 5), np.nan, dtype=np.float32)
    for index in range(n):
        cases = synth_cases(o, k, driver, index, truths)
        for j, name in enumerate(truths):
            y, f1, f2, lag1, lag2 = cases[name]
            lags[j] = lag1, lag2
            if not (np.all(f1 > 0) and np.all(f2 > 0)):
                continue
            a, b, off, fa, fb = fit_pair(o, prep, y, f1, f2)
            fits[j, index, 0], fits[j, index, 1], offs[j, index] = a, b, off
            if kind == "ident" and name == "A_lyu" and driver == "fast" and index < 10:
                jt = d10.fit_joint(prep, y, f1, o["s1"], f2, o["s2"], o["r"], fa, fb)
                joint[index] = jt["dtau_q16"], jt["dtau_q50"], jt["dtau_q84"], jt["dtau_q025"], jt["dtau_q975"]
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, fits=fits, offs=offs, lags=lags, joint=joint, truths=np.array(truths))
    return path.name


def run_samples(kind: str, truths: list[str], workers: int, n: int = N_IDENT) -> None:
    objs = load()
    todo = [(kind, k, driver, truths, n) for driver in DRIVERS for k, o in enumerate(objs)
            if not (CKPT / kind / f"{driver}_{o['name']}.npz").exists()]
    print(f"{kind}: {len(todo)} object files to make ({len(objs)} objects, {len(truths)} truth models, {n} realisations)", flush=True)
    t0 = time.time()
    with multiprocessing.get_context("fork").Pool(workers, maxtasksperchild=4) as pool:
        for done, name in enumerate(pool.imap_unordered(object_samples, todo), start=1):
            if done % 5 == 0 or done == len(todo):
                print(f"{time.time() - t0:7.0f} s  {done}/{len(todo)}  {name}", flush=True)


def read_samples(kind: str, driver: str) -> dict:
    """Arrays over (object, truth, realisation, band, key)."""
    objs = load()
    packs = [np.load(CKPT / kind / f"{driver}_{o['name']}.npz") for o in objs]
    return {"truths": [str(x) for x in packs[0]["truths"]], "fits": np.stack([p["fits"] for p in packs]),
            "offs": np.stack([p["offs"] for p in packs]), "lags": np.stack([p["lags"] for p in packs]),
            "joint": np.stack([p["joint"] for p in packs])}


def main() -> None:
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    workers = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 2
    if cmd == "build":
        build()
    elif cmd == "flags":
        update_flags()
    elif cmd == "ident":
        run_samples("ident", list(TRUTHS), workers)
    else:
        import stage10_analysis as an

        an.dispatch(cmd, workers)


if __name__ == "__main__":
    main()
