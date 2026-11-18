"""Stage 7b: post hoc supplement to the X-ray analysis (docs/XRAY_PLAN.md section 9).

Everything here was added AFTER the stage 7 results were seen.

    python scripts/xray_7b.py light     -> results/s7b_light.json   (no long fits)
    python scripts/xray_7b.py fits      -> results/s7b_fits.json    (checkpointed, one step at a time)

Run inside the HEASoft environment (source scripts/heasoft_env.sh).
Logs and .xcm files of the final fits go to results/s7b_xcm/.
"""

from __future__ import annotations

import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from astropy.io import fits as pf

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import xray_fit as xf  # noqa: E402
from irecho import config, xray  # noqa: E402

xspec = xf.xspec
RES = config.RESULTS_DIR
WORK = xf.WORK
FITDIR = xf.FITDIR
XCM = RES / "s7b_xcm"
OBSID = {"xrt10": "00040946001", "xrt22": "00081051001"}
TIME_LIMIT_S = 3600.0


# ----------------------------------------------------------------------------- A1: BAT common band
def bat_band(nh_gal: float) -> dict:
    with pf.open(FITDIR / "bat.pha") as h:
        eb = h["EBOUNDS"].data
        e_lo, e_hi = [float(x) for x in eb["E_MIN"]], [float(x) for x in eb["E_MAX"]]
    chans = xray.bat_channels_in_band(e_lo, e_hi, 14.0, 50.0)
    out = {"ebounds_kev": [[a, b] for a, b in zip(e_lo, e_hi)], "channels_used": chans,
           "channel_bounds_used_kev": [[e_lo[c - 1], e_hi[c - 1]] for c in chans]}

    def bat_fit(gamma_fixed: float | None, channels: list[int]) -> dict:
        xf.load(["bat"])
        s = xspec.AllData(1)
        drop = [c for c in range(1, len(e_lo) + 1) if c not in channels]
        if drop:
            s.ignore(",".join(str(c) for c in drop))        # integers are channel numbers
        xspec.AllModels.clear()
        m = xspec.Model("cflux*powerlaw")
        m.cflux.Emin.values, m.cflux.Emax.values = [14.0, -1], [50.0, -1]
        m.cflux.lg10Flux.values = [-11.4, 0.01, -14, -14, -9, -9]
        m.powerlaw.PhoIndex.values = [1.9, 0.01, -1.0, -1.0, 5.0, 5.0] if gamma_fixed is None else [gamma_fixed, -1]
        m.powerlaw.norm.values = [1.0, -1]
        stat = xf.fit(3)
        rec = {"channels_noticed": [int(c) for c in s.noticed], "n_free": 2 if gamma_fixed is None else 1,
               "chi2": stat, "dof": int(xspec.Fit.dof), "gamma": m.powerlaw.PhoIndex.values[0],
               "gamma_fixed": gamma_fixed is not None, "log_flux_14_50": m.cflux.lg10Flux.values[0],
               "log_flux_14_50_90": xf.interval_of(m.cflux.lg10Flux)}
        if gamma_fixed is None:
            rec["gamma_90"] = xf.interval_of(m.powerlaw.PhoIndex)
        return rec

    # NuSTAR 14-50 keV (as in stage 7; the two channels that contain 14.0 and 50.0 keV are dropped, 40 eV each)
    over = {k: {"ignore": "**-14.0 50.0-**"} for k in ("fpma", "fpmb")}
    xf.load(["fpma", "fpmb"], over)
    xspec.AllModels.clear()
    m = xspec.Model("constant*cflux*powerlaw")
    m.constant.factor.values = [1.0, -1]
    m.cflux.Emin.values, m.cflux.Emax.values = [14.0, -1], [50.0, -1]
    m.cflux.lg10Flux.values = [-10.8, 0.01, -14, -14, -9, -9]
    m.powerlaw.PhoIndex.values = [1.6, 0.01, 0.5, 0.5, 3.5, 3.5]
    m.powerlaw.norm.values = [1.0, -1]
    c = xspec.AllModels(2).constant.factor
    c.link = ""
    c.values = [1.0, 0.01, 0.5, 0.5, 2.0, 2.0]
    stat = xf.fit(3)
    nu = {"cstat": stat, "n_bins": int(sum(len(xspec.AllData(i).noticed) for i in (1, 2))), "n_free": 3,
          "gamma": m.powerlaw.PhoIndex.values[0], "gamma_90": xf.interval_of(m.powerlaw.PhoIndex),
          "log_flux_14_50": m.cflux.lg10Flux.values[0], "log_flux_14_50_90": xf.interval_of(m.cflux.lg10Flux),
          "const_fpmb": c.values[0]}
    out["nustar_14_50"] = nu

    def ratio(bat: dict) -> dict:
        d = nu["log_flux_14_50"] - bat["log_flux_14_50"]
        lo = math.hypot(nu["log_flux_14_50"] - nu["log_flux_14_50_90"][0], bat["log_flux_14_50_90"][1] - bat["log_flux_14_50"])
        hi = math.hypot(nu["log_flux_14_50_90"][1] - nu["log_flux_14_50"], bat["log_flux_14_50"] - bat["log_flux_14_50_90"][0])
        return {"log_ratio": d, "lo": d - lo, "hi": d + hi, "ratio": 10 ** d,
                "note": "90% intervals of the two fluxes combined in quadrature; statistical only"}

    out["bat_ch1_4_free_gamma"] = bat_fit(None, chans)
    out["bat_ch1_4_gamma_fixed_to_nustar"] = bat_fit(nu["gamma"], chans)
    out["bat_all_channels_integrated_14_50"] = bat_fit(None, list(range(1, len(e_lo) + 1)))   # what stage 7 did
    for key in ("bat_ch1_4_free_gamma", "bat_ch1_4_gamma_fixed_to_nustar", "bat_all_channels_integrated_14_50"):
        out[key]["nustar_over_bat"] = ratio(out[key])
    # joint fit: one power law shape for NuSTAR (14-50 keV) and BAT channels 1-4, one constant on BAT
    xf.load(["fpma", "fpmb", "bat"], over)
    s = xspec.AllData(3)
    s.ignore(",".join(str(c_) for c_ in range(1, len(e_lo) + 1) if c_ not in chans))
    xspec.AllModels.clear()
    m = xspec.Model("constant*cflux*powerlaw")
    m.constant.factor.values = [1.0, -1]
    m.cflux.Emin.values, m.cflux.Emax.values = [14.0, -1], [50.0, -1]
    m.cflux.lg10Flux.values = [nu["log_flux_14_50"], 0.01, -14, -14, -9, -9]
    m.powerlaw.PhoIndex.values = [nu["gamma"], 0.01, 0.5, 0.5, 3.5, 3.5]
    m.powerlaw.norm.values = [1.0, -1]
    for g, start in ((2, 1.0), (3, 0.3)):
        cc = xspec.AllModels(g).constant.factor
        cc.link = ""
        cc.values = [start, 0.01, 0.01, 0.01, 20.0, 20.0]
    stat = xf.fit(3)
    k = xspec.AllModels(3).constant.factor
    idx = xspec.AllModels(3).startParIndex + k.index - 1
    xspec.Fit.error(f"maximum 50 {xf.DELTA_90} {idx}")
    lo, hi, code = k.error
    out["joint_same_shape_14_50"] = {
        "statistic": stat, "bat_channels_noticed": [int(c_) for c_ in s.noticed], "gamma": m.powerlaw.PhoIndex.values[0],
        "k_bat": k.values[0], "k_bat_90": [float(lo), float(hi), code], "log_ratio": -math.log10(k.values[0]),
        "lo": -math.log10(hi), "hi": -math.log10(lo),
        "note": "NuSTAR 14-50 keV and BAT channels 1-4 with one photon index; assumes the same shape in this band"}
    return out


# ----------------------------------------------------------------------------- A2: XRT scale
def run_tool(args: list[str], cwd: Path, log: str) -> str:
    res = subprocess.run(args, cwd=cwd, capture_output=True, text=True)
    (cwd / log).write_text(res.stdout + res.stderr)
    return res.stdout + res.stderr


def xrt_corrections() -> dict:
    """What xrtmkarf applied: ARF with PSF and exposure-map corrections against ARFs without them."""
    reg = json.loads((RES / "s7_regions.json").read_text())
    out = {}
    for key, obsid in OBSID.items():
        d = WORK / key / "prod"
        expo = f"../pipe/sw{obsid}xpcw3po_ex.img"
        variants = {"used_psf_and_expomap": ["psfflag=yes", f"expofile={expo}"],
                    "psf_only": ["psfflag=yes", "expofile=NONE"],
                    "no_psf_no_expomap": ["psfflag=no", "expofile=NONE"]}
        arfs = {}
        logs = {}
        for name, opts in variants.items():
            outfile = f"arf7b_{name}.arf"
            logs[name] = run_tool(["xrtmkarf", "phafile=src.pha", "srcx=-1", "srcy=-1", f"outfile={outfile}",
                                   "clobber=yes", *opts], d, f"arf7b_{name}.log")
            with pf.open(d / outfile) as h:
                t = h["SPECRESP"].data
                arfs[name] = (0.5 * (t["ENERG_LO"] + t["ENERG_HI"]), np.array(t["SPECRESP"], dtype=float))
        with pf.open(d / "src.arf") as h:
            used = np.array(h["SPECRESP"].data["SPECRESP"], dtype=float)
        e, a_used = arfs["used_psf_and_expomap"]
        rec = {"reproduces_the_arf_used_in_stage_7": bool(np.allclose(a_used, used, rtol=1e-5))}
        m = re.search(r"Off-axis Angle \(arcmin\):\s*([0-9.]+)", logs["used_psf_and_expomap"])
        rec["off_axis_arcmin"] = float(m.group(1)) if m else None
        energies = (0.5, 1.0, 1.5, 3.0, 5.0, 7.0, 9.0)
        rec["energies_kev"] = list(energies)
        rec["factor_psf_and_expomap"] = [float(np.interp(x, e, a_used / arfs["no_psf_no_expomap"][1])) for x in energies]
        rec["factor_psf_only"] = [float(np.interp(x, e, arfs["psf_only"][1] / arfs["no_psf_no_expomap"][1])) for x in energies]
        rec["factor_expomap_only"] = [float(np.interp(x, e, a_used / arfs["psf_only"][1])) for x in energies]
        rec["area_cm2_no_correction"] = [float(np.interp(x, e, arfs["no_psf_no_expomap"][1])) for x in energies]
        # exposure map at the source position (sky image; values in seconds)
        with pf.open(WORK / key / "pipe" / f"sw{obsid}xpcw3po_ex.img") as h:
            img = np.array(h[0].data, dtype=float)
        cx, cy = reg[key]["centroid_xy"]
        yy, xx = np.indices(img.shape)
        inreg = np.hypot(xx + 1 - cx, yy + 1 - cy) <= 20.0
        rec["exposure_map"] = {"max_s": float(img.max()), "at_source_pixel_s": float(img[int(round(cy)) - 1, int(round(cx)) - 1]),
                               "mean_in_region_s": float(img[inreg].mean()), "min_in_region_s": float(img[inreg].min()),
                               "mean_in_region_over_max": float(img[inreg].mean() / img.max())}
        # bad columns: detector position of the source in every good-time interval
        with pf.open(WORK / key / "pipe" / f"sw{obsid}xpcw3po_cl.evt") as h:
            ev, gti, bp = h["EVENTS"].data, h["GTI"].data, h["BADPIX"].data
            src = np.hypot(ev["X"] - cx, ev["Y"] - cy) < 20.0
            cols = sorted({int(x) for x in bp["RAWX"][bp["YEXTENT"] > 100]})
            rows = []
            for a, b in zip(gti["START"], gti["STOP"]):
                m_ = src & (ev["TIME"] >= a) & (ev["TIME"] <= b)
                if m_.sum() < 3:
                    continue
                rx = float(np.median(ev["RAWX"][m_]))
                rows.append({"duration_s": float(b - a), "n_source_events": int(m_.sum()), "rawx_median": rx,
                             "rawy_median": float(np.median(ev["RAWY"][m_])),
                             "nearest_bad_column_pix": float(min(abs(rx - c) for c in cols))})
        rec["bad_columns_rawx"] = cols
        rec["per_gti"] = rows
        rec["nearest_bad_column_pix_min"] = min(r["nearest_bad_column_pix"] for r in rows)
        rec["n_gti_with_bad_column_inside_region"] = int(sum(r["nearest_bad_column_pix"] <= 20 for r in rows))
        # count-weighted correction factor in the two bands (source region events)
        with pf.open(WORK / key / "pipe" / f"sw{obsid}xpcw3po_cl.evt") as h:
            ev = h["EVENTS"].data
            src = np.hypot(ev["X"] - cx, ev["Y"] - cy) < 20.0
            en = ev["PI"][src] * 0.01
        f_all = a_used / arfs["no_psf_no_expomap"][1]
        rec["count_weighted_factor"] = {}
        for band, (lo, hi) in (("0.3-10", (0.3, 10.0)), ("2-10", (2.0, 10.0))):
            ee = en[(en >= lo) & (en <= hi)]
            rec["count_weighted_factor"][band] = float(ee.size / np.sum(1.0 / np.interp(ee, e, f_all)))
        out[key] = rec
    # comparison with 2SXPS (2010 observation only)
    rates = json.loads((RES / "s7_summary.json").read_text())["xrt_count_rates"]["xrt10"]
    cat = {"source": "2SXPS J092937.8+623239 (VizieR IX/58, queried 2026-11-17)", "n_obs": 1, "exposure_s": 8933,
           "first_obs": "2010-12-29T00:20:27", "last_obs": "2010-12-29T16:23:50", "off_axis_arcmin": 1.32,
           "cr_0p3_10": [8.93e-3, 1.22e-3], "cr_0p3_1": [8.74e-4, 4.00e-4], "cr_1_2": [5.89e-4, 3.16e-4],
           "cr_2_10": [7.20e-3, 1.09e-3]}
    cmp_ = {"catalogue": cat}
    for band, key_rate, cat_key in (("0.3-10", None, "cr_0p3_10"), ("2-10", "net_rate_2_10", "cr_2_10")):
        if key_rate is None:
            raw = rates["net_counts_0p3_10"] / rates["exposure_s"]
            err = math.hypot(rates["net_rate_0p3_2_err"], rates["net_rate_2_10_err"])
        else:
            raw, err = rates[key_rate], rates[key_rate + "_err"]
        f = out["xrt10"]["count_weighted_factor"][band]
        corr, corr_err = raw / f, err / f
        cv, ce = cat[cat_key]
        cmp_[band] = {"own_net_rate_raw": raw, "own_net_rate_raw_err": err, "correction_factor": f,
                      "own_rate_corrected": corr, "own_rate_corrected_err": corr_err, "catalogue_rate": cv,
                      "catalogue_err": ce, "difference_sigma": (corr - cv) / math.hypot(corr_err, ce)}
    cmp_["exposure_ratio_catalogue_over_nominal"] = cat["exposure_s"] / rates["exposure_s"]
    out["comparison_with_2sxps"] = cmp_
    return out


def gti_overlap() -> dict:
    """Seconds of the 2022 XRT good-time intervals that fall inside NuSTAR good-time intervals."""
    def gtis(path):
        with pf.open(path) as h:
            hd = h["GTI"].header
            ref = hd["MJDREFI"] + hd["MJDREFF"]
            tz = float(hd.get("TIMEZERO", 0.0))
            return [((ref + (a + tz) / 86400.0), (ref + (b + tz) / 86400.0)) for a, b in zip(h["GTI"].data["START"], h["GTI"].data["STOP"])]

    xrt = gtis(WORK / "xrt22" / "pipe" / "sw00081051001xpcw3po_cl.evt")
    out = {"xrt22_gti_total_s": float(sum(b - a for a, b in xrt) * 86400.0), "n_xrt_gti": len(xrt)}
    for mod in ("A", "B"):
        nu = gtis(WORK / "nustar" / "pipe" / f"nu60160365002{mod}01_cl.evt")
        out[f"fpm{mod.lower()}_gti_total_s"] = float(sum(b - a for a, b in nu) * 86400.0)
        out[f"overlap_with_fpm{mod.lower()}_s"] = xray.gti_overlap_seconds(xrt, nu) * 86400.0
    out["xrt_span_mjd"] = [xrt[0][0], xrt[-1][1]]
    out["note"] = "MJD from each file's MJDREF; the Swift and NuSTAR clocks are compared at the level of seconds only"
    return out


# ----------------------------------------------------------------------------- A3: chain, sims, hardness
def mcmc_diagnostics() -> dict:
    path = FITDIR / "chain_c2.fits"
    if not path.exists():
        return {"status": "chain file not present; not checked"}
    with pf.open(path) as h:
        d = h[1].data
        nw = int(h[1].header.get("NWALKERS", 20))
        cols = h[1].columns.names
        arr = {c: np.array(d[c], dtype=float) for c in cols}
    n = (len(arr[cols[0]]) // nw) * nw
    d2 = (arr["lg10Lum__8"] - arr["lg10Lum__41"])[:n].reshape(-1, nw)
    dn = np.log10(arr["nH__3"] / arr["nH__36"])[:n].reshape(-1, nw)
    out = {"file": "data/cache/s7/work/fit/chain_c2.fits", "n_rows": int(n), "n_walkers": nw, "n_steps_per_walker": int(n // nw),
           "layout_assumed": "rows cycle through the walkers (row i belongs to walker i mod n_walkers)"}
    for name, x in (("delta2", d2), ("delta_nh", dn)):
        tau = xray.autocorr_time(x)
        out[name] = {"autocorr_time_steps": tau, "effective_sample_size": float(n / tau), "gelman_rubin_split": xray.gelman_rubin(x),
                     "walker_means_min_max": [float(x.mean(axis=0).min()), float(x.mean(axis=0).max())],
                     "p05_p50_p95": [float(v) for v in np.percentile(x, [5, 50, 95])],
                     "first_half_p05_p95": [float(v) for v in np.percentile(x[: x.shape[0] // 2], [5, 95])],
                     "second_half_p05_p95": [float(v) for v in np.percentile(x[x.shape[0] // 2:], [5, 95])]}
    # the alternative layout (walker blocks) as a sanity check on the assumption
    alt = d2.reshape(-1).reshape(nw, -1).T
    out["delta2_autocorr_time_if_rows_were_walker_blocks"] = xray.autocorr_time(alt)
    return out


def sims_record() -> dict:
    s = json.loads((RES / "s7_sim_validation.json").read_text())
    rows = s["rows"]
    out = {"n_sims": s["n_sims"], "scope": "spectra made with M0 and fitted with M0; XRT relative constant fixed at 1; "
           "2022 XRT has about 150 net counts in the synthetic spectra (107 in the data)",
           "coverage_90": s["coverage_90"], "binomial_se_at_90": xray.binomial_se(0.9, s["n_sims"])}
    out["coverage_minus_nominal_in_sigma"] = {k: (v - 0.9) / out["binomial_se_at_90"] for k, v in s["coverage_90"].items()}
    out["n_nh10_at_limit"] = int(sum(r["nh10"] <= 0.011 or r["nh10"] >= 499 for r in rows))
    out["n_nh22_at_limit"] = int(sum(r["nh"] <= 0.011 or r["nh"] >= 499 for r in rows))
    out["n_gamma_at_limit"] = int(sum(r["gamma"] <= 1.001 or r["gamma"] >= 2.999 for r in rows))
    out["n_negative_dstat_below_minus_0p05"] = {k: int(sum(r[f"dstat_{k}"] < -0.05 for r in rows))
                                                for k in ("nh", "gamma", "delta2", "delta_nh")}
    out["n_non_finite"] = int(sum(not all(math.isfinite(r[k]) for k in ("stat", "nh", "gamma", "delta2", "delta_nh")) for r in rows))
    out["xrt22_net_counts_median"] = float(np.median([r["net_counts"][0] for r in rows]))
    out["xrt10_net_counts_median"] = float(np.median([r["net_counts"][3] for r in rows]))
    out["not_recorded"] = "the 90% intervals themselves were not computed per realisation (coverage was tested by fixing the truth)"
    return out


def hardness() -> dict:
    reg = json.loads((RES / "s7_regions.json").read_text())
    arf = json.loads((RES / "s7_arf.json").read_text())
    rng = np.random.default_rng([config.SEED, 72])
    out = {"method": "posterior samples of the net counts with Poisson source and background (flat priors, net >= 0); "
                     "HR = (H - S) / (H + S) from the samples"}
    for key in ("xrt10", "xrt22"):
        ratio = arf[key]["backscal_source"] / arf[key]["backscal_background"]
        s = xray.poisson_net_samples(reg[key]["source_counts_0p3_2"], reg[key]["background_counts_0p3_2"], ratio, 400000, rng)
        h = xray.poisson_net_samples(reg[key]["source_counts_2_10"], reg[key]["background_counts_2_10"], ratio, 400000, rng)
        hr = (h - s) / (h + s)
        out[key] = {"soft_counts_source_region": reg[key]["source_counts_0p3_2"],
                    "soft_background_expected": reg[key]["background_counts_0p3_2"] * ratio,
                    "hard_counts_source_region": reg[key]["source_counts_2_10"],
                    "hard_background_expected": reg[key]["background_counts_2_10"] * ratio,
                    "soft_net_p05_p50_p95": [float(v) for v in np.percentile(s, [5, 50, 95])],
                    "soft_net_upper_limit_90": float(np.percentile(s, 90)),
                    "hr_p05_p50_p95": [float(v) for v in np.percentile(hr, [5, 50, 95])],
                    "hr_lower_limit_90": float(np.percentile(hr, 10))}
    a, b = out["xrt10"], out["xrt22"]
    out["note"] = "2022: the soft band is not detected, so only a lower limit on HR is given"
    out["background_fraction_in_source_region"] = {
        "fpma": reg["fpma"]["expected_background_in_source"] / reg["fpma"]["source_counts"],
        "fpmb": reg["fpmb"]["expected_background_in_source"] / reg["fpmb"]["source_counts"],
        "xrt10_0p3_10": (a["soft_background_expected"] + a["hard_background_expected"]) / (a["soft_counts_source_region"] + a["hard_counts_source_region"]),
        "xrt22_0p3_10": (b["soft_background_expected"] + b["hard_background_expected"]) / (b["soft_counts_source_region"] + b["hard_counts_source_region"])}
    return out


def run_light() -> None:
    xf.setup()
    xf.prepare()
    nh_gal = xf.galactic_nh()
    out = {"note": "stage 7b: computed after the stage 7 results were seen (docs/XRAY_PLAN.md section 9)"}
    out["bat_common_band"] = bat_band(nh_gal)
    print("A1 done", flush=True)
    out["xrt_corrections"] = xrt_corrections()
    out["gti_overlap"] = gti_overlap()
    out["mcmc_diagnostics"] = mcmc_diagnostics()
    out["synthetic_validation_record"] = sims_record()
    out["hardness_ratio"] = hardness()
    (RES / "s7b_light.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    print("written s7b_light.json", flush=True)


# ----------------------------------------------------------------------------- fits (checkpointed)
def start_from(values: dict) -> dict:
    return {"nh": values["nh_1e22"], "gamma": values["gamma"], "lum": values["log_lint_2_10"],
            "lum_s": values.get("log_lscat_2_10", 40.5), "refl_norm": values.get("refl_norm", 1e-3),
            "line_norm": values.get("line_norm", 3e-6), "nh10": values.get("nh10_1e22", values["nh_1e22"]),
            "lum10": values.get("log_lint10_2_10", values["log_lint_2_10"])}


def pegged(model: str, keys) -> list[str]:
    """Free parameters that sit on a hard limit."""
    hits = []
    for g in range(1, len(keys) + 1):
        m = xspec.AllModels(g)
        for i in range(1, m.nParameters + 1):
            p = m(i)
            if p.frozen or p.link:
                continue
            v, _, lo, _, _, hi = p.values
            if abs(v - lo) <= 1e-6 * max(1.0, abs(lo)) or abs(v - hi) <= 1e-6 * max(1.0, abs(hi)):
                hits.append(f"group {g} {p.name} = {v:g}")
    return hits


def save_xcm(tag: str, keys) -> dict:
    """Save the current fit as .xcm and a 'show all' log; return what was saved."""
    XCM.mkdir(parents=True, exist_ok=True)
    xcm = FITDIR / f"{tag}.xcm"
    if xcm.exists():
        xcm.unlink()
    xspec.Xset.save(str(xcm), info="a")
    log = FITDIR / f"{tag}.log"
    if log.exists():
        log.unlink()
    xspec.Xset.openLog(str(log))
    xspec.Xset.logChatter = 10
    xspec.AllData.show()
    xspec.AllModels.show()
    xspec.Fit.show()
    xspec.Xset.logChatter = 0
    xspec.Xset.closeLog()
    text = xcm.read_text()
    text = "\n".join(ln for ln in text.splitlines() if not ln.startswith("cd ")) + "\n"
    if "energies" not in text:
        text = text.replace("method ", f"energies {xf.ENERGY_GRID}\nmethod ", 1)
    header = (f"# {tag}: saved by scripts/xray_7b.py. Run from a directory that holds the spectra and responses\n"
              f"# (data/cache/s7/work/fit/, made by scripts/xray_fit.py prepare). Model energies: {xf.ENERGY_GRID}.\n")
    (XCM / f"{tag}.xcm").write_text(header + text)
    logtext = log.read_text().replace(str(FITDIR) + "/", "").replace(str(config.PROJECT_ROOT), "<project>")
    (XCM / f"{tag}.log").write_text(logtext)
    return {"xcm": f"results/s7b_xcm/{tag}.xcm", "log": f"results/s7b_xcm/{tag}.log", "statistic": float(xspec.Fit.statistic),
            "spectra": list(keys), "statistics": [xf.SPECTRA[k]["stat"] for k in keys],
            "bands": [xf.SPECTRA[k]["ignore"] for k in keys], "files": [xf.SPECTRA[k]["pha"] for k in keys]}


def reload_in_new_session(tag: str) -> dict:
    """Start the xspec program in a new process, load the .xcm, refit and read the result back."""
    script = FITDIR / f"reload_{tag}.tcl"
    script.write_text(f"@{tag}.xcm\nquery yes\nfit 100\nlog reload_{tag}.log\nshow fit\nshow free\nlog none\nexit\n")
    shutil.copyfile(XCM / f"{tag}.xcm", FITDIR / f"{tag}.xcm")
    res = subprocess.run(["xspec", "-", script.name], cwd=FITDIR, capture_output=True, text=True, timeout=1200)
    text = (FITDIR / f"reload_{tag}.log").read_text() if (FITDIR / f"reload_{tag}.log").exists() else res.stdout
    (XCM / f"reload_{tag}.log").write_text(text.replace(str(FITDIR) + "/", ""))
    m = re.search(r"Total fit statistic\s+([-+0-9.eE]+)", text)
    pars = {}
    for line in text.splitlines():
        mm = re.match(r"#?\s*(\d+)\s+(\d+)\s+(\w+)\s+(\w+)\s+(?:\S+\s+)?([-+0-9.eE]+)\s+\+/-", line)
        if mm:
            pars[f"{mm.group(3)}.{mm.group(4)}#{mm.group(1)}"] = float(mm.group(5))
    return {"statistic": float(m.group(1)) if m else None, "free_parameters": pars, "returncode": res.returncode}


def c2_with_xrt_constant(mode: str, nh_gal: float, start: dict) -> dict:
    """M0, 2022 against 2010, with a free XRT constant.

    mode 'both': one constant shared by the two XRT spectra (an instrument offset).
    mode '2022': a constant on the 2022 XRT spectrum only (a problem of that observation).
    """
    keys = ["xrt22", "fpma", "fpmb", "xrt10"]
    xf.load(keys)
    xf.build(keys, "M0", nh_gal, start)
    c22 = xf.par(1, "const")
    c22.values = [0.7, 0.01, 0.2, 0.2, 3.0, 3.0]
    c22.frozen = False
    xf.par(2, "const").link = ""
    xf.par(2, "const").values = [1.0, -1]              # FPMA is the reference
    c10 = xf.par(4, "const")
    if mode == "both":
        c10.link = f"p{c22.index}"
    else:
        c10.link = ""
        c10.values = [1.0, -1]
    best = xf.fit(3)
    v = xf.read_values(keys, "M0")
    res = {"mode": mode, "statistic": best, "n_free": xf.n_free(), "values": v, "const_xrt22": c22.values[0],
           "const_xrt22_90": xf.interval(1, "const"), "nh22_90": xf.interval(1, "nh"), "nh10_90": xf.interval(4, "nh"),
           "lum22_90": xf.interval(1, "lum"), "lum10_90": xf.interval(4, "lum"), "by_spectrum": xf.stat_by_spectrum(keys),
           "pegged": pegged("M0", keys)}
    p_l22, p_l10, p_n22, p_n10 = xf.par(1, "lum"), xf.par(4, "lum"), xf.par(1, "nh"), xf.par(4, "nh")
    d2 = v["log_lint_2_10"] - v["log_lint10_2_10"]
    dn = math.log10(v["nh_1e22"] / v["nh10_1e22"])
    cval = c22.values[0]

    def restore():
        p_l10.link = ""
        p_n10.link = ""
        p_l10.values = [v["log_lint10_2_10"], 0.01, 38, 38, 46, 46]
        p_n10.values = [v["nh10_1e22"], 0.1, 0.01, 0.01, 500, 500]
        p_n22.values = [v["nh_1e22"], 0.1, 0.01, 0.01, 500, 500]
        p_l22.values = [v["log_lint_2_10"], 0.01, 38, 38, 46, 46]
        xf.par(1, "gamma").values = [v["gamma"], 0.01, 1.0, 1.0, 3.0, 3.0]
        c22.values = [cval, 0.01, 0.2, 0.2, 3.0, 3.0]
        xf.fit(1)

    res["delta2"] = dict(xf.profile(lambda x: setattr(p_l10, "link", f"p{p_l22.index} - ({x:.6f})"), restore, best, d2,
                                    0.04, -2.0, 2.0, label=f"{mode} d2"), value=d2)
    res["delta_nh"] = dict(xf.profile(lambda x: setattr(p_n10, "link", f"p{p_n22.index} * {10 ** (-x):.8f}"), restore, best,
                                      dn, 0.04, -2.5, 2.5, label=f"{mode} dN"), value=dn)
    return res


def run_fits() -> None:
    xf.setup()
    xf.prepare()
    nh_gal = xf.galactic_nh()
    s7 = json.loads((RES / "s7_fits.json").read_text())
    ckpt_path = FITDIR / "s7b_checkpoint.json"
    ckpt = json.loads(ckpt_path.read_text()) if ckpt_path.exists() else {}

    def step(name, fn):
        if name not in ckpt:
            t0 = time.time()
            try:
                ckpt[name] = fn()
            except Exception as exc:  # noqa: BLE001
                ckpt[name] = {"status": f"failed: {type(exc).__name__}: {str(exc)[:200]}"}
            if isinstance(ckpt[name], dict):
                ckpt[name]["seconds"] = time.time() - t0
                if ckpt[name]["seconds"] > TIME_LIMIT_S:
                    ckpt[name]["over_time_limit"] = True
            ckpt_path.write_text(json.dumps(ckpt) + "\n")
            print(f"done {name} in {time.time() - t0:.0f} s", flush=True)
        return ckpt[name]

    keys22 = ["xrt22", "fpma", "fpmb"]

    # --- A3: .xcm of the final fits and a reload in a new session
    def xcm_all():
        out = {}
        for model in ("M0", "M2"):
            xf.load(keys22)
            xf.build(keys22, model, nh_gal, start_from(s7["fits_2022"][model]["values"]))
            stat = xf.fit(2)
            out[f"fit2022_{model}"] = dict(save_xcm(f"fit2022_{model}", keys22), stage7_statistic=s7["fits_2022"][model]["statistic"],
                                           refit_statistic=stat, values=xf.read_values(keys22, model))
            k1 = ["xrt22", "fpma", "fpmb", "bat"]
            xf.load(k1)
            xf.build(k1, model, nh_gal, start_from(s7["c1"][model]["values"]))
            xf.par(4, "const").values = [s7["c1"][model]["values"]["const_bat"], 0.01, 0.05, 0.05, 20.0, 20.0]
            stat = xf.fit(2)
            out[f"c1_{model}"] = dict(save_xcm(f"c1_{model}", k1), stage7_statistic=s7["c1"][model]["statistic"], refit_statistic=stat)
            k2 = ["xrt22", "fpma", "fpmb", "xrt10"]
            xf.load(k2)
            xf.build(k2, model, nh_gal, start_from(s7["c2"][model]["values"]))
            stat = xf.fit(2)
            out[f"c2_{model}"] = dict(save_xcm(f"c2_{model}", k2), stage7_statistic=s7["c2"][model]["statistic"], refit_statistic=stat)
            print("xcm", model, flush=True)
        out["grouping"] = {"xrt and nustar": "ftgrouppha grouptype=min groupscale=1 (source spectrum); background ungrouped, read through BACKFILE",
                           "bat": "8 channels, not grouped, chi2"}
        return out

    step("xcm", xcm_all)

    def reload_m0():
        r = reload_in_new_session("fit2022_M0")
        r["stage7_statistic"] = s7["fits_2022"]["M0"]["statistic"]
        r["stage7_values"] = {k: s7["fits_2022"]["M0"]["values"][k] for k in ("nh_1e22", "gamma", "log_lint_2_10", "const_fpmb")}
        return r

    step("reload_m0", reload_m0)

    # --- A3: model energy grid 3000 against 1500, M2 on the 2022 data
    def grid_3000():
        xspec.AllModels.setEnergies("0.1 300. 3000 log")
        try:
            xf.load(keys22)
            v0 = s7["fits_2022"]["M0"]["values"]
            xf.build(keys22, "M2", nh_gal, {"nh": v0["nh_1e22"], "gamma": v0["gamma"], "lum": v0["log_lint_2_10"]})
            stat = xf.fit(3)
            vals = xf.read_values(keys22, "M2")
        finally:
            xspec.AllModels.setEnergies(xf.ENERGY_GRID)
        ref = s7["fits_2022"]["M2"]
        return {"grid": "0.1 300. 3000 log", "statistic": stat, "values": vals, "pegged": pegged("M2", keys22),
                "reference_1500": {"statistic": ref["statistic"], "values": ref["values"]},
                "difference": {"statistic": stat - ref["statistic"], "gamma": vals["gamma"] - ref["values"]["gamma"],
                               "nh_fraction": vals["nh_1e22"] / ref["values"]["nh_1e22"] - 1,
                               "refl_norm_fraction": vals["refl_norm"] / ref["values"]["refl_norm"] - 1,
                               "line_norm_fraction": vals["line_norm"] / ref["values"]["line_norm"] - 1,
                               "log_lint": vals["log_lint_2_10"] - ref["values"]["log_lint_2_10"]}}

    step("grid_3000_m2", grid_3000)

    # --- A2: M2 on the 2022 data with a free XRT constant
    def m2_xrt_free():
        xf.load(keys22)
        xf.build(keys22, "M2", nh_gal, start_from(s7["fits_2022"]["M2"]["values"]))
        c = xf.par(1, "const")
        c.values = [0.8, 0.01, 0.2, 0.2, 3.0, 3.0]
        c.frozen = False
        xf.par(2, "const").link = ""
        xf.par(2, "const").values = [1.0, -1]
        stat = xf.fit(3)
        res = {"statistic": stat, "n_free": xf.n_free(), "values": xf.read_values(keys22, "M2"), "const_xrt22": c.values[0],
               "by_spectrum": xf.stat_by_spectrum(keys22), "pegged": pegged("M2", keys22),
               "reference_const_fixed": {"statistic": s7["fits_2022"]["M2"]["statistic"],
                                         "by_spectrum": s7["fits_2022"]["M2"]["by_spectrum"]}}
        res["const_xrt22_90"] = xf.interval(1, "const")
        for nme in ("nh", "gamma", "lum", "refl_norm", "line_norm"):
            res[f"{nme}_90"] = xf.interval(1, nme)
        xf.fit(1)
        res["statistic"] = float(xspec.Fit.statistic)
        res["plot"] = {k: {kk: vv for kk, vv in v.items() if kk in ("energy", "energy_err", "ratio", "ratio_err")}
                       for k, v in xf.plot_arrays(keys22).items()}
        res.update(save_xcm("fit2022_M2_xrt_constant_free", keys22))
        return res

    step("m2_xrt_free", m2_xrt_free)

    # --- A2: Delta2 and DeltaN (M0) under the two readings of the XRT constant
    st0 = start_from(s7["c2"]["M0"]["values"])
    step("c2_m0_const_both_xrt", lambda: c2_with_xrt_constant("both", nh_gal, st0))
    step("c2_m0_const_2022_only", lambda: c2_with_xrt_constant("2022", nh_gal, st0))

    # --- A3: NuSTAR spectra from the default SAA setting, M0
    def saa_default():
        d = WORK / "nustar"
        if not (d / "prod_none_fpmb" / "src_min1.pha").exists():
            return {"status": "spectra of the default SAA setting not extracted (scripts/xray_process.sh products_saa_none)"}
        over = {}     # file names must not start with 'none': XSPEC reads that as 'no response'
        for k, mod in (("fpma", "A"), ("fpmb", "B")):
            src = d / f"prod_none_{k}"
            for name, target in ((f"nu60160365002{mod}01_bk.pha", f"saadef_{k}_bkg.pha"), (f"nu60160365002{mod}01_sr.rmf", f"saadef_{k}.rmf"),
                                 (f"nu60160365002{mod}01_sr.arf", f"saadef_{k}.arf")):
                shutil.copyfile(src / name, FITDIR / target)
            shutil.copyfile(src / "src_min1.pha", FITDIR / f"saadef_{k}_min1.pha")
            with pf.open(FITDIR / f"saadef_{k}_min1.pha", mode="update") as h:
                hd = h["SPECTRUM"].header
                hd["BACKFILE"], hd["RESPFILE"], hd["ANCRFILE"] = f"saadef_{k}_bkg.pha", f"saadef_{k}.rmf", f"saadef_{k}.arf"
            over[k] = {"pha": f"saadef_{k}_min1.pha"}
        r = xf.fit_model(keys22, "M0", nh_gal, start_from(s7["fits_2022"]["M0"]["values"]), errors=["nh", "gamma", "lum"],
                         overrides=over)
        r["fluxes_observed"] = xf.observed_fluxes(keys22)
        ref = s7["fits_2022"]["M0"]
        r["reference_optimized"] = {"values": ref["values"], "fluxes_observed": ref["fluxes_observed"],
                                    "exposure_s": {k: ref["by_spectrum"][k]["exposure_s"] for k in ("fpma", "fpmb")}}
        r["difference"] = {"gamma": r["values"]["gamma"] - ref["values"]["gamma"],
                           "nh_fraction": r["values"]["nh_1e22"] / ref["values"]["nh_1e22"] - 1,
                           "log_lint": r["values"]["log_lint_2_10"] - ref["values"]["log_lint_2_10"],
                           "flux_3_79_fraction": r["fluxes_observed"]["fpma"]["3-79"] / ref["fluxes_observed"]["fpma"]["3-79"] - 1,
                           "flux_10_40_fraction": r["fluxes_observed"]["fpma"]["10-40"] / ref["fluxes_observed"]["fpma"]["10-40"] - 1}
        return r

    step("saa_default_m0", saa_default)
    out = {"note": "stage 7b: computed after the stage 7 results were seen (docs/XRAY_PLAN.md section 9)",
           "time_limit_s": TIME_LIMIT_S, **ckpt}
    (RES / "s7b_fits.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    print("written s7b_fits.json", flush=True)


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else "light"
    os.chdir(config.PROJECT_ROOT)
    if stage == "light":
        run_light()
    elif stage == "fits":
        run_fits()
    else:
        raise SystemExit(__doc__)
