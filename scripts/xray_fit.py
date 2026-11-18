"""Stage 7: spectral fits of CGCG 312-012 with PyXspec (docs/XRAY_PLAN.md section 4).

    python scripts/xray_fit.py fits        -> results/s7_fits.json
    python scripts/xray_fit.py plotdata    -> results/s7_spectra_plot.json
    python scripts/xray_fit.py arf         -> results/s7_arf.json
    python scripts/xray_fit.py lumcheck    -> results/s7_lum_check.json
    python scripts/xray_fit.py posthoc     -> results/s7_posthoc.json (not registered)
    python scripts/xray_fit.py sims [N]    -> results/s7_sim_validation.json
    python scripts/xray_fit.py mcmc        -> results/s7_mcmc.json

Run inside the HEASoft environment (source scripts/heasoft_env.sh). Spectra
come from scripts/xray_process.sh (data/cache/s7/work/).

One model expression is used for M0, M1 and M2; components that a model
does not have are switched off:

    constant * TBabs * ( zTBabs * clumin * zpowerlw      primary (absorbed)
                         + clumin * zpowerlw             scattered (same photon index)
                         + pexrav                        reflection only (rel_refl = -1)
                         + zgauss )                      6.4 keV line

clumin makes the unabsorbed rest-frame 2-10 keV luminosity of a power law a
fit parameter (lg10Lum), so its interval comes from the fit statistic like
that of N_H and the photon index.
"""

from __future__ import annotations

import json
import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from irecho import config, physics, xray  # noqa: E402

import xspec  # noqa: E402

WORK = xray.S7_DIR / "work"
Z = xray.TARGET["z"]
DELTA_90 = 2.706
EXPR = "constant*TBabs*(zTBabs*clumin*zpowerlw + clumin*zpowerlw + pexrav + zgauss)"
# the same expression with the unused components left out (pexrav is slow to evaluate)
EXPRS = {"M0": "constant*TBabs*zTBabs*clumin*zpowerlw",
         "M1": "constant*TBabs*(zTBabs*clumin*zpowerlw + clumin*zpowerlw)",
         "M2": EXPR}
ENERGY_GRID = "0.1 300. 1500 log"
OFF_LUM = 30.0   # lg10Lum of a switched-off power law (erg/s)

# component.parameter names used below
P = {"const": ("constant", "factor"), "nh_gal": ("TBabs", "nH"), "nh": ("zTBabs", "nH"),
     "lum": ("clumin", "lg10Lum"), "gamma": ("zpowerlw", "PhoIndex"),
     "lum_s": ("clumin_6", "lg10Lum"), "gamma_s": ("zpowerlw_7", "PhoIndex"),
     "refl_gamma": ("pexrav", "PhoIndex"), "refl_norm": ("pexrav", "norm"), "line_norm": ("zgauss", "norm")}

FITDIR = WORK / "fit"
SPECTRA = {
    "xrt22": {"pha": "xrt22_min1.pha", "ignore": "**-0.3 10.0-**", "stat": "cstat"},
    "fpma": {"pha": "fpma_min1.pha", "ignore": "**-3.0 79.0-**", "stat": "cstat"},
    "fpmb": {"pha": "fpmb_min1.pha", "ignore": "**-3.0 79.0-**", "stat": "cstat"},
    "xrt10": {"pha": "xrt10_min1.pha", "ignore": "**-0.3 10.0-**", "stat": "cstat"},
    "bat": {"pha": "bat.pha", "ignore": None, "stat": "chi"},   # all 8 channels (14-195 keV)
}
SOURCES = {  # key: (directory, ungrouped, grouped files, background, rmf, arf)
    "xrt22": (WORK / "xrt22" / "prod", "src.pha", {"min1": "src_min1.pha"}, "bkg.pha", "src.rmf", "src.arf"),
    "xrt10": (WORK / "xrt10" / "prod", "src.pha", {"min1": "src_min1.pha"}, "bkg.pha", "src.rmf", "src.arf"),
    "fpma": (WORK / "nustar" / "prod_fpma", "nu60160365002A01_sr.pha",
             {"min1": "src_min1.pha", "min20": "src_min20.pha"}, "nu60160365002A01_bk.pha",
             "nu60160365002A01_sr.rmf", "nu60160365002A01_sr.arf"),
    "fpmb": (WORK / "nustar" / "prod_fpmb", "nu60160365002B01_sr.pha",
             {"min1": "src_min1.pha", "min20": "src_min20.pha"}, "nu60160365002B01_bk.pha",
             "nu60160365002B01_sr.rmf", "nu60160365002B01_sr.arf"),
}


def prepare() -> None:
    """One flat directory with short names; the PHA headers point at background and responses.

    XSPEC then loads everything from the keywords. All later work runs with
    this directory as the working directory.
    """
    import shutil

    from astropy.io import fits as pf

    FITDIR.mkdir(parents=True, exist_ok=True)
    for key, (d, raw, grouped, bkg, rmf, arf) in SOURCES.items():
        for name, src in ((f"{key}_bkg.pha", bkg), (f"{key}.rmf", rmf), (f"{key}.arf", arf)):
            shutil.copyfile(d / src, FITDIR / name)
        for tag, src in dict(grouped, src=raw).items():
            out = FITDIR / f"{key}_{tag}.pha"
            shutil.copyfile(d / src, out)
            with pf.open(out, mode="update") as h:
                hd = h["SPECTRUM"].header
                hd["BACKFILE"], hd["RESPFILE"], hd["ANCRFILE"] = f"{key}_bkg.pha", f"{key}.rmf", f"{key}.arf"
    shutil.copyfile(WORK / "bat" / "bat157.rsp", FITDIR / "bat.rsp")
    shutil.copyfile(WORK / "bat" / "bat157.pha", FITDIR / "bat.pha")
    with pf.open(FITDIR / "bat.pha", mode="update") as h:
        h[1].header["RESPFILE"] = "bat.rsp"
    os.chdir(FITDIR)


def galactic_nh() -> float:
    """Weighted HI4PI column from the log of the HEASoft nh tool, in 1e22 cm^-2."""
    for line in (WORK / "nh.log").read_text().splitlines():
        if "Weighted average nH" in line:
            return float(line.split()[-1]) / 1e22
    raise SystemExit("no nH in nh.log")


def setup() -> None:
    xspec.Xset.chatter = 0
    xspec.Xset.allowPrompting = False
    xspec.Xset.logChatter = 0
    xspec.Xset.abund = "wilm"
    xspec.Xset.xsect = "vern"
    xspec.Xset.cosmo = f"{physics.H0_KM_S_MPC} 0 {1 - physics.OMEGA_M}"
    xspec.Fit.query = "yes"
    xspec.Fit.nIterations = 200
    xspec.Fit.criticalDelta = 1e-4
    xspec.Plot.xAxis = "keV"
    # Models are computed on one logarithmic grid (0.53% per bin) and interpolated onto the
    # response energies. With the native grids (about 10,000 energies) one fit with pexrav
    # took two minutes.
    xspec.AllModels.setEnergies(ENERGY_GRID)
    xspec.Xset.parallel.leven = 3
    xspec.Xset.parallel.error = 3


def load(keys, overrides: dict | None = None) -> None:
    """Load spectra (working directory FITDIR), each in its own data group, and set the statistics."""
    xspec.AllData.clear()
    xspec.AllModels.clear()
    cfgs = [dict(SPECTRA[key], **(overrides or {}).get(key, {})) for key in keys]
    xspec.AllData(" ".join(f"{i}:{i} {cfg['pha']}" for i, cfg in enumerate(cfgs, start=1)))
    for i, cfg in enumerate(cfgs, start=1):
        if cfg["ignore"]:   # XSPEC also drops the channel that contains the limit itself
            xspec.AllData(i).ignore(cfg["ignore"])
    xspec.AllData.ignore("bad")
    xspec.Fit.statMethod = "cstat"
    for i, cfg in enumerate(cfgs, start=1):
        if cfg["stat"] != "cstat":
            xspec.Fit.statMethod = f"{cfg['stat']} {i}"


def par(group: int, name: str):
    comp, pname = P[name]
    return getattr(getattr(xspec.AllModels(group), comp), pname)


def build(keys, model: str, nh_gal: float, start: dict | None = None) -> None:
    """Model for the loaded groups. Group 1 carries the 2022 parameters."""
    xspec.AllModels.clear()
    m = xspec.Model(EXPRS[model])
    start = dict({"nh": 15.0, "gamma": 1.8, "lum": 42.5, "lum_s": 40.5, "refl_norm": 1e-3, "line_norm": 3e-6},
                 **(start or {}))
    m.constant.factor.values = [1.0, -1]
    m.TBabs.nH.values = [nh_gal, -1]
    m.zTBabs.nH.values = [start["nh"], 0.1, 0.01, 0.01, 500, 500]
    m.zTBabs.Redshift.values = [Z, -1]
    m.clumin.Emin.values, m.clumin.Emax.values, m.clumin.Redshift.values = [2.0, -1], [10.0, -1], [Z, -1]
    m.clumin.lg10Lum.values = [start["lum"], 0.01, 38, 38, 46, 46]
    m.zpowerlw.PhoIndex.values = [start["gamma"], 0.01, 1.0, 1.0, 3.0, 3.0]
    m.zpowerlw.Redshift.values, m.zpowerlw.norm.values = [Z, -1], [1.0, -1]
    if model in ("M1", "M2"):
        c = m.clumin_6
        c.Emin.values, c.Emax.values, c.Redshift.values = [2.0, -1], [10.0, -1], [Z, -1]
        c.lg10Lum.values = [start["lum_s"], 0.01, 36, 36, 46, 46]
        m.zpowerlw_7.PhoIndex.link = f"p{m.zpowerlw.PhoIndex.index}"
        m.zpowerlw_7.Redshift.values, m.zpowerlw_7.norm.values = [Z, -1], [1.0, -1]
    if model == "M2":
        m.pexrav.PhoIndex.link = f"p{m.zpowerlw.PhoIndex.index}"
        # no cutoff: the hard upper limit of foldE (1e6 keV); 0 is not accepted by this XSPEC version
        # rel_refl = -1: reflection only (the default hard limits 0-1e6 have to be widened for this)
        m.pexrav.foldE.values, m.pexrav.rel_refl.values = [1e6, -1], [-1.0, -1, -10, -10, 10, 10]
        m.pexrav.Redshift.values, m.pexrav.abund.values, m.pexrav.Fe_abund.values = [Z, -1], [1.0, -1], [1.0, -1]
        m.pexrav.cosIncl.values = [0.45, -1]
        m.zgauss.LineE.values, m.zgauss.Sigma.values, m.zgauss.Redshift.values = [6.4, -1], [0.01, -1], [Z, -1]
        m.pexrav.norm.values = [start["refl_norm"], 1e-5, 0, 0, 1, 1]
        m.zgauss.norm.values = [start["line_norm"], 1e-7, 0, 0, 1e-2, 1e-2]
    for g, key in enumerate(keys, start=1):
        if key == "fpmb":
            par(g, "const").link = ""
            par(g, "const").values = [1.0, 0.01, 0.5, 0.5, 2.0, 2.0]
        elif key == "bat":
            par(g, "const").link = ""
            par(g, "const").values = [1.0, 0.01, 0.05, 0.05, 20.0, 20.0]
        elif key == "xrt10":
            par(g, "nh").link = ""
            par(g, "nh").values = [start.get("nh10", start["nh"]), 0.1, 0.01, 0.01, 500, 500]
            par(g, "lum").link = ""
            par(g, "lum").values = [start.get("lum10", start["lum"]), 0.01, 38, 38, 46, 46]


def fit(n: int = 3) -> float:
    for _ in range(n):
        try:
            xspec.Fit.perform()
        except Exception:  # noqa: BLE001
            pass
    return float(xspec.Fit.statistic)


def n_free() -> int:
    count = 0
    for g in range(1, xspec.AllData.nGroups + 1):
        m = xspec.AllModels(g)
        for i in range(1, m.nParameters + 1):
            p = m(i)
            if not p.frozen and not p.link:
                count += 1
    return count


def global_index(group: int, name: str) -> int:
    """Parameter number as the XSPEC commands count it (p.index is the number inside its data group)."""
    return xspec.AllModels(group).startParIndex + par(group, name).index - 1


def interval(group: int, name: str) -> list:
    """90% interval of a free parameter with Fit.error (refits when a lower statistic is found)."""
    p = par(group, name)
    idx = global_index(group, name)
    for _ in range(3):
        try:
            xspec.Fit.error(f"maximum 50 {DELTA_90} {idx}")
        except Exception:  # noqa: BLE001
            fit(1)
            continue
        lo, hi, code = p.error
        if code[1] == "T":      # new minimum found
            fit(2)
            continue
        return [float(lo), float(hi), code]
    lo, hi, code = p.error
    return [float(lo), float(hi), code]


def interval_of(p) -> list:
    """Fit.error for a parameter of the first data group of an ad-hoc model."""
    for _ in range(3):
        try:
            xspec.Fit.error(f"maximum 50 {DELTA_90} {p.index}")
        except Exception:  # noqa: BLE001
            fit(1)
            continue
        lo, hi, code = p.error
        if code[1] == "T":
            fit(2)
            continue
        return [float(lo), float(hi), code]
    lo, hi, code = p.error
    return [float(lo), float(hi), code]


def profile(setter, release, best_stat: float, centre: float, step: float, lo_lim: float, hi_lim: float,
            max_steps: int = 30, label: str = "") -> dict:
    """Profile of the fit statistic over a quantity that is not a single parameter.

    setter(value) fixes the quantity, the rest is refitted. Walks outward from
    the best value until the statistic has risen by more than DELTA_90 and
    interpolates the crossing. release() restores the free state. A side that
    does not cross within max_steps is reported as None (not bounded).
    """
    out = {"grid": [], "delta_stat": []}
    bounds = {}
    for sign in (-1, +1):
        prev_v, prev_d, crossing = centre, 0.0, None
        for k in range(1, max_steps + 1):
            v = centre + sign * k * step
            if v < lo_lim or v > hi_lim:
                break
            setter(v)
            d = fit(1) - best_stat
            out["grid"].append(v)
            out["delta_stat"].append(d)
            print(f"   profile {label} {v:+.3f} dstat {d:.3f}", flush=True)
            if d > DELTA_90:
                crossing = prev_v + (DELTA_90 - prev_d) / (d - prev_d) * (v - prev_v)
                break
            prev_v, prev_d = v, d
        bounds[sign] = crossing
        release()
    order = np.argsort(out["grid"])
    out["grid"] = [float(out["grid"][i]) for i in order]
    out["delta_stat"] = [float(out["delta_stat"][i]) for i in order]
    out["lo"], out["hi"] = bounds[-1], bounds[+1]
    out["step"] = step
    out["min_delta_stat"] = float(min(out["delta_stat"])) if out["delta_stat"] else 0.0
    return out


def read_values(keys, model: str) -> dict:
    out = {}
    g1 = 1
    out["nh_1e22"] = par(g1, "nh").values[0]
    out["gamma"] = par(g1, "gamma").values[0]
    out["log_lint_2_10"] = par(g1, "lum").values[0]
    if model in ("M1", "M2"):
        out["log_lscat_2_10"] = par(g1, "lum_s").values[0]
        out["scattered_fraction"] = 10 ** (out["log_lscat_2_10"] - out["log_lint_2_10"])
    if model == "M2":
        out["refl_norm"] = par(g1, "refl_norm").values[0]
        out["line_norm"] = par(g1, "line_norm").values[0]
        d_l = physics.luminosity_distance_mpc(Z)
        k_prim = 10 ** out["log_lint_2_10"] / xray.zpowerlw_rest_luminosity(1.0, out["gamma"], Z, d_l)
        out["primary_norm_1kev"] = k_prim
        out["reflection_r"] = out["refl_norm"] / (k_prim * (1 + Z) ** (-out["gamma"]))
    for g, key in enumerate(keys, start=1):
        if key in ("fpmb", "bat"):
            out[f"const_{key}"] = par(g, "const").values[0]
        if key == "xrt10":
            out["nh10_1e22"] = par(g, "nh").values[0]
            out["log_lint10_2_10"] = par(g, "lum").values[0]
    return out


def stat_by_spectrum(keys) -> dict:
    """Statistic of each spectrum, read from the text of the XSPEC 'show fit' command."""
    import re

    log = FITDIR / "show_fit.log"
    if log.exists():
        log.unlink()
    xspec.Xset.openLog(str(log))
    xspec.Xset.logChatter = 10
    xspec.Fit.show()
    xspec.Xset.logChatter = 0
    xspec.Xset.closeLog()
    parts = {int(m.group(3)): float(m.group(1)) for m in re.finditer(
        r"([-+0-9.eE]+)\s+using (\d+) bins, spectrum (\d+)", log.read_text())}
    out = {}
    for i, key in enumerate(keys, start=1):
        s = xspec.AllData(i)
        out[key] = {"statistic": parts.get(i), "n_bins": len(s.noticed), "stat_method": SPECTRA[key]["stat"],
                    "exposure_s": float(s.exposure), "net_rate": float(s.rate[0]), "net_rate_err": float(s.rate[1])}
    return out


def fit_model(keys, model: str, nh_gal: float, start: dict | None = None, errors=True,
              overrides: dict | None = None) -> dict:
    """Fit one model to the listed spectra. errors: True (all free parameters), False, or a list of names."""
    load(keys, overrides)
    build(keys, model, nh_gal, start)
    stat = fit(3)
    res = {"model": model, "spectra": list(keys), "statistic": stat, "n_free": n_free(),
           "n_bins": int(sum(len(xspec.AllData(i).noticed) for i in range(1, len(keys) + 1)))}
    res["aic"] = stat + 2 * res["n_free"]
    if errors:
        names = ["nh", "gamma", "lum"] + (["lum_s"] if model != "M0" else []) + \
                (["refl_norm", "line_norm"] if model == "M2" else [])
        extra = ["const_fpmb", "const_bat", "nh10", "lum10"]
        if errors is not True:
            names = [n for n in names if n in errors]
            extra = [n for n in extra if n in errors]
        res["intervals_90"] = {}
        for nme in names:
            res["intervals_90"][nme] = interval(1, nme)
            print(f"   interval {model} {nme} {res['intervals_90'][nme]}", flush=True)
        for g, key in enumerate(keys, start=1):
            if key in ("fpmb", "bat") and f"const_{key}" in extra:
                res["intervals_90"][f"const_{key}"] = interval(g, "const")
            if key == "xrt10":
                for nme, short in (("nh10", "nh"), ("lum10", "lum")):
                    if nme in extra:
                        res["intervals_90"][nme] = interval(g, short)
                        print(f"   interval {model} {nme} {res['intervals_90'][nme]}", flush=True)
        stat = fit(1)
        res["statistic"] = stat
        res["aic"] = stat + 2 * res["n_free"]
    res["values"] = read_values(keys, model)
    res["by_spectrum"] = stat_by_spectrum(keys)
    return res


def observed_fluxes(keys) -> dict:
    """Model fluxes (absorbed, observed frame) per spectrum in the bands of the plan."""
    bands = {"xrt22": [(0.3, 2.0), (2.0, 10.0)], "xrt10": [(0.3, 2.0), (2.0, 10.0)],
             "fpma": [(3.0, 10.0), (10.0, 40.0), (3.0, 79.0)], "fpmb": [(3.0, 10.0), (10.0, 40.0), (3.0, 79.0)],
             "bat": [(14.0, 50.0), (14.0, 195.0)]}
    out = {}
    for i, key in enumerate(keys, start=1):
        out[key] = {}
        for lo, hi in bands[key]:
            xspec.AllModels.calcFlux(f"{lo} {hi}")
            out[key][f"{lo:g}-{hi:g}"] = float(xspec.AllData(i).flux[0])
    return out


def plot_arrays(keys, rebin=(3.0, 30)) -> dict:
    """Data, model and ratio per spectrum for the figure (rebinned for display only)."""
    xspec.Plot.device = "/null"
    xspec.Plot.setRebin(*rebin)
    xspec.Plot.background = False
    out = {}
    xspec.Plot("ldata")
    for i, key in enumerate(keys, start=1):
        out[key] = {"energy": list(map(float, xspec.Plot.x(i))), "energy_err": list(map(float, xspec.Plot.xErr(i))),
                    "rate": list(map(float, xspec.Plot.y(i))), "rate_err": list(map(float, xspec.Plot.yErr(i))),
                    "model": list(map(float, xspec.Plot.model(i)))}
    xspec.Plot("ratio")
    for i, key in enumerate(keys, start=1):
        out[key]["ratio"] = list(map(float, xspec.Plot.y(i)))
        out[key]["ratio_err"] = list(map(float, xspec.Plot.yErr(i)))
    xspec.Plot("eeufspec")
    for i, key in enumerate(keys, start=1):
        out[key]["eeuf"] = list(map(float, xspec.Plot.y(i)))
        out[key]["eeuf_err"] = list(map(float, xspec.Plot.yErr(i)))
        out[key]["eeuf_model"] = list(map(float, xspec.Plot.model(i)))
    xspec.Plot.setRebin(0, 1)
    return out


def comparison_c2(model: str, nh_gal: float, start: dict) -> dict:
    """2022 against 2010: joint fit and profiles of Delta2 and DeltaN (section 4.3)."""
    keys = ["xrt22", "fpma", "fpmb", "xrt10"]
    res = fit_model(keys, model, nh_gal, start, errors=["nh", "gamma", "lum", "nh10", "lum10"])
    best = res["statistic"]
    v = res["values"]
    g10 = keys.index("xrt10") + 1
    p_l22, p_l10 = par(1, "lum"), par(g10, "lum")
    p_n22, p_n10 = par(1, "nh"), par(g10, "nh")
    d2 = v["log_lint_2_10"] - v["log_lint10_2_10"]
    dn = math.log10(v["nh_1e22"] / v["nh10_1e22"])

    def restore():
        p_l10.link = ""
        p_n10.link = ""
        p_l10.values = [v["log_lint10_2_10"], 0.01, 38, 38, 46, 46]
        p_n10.values = [v["nh10_1e22"], 0.1, 0.01, 0.01, 500, 500]
        p_n22.values = [v["nh_1e22"], 0.1, 0.01, 0.01, 500, 500]
        p_l22.values = [v["log_lint_2_10"], 0.01, 38, 38, 46, 46]
        par(1, "gamma").values = [v["gamma"], 0.01, 1.0, 1.0, 3.0, 3.0]
        fit(1)

    def set_d2(x):
        p_l10.link = f"p{p_l22.index} - ({x:.6f})"

    def set_dn(x):
        p_n10.link = f"p{p_n22.index} * {10 ** (-x):.8f}"

    res["delta2"] = dict(profile(set_d2, restore, best, d2, 0.04, -2.0, 2.0, label=f"{model} d2"), value=d2)
    res["delta_nh"] = dict(profile(set_dn, restore, best, dn, 0.04, -2.5, 2.5, label=f"{model} dN"), value=dn)
    # likelihood-ratio checks: N_H, then luminosity, tied between the epochs
    p_n10.link = f"p{p_n22.index}"
    res["stat_nh_tied"] = fit(2)
    res["delta_stat_nh_tied"] = res["stat_nh_tied"] - best
    restore()
    p_l10.link = f"p{p_l22.index}"
    res["stat_lum_tied"] = fit(2)
    res["delta_stat_lum_tied"] = res["stat_lum_tied"] - best
    restore()
    res["fluxes_observed"] = observed_fluxes(keys)
    return res


def comparison_c1(model: str, nh_gal: float, start: dict) -> dict:
    """2022 against the BAT 157-month average: one free constant on BAT (section 4.3)."""
    keys = ["xrt22", "fpma", "fpmb", "bat"]
    res = fit_model(keys, model, nh_gal, start, errors=["nh", "gamma", "lum", "const_bat"])
    k = res["values"]["const_bat"]
    lo, hi, code = res["intervals_90"]["const_bat"]
    res["delta1"] = {"value": -math.log10(k), "lo": -math.log10(hi) if hi > 0 else None,
                     "hi": -math.log10(lo) if lo > 0 else None, "k_bat": k, "k_bat_interval": [lo, hi], "code": code}
    res["fluxes_observed"] = observed_fluxes(keys)
    return res


def simple_checks(nh_gal: float) -> dict:
    """Checks that lean less on the model (section 4.4, item 3a) and BAT alone.

    Note (stage 7b): the BAT power law here is fitted over all 8 channels
    (14-195 keV) and only the flux is integrated over 14-50 keV. The fit of
    the 14-50 keV channels alone is in scripts/xray_7b.py (bat_band).
    """
    out = {}
    # BAT alone, power law
    load(["bat"])
    xspec.AllModels.clear()
    m = xspec.Model("cflux*powerlaw")
    m.cflux.Emin.values, m.cflux.Emax.values = [14.0, -1], [50.0, -1]
    m.cflux.lg10Flux.values = [-11.5, 0.01, -14, -14, -9, -9]
    m.powerlaw.PhoIndex.values = [1.9, 0.01, 0.5, 0.5, 3.5, 3.5]
    m.powerlaw.norm.values = [1.0, -1]
    stat = fit(3)
    out["bat_powerlaw"] = {"chi2": stat, "dof": int(xspec.Fit.dof), "gamma": m.powerlaw.PhoIndex.values[0],
                           "gamma_90": interval_of(m.powerlaw.PhoIndex), "log_flux_14_50": m.cflux.lg10Flux.values[0],
                           "log_flux_14_50_90": interval_of(m.cflux.lg10Flux)}
    m.cflux.Emin.values, m.cflux.Emax.values = [14.0, -1], [195.0, -1]
    fit(3)
    out["bat_powerlaw"]["log_flux_14_195"] = m.cflux.lg10Flux.values[0]
    out["bat_powerlaw"]["log_flux_14_195_90"] = interval_of(m.cflux.lg10Flux)
    # NuSTAR 14-50 keV, power law
    over = {k: {"ignore": "**-14.0 50.0-**"} for k in ("fpma", "fpmb")}
    load(["fpma", "fpmb"], over)
    xspec.AllModels.clear()
    m = xspec.Model("constant*cflux*powerlaw")
    m.constant.factor.values = [1.0, -1]
    m.cflux.Emin.values, m.cflux.Emax.values = [14.0, -1], [50.0, -1]
    m.cflux.lg10Flux.values = [-11.3, 0.01, -14, -14, -9, -9]
    m.powerlaw.PhoIndex.values = [1.8, 0.01, 0.5, 0.5, 3.5, 3.5]
    m.powerlaw.norm.values = [1.0, -1]
    m2 = xspec.AllModels(2)
    m2.constant.factor.link = ""
    m2.constant.factor.values = [1.0, 0.01, 0.5, 0.5, 2.0, 2.0]
    stat = fit(3)
    out["nustar_14_50_powerlaw"] = {"cstat": stat, "n_bins": int(sum(len(xspec.AllData(i).noticed) for i in (1, 2))),
                                    "gamma": m.powerlaw.PhoIndex.values[0], "gamma_90": interval_of(m.powerlaw.PhoIndex),
                                    "log_flux_14_50": m.cflux.lg10Flux.values[0],
                                    "log_flux_14_50_90": interval_of(m.cflux.lg10Flux),
                                    "const_fpmb": m2.constant.factor.values[0]}
    a, b = out["nustar_14_50_powerlaw"], out["bat_powerlaw"]
    d = a["log_flux_14_50"] - b["log_flux_14_50"]
    err_lo = math.hypot(a["log_flux_14_50"] - a["log_flux_14_50_90"][0], b["log_flux_14_50_90"][1] - b["log_flux_14_50"])
    err_hi = math.hypot(a["log_flux_14_50_90"][1] - a["log_flux_14_50"], b["log_flux_14_50"] - b["log_flux_14_50_90"][0])
    out["delta_14_50_powerlaw"] = {"value": d, "lo": d - err_lo, "hi": d + err_hi,
                                   "note": "log10 flux ratio NuSTAR 2022 / BAT average, 14-50 keV, separate power laws "
                                           "(BAT fitted over 14-195 keV, flux integrated over 14-50 keV); "
                                           "errors combined in quadrature, statistical only"}
    return out


def run_fits() -> None:
    """All registered fits. Every step is saved to a checkpoint and skipped on a restart."""
    setup()
    prepare()
    nh_gal = galactic_nh()
    ckpt_path = FITDIR / "fits_checkpoint.json"
    ckpt = json.loads(ckpt_path.read_text()) if ckpt_path.exists() else {}

    def step(name, fn):
        if name not in ckpt:
            import time
            t0 = time.time()
            ckpt[name] = fn()
            ckpt_path.write_text(json.dumps(ckpt) + "\n")
            print(f"done {name} in {time.time() - t0:.0f} s", flush=True)
        return ckpt[name]

    def start_of(model):
        vv = ckpt[f"fits_2022/{model}"]["values"]
        return {"nh": vv["nh_1e22"], "gamma": vv["gamma"], "lum": vv["log_lint_2_10"],
                "lum_s": vv.get("log_lscat_2_10", 40.5), "refl_norm": vv.get("refl_norm", 1e-3),
                "line_norm": vv.get("line_norm", 3e-6)}

    out = {"xspec_version": list(xspec.Xset.version), "galactic_nh_1e22": nh_gal, "expression": EXPR,
           "expressions": EXPRS, "energy_grid": ENERGY_GRID, "abund": "wilm", "xsect": "vern",
           "cosmology": {"H0": physics.H0_KM_S_MPC, "omega_m": physics.OMEGA_M,
                         "d_l_mpc": physics.luminosity_distance_mpc(Z)}, "z": Z}
    # 1. model choice on the 2022 data
    keys22 = ["xrt22", "fpma", "fpmb"]

    def fit22(model):
        st = start_of("M0") if model != "M0" else {}
        r = fit_model(keys22, model, nh_gal, {k: st[k] for k in ("nh", "gamma", "lum")} if st else {})
        r["fluxes_observed"] = observed_fluxes(keys22)
        return r

    out["fits_2022"] = {}
    for model in ("M0", "M1", "M2"):
        r = step(f"fits_2022/{model}", lambda m=model: fit22(m))
        out["fits_2022"][model] = r
        print(model, round(r["statistic"], 2), r["n_free"], {k: round(v, 4) for k, v in r["values"].items()}, flush=True)
    aic = {m: out["fits_2022"][m]["aic"] for m in out["fits_2022"]}
    chosen = xray.choose_model(aic)
    out["aic_2022"] = aic
    out["chosen_model"] = chosen
    print("chosen", chosen, aic, flush=True)
    # 2. comparisons with every model (the registered verdict uses the chosen one)
    out["c1"], out["c2"], out["c3"] = {}, {}, {}
    for model in ("M0", "M1", "M2"):
        c1 = step(f"c1/{model}", lambda m=model: comparison_c1(m, nh_gal, start_of(m)))
        out["c1"][model] = c1
        print("C1", model, c1["delta1"], flush=True)
    for model in ("M0", "M1", "M2"):
        c2 = step(f"c2/{model}", lambda m=model: comparison_c2(m, nh_gal, start_of(m)))
        out["c2"][model] = c2
        print("C2", model, {k: c2[k]["value"] for k in ("delta2", "delta_nh")},
              [c2["delta2"]["lo"], c2["delta2"]["hi"]], [c2["delta_nh"]["lo"], c2["delta_nh"]["hi"]],
              "tied NH dstat", round(c2["delta_stat_nh_tied"], 2), flush=True)

    # 3. reference only: 2010 XRT with the BAT average, constant fixed at 1 (C3)
    def fit_c3(model):
        keys = ["xrt10", "bat"]
        load(keys)
        build(["xrt10_plain", "bat_plain"], model, nh_gal, start_of(model))   # no special groups: all tied to group 1
        if model == "M2":   # 59 XRT counts and 8 BAT channels cannot fix a reflection component and a line
            par(1, "line_norm").values = [0.0, -1]
        stat = fit(3)
        names = ["nh", "gamma", "lum"]
        return {"statistic": stat, "n_free": n_free(), "values": read_values([], model),
                "intervals_90": {n_: interval(1, n_) for n_ in names},
                "by_spectrum": stat_by_spectrum(keys), "note": "BAT constant fixed at 1; reference only"}

    for model in ("M0", "M1", "M2"):
        r = step(f"c3/{model}", lambda m=model: fit_c3(m))
        out["c3"][model] = r
        print("C3", model, round(r["statistic"], 2), {k: round(x, 4) for k, x in r["values"].items()}, flush=True)

    # 4. sensitivity: XRT constant free in the 2022 fit; NuSTAR grouped to 20 counts with chi2
    def xrt_free():
        load(keys22)
        build(keys22, chosen, nh_gal, start_of(chosen))
        c = par(1, "const")
        c.values = [1.0, 0.01, 0.3, 0.3, 3.0, 3.0]
        c.frozen = False
        par(2, "const").link = ""          # FPMA stays the reference
        par(2, "const").values = [1.0, -1]
        stat = fit(3)
        return {"statistic": stat, "values": read_values(keys22, chosen), "const_xrt22": c.values[0],
                "const_xrt22_90": interval(1, "const"), "nh_90": interval(1, "nh"), "lum_90": interval(1, "lum")}

    out["sensitivity_xrt_constant_free"] = step("sens/xrt_free", xrt_free)

    def nustar_chi2():
        over = {k: {"pha": f"{k}_min20.pha", "stat": "chi"} for k in ("fpma", "fpmb")}
        res = {}
        for model in ("M0", "M2"):
            r = fit_model(keys22, model, nh_gal, start_of(model), errors=["nh", "gamma", "lum"], overrides=over)
            nu = sum(r["by_spectrum"][k]["statistic"] for k in ("fpma", "fpmb"))
            nb = sum(r["by_spectrum"][k]["n_bins"] for k in ("fpma", "fpmb"))
            r["nustar_chi2"], r["nustar_bins"] = nu, nb
            res[model] = r
        res["note"] = "FPMA and FPMB grouped to >= 20 counts and fitted with chi2; XRT with cstat"
        return res

    out["sensitivity_nustar_chi2"] = step("sens/nustar_chi2", nustar_chi2)
    out["simple_checks"] = step("simple_checks", lambda: simple_checks(nh_gal))
    (config.RESULTS_DIR / "s7_fits.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    print("written", flush=True)


# ----------------------------------------------------------------------------------------------
def run_sims(n_sims: int) -> None:
    """Section 4.4 item 4: fake spectra with the real responses, exposures and backgrounds."""
    setup()
    fits = json.loads((config.RESULTS_DIR / "s7_fits.json").read_text())
    nh_gal, model = fits["galactic_nh_1e22"], fits["chosen_model"]
    truth_fit = fits["c2"][model]["values"]
    keys = ["xrt22", "fpma", "fpmb", "xrt10"]
    start = {"nh": truth_fit["nh_1e22"], "gamma": truth_fit["gamma"], "lum": truth_fit["log_lint_2_10"],
             "lum_s": truth_fit.get("log_lscat_2_10", 40.5), "refl_norm": truth_fit.get("refl_norm", 1e-3),
             "line_norm": truth_fit.get("line_norm", 3e-6), "nh10": truth_fit["nh10_1e22"],
             "lum10": truth_fit["log_lint10_2_10"]}
    truth = {"nh": start["nh"], "gamma": start["gamma"], "lum": start["lum"], "nh10": start["nh10"],
             "lum10": start["lum10"], "delta2": start["lum"] - start["lum10"],
             "delta_nh": math.log10(start["nh"] / start["nh10"])}
    prepare()
    ckpt = FITDIR / "sim_checkpoint.json"
    rows = json.loads(ckpt.read_text()) if ckpt.exists() else []
    over_un = {k: {"pha": f"{k}_src.pha"} for k in keys}
    over = {k: {"pha": f"sim_{k}_g.pha"} for k in keys}
    for i in range(len(rows), n_sims):
        xspec.Xset.seed = 7000 + i
        load(keys, over_un)
        build(keys, model, nh_gal, start)
        for old_file in FITDIR.glob("sim_*.pha"):     # fakeit does not overwrite existing files
            old_file.unlink()
        settings = [xspec.FakeitSettings(fileName=f"sim_{k}.pha") for k in keys]
        xspec.AllData.fakeit(len(keys), settings, applyStats=True, filePrefix="", noWrite=False)
        xspec.AllData.clear()
        for k in keys:
            subprocess.run(["ftgrouppha", f"infile=sim_{k}.pha", f"outfile=sim_{k}_g.pha", "grouptype=min",
                            "groupscale=1", "clobber=yes"], check=True, capture_output=True)
        load(keys, over)
        build(keys, model, nh_gal, start)
        counts = [float(xspec.AllData(j).rate[0] * xspec.AllData(j).exposure) for j in range(1, len(keys) + 1)]
        best = fit(2)
        v = read_values(keys, model)
        row = {"i": i, "stat": best, "net_counts": counts, "nh": v["nh_1e22"], "gamma": v["gamma"], "lum": v["log_lint_2_10"],
               "nh10": v["nh10_1e22"], "lum10": v["log_lint10_2_10"]}
        row["delta2"] = row["lum"] - row["lum10"]
        row["delta_nh"] = math.log10(row["nh"] / row["nh10"])
        g10 = keys.index("xrt10") + 1
        # coverage: the 90% interval contains the truth iff fixing the truth raises the statistic by < 2.706
        tests = {"nh": lambda: setattr(par(1, "nh"), "frozen", True) or par(1, "nh").__setattr__("values", [truth["nh"], -1]),
                 "gamma": lambda: par(1, "gamma").__setattr__("values", [truth["gamma"], -1]),
                 "delta2": lambda: par(g10, "lum").__setattr__("link", f"p{par(1, 'lum').index} - ({truth['delta2']:.6f})"),
                 "delta_nh": lambda: par(g10, "nh").__setattr__(
                     "link", f"p{par(1, 'nh').index} * {10 ** (-truth['delta_nh']):.8f}")}
        for name, constrain in tests.items():
            load(keys, over)
            build(keys, model, nh_gal, dict(start, nh=v["nh_1e22"], gamma=v["gamma"], lum=v["log_lint_2_10"],
                                            nh10=v["nh10_1e22"], lum10=v["log_lint10_2_10"]))
            constrain()
            row[f"dstat_{name}"] = fit(1) - best
        rows.append(row)
        ckpt.write_text(json.dumps(rows) + "\n")
        if i % 10 == 0:
            print(i, {k: round(x, 3) for k, x in row.items() if k not in ("i", "net_counts")},
                  [round(c) for c in counts], flush=True)
    arr = {k: np.array([r[k] for r in rows]) for k in rows[0] if k not in ("i", "net_counts")}
    assert len({round(r["stat"], 3) for r in rows}) > 0.9 * len(rows), "synthetic spectra are not independent"
    summ = {"n_sims": len(rows), "model": model, "truth": truth, "recovered_median": {}, "coverage_90": {}, "p05_p95": {}}
    for k in ("nh", "gamma", "lum", "nh10", "lum10", "delta2", "delta_nh"):
        summ["recovered_median"][k] = float(np.median(arr[k]))
        summ["p05_p95"][k] = [float(x) for x in np.percentile(arr[k], [5, 95])]
    for k in ("nh", "gamma", "delta2", "delta_nh"):
        # a negative difference means the 'best' fit was not the minimum; counted as covered
        summ["coverage_90"][k] = float(np.mean(arr[f"dstat_{k}"] <= DELTA_90))
        summ.setdefault("n_negative_dstat", {})[k] = int(np.sum(arr[f"dstat_{k}"] < -0.05))
    bias = {"nh_fraction": summ["recovered_median"]["nh"] / truth["nh"] - 1,
            "gamma": summ["recovered_median"]["gamma"] - truth["gamma"],
            "delta2_dex": summ["recovered_median"]["delta2"] - truth["delta2"],
            "delta_nh_dex": summ["recovered_median"]["delta_nh"] - truth["delta_nh"]}
    summ["bias"] = bias
    cov = summ["coverage_90"]
    summ["pass"] = {"nh": bool(abs(bias["nh_fraction"]) <= 0.10 and 0.80 <= cov["nh"] <= 0.97),
                    "gamma": bool(abs(bias["gamma"]) <= 0.1 and 0.80 <= cov["gamma"] <= 0.97),
                    "delta2": bool(abs(bias["delta2_dex"]) <= 0.05 and 0.80 <= cov["delta2"] <= 0.97),
                    "delta_nh_coverage_only": bool(0.80 <= cov["delta_nh"] <= 0.97)}
    summ["rows"] = rows
    (config.RESULTS_DIR / "s7_sim_validation.json").write_text(json.dumps(summ, indent=1) + "\n")
    print(json.dumps({k: summ[k] for k in ("n_sims", "truth", "recovered_median", "coverage_90", "bias", "pass")}, indent=1))


def run_mcmc() -> None:
    """Section 4.4 item 2: intervals of Delta2 and DeltaN from an XSPEC chain."""
    from astropy.io import fits as pf

    setup()
    fits = json.loads((config.RESULTS_DIR / "s7_fits.json").read_text())
    nh_gal, model = fits["galactic_nh_1e22"], fits["chosen_model"]
    v = fits["c2"][model]["values"]
    keys = ["xrt22", "fpma", "fpmb", "xrt10"]
    start = {"nh": v["nh_1e22"], "gamma": v["gamma"], "lum": v["log_lint_2_10"], "lum_s": v.get("log_lscat_2_10", 40.5),
             "refl_norm": v.get("refl_norm", 1e-3), "line_norm": v.get("line_norm", 3e-6), "nh10": v["nh10_1e22"],
             "lum10": v["log_lint10_2_10"]}
    prepare()
    load(keys)
    build(keys, model, nh_gal, start)
    fit(4)
    chain_file = FITDIR / "chain_c2.fits"
    if chain_file.exists():
        chain_file.unlink()
    xspec.Xset.seed = 77
    xspec.AllChains.defBurn, xspec.AllChains.defLength, xspec.AllChains.defWalkers = 20000, 100000, 20
    xspec.AllChains.defAlgorithm = "gw"
    c = xspec.Chain(str(chain_file))
    with pf.open(chain_file) as h:
        d = h[1].data
        cols = h[1].columns.names
        g10 = keys.index("xrt10") + 1

        def col(group, name):
            want = f"{par(group, name).name}__{global_index(group, name)}"
            return np.array(d[want], dtype=float)

        l22, l10, n22, n10 = col(1, "lum"), col(g10, "lum"), col(1, "nh"), col(g10, "nh")
    d2, dn = l22 - l10, np.log10(n22 / n10)
    out = {"model": model, "algorithm": "gw", "walkers": 20, "burn": 20000, "length": int(l22.size), "columns": cols,
           "delta2": {"median": float(np.median(d2)), "p05": float(np.percentile(d2, 5)), "p95": float(np.percentile(d2, 95))},
           "delta_nh": {"median": float(np.median(dn)), "p05": float(np.percentile(dn, 5)), "p95": float(np.percentile(dn, 95))},
           "nh22": [float(x) for x in np.percentile(n22, [5, 50, 95])],
           "nh10": [float(x) for x in np.percentile(n10, [5, 50, 95])],
           "lum22": [float(x) for x in np.percentile(l22, [5, 50, 95])],
           "lum10": [float(x) for x in np.percentile(l10, [5, 50, 95])],
           "note": "flat priors in the fit parameters (N_H linear, log L); 5-95% of the chain"}
    (config.RESULTS_DIR / "s7_mcmc.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: out[k] for k in ("delta2", "delta_nh", "nh22", "nh10")}, indent=1))


def run_arf() -> None:
    """Mean effective area of the two XRT extractions (PSF, vignetting and exposure-map losses included)."""
    from astropy.io import fits as pf

    out = {}
    for key in ("xrt10", "xrt22"):
        d = SOURCES[key][0]
        with pf.open(d / SOURCES[key][5]) as h:
            t = h["SPECRESP"].data
            e = 0.5 * (t["ENERG_LO"] + t["ENERG_HI"])
            a = np.array(t["SPECRESP"], dtype=float)
        with pf.open(d / SOURCES[key][1]) as h:
            hd = h["SPECTRUM"].header
            exposure, backscal = float(hd["EXPOSURE"]), float(hd["BACKSCAL"])
        with pf.open(d / SOURCES[key][3]) as h:
            bk = h["SPECTRUM"].header
            bkg_backscal, bkg_exposure = float(bk["BACKSCAL"]), float(bk["EXPOSURE"])
        out[key] = {"exposure_s": exposure, "backscal_source": backscal, "backscal_background": bkg_backscal,
                    "background_exposure_s": bkg_exposure,
                    "mean_area_cm2": {"0.3-2": float(a[(e >= 0.3) & (e < 2)].mean()),
                                      "2-10": float(a[(e >= 2) & (e <= 10)].mean())},
                    "area_cm2_at": {f"{x:g}": float(np.interp(x, e, a)) for x in (1.0, 1.5, 3.0, 5.0, 7.0)},
                    "rmf": (d / "rmf_name.txt").read_text().strip()}
    (config.RESULTS_DIR / "s7_arf.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


def run_lumcheck() -> None:
    """Independent check of the luminosity definition: M0 with a plain zpowerlw normalisation.

    The fit with clumin gives lg10Lum directly (XSPEC cosmology). Here the same
    data are fitted with the norm of zpowerlw as the parameter; the test suite
    converts norm and photon index to a luminosity with the project's own
    distance and compares.
    """
    setup()
    prepare()
    nh_gal = galactic_nh()
    keys = ["xrt22", "fpma", "fpmb"]
    load(keys)
    xspec.AllModels.clear()
    m = xspec.Model("constant*TBabs*zTBabs*zpowerlw")
    m.constant.factor.values = [1.0, -1]
    m.TBabs.nH.values = [nh_gal, -1]
    m.zTBabs.nH.values = [20.0, 0.1, 0.01, 0.01, 500, 500]
    m.zTBabs.Redshift.values = [Z, -1]
    m.zpowerlw.PhoIndex.values = [1.6, 0.01, 1.0, 1.0, 3.0, 3.0]
    m.zpowerlw.Redshift.values = [Z, -1]
    m.zpowerlw.norm.values = [1e-3, 1e-5, 1e-7, 1e-7, 1.0, 1.0]
    c = xspec.AllModels(3).constant.factor
    c.link = ""
    c.values = [1.0, 0.01, 0.5, 0.5, 2.0, 2.0]
    stat = fit(4)
    out = {"expression": "constant*TBabs*zTBabs*zpowerlw", "statistic": stat, "nh_1e22": m.zTBabs.nH.values[0],
           "gamma": m.zpowerlw.PhoIndex.values[0], "norm": m.zpowerlw.norm.values[0], "const_fpmb": c.values[0],
           "z": Z, "d_l_mpc_project": physics.luminosity_distance_mpc(Z)}
    out["log_lint_2_10_from_norm"] = math.log10(xray.zpowerlw_rest_luminosity(
        out["norm"], out["gamma"], Z, out["d_l_mpc_project"]))
    (config.RESULTS_DIR / "s7_lum_check.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


def run_posthoc() -> None:
    """Checks added AFTER the registered results were seen (not in docs/XRAY_PLAN.md).

    1. Is the curvature that M2 fits with reflection also fitted by a high-energy
       cutoff? M0 with zcutoffpl instead of zpowerlw, 2022 data.
    2. Model choice without the sequential step: every model against M0.
    3. Reflection only (M0 + pexrav, no scattered component, no line), 2022 data.
    """
    setup()
    prepare()
    nh_gal = galactic_nh()
    fits = json.loads((config.RESULTS_DIR / "s7_fits.json").read_text())
    keys = ["xrt22", "fpma", "fpmb"]
    out = {"note": "computed after the registered results were seen; not part of docs/XRAY_PLAN.md"}
    v0 = fits["fits_2022"]["M0"]["values"]
    load(keys)
    xspec.AllModels.clear()
    m = xspec.Model("constant*TBabs*zTBabs*clumin*zcutoffpl")
    m.constant.factor.values = [1.0, -1]
    m.TBabs.nH.values = [nh_gal, -1]
    m.zTBabs.nH.values = [v0["nh_1e22"], 0.1, 0.01, 0.01, 500, 500]
    m.zTBabs.Redshift.values = [Z, -1]
    m.clumin.Emin.values, m.clumin.Emax.values, m.clumin.Redshift.values = [2.0, -1], [10.0, -1], [Z, -1]
    m.clumin.lg10Lum.values = [v0["log_lint_2_10"], 0.01, 38, 38, 46, 46]
    m.zcutoffpl.PhoIndex.values = [v0["gamma"], 0.01, 0.5, 0.5, 3.0, 3.0]
    m.zcutoffpl.HighECut.values = [100.0, 1.0, 5.0, 5.0, 1000.0, 1000.0]
    m.zcutoffpl.Redshift.values, m.zcutoffpl.norm.values = [Z, -1], [1.0, -1]
    c = xspec.AllModels(3).constant.factor
    c.link = ""
    c.values = [1.0, 0.01, 0.5, 0.5, 2.0, 2.0]
    stat = fit(4)
    rec = {"expression": "constant*TBabs*zTBabs*clumin*zcutoffpl", "statistic": stat, "n_free": 5, "aic": stat + 10,
           "nh_1e22": m.zTBabs.nH.values[0], "gamma": m.zcutoffpl.PhoIndex.values[0],
           "e_cut_kev": m.zcutoffpl.HighECut.values[0], "log_lint_2_10": m.clumin.lg10Lum.values[0],
           "nh_90": interval_of(m.zTBabs.nH), "gamma_90": interval_of(m.zcutoffpl.PhoIndex),
           "e_cut_90": interval_of(m.zcutoffpl.HighECut), "log_lint_90": interval_of(m.clumin.lg10Lum)}
    out["cutoff_power_law_2022"] = rec
    print("cutoff", {k: rec[k] for k in ("statistic", "nh_1e22", "gamma", "e_cut_kev", "log_lint_2_10", "e_cut_90")}, flush=True)
    # reflection only
    load(keys)
    build(keys, "M2", nh_gal, {"nh": v0["nh_1e22"], "gamma": 1.8, "lum": v0["log_lint_2_10"]})
    par(1, "lum_s").values = [OFF_LUM, -1, 29, 29, 46, 46]
    par(1, "line_norm").values = [0.0, -1]
    stat = fit(4)
    vals = read_values(keys, "M2")
    out["reflection_only_2022"] = {"statistic": stat, "n_free": n_free(), "aic": stat + 2 * n_free(), "values": vals,
                                   "refl_norm_90": interval(1, "refl_norm"), "gamma_90": interval(1, "gamma"),
                                   "nh_90": interval(1, "nh"), "lum_90": interval(1, "lum")}
    print("reflection only", stat, vals, flush=True)
    # 4. M2, 2022 against 2010: the registered profile started outside its range because the best
    #    fit puts the 2010 primary at its lower limit. Fixed grid of Delta2 instead.
    keys4 = ["xrt22", "fpma", "fpmb", "xrt10"]
    v2 = fits["c2"]["M2"]["values"]
    st = {"nh": v2["nh_1e22"], "gamma": v2["gamma"], "lum": v2["log_lint_2_10"], "lum_s": v2["log_lscat_2_10"],
          "refl_norm": v2["refl_norm"], "line_norm": v2["line_norm"], "nh10": 20.0, "lum10": v2["log_lint_2_10"] - 1.0}
    best = fits["c2"]["M2"]["statistic"]
    grid = []
    for d2 in (2.0, 1.6, 1.2, 1.0, 0.8, 0.6, 0.4, 0.2, 0.0):
        load(keys4)
        build(keys4, "M2", nh_gal, dict(st, lum10=st["lum"] - d2))
        par(4, "lum").link = f"p{par(1, 'lum').index} - ({d2:.6f})"
        stat = fit(2)
        grid.append({"delta2": d2, "delta_stat": stat - best, "nh10_1e22": par(4, "nh").values[0]})
        print("M2 grid", grid[-1], flush=True)
    out["m2_delta2_grid"] = {"best_statistic": best, "best_delta2": fits["c2"]["M2"]["delta2"]["value"], "grid": grid,
                             "note": "N_H(2010) free at every point; reflection, line and scattered component tied "
                                     "to 2022 in absolute terms"}
    aic = dict(fits["aic_2022"], cutoff=rec["aic"], reflection_only=out["reflection_only_2022"]["aic"])
    out["aic_2022_all"] = aic
    out["aic_relative_to_m0"] = {k: a - aic["M0"] for k, a in aic.items()}
    out["lowest_aic"] = min(aic, key=aic.get)
    (config.RESULTS_DIR / "s7_posthoc.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out["aic_relative_to_m0"]), flush=True)


def run_plotdata() -> None:
    """Plot arrays for every model (refits from the stored best values)."""
    setup()
    prepare()
    fits = json.loads((config.RESULTS_DIR / "s7_fits.json").read_text())
    nh_gal = fits["galactic_nh_1e22"]
    out = {}
    for model in ("M0", "M1", "M2"):
        v = fits["c2"][model]["values"]
        start = {"nh": v["nh_1e22"], "gamma": v["gamma"], "lum": v["log_lint_2_10"],
                 "lum_s": v.get("log_lscat_2_10", 40.5), "refl_norm": v.get("refl_norm", 1e-3),
                 "line_norm": v.get("line_norm", 3e-6), "nh10": v["nh10_1e22"], "lum10": v["log_lint10_2_10"]}
        keys = ["xrt22", "fpma", "fpmb", "xrt10"]
        load(keys)
        build(keys, model, nh_gal, start)
        stat = fit(1)
        rec = {"c2": plot_arrays(keys), "c2_statistic": stat}
        v = fits["c1"][model]["values"]
        start = {"nh": v["nh_1e22"], "gamma": v["gamma"], "lum": v["log_lint_2_10"],
                 "lum_s": v.get("log_lscat_2_10", 40.5), "refl_norm": v.get("refl_norm", 1e-3),
                 "line_norm": v.get("line_norm", 3e-6)}
        keys = ["xrt22", "fpma", "fpmb", "bat"]
        load(keys)
        build(keys, model, nh_gal, start)
        par(4, "const").values = [v["const_bat"], 0.01, 0.05, 0.05, 20.0, 20.0]
        rec["c1_statistic"] = fit(1)
        rec["bat"] = plot_arrays(keys, rebin=(0, 1))["bat"]
        out[model] = rec
        print(model, rec["c2_statistic"], fits["c2"][model]["statistic"], rec["c1_statistic"],
              fits["c1"][model]["statistic"], flush=True)
    (config.RESULTS_DIR / "s7_spectra_plot.json").write_text(json.dumps(out) + "\n")


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else "fits"
    if stage == "fits":
        run_fits()
    elif stage == "arf":
        run_arf()
    elif stage == "plotdata":
        run_plotdata()
    elif stage == "lumcheck":
        run_lumcheck()
    elif stage == "posthoc":
        run_posthoc()
    elif stage == "sims":
        run_sims(int(sys.argv[2]) if len(sys.argv) > 2 else 200)
    elif stage == "mcmc":
        run_mcmc()
    else:
        raise SystemExit(__doc__)
