"""Stage 0 (corrections), stage 3 (physical quantities) and stage 4 (detectors,
injection tests) analyses. Called by scripts/run_stage34.py; every function
works from the local cache only."""

from __future__ import annotations

import json

import numpy as np

from . import analysis, config, differential, injection, mlcompare, physics, theory, untimely
from .analysis import mjd_to_date, visit_arrays
from .config import BANDS, FLARE, SLOW_RISE
from .epochs import MEDIAN_SE_FACTOR
from .flare import detect_flares
from .lightcurve import LightCurve, load_target, mjd_to_year
from .slowrise import detect_slow_rise
from .stats import constant_chi2, pooled_error_floor

# Visit error used from stage 3 on: the in-visit scatter estimate plus the
# uncertainty of the field offset (see README section 5.1 for why).
ERR_KEY = "err_scatter"
N_DRAWS = 4000

# AGN bolometric luminosities of NGC 3786 from the literature, log10(erg/s).
LITERATURE_LBOL = (
    {"label": "Son+2022, before the flare (W1-based)", "log_lbol": 42.89, "date": None,
     "source": "Son et al. 2022, ApJ 937, 3, Sec. 3.3"},
    {"label": "Son+2022, after the flare (W1-based)", "log_lbol": 43.35, "date": None,
     "source": "Son et al. 2022, ApJ 937, 3, Sec. 3.3"},
    {"label": "Son+2026, broad H-alpha, 2022-02-02", "log_lbol": 43.31, "date": "2022-02-02",
     "source": "Son et al. 2026, arXiv:2604.17173, Sec. 4"},
    {"label": "Son+2026, broad H-alpha, 2023-06-10", "log_lbol": 43.24, "date": "2023-06-10",
     "source": "Son et al. 2026, arXiv:2604.17173, Sec. 4"},
    {"label": "Son+2026, broad H-alpha, 2023-12-23", "log_lbol": 43.34, "date": "2023-12-23",
     "source": "Son et al. 2026, arXiv:2604.17173, Sec. 4"},
    {"label": "Son+2026, extinction-corrected, 2022-02-02", "log_lbol": 44.05,
     "date": "2022-02-02", "source": "Son et al. 2026, arXiv:2604.17173, Sec. 4"},
    {"label": "Son+2026, extinction-corrected, 2023-06-10", "log_lbol": 44.39,
     "date": "2023-06-10", "source": "Son et al. 2026, arXiv:2604.17173, Sec. 4"},
    {"label": "Son+2026, extinction-corrected, 2023-12-23", "log_lbol": 44.64,
     "date": "2023-12-23", "source": "Son et al. 2026, arXiv:2604.17173, Sec. 4"},
)


# --- loading ------------------------------------------------------------------
def load_fields(allow_partial: bool = False) -> dict:
    """Light curves grouped by field: NGC 3786 with its 30 reference stars and
    each control galaxy with its comparison stars. ``allow_partial`` (for
    development only) skips comparison stars that are not in the cache yet."""
    from . import irsa
    t = analysis.load_ngc3786()
    fields = {"ngc3786": {
        "target": load_target(t["name"], t["ra"], t["dec"]),
        "stars": [load_target(s["name"], s["ra"], s["dec"])
                  for s in analysis.load_reference_stars()],
    }}
    field_stars = analysis.load_field_stars()
    for c in analysis.load_controls():
        fields[c["name"]] = {
            "target": load_target(c["name"], c["ra"], c["dec"]),
            "stars": [load_target(s["name"], s["ra"], s["dec"], surveys=("neowise",))
                      for s in field_stars if s["field"] == c["name"]
                      and (not allow_partial or irsa.cache_paths(s["name"])["neowise"].exists())],
        }
        if allow_partial and len(fields[c["name"]]["stars"]) < differential.MIN_STARS + 1:
            fields[c["name"]]["stars"] = []
    return fields


def star_visits(stars: list[LightCurve]) -> list[dict]:
    return [visit_arrays(lc, "neowise") for lc in stars]


def corrected_target(field: dict, floor_mag: dict | None = None) -> dict:
    offsets = differential.field_offsets(star_visits(field["stars"]))
    return differential.apply_offsets(visit_arrays(field["target"], "neowise"), offsets,
                                      err_key=ERR_KEY, floor_mag=floor_mag)


def corrected_stars(field: dict) -> list[dict]:
    """Leave-one-out corrected star light curves of one field."""
    return differential.leave_one_out(star_visits(field["stars"]), err_key=ERR_KEY)


# --- stage 0a: differential photometry -----------------------------------------
def _chi2_table(before: list[dict], after: list[dict], model: str) -> dict:
    out = {}
    for b in BANDS:
        key = f"{b}_{model}" if model != "err" else f"{b}_err"
        c0 = [constant_chi2(v[f"{b}_mag"], v[key]) for v in before]
        c1 = [constant_chi2(v[f"{b}_mag"], v[f"{b}_err"]) for v in after]
        out[b] = {
            "chi2_red_median_before": float(np.median([c.chi2_red for c in c0])),
            "chi2_red_median_after": float(np.median([c.chi2_red for c in c1])),
            "chi2_red_mean_before": float(np.mean([c.chi2_red for c in c0])),
            "chi2_red_mean_after": float(np.mean([c.chi2_red for c in c1])),
            "n_p_below_0.01_before": int(sum(c.p_value < 0.01 for c in c0)),
            "n_p_below_0.01_after": int(sum(c.p_value < 0.01 for c in c1)),
            "pooled_floor_mag_before": pooled_error_floor(
                [v[f"{b}_mag"] for v in before], [v[key] for v in before]),
            "pooled_floor_mag_after": pooled_error_floor(
                [v[f"{b}_mag"] for v in after], [v[f"{b}_err"] for v in after]),
            "rms_median_mag_before": float(np.median([c.rms for c in c0])),
            "rms_median_mag_after": float(np.median([c.rms for c in c1])),
        }
    return out


def differential_validation(stars: list[LightCurve], seed: int = config.SEED) -> dict:
    """Does subtracting the field offset bring the reference stars' reduced
    chi-square to 1? Leave-one-out, for the three visit-error models."""
    sv = star_visits(stars)
    out: dict = {"n_stars": len(sv), "models": {}}
    for model in analysis.ERROR_MODELS:
        loo = differential.leave_one_out(sv, err_key=model)
        out["models"][model] = _chi2_table(sv, loo, model)
    offs = differential.field_offsets(sv)
    out["offsets"] = {
        "visit_date": [mjd_to_date(m) for m in offs.mjd],
        **{f"{b}_offset_mmag": (1e3 * offs.offset[b]).tolist() for b in BANDS},
        **{f"{b}_offset_err_mmag": (1e3 * offs.offset_err[b]).tolist() for b in BANDS},
    }
    # Do stars of different brightness see the same offset?
    w1 = np.array([np.median(v["w1_mag"]) for v in sv])
    order = np.argsort(w1)
    thirds = np.array_split(order, 3)
    per_bin = [differential.field_offsets([sv[i] for i in idx]) for idx in thirds]
    out["brightness_dependence"] = {
        "bins_w1_mag": [[float(w1[idx].min()), float(w1[idx].max())] for idx in thirds]}
    for b in BANDS:
        diff = per_bin[0].offset[b] - per_bin[2].offset[b]
        err = np.hypot(per_bin[0].offset_err[b], per_bin[2].offset_err[b])
        out["brightness_dependence"][b] = {
            "bright_minus_faint_median_mmag": float(1e3 * np.median(diff)),
            "bright_minus_faint_rms_mmag": float(1e3 * np.sqrt(np.mean(diff**2))),
            "expected_rms_from_errors_mmag": float(1e3 * np.sqrt(np.mean(err**2))),
        }
    # Ten comparison stars (as in the control fields) applied to the other twenty.
    rng = np.random.default_rng(seed)
    pick = rng.permutation(len(sv))
    ten, rest = [sv[i] for i in pick[:10]], [sv[i] for i in pick[10:]]
    offs10 = differential.field_offsets(ten)
    after = [differential.apply_offsets(v, offs10, err_key=ERR_KEY) for v in rest]
    out["ten_star_test"] = _chi2_table(rest, after, ERR_KEY)
    out["adopted_error_model"] = ERR_KEY
    return out


def field_star_validation(fields: dict) -> dict:
    """Leave-one-out chi-square of the comparison stars in every control field."""
    out = {}
    for name, field in fields.items():
        if name == "ngc3786" or not field["stars"]:
            continue
        sv = star_visits(field["stars"])
        loo = differential.leave_one_out(sv, err_key=ERR_KEY)
        offs = differential.field_offsets(sv)
        out[name] = {
            "n_stars": len(sv), "n_field_visits": int(offs.mjd.size),
            **{b: {**_chi2_table(sv, loo, ERR_KEY)[b],
                   "offset_std_mmag": float(1e3 * np.nanstd(offs.offset[b], ddof=1)),
                   "offset_max_abs_mmag": float(1e3 * np.nanmax(np.abs(offs.offset[b]))),
                   "offset_err_median_mmag": float(1e3 * np.nanmedian(offs.offset_err[b]))}
               for b in BANDS},
        }
    return out


# --- stage 0b: unTimely ---------------------------------------------------------
def untimely_comparison(fields: dict) -> dict:
    """Compare the profile-fit visit medians with unTimely coadd photometry."""
    out = {}
    for name, field in fields.items():
        if not untimely.cache_path(name).exists():
            out[name] = {"status": "no unTimely file"}
            continue
        u = untimely.load(name)
        v = visit_arrays(field["target"], "neowise")
        pairs = untimely.match_visits(v["mjd"], u["mjd"])
        i = np.array([p[0] for p in pairs], dtype=int)
        j = np.array([p[1] for p in pairs], dtype=int)
        row: dict = {"status": "ok", "n_untimely_epochs": int(u["mjd"].size),
                     "n_matched_neowise_visits": len(pairs),
                     "untimely_last_date": mjd_to_date(np.nanmax(u["mjd"])),
                     "max_detections_within_5arcsec": int(max(u["w1_n_det"].max(),
                                                              u["w2_n_det"].max())),
                     "flags_info_values": sorted({int(x) for b in BANDS
                                                  for x in u[f"{b}_flags_info"]})}
        for b in BANDS:
            mp, mu = v[f"{b}_mag"][i], u[f"{b}_mag"][j]
            ok = np.isfinite(mp) & np.isfinite(mu)
            d = (mu - mp)[ok]
            dp, du = mp[ok] - np.median(mp[ok]), mu[ok] - np.median(mu[ok])
            slope = float(np.sum(dp * du) / np.sum(dp**2)) if np.sum(dp**2) > 0 else np.nan
            row[b] = {
                "median_untimely_minus_pipeline_mag": float(np.median(d)),
                "rms_about_median_mag": float(np.std(d, ddof=1)),
                "pipeline_rms_mag": float(np.std(mp[ok], ddof=1)),
                "untimely_rms_mag": float(np.std(mu[ok], ddof=1)),
                "slope_untimely_vs_pipeline": slope,
            }
        out[name] = row
    return out


def untimely_flare(field: dict, res_a) -> dict:
    """NGC 3786: amplitude and excess flux at the flare onset, measured the
    same way on both photometries (baseline = the same pre-flare visits)."""
    v = visit_arrays(field["target"], "neowise")
    u = untimely.load("ngc3786")
    onset = res_a.onset_index
    base_idx = np.flatnonzero(~res_a.flagged & (np.arange(v["mjd"].size) < onset))
    pairs = dict(untimely.match_visits(v["mjd"], u["mjd"]))
    base_pairs = [(i, pairs[i]) for i in base_idx if i in pairs]
    out: dict = {"n_baseline_visits_matched": len(base_pairs),
                 "onset_matched": onset in pairs,
                 "onset_date_pipeline": mjd_to_date(v["mjd"][onset])}
    if onset not in pairs:
        return out
    out["onset_date_untimely"] = mjd_to_date(u["mjd"][pairs[onset]])
    excess = {}
    for b in BANDS:
        ip = np.array([p[0] for p in base_pairs])
        iu = np.array([p[1] for p in base_pairs])
        row = {}
        for label, mag, flux, k, idx in (
            ("pipeline", v[f"{b}_mag"], v[f"{b}_flux_jy"], onset, ip),
            ("untimely", u[f"{b}_mag"], u[f"{b}_flux_jy"], pairs[onset], iu),
        ):
            row[f"{label}_baseline_mag"] = float(np.median(mag[idx]))
            row[f"{label}_onset_mag"] = float(mag[k])
            row[f"{label}_dmag"] = float(np.median(mag[idx]) - mag[k])
            row[f"{label}_excess_flux_mjy"] = float(1e3 * (flux[k] - np.median(flux[idx])))
        excess[b] = row
        out[b] = row
    for label in ("pipeline", "untimely"):
        r = excess["w1"][f"{label}_excess_flux_mjy"] / excess["w2"][f"{label}_excess_flux_mjy"]
        from .photometry import blackbody_temperature
        out[f"{label}_excess_ratio_w1_w2"] = float(r)
        out[f"{label}_excess_temperature_k"] = blackbody_temperature(r)
    return out


def followup_comparisons(field: dict) -> dict:
    """Two checks on NGC 3786 that became possible in stage 3-4.

    * Son et al. (2026, Sec. 2) report a mean brightening of 0.68 mag (W1) and
      1.02 mag (W2) relative to the pre-outburst state (MJD < 59000) from
      6-month bins. The same quantity from this pipeline, with a mean and with
      a median reference.
    * Stage 2 found the 2010 (AllWISE) colour redder than the 2014-2020
      (NEOWISE) baseline and could not tell a real change from a difference
      between the two tables. unTimely measures both eras the same way.
    """
    raw = visit_arrays(field["target"], "neowise")
    cor = corrected_target(field)
    out: dict = {"son2026_reported_mean_brightening_mag": {"w1": 0.68, "w2": 1.02},
                 "pre_outburst_mjd_limit": 59000.0}
    for label, v in (("uncorrected", raw), ("corrected", cor)):
        pre, post = v["mjd"] < 59000.0, v["mjd"] >= 59000.0
        out[label] = {"n_pre": int(pre.sum()), "n_post": int(post.sum())}
        for b in BANDS:
            m = v[f"{b}_mag"]
            out[label][f"{b}_mean_minus_mean"] = float(m[pre].mean() - m[post].mean())
            out[label][f"{b}_median_minus_mean"] = float(np.median(m[pre]) - m[post].mean())
    e = field["target"].epochs
    aw = field["target"].select("allwise")
    u = untimely.load("ngc3786")
    early, quiet = u["mjd"] < 56000.0, (u["mjd"] > 56000.0) & (u["mjd"] < 59000.0)
    u_color = u["w1_mag"] - u["w2_mag"]
    out["color_2010_vs_2014_2020"] = {
        "pipeline_allwise_2010_colors": [float(c) for c in e["color"][aw]],
        "pipeline_neowise_pre_outburst_median_color": float(
            np.median(raw["color"][raw["mjd"] < 59000.0])),
        "untimely_2010_colors": [float(c) for c in u_color[early]],
        "untimely_2014_2020_pre_outburst_median_color": float(np.median(u_color[quiet])),
        "untimely_2010_dates": [mjd_to_date(m) for m in u["mjd"][early]],
    }
    return out


# --- stage 3: physical quantities ------------------------------------------------
def _date_to_mjd(date: str) -> float:
    from datetime import datetime
    return (datetime.strptime(date, "%Y-%m-%d") - datetime(1858, 11, 17)).days


def physics_analysis(v: dict, z: float, seed: int = config.SEED) -> dict:
    """T, L_dust and R_bb of the variable component at every flagged visit of
    the (differentially corrected) NGC 3786 light curve."""
    res = detect_flares(v["mjd"], v["w1_mag"], v["w1_err"], v["w2_mag"], v["w2_err"], FLARE)
    d_l = physics.luminosity_distance_mpc(z)
    n = v["mjd"].size
    base = ~res.flagged & (np.arange(n) < res.onset_index)
    rng = np.random.default_rng(seed)

    def baseline(mask: np.ndarray) -> tuple[dict, dict]:
        val = {b: float(np.median(v[f"{b}_flux_jy"][mask])) for b in BANDS}
        if mask.sum() > 1:
            err = {b: MEDIAN_SE_FACTOR * float(np.std(v[f"{b}_flux_jy"][mask], ddof=1))
                   / np.sqrt(mask.sum()) for b in BANDS}
        else:
            err = {b: float(v[f"{b}_flux_err_jy"][mask][0]) for b in BANDS}
        return val, err

    base_val, base_err = baseline(base)
    base_draws = {b: rng.normal(base_val[b], base_err[b], N_DRAWS) for b in BANDS}
    # Alternative baselines (systematic check): the faint half of the pre-flare
    # visits, and the single last pre-flare visit.
    faint = base & (v["w2_mag"] >= np.median(v["w2_mag"][base]))
    last = np.zeros(n, dtype=bool)
    last[np.flatnonzero(base)[-1]] = True
    alt = {"faint_half": baseline(faint)[0], "last_pre_flare_visit": baseline(last)[0]}

    visits = []
    for k in np.flatnonzero(res.flagged):
        flux = {b: float(v[f"{b}_flux_jy"][k]) for b in BANDS}
        ferr = {b: float(v[f"{b}_flux_err_jy"][k]) for b in BANDS}
        d = {b: flux[b] - base_val[b] for b in BANDS}
        point = physics.blackbody_quantities(d["w1"], d["w2"], z, d_l)
        cc = physics.blackbody_quantities(d["w1"], d["w2"], z, d_l, color_correct=True)
        mc = physics.monte_carlo(flux, ferr, base_val, base_err, z, d_l, N_DRAWS, rng, base_draws)
        row = {
            "date": mjd_to_date(v["mjd"][k]), "mjd": float(v["mjd"][k]),
            "year": float(mjd_to_year(v["mjd"][k])),
            "w1_mag": float(v["w1_mag"][k]), "w2_mag": float(v["w2_mag"][k]),
            "w1_excess_mjy": 1e3 * d["w1"], "w2_excess_mjy": 1e3 * d["w2"],
            "w1_excess_err_mjy": 1e3 * float(np.hypot(ferr["w1"], base_err["w1"])),
            "w2_excess_err_mjy": 1e3 * float(np.hypot(ferr["w2"], base_err["w2"])),
            **point, **mc,
            "color_corrected": {k2: cc[k2] for k2 in ("t_rest_k", "l_dust_erg_s",
                                                      "r_bb_light_days")},
            "alt_baseline": {},
        }
        for label, bv in alt.items():
            q = physics.blackbody_quantities(flux["w1"] - bv["w1"], flux["w2"] - bv["w2"], z, d_l)
            row["alt_baseline"][label] = {k2: q[k2] for k2 in ("t_rest_k", "l_dust_erg_s",
                                                               "r_bb_light_days")}
        visits.append(row)

    # Sublimation radii for the literature luminosities, and the nearest visit.
    rsub = []
    for lit in LITERATURE_LBOL:
        l_agn = 10.0 ** lit["log_lbol"]
        row = {**lit}
        for grain in ("graphite", "silicate"):
            row[f"r_sub_{grain}_light_days"] = theory.light_travel_days(
                theory.sublimation_radius_pc(l_agn, grain))
        if lit["date"] is not None:
            near = min(visits, key=lambda r: abs(r["mjd"] - _date_to_mjd(lit["date"])))
            row["nearest_visit"] = near["date"]
            row["nearest_visit_gap_days"] = float(abs(near["mjd"] - _date_to_mjd(lit["date"])))
            row["r_bb_light_days"] = near["r_bb_light_days"]
            row["r_bb_over_r_sub_graphite"] = near["r_bb_light_days"] / row["r_sub_graphite_light_days"]
            row["r_bb_over_r_sub_silicate"] = near["r_bb_light_days"] / row["r_sub_silicate_light_days"]
            row["l_dust_over_l_bol"] = near["l_dust_erg_s"] / l_agn
        rsub.append(row)

    t = np.array([r["t_rest_k"] for r in visits])
    lum = np.array([r["l_dust_erg_s"] for r in visits])
    rad = np.array([r["r_bb_light_days"] for r in visits])
    slope = float(np.polyfit(np.log10(lum), np.log10(rad), 1)[0]) if len(visits) > 2 else np.nan
    onset_prev = float(v["mjd"][res.onset_index] - v["mjd"][res.onset_index - 1])
    return {
        "z": z, "luminosity_distance_mpc": d_l,
        "cosmology": {"h0_km_s_mpc": physics.H0_KM_S_MPC, "omega_m": physics.OMEGA_M,
                      "note": "flat LCDM as in Son et al. (2022); no peculiar-velocity correction"},
        "n_draws": N_DRAWS, "seed": seed, "error_model": ERR_KEY,
        "baseline": {
            "n_visits": int(base.sum()),
            "span": [mjd_to_date(v["mjd"][base].min()), mjd_to_date(v["mjd"][base].max())],
            **{f"{b}_flux_mjy": 1e3 * base_val[b] for b in BANDS},
            **{f"{b}_flux_err_mjy": 1e3 * base_err[b] for b in BANDS},
            "alt_faint_half_n_visits": int(faint.sum()),
        },
        "flagged_visits": visits,
        "summary": {
            "t_rest_k_min": float(t.min()), "t_rest_k_max": float(t.max()),
            "t_rest_k_mean": float(t.mean()), "t_rest_k_std": float(t.std(ddof=1)),
            "l_dust_ratio_last_over_first": float(lum[-1] / lum[0]),
            "r_bb_ratio_last_over_first": float(rad[-1] / rad[0]),
            "dlog_r_dlog_l": slope,
            "onset_light_travel_limit_light_days": onset_prev,
        },
        "sublimation_radius": rsub,
    }


# --- stage 3: hypotheses of docs/PREDICTIONS.md -----------------------------------
SIGMA_SYS_K = 70.0  # pipeline vs. unTimely temperature at the onset visit
H1_MIN_DECLINE_K = 100.0
H2_RANGE_K = (1000.0, 1200.0)
H3_LIMIT_K = 1500.0


def hypothesis_tests(phys: dict) -> dict:
    """Apply the decision rules of docs/PREDICTIONS.md section 2 to the
    temperatures of the flagged visits. The rules were committed before this
    function was written (commit 0aefbe4); the temperatures themselves had
    been seen earlier."""
    from scipy.stats import chi2 as chi2_dist
    from scipy.stats import t as t_dist

    from .slowrise import mann_kendall_z

    rows = phys["flagged_visits"]
    year = np.array([r["year"] for r in rows])
    temp = np.array([r["t_rest_k"] for r in rows])
    s_stat = np.array([(r["t_rest_k_p84"] - r["t_rest_k_p16"]) / 2.0 for r in rows])
    s_tot = np.hypot(s_stat, SIGMA_SYS_K)
    n = temp.size
    span = float(year[-1] - year[0])

    # H1: T decreases with time.
    x = year - year.mean()
    b = float((x * (temp - temp.mean())).sum() / (x**2).sum())
    resid = temp - temp.mean() - b * x
    se = float(np.sqrt((resid**2).sum() / (n - 2) / (x**2).sum()))
    tcrit = float(t_dist.ppf(0.975, n - 2))
    b_lo, b_hi = b - tcrit * se, b + tcrit * se
    dt_lo, dt_hi = b_lo * span, b_hi * span
    if b_hi < 0:
        v1 = "supported"
    elif dt_lo > -H1_MIN_DECLINE_K:
        v1 = "rejected"
    else:
        v1 = "undetermined"
    h1 = {"slope_k_per_yr": b, "slope_se_k_per_yr": se, "slope_ci95_k_per_yr": [b_lo, b_hi],
          "span_yr": span, "delta_t_k": b * span, "delta_t_ci95_k": [dt_lo, dt_hi],
          "mann_kendall_z": mann_kendall_z(temp), "residual_rms_k": float(np.sqrt((resid**2).mean())),
          "verdict": v1}

    # H2: T nearly constant near 1000-1200 K.
    def constant(sig: np.ndarray) -> dict:
        w = 1.0 / sig**2
        mean = float((w * temp).sum() / w.sum())
        chi2 = float((w * (temp - mean) ** 2).sum())
        return {"weighted_mean_k": mean, "mean_err_k": float(1.0 / np.sqrt(w.sum())),
                "chi2": chi2, "dof": n - 1, "p_value": float(chi2_dist.sf(chi2, n - 1))}

    c_stat, c_tot = constant(s_stat), constant(s_tot)
    mean, mean_err = c_stat["weighted_mean_k"], c_stat["mean_err_k"]
    inside = H2_RANGE_K[0] <= mean <= H2_RANGE_K[1]
    far_outside = (mean < H2_RANGE_K[0] - 2 * mean_err) or (mean > H2_RANGE_K[1] + 2 * mean_err)
    if inside and c_stat["p_value"] > 0.05:
        v2 = "supported"
    elif far_outside or c_tot["p_value"] < 0.01:
        v2 = "rejected"
    else:
        v2 = "undetermined"
    h2 = {"statistical_errors_only": c_stat, "with_systematic_70k": c_tot,
          "mean_inside_1000_1200": bool(inside),
          "t_min_k": float(temp.min()), "t_max_k": float(temp.max()),
          "t_std_k": float(temp.std(ddof=1)),
          "intrinsic_scatter_k": float(np.sqrt(max(temp.var(ddof=1) - (s_stat**2).mean(), 0.0))),
          "verdict": v2}

    # H3: early T does not exceed ~1500 K; later re-rise possible.
    early = slice(0, 2)
    upper = temp[early] + 2 * s_tot[early]
    lower = temp[early] - 2 * s_tot[early]
    if np.all(upper < H3_LIMIT_K):
        v3 = "supported"
    elif np.any(lower > H3_LIMIT_K):
        v3 = "rejected"
    else:
        v3 = "undetermined"
    k_min = int(np.argmin(temp[early]))
    later = np.arange(2, n)
    excess_stat = (temp[later] - temp[k_min]) / np.hypot(s_stat[later], s_stat[k_min])
    excess_tot = (temp[later] - temp[k_min]) / np.hypot(s_tot[later], s_tot[k_min])
    h3 = {"early_dates": [rows[i]["date"] for i in range(2)],
          "early_t_k": temp[early].tolist(), "early_upper_2sigma_k": upper.tolist(),
          "early_lower_2sigma_k": lower.tolist(),
          "rerise_dates_2sigma_statistical": [rows[j]["date"] for j, e in zip(later, excess_stat)
                                              if e >= 2],
          "rerise_dates_2sigma_with_systematic": [rows[j]["date"] for j, e in zip(later, excess_tot)
                                                  if e >= 2],
          "max_later_t_k": float(temp[later].max()),
          "verdict": v3,
          "note": "the first clause is also satisfied under H1 and H2; the first visit is up "
                  "to 206 days after the onset"}
    return {"rules": "docs/PREDICTIONS.md section 2 (commit 0aefbe4)",
            "dates": [r["date"] for r in rows], "t_rest_k": temp.tolist(),
            "sigma_stat_k": s_stat.tolist(), "sigma_sys_k": SIGMA_SYS_K,
            "H1_temperature_mixing": h1, "H2_silicate": h2, "H3_dust_destruction": h3}


def optical_comparison(v: dict, base_w1_flux_mjy: float, name: str = "ngc3786") -> dict:
    """Alternative of docs/PREDICTIONS.md: did the optical (ZTF g, r) brighten
    with the mid-IR? ``v`` is the corrected visit table of the galaxy and
    ``base_w1_flux_mjy`` its pre-outburst W1 flux."""
    from scipy.stats import spearmanr

    from . import ztf

    out: dict = {"rules": "docs/PREDICTIONS.md section 2 (commit 0aefbe4)", "bands": {}}
    excess_w1 = 1e3 * v["w1_flux_jy"] - base_w1_flux_mjy
    post = v["mjd"] >= 59000.0
    passes = True
    for band in ("g", "r"):
        path = ztf.cache_path(name, band)
        if not path.exists():
            out["bands"][band] = {"status": "no ZTF file"}
            out["verdict"] = "unconfirmed (no data)"
            return out
        lc = ztf.read_clean(path)
        bins = ztf.bin_to_visits(lc["mjd"], lc["mag"], v["mjd"])
        ok = np.isfinite(bins["mag"])
        pre_bins, post_bins = bins["mag"][ok & ~post], bins["mag"][ok & post]
        brightening = float(np.median(pre_bins) - np.median(post_bins))
        pre_scatter = float(np.std(pre_bins, ddof=1)) if pre_bins.size > 1 else float("nan")
        near = ok & post & (np.abs(bins["mjd"] - v["mjd"]) <= 60.0)
        if near.sum() >= 4:
            rho, p = spearmanr(-bins["mag"][near], excess_w1[near])
            rho, p = float(rho), float(p)
        else:
            rho, p = float("nan"), float("nan")
        band_pass = (brightening >= 3 * pre_scatter) and (rho > 0) and (p < 0.05)
        passes &= bool(band_pass)
        out["bands"][band] = {
            "status": "ok", "oid": lc["oid"], "n_rows_total": lc["n_rows_total"],
            "n_rows_clean": lc["n_rows_clean"], "n_rows_used": lc["n_rows_used"],
            "mjd_range": [float(lc["mjd"].min()), float(lc["mjd"].max())],
            "date_range": [mjd_to_date(lc["mjd"].min()), mjd_to_date(lc["mjd"].max())],
            "bins": [{"wise_visit": mjd_to_date(v["mjd"][k]), "n": int(bins["n"][k]),
                      "mag": float(bins["mag"][k]), "err": float(bins["err"][k]),
                      "offset_days": float(bins["mjd"][k] - v["mjd"][k]),
                      "w1_excess_mjy": float(excess_w1[k]), "post_outburst": bool(post[k])}
                     for k in range(v["mjd"].size) if ok[k]],
            "n_pre_bins": int(pre_bins.size), "n_post_bins": int(post_bins.size),
            "brightening_mag": brightening, "pre_outburst_bin_scatter_mag": pre_scatter,
            "brightening_over_scatter": brightening / pre_scatter if pre_scatter > 0 else None,
            "post_bin_range_mag": float(post_bins.max() - post_bins.min()),
            "n_pairs_for_correlation": int(near.sum()),
            "spearman_rho": rho, "spearman_p": p, "passes": bool(band_pass),
        }
    out["verdict"] = "supported" if passes else "undetermined"
    if name == "ngc3786":
        out["son2026_reported_mean_brightening_mag"] = {"g": 0.11, "r": 0.14}
    return out


def exploratory_temperature_checks(phys: dict, optical: dict) -> dict:
    """NOT part of docs/PREDICTIONS.md; added after the pre-specified tests
    were run. Does the temperature of the variable component follow the
    brightness? T and L_dust come from the same two fluxes, so the optical
    magnitude (an independent instrument) is the cleaner comparison."""
    from scipy.stats import spearmanr

    rows = phys["flagged_visits"]
    temp = np.array([r["t_rest_k"] for r in rows])
    out: dict = {"status": "exploratory, defined after seeing the results"}
    for key in ("w1_excess_mjy", "w2_excess_mjy", "l_dust_erg_s"):
        rho, p = spearmanr(temp, [r[key] for r in rows])
        out[f"spearman_t_vs_{key}"] = {"rho": float(rho), "p": float(p), "n": len(rows)}
    lum = np.array([r["l_dust_erg_s"] for r in rows])
    out["dlog_t_dlog_l_dust"] = float(np.polyfit(np.log10(lum), np.log10(temp), 1)[0])
    for band, info in optical.get("bands", {}).items():
        if info.get("status") != "ok":
            continue
        by_date = {b["wise_visit"]: b["mag"] for b in info["bins"]}
        pairs = [(t, by_date[r["date"]]) for t, r in zip(temp, rows) if r["date"] in by_date]
        if len(pairs) >= 4:
            rho, p = spearmanr([a for a, _ in pairs], [-m for _, m in pairs])
            out[f"spearman_t_vs_ztf_{band}_brightness"] = {"rho": float(rho), "p": float(p),
                                                          "n": len(pairs)}
    return out


# --- stage 4: detectors -----------------------------------------------------------
def run_both(v: dict) -> dict:
    res_a = detect_flares(v["mjd"], v["w1_mag"], v["w1_err"], v["w2_mag"], v["w2_err"], FLARE)
    res_b = detect_slow_rise(v["mjd"], v["w1_mag"], v["w1_err"], v["w2_mag"], v["w2_err"],
                             SLOW_RISE)
    b = res_b.as_dict()
    if np.isfinite(b["onset_mjd"]):
        b["onset_date"] = mjd_to_date(b["onset_mjd"])
    return {
        "n_visits": int(v["mjd"].size),
        "a_n_flagged": int(res_a.flagged.sum()),
        "a_flagged_dates": [mjd_to_date(m) for m in res_a.mjd[res_a.flagged]],
        "a_w1_max_significance": res_a.max_significance("w1"),
        "a_w2_max_significance": res_a.max_significance("w2"),
        "b": b,
    }


def detector_analysis(fields: dict) -> dict:
    """Detectors (a) and (b) on the real, differentially corrected light curves."""
    out: dict = {"criteria": {
        "a": {"n_sigma": FLARE.n_sigma, "min_amp_mag": FLARE.min_amp_mag,
              "min_baseline_epochs": FLARE.min_baseline_epochs},
        "b": {"p_max": SLOW_RISE.p_max, "min_cum_mag": SLOW_RISE.min_cum_mag,
              "min_visits": SLOW_RISE.min_visits,
              "min_baseline_visits": SLOW_RISE.min_baseline_visits,
              "min_ramp_visits": SLOW_RISE.min_ramp_visits,
              "tau_grid_yr": list(SLOW_RISE.tau_grid_yr)}},
        "galaxies": {}, "galaxies_uncorrected": {}}
    star_rows = []
    for name, field in fields.items():
        raw = visit_arrays(field["target"], "neowise")
        if field["stars"]:
            out["galaxies"][name] = run_both(corrected_target(field))
        out["galaxies_uncorrected"][name] = run_both(raw)
        if name == "ngc3786" and field["stars"]:
            v = corrected_target(field)
            cut = v["mjd"] < analysis.PAPER_DATA_END_MJD
            out["ngc3786_paper_range"] = run_both({k: val[cut] for k, val in v.items()})
        raw_stars = star_visits(field["stars"])
        for lc, v_raw, v_cor in zip(field["stars"], raw_stars, corrected_stars(field)):
            r_cor, r_raw = run_both(v_cor), run_both(v_raw)
            star_rows.append({
                "name": lc.name, "field": name, "n_visits": r_cor["n_visits"],
                "a_flagged": r_cor["a_n_flagged"] > 0, "b_flagged": r_cor["b"]["flagged"],
                "b_p_w1": r_cor["b"].get("w1_p"), "b_p_w2": r_cor["b"].get("w2_p"),
                "b_dm_cum_w1": r_cor["b"].get("w1_dm_cum"), "b_dm_cum_w2": r_cor["b"].get("w2_dm_cum"),
                "raw_a_flagged": r_raw["a_n_flagged"] > 0, "raw_b_flagged": r_raw["b"]["flagged"],
                "raw_b_p_w1": r_raw["b"].get("w1_p"), "raw_b_p_w2": r_raw["b"].get("w2_p"),
            })
    out["stars"] = {
        "n": len(star_rows),
        "n_flagged_a": int(sum(r["a_flagged"] for r in star_rows)),
        "n_flagged_b": int(sum(r["b_flagged"] for r in star_rows)),
        # How often would the F-test alone (no amplitude floor) fire in both bands?
        "n_p_both_below_p_max_corrected": int(sum(
            (r["b_p_w1"] or 1) <= SLOW_RISE.p_max and (r["b_p_w2"] or 1) <= SLOW_RISE.p_max
            for r in star_rows)),
        "n_p_both_below_p_max_uncorrected": int(sum(
            (r["raw_b_p_w1"] or 1) <= SLOW_RISE.p_max and (r["raw_b_p_w2"] or 1) <= SLOW_RISE.p_max
            for r in star_rows)),
        "max_abs_cum_change_mag": float(max(
            max(abs(r["b_dm_cum_w1"] or 0), abs(r["b_dm_cum_w2"] or 0)) for r in star_rows))
        if star_rows else None,
        "per_star": star_rows,
    }
    return out


# --- stage 4: injection grid ---------------------------------------------------------
def injection_hosts(fields: dict) -> list[dict]:
    hosts = []
    for name, field in fields.items():
        if name != "ngc3786" and field["stars"]:
            hosts.append({"name": name, "kind": "control", "visits": corrected_target(field)})
        for lc, v in zip(field["stars"], corrected_stars(field)):
            hosts.append({"name": lc.name, "kind": "star", "visits": v})
    return hosts


METHODS = ("a", "b", "a_or_b", "ml_matched", "ml_far01", "ml_far05")


def _rate(rows: list[dict]) -> dict:
    n = len(rows)
    out = {"n": n}
    for m in METHODS:
        hits = sum((r["a"] or r["b"]) if m == "a_or_b" else r[m] for r in rows)
        out[m] = int(hits)
    return out


def injection_analysis(hosts: list[dict]) -> dict:
    rows = mlcompare.evaluate(injection.run_grid(hosts))
    grid = config.INJECTION
    recovery, false_alarm = [], []
    for kind in ("control", "star"):
        sub = [r for r in rows if r["kind"] == kind]
        for n_base in (*grid.baseline_visits, "all"):
            sel = [r for r in sub if n_base == "all" or r["baseline_visits"] == n_base]
            null = [r for r in sel if r["amp"] == 0.0]
            false_alarm.append({"kind": kind, "baseline_visits": n_base, **_rate(null),
                                "n_hosts": len({r["host"] for r in null}),
                                "hosts_flagged_a_or_b": sorted(
                                    {r["host"] for r in null if r["a"] or r["b"]})})
            for rise in grid.rise_years:
                for amp in grid.amplitudes_mag:
                    inj = [r for r in sel if r["amp"] == amp and r["rise"] == rise]
                    recovery.append({"kind": kind, "baseline_visits": n_base, "rise_years": rise,
                                     "amp_mag": amp, **_rate(inj)})
    return {
        "grid": {"amplitudes_mag": list(grid.amplitudes_mag), "rise_years": list(grid.rise_years),
                 "baseline_visits": list(grid.baseline_visits), "post_visits": grid.post_visits},
        "n_hosts": {k: len({r["host"] for r in rows if r["kind"] == k}) for k in ("control", "star")},
        "n_trials": len(rows),
        "ml": {"model": "IsolationForest", "n_trees": mlcompare.N_TREES,
               "n_folds": mlcompare.N_FOLDS, "features": list(injection.FEATURE_NAMES),
               "seed": config.SEED},
        "false_alarm": false_alarm,
        # The control galaxies are not guaranteed to be constant (NGC 2273 is
        # itself flagged by detector (b)), so the null rate is also listed per
        # control host.
        "false_alarm_by_control_host": _false_alarm_by_host(rows),
        "recovery": recovery,
    }


def _false_alarm_by_host(rows: list[dict]) -> dict:
    null = [r for r in rows if r["kind"] == "control" and r["amp"] == 0.0]
    hosts = sorted({r["host"] for r in null})
    per_host = {h: _rate([r for r in null if r["host"] == h]) for h in hosts}
    return {"per_host": per_host}


def dimming_check(hosts: list[dict], amp_mag: float = 0.3) -> dict:
    """Specificity: an instantaneous *dimming* is not a flare. How often does
    each method flag it? (The isolation forest is trained as in the main test.)"""
    grid = config.INJECTION
    rows = []
    for host in hosts:
        v = host["visits"]
        for n_base in grid.baseline_visits:
            for onset in injection.segment_onsets(v["mjd"].size, n_base, grid.post_visits):
                seg = injection.cut_segment(v, onset, n_base, grid.post_visits)
                common = {"host": host["name"], "kind": host["kind"], "baseline_visits": n_base,
                          "onset": onset}
                rows.append({**common, "amp": 0.0, "rise": -1.0,
                             **injection.run_detectors(seg, None),
                             "features": injection.features(seg)})
                dim = injection.inject(seg, n_base, -amp_mag, 0.0)
                rows.append({**common, "amp": -amp_mag, "rise": 0.0,
                             **injection.run_detectors(dim, n_base),
                             "features": injection.features(dim)})
    rows = mlcompare.evaluate(rows)
    out = {"amp_mag": -amp_mag}
    for kind in ("control", "star"):
        out[kind] = _rate([r for r in rows if r["kind"] == kind and r["amp"] < 0])
    return out


def to_jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_jsonable(v) for v in obj]
    if isinstance(obj, (np.floating, float)):
        return None if not np.isfinite(obj) else float(f"{float(obj):.6g}")
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, np.ndarray):
        return to_jsonable(obj.tolist())
    return obj


def write_json(name: str, data: dict) -> None:
    path = config.RESULTS_DIR / name
    path.write_text(json.dumps(to_jsonable(data), indent=1) + "\n")
