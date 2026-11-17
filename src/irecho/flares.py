"""Stage 6A: properties of infrared flare candidates and the pre-registered
consistency table (docs/FLARES_PLAN.md, committed before these quantities were
computed). Everything refers to the variable component: visit flux minus the
median flux of the baseline visits.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

from . import physics
from .config import BANDS, LAMBDA_ISO_UM, SEED
from .epochs import MEDIAN_SE_FACTOR
from .photometry import C_LIGHT

DAYS_PER_YEAR = 365.25
N_DRAWS = 4000
JY_CGS = 1e-23  # erg s^-1 cm^-2 Hz^-1
MPC_CM = 3.0856775814913673e24
# FLARES_PLAN.md section 4
EVENT_MIN_MAG = 0.1
EVENT_MIN_SIGMA = 3.0
SN_MAX_NULNU_W2 = 1e42  # erg/s
SN_MAX_ENERGY = 1e51  # erg
TDE_MAX_ENERGY = 9e52  # erg
TDE_MAX_RISE_YR = 1.5
NUCLEAR_MAX_ARCSEC = 1.0
ARTEFACT_MAX_ARCSEC = 3.0
BOTH_BANDS_MIN_FRACTION = 0.8
UNTIMELY_MAX_DIFF_MAG = 0.1
REBRIGHTEN_MAG = 0.1


# --- events -----------------------------------------------------------------------------
def baseline(v: dict, n_base: int) -> tuple[dict, dict]:
    """Median flux of the first ``n_base`` visits and its standard error (as
    stage 3: 1.2533 x std / sqrt(n))."""
    val = {b: float(np.median(v[f"{b}_flux_jy"][:n_base])) for b in BANDS}
    err = {b: MEDIAN_SE_FACTOR * float(np.std(v[f"{b}_flux_jy"][:n_base], ddof=1)) / np.sqrt(n_base)
           for b in BANDS}
    return val, err


def excess(v: dict, base: dict, base_err: dict) -> dict:
    out = {}
    for b in BANDS:
        out[b] = v[f"{b}_flux_jy"] - base[b]
        out[f"{b}_err"] = np.hypot(v[f"{b}_flux_err_jy"], base_err[b])
        out[f"{b}_dmag"] = -2.5 * np.log10(v[f"{b}_flux_jy"] / base[b])  # negative = brighter
    return out


def episodes(ex: dict, n_base: int) -> list[tuple[int, int]]:
    """Runs of visits (index ranges, inclusive) after the baseline whose W2
    excess is >= 3 sigma and >= 0.1 mag."""
    on = (ex["w2"] / ex["w2_err"] >= EVENT_MIN_SIGMA) & (-ex["w2_dmag"] >= EVENT_MIN_MAG)
    on[:n_base] = False
    runs, start = [], None
    for k, flag in enumerate(on):
        if flag and start is None:
            start = k
        if not flag and start is not None:
            runs.append((start, k - 1))
            start = None
    if start is not None:
        runs.append((start, on.size - 1))
    return runs


def main_episode(ex: dict, runs: list[tuple[int, int]]) -> tuple[int, int, int] | None:
    """(start, end, peak index) of the run with the largest W2 excess."""
    if not runs:
        return None
    best = max(runs, key=lambda r: ex["w2"][r[0]: r[1] + 1].max())
    peak = best[0] + int(np.argmax(ex["w2"][best[0]: best[1] + 1]))
    return best[0], best[1], peak


# --- time scales and decay shape ------------------------------------------------------------
def _fit(model, p0, t, y, err, bounds):
    sol = least_squares(lambda p: (model(t, *p) - y) / err, p0, bounds=bounds)
    return sol.x, float(2.0 * sol.cost)


def decay_fits(t_yr, y, err) -> dict:
    """Exponential, (t - t0)^(-5/3) and free power law fitted to the excess
    after the peak (``t_yr`` measured from the peak, first point = peak).
    Descriptive only (FLARES_PLAN.md section 3)."""
    t, y, err = (np.asarray(x, dtype=float) for x in (t_yr, y, err))
    out = {"n_points": int(t.size)}
    if t.size < 4:
        return out
    a0 = float(y[0])
    p, chi2 = _fit(lambda t, a, tau: a * np.exp(-t / tau), [a0, 2.0], t, y, err,
                   ([0.0, 0.05], [np.inf, 100.0]))
    out["exponential"] = {"tau_yr": float(p[1]), "chi2": chi2, "dof": int(t.size - 2)}
    # power laws: t0 between 3 yr before the peak and just before it
    p, chi2 = _fit(lambda t, a, t0: a * ((t + t0) / t0) ** (-5.0 / 3.0), [a0, 1.0], t, y, err,
                   ([0.0, 0.02], [np.inf, 3.0]))
    out["power_law_5_3"] = {"t0_before_peak_yr": float(p[1]), "chi2": chi2, "dof": int(t.size - 2)}
    if t.size >= 5:
        p, chi2 = _fit(lambda t, a, t0, n: a * ((t + t0) / t0) ** (-n), [a0, 1.0, 1.5], t, y, err,
                       ([0.0, 0.02, 0.1], [np.inf, 3.0, 10.0]))
        out["power_law_free"] = {"t0_before_peak_yr": float(p[1]), "index": float(p[2]),
                                 "chi2": chi2, "dof": int(t.size - 3)}
    return out


def time_scales(v: dict, ex: dict, n_base: int, ep: tuple[int, int, int]) -> dict:
    start, end, peak = ep
    yr = v["mjd"] / DAYS_PER_YEAR
    quiet = np.flatnonzero(ex["w2"][:peak] / ex["w2_err"][:peak] < EVENT_MIN_SIGMA)
    t_rise = float(yr[peak] - yr[quiet[-1]]) if quiet.size else float("nan")
    post = np.arange(peak, v["mjd"].size)
    fits = decay_fits(yr[post] - yr[peak], ex["w2"][post], ex["w2_err"][post])
    t_decay = fits.get("exponential", {}).get("tau_yr", float("nan"))
    # re-brightening after the peak: a later visit brighter than an earlier post-peak visit
    # by more than 0.1 mag and 3 sigma
    rebright = False
    for i in range(peak + 1, v["mjd"].size):
        for j in range(i + 1, v["mjd"].size):
            dm = v["w2_mag"][i] - v["w2_mag"][j]
            if dm >= REBRIGHTEN_MAG and dm / np.hypot(v["w2_err"][i], v["w2_err"][j]) >= 3.0:
                rebright = True
    return {"t_rise_yr": t_rise, "t_decay_efold_yr": float(t_decay), "decay_fits": fits,
            "rebrightening_after_peak": bool(rebright),
            "n_visits_after_peak": int(post.size - 1), "ended": bool(end < v["mjd"].size - 1)}


# --- physical quantities ------------------------------------------------------------------------
def nu_l_nu(excess_jy: float, band: str, d_l_mpc: float) -> float:
    nu = C_LIGHT / (LAMBDA_ISO_UM[band] * 1e-6)
    return float(4.0 * np.pi * (d_l_mpc * MPC_CM) ** 2 * nu * excess_jy * JY_CGS)


def visit_physics(v: dict, n_base: int, z: float, seed_tag: int) -> dict:
    """T, L, R_bb of the variable component per visit, the energy integral and
    their Monte Carlo intervals. One set of baseline draws is shared by all
    visits."""
    rng = np.random.default_rng([SEED, 61, seed_tag])
    base, base_err = baseline(v, n_base)
    d_l = physics.luminosity_distance_mpc(z)
    base_draws = {b: rng.normal(base[b], base_err[b], N_DRAWS) for b in BANDS}
    n = v["mjd"].size
    rows = []
    l_draws = np.full((n, N_DRAWS), np.nan)
    for k in range(n_base, n):
        flux = {b: float(v[f"{b}_flux_jy"][k]) for b in BANDS}
        ferr = {b: float(v[f"{b}_flux_err_jy"][k]) for b in BANDS}
        f = {b: rng.normal(flux[b], ferr[b], N_DRAWS) for b in BANDS}
        keys = ("t_rest_k", "l_dust_erg_s", "r_bb_light_days")
        draws = {key: np.full(N_DRAWS, np.nan) for key in keys}
        for i in range(N_DRAWS):
            q = physics.blackbody_quantities(f["w1"][i] - base_draws["w1"][i],
                                             f["w2"][i] - base_draws["w2"][i], z, d_l)
            for key in keys:
                draws[key][i] = q[key]
        l_draws[k] = draws["l_dust_erg_s"]
        best = physics.blackbody_quantities(flux["w1"] - base["w1"], flux["w2"] - base["w2"], z, d_l)
        row = {"index": k, "mjd": float(v["mjd"][k]),
               "valid_draw_fraction": float(np.isfinite(draws["t_rest_k"]).mean())}
        for key in keys:
            good = draws[key][np.isfinite(draws[key])]
            lo, med, hi = np.percentile(good, [16, 50, 84]) if good.size else (np.nan,) * 3
            row.update({key: best[key], f"{key}_p16": float(lo), f"{key}_p50": float(med),
                        f"{key}_p84": float(hi)})
        row["nu_l_nu_w2_erg_s"] = nu_l_nu(flux["w2"] - base["w2"], "w2", d_l)
        row["nu_l_nu_w1_erg_s"] = nu_l_nu(flux["w1"] - base["w1"], "w1", d_l)
        rows.append(row)
    return {"baseline_flux_jy": base, "baseline_flux_err_jy": base_err, "d_l_mpc": float(d_l),
            "visits": rows, "_l_draws": l_draws}


def energy(v: dict, l_draws: np.ndarray, start: int, end: int, z: float) -> dict:
    """Trapezoid integral of L_dust over the visits of one episode, rest frame.
    Draws with an undefined luminosity at a visit count as zero there."""
    idx = np.arange(start, end + 1)
    t = v["mjd"][idx] * 86400.0 / (1.0 + z)
    lum = np.nan_to_num(l_draws[idx], nan=0.0)
    if idx.size < 2:
        return {"e_dust_erg": float("nan"), "e_dust_erg_p16": float("nan"),
                "e_dust_erg_p84": float("nan"), "n_visits": int(idx.size)}
    e = np.trapezoid(lum, t, axis=0)
    lo, med, hi = np.percentile(e, [16, 50, 84])
    return {"e_dust_erg": float(med), "e_dust_erg_p16": float(lo), "e_dust_erg_p84": float(hi),
            "n_visits": int(idx.size)}


# --- is the signal real (A1) ------------------------------------------------------------------------
def exposure_checks(lc, v: dict, base_idx, event_idx, base_mag: dict) -> dict:
    """Single-exposure level: both bands brighter in the same exposures, the
    position of the source, and the fraction of exposures that pass the cuts."""
    obs, eid = lc.obs, lc.epoch_id
    visit_mjd = np.array([obs["mjd"][eid == k].mean() for k in range(int(eid.max()) + 1)])

    def rows_of(indices):
        sel = np.zeros(obs["mjd"].size, dtype=bool)
        for k in indices:
            j = int(np.argmin(np.abs(visit_mjd - v["mjd"][k])))
            if abs(visit_mjd[j] - v["mjd"][k]) < 15.0:
                sel |= eid == j
        return sel

    good = lc.masks["w1"] & lc.masks["w2"]
    ev, ba = rows_of(event_idx), rows_of(base_idx)
    both = ((obs["w1mag"] < base_mag["w1"]) & (obs["w2mag"] < base_mag["w2"]))[ev & good]
    cosd = np.cos(np.deg2rad(lc.dec))

    def centroid(sel):
        ra, dec = obs["ra"][sel & good], obs["dec"][sel & good]
        n = ra.size
        return (float(np.median(ra)), float(np.median(dec)),
                MEDIAN_SE_FACTOR * float(np.std(ra, ddof=1)) * cosd * 3600.0 / np.sqrt(n),
                MEDIAN_SE_FACTOR * float(np.std(dec, ddof=1)) * 3600.0 / np.sqrt(n), n)

    ra_e, dec_e, sra_e, sdec_e, n_e = centroid(ev)
    ra_b, dec_b, sra_b, sdec_b, n_b = centroid(ba)
    d_ra = (ra_e - ra_b) * cosd * 3600.0
    d_dec = (dec_e - dec_b) * 3600.0
    sep = float(np.hypot(d_ra, d_dec))
    sep_err = float(np.hypot(d_ra * np.hypot(sra_e, sra_b), d_dec * np.hypot(sdec_e, sdec_b))
                    / sep) if sep > 0 else float("nan")
    return {
        "n_event_exposures": int((ev & good).sum()), "n_baseline_exposures": int((ba & good).sum()),
        "fraction_event_exposures_brighter_in_both_bands": float(both.mean()) if both.size else None,
        "position_shift_arcsec": sep, "position_shift_err_arcsec": sep_err,
        "position_shift_ra_arcsec": float(d_ra), "position_shift_dec_arcsec": float(d_dec),
        "kept_fraction_event": float((ev & good).sum() / max(ev.sum(), 1)),
        "kept_fraction_baseline": float((ba & good).sum() / max(ba.sum(), 1)),
    }


def untimely_check(v: dict, u: dict, n_base: int, event_idx, pairs) -> dict:
    """Brightening relative to the baseline visits in the pipeline photometry
    and in unTimely, for the event visits that unTimely covers."""
    pairs = [(i, j) for i, j in pairs if np.isfinite(u["w1_mag"][j]) and np.isfinite(u["w2_mag"][j])]
    base = [(i, j) for i, j in pairs if i < n_base]
    event = [(i, j) for i, j in pairs if i in set(int(k) for k in event_idx)]
    out = {"n_baseline_pairs": len(base), "n_event_pairs": len(event)}
    if len(base) < 2 or not event:
        out["status"] = "not covered"
        return out
    diffs = {}
    for b in BANDS:
        dm_pipe = np.median([v[f"{b}_mag"][i] for i, _ in base]) - np.array(
            [v[f"{b}_mag"][i] for i, _ in event])
        dm_unt = np.median([u[f"{b}_mag"][j] for _, j in base]) - np.array(
            [u[f"{b}_mag"][j] for _, j in event])
        diffs[b] = dm_pipe - dm_unt
        out[f"{b}_brightening_pipeline_mag"] = [round(float(x), 3) for x in dm_pipe]
        out[f"{b}_brightening_untimely_mag"] = [round(float(x), 3) for x in dm_unt]
    out["max_abs_difference_mag"] = float(max(np.max(np.abs(d)) for d in diffs.values()))
    out["median_abs_difference_mag"] = float(np.median(np.abs(np.concatenate(list(diffs.values())))))
    out["status"] = "covered"
    return out


# --- consistency table (FLARES_PLAN.md section 4.2) ------------------------------------------------------
def classify(m: dict) -> dict:
    """``m``: measured values of one object. Returns the verdict per
    interpretation with the reasons."""
    out = {}
    # artefact
    fails = [name for name, ok in m["o1_checks"].items() if ok is False]
    unknown = [name for name, ok in m["o1_checks"].items() if ok is None]
    out["artefact"] = {"verdict": "부합" if fails else "불일치", "failed_checks": fails,
                       "not_checkable": unknown}
    # supernova
    reasons = []
    if m["nu_l_nu_w2_peak"] >= SN_MAX_NULNU_W2:
        reasons.append("O2: peak nuLnu(4.6um) >= 1e42 erg/s")
    if m["e_dust"] > SN_MAX_ENERGY:
        reasons.append("O3: E_dust > 1e51 erg")
    if m["n_episodes"] >= 2:
        reasons.append("O5: more than one episode")
    out["supernova"] = {"verdict": "불일치" if reasons else "부합", "reasons": reasons}
    # obscuration change
    diff, err = m["dm_w2_minus_dm_w1"], m["dm_diff_err"]
    if diff > 2.0 * err:
        out["obscuration_change"] = {"verdict": "불일치",
                                     "reasons": ["O7: W2 brightened more than W1 by > 2 sigma"]}
    else:
        out["obscuration_change"] = {"verdict": "부합", "reasons": []}
    # TDE inside the AGN
    nuclear = (m["position_shift_arcsec"] < NUCLEAR_MAX_ARCSEC
               and m["position_shift_arcsec"] < 3.0 * m["position_shift_err_arcsec"])
    off_nucleus = (m["position_shift_arcsec"] >= NUCLEAR_MAX_ARCSEC
                   and m["position_shift_arcsec"] >= 3.0 * m["position_shift_err_arcsec"])
    reasons, missing = [], []
    if off_nucleus:
        reasons.append("O4: not at the nucleus")
    if m["n_episodes"] >= 2:
        reasons.append("O5: more than one episode")
    if m["rebrightening_after_peak"]:
        reasons.append("O5: re-brightening after the peak")
    if np.isfinite(m["t_rise_yr"]) and m["t_rise_yr"] > TDE_MAX_RISE_YR:
        reasons.append("O5: rise longer than 1.5 yr")
    if m["e_dust"] > TDE_MAX_ENERGY:
        reasons.append("O3: E_dust > 9e52 erg")
    if not np.isfinite(m["t_rise_yr"]):
        missing.append("t_rise")
    if not np.isfinite(m["t_decay_efold_yr"]):
        missing.append("t_decay")
    if reasons:
        verdict = "불일치"
    elif (nuclear and not missing and m["t_decay_efold_yr"] > m["t_rise_yr"]):
        verdict = "부합"
    else:
        verdict = "판정불가"
    out["tde_in_agn"] = {"verdict": verdict, "reasons": reasons, "not_measured": missing}
    # AGN flare / changing look
    out["agn_flare"] = {"verdict": "불일치" if off_nucleus else ("부합" if nuclear else "판정불가"),
                        "reasons": ["O4: not at the nucleus"] if off_nucleus else []}
    return out
