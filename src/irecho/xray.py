"""Stage 7: X-ray analysis of CGCG 312-012 (docs/XRAY_PLAN.md).

Constants of the data set, the decision rules of the plan (section 6) and
small numerical helpers. Nothing here needs HEASoft; the spectral fits are
run by scripts/xray_fit.py inside the HEASoft environment and write JSON
that the functions below turn into labels.
"""

from __future__ import annotations

import csv
import math

import numpy as np

from . import config

S7_DIR = config.CACHE_DIR / "s7"

TARGET = {"name": "bat0463", "counterpart": "CGCG 312-012", "bat_name": "SWIFT J0929.7+6232",
          "ra_deg": 142.4078319, "dec_deg": 62.5440626, "z": 0.0256}
OBS = {"nustar22": "60160365002", "xrt10": "00040946001", "xrt22": "00081051001"}

_HEASARC = "https://heasarc.gsfc.nasa.gov/FTP"
ARCHIVE_DIRS = (
    ("nustar/60160365002", f"{_HEASARC}/nustar/data/obs/01/6/60160365002/", ("auxil", "hk", "event_uf", "event_cl")),
    ("swift/00040946001", f"{_HEASARC}/swift/data/obs/2010_12/00040946001/", ("auxil", "xrt")),
    ("swift/00081051001", f"{_HEASARC}/swift/data/obs/2022_05/00081051001/", ("auxil", "xrt")),
)
_BAT157 = "https://swift.gsfc.nasa.gov/bat_survey/bs157mon"
BAT157_FILES = {
    "SWIFT_J0929.7+6232.pha": f"{_BAT157}/spectra/SWIFT_J0929.7+6232.pha",
    "BAT_157m_swiftbat_survey_full_157m.rsp": f"{_BAT157}/data/BAT_157m_swiftbat_survey_full_157m.rsp",
    "BAT_157m_catalog_20250525.txt": f"{_BAT157}/data/BAT_157m_catalog_20250525.txt",
}

# --- numbers fixed in docs/XRAY_PLAN.md before the data were opened -----------------
DELTA_L_MIN = 0.18          # dex, factor 1.5 (section 6.1)
DELTA_L_SAME_BOUND = 0.30   # dex, the 90% interval must lie inside for 'same'
DELTA_NH_MIN = 0.10         # dex (section 6.2)
DELTA_NH_SAME_BOUND = 0.20  # dex
SYS_NUSTAR_BAT_DEX = 0.06   # cross-calibration, added in quadrature to Delta1
SYS_XRT_NUSTAR_DEX = 0.04   # cross-calibration, added in quadrature to Delta2
AIC_STEP = 4.0              # a more complex model is adopted only if AIC drops by this much
MIN_XRT_COUNTS = 25         # below this no N_H is measured from an XRT observation
MAX_NH_INTERVAL_DEX = 1.0   # wider 90% interval: 'not constrained'
LOG_L_DUST_PEAK = 43.11     # stage 6A, results/s6a_flares.json


def _bass_log_lbol(name: str = "bat0463") -> float | None:
    """BASS DR2 log L_bol of the target, read from targets/sample_v1_props.csv.

    None when the column is blank (the public copy of the repository leaves the
    catalogue values out; scripts/restore_catalog_columns.py fills them in).
    """
    path = config.TARGETS_DIR / "sample_v1_props.csv"
    if not path.exists():
        return None
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            if row["name"] == name:
                return float(row["dr2_loglbol"]) if row["dr2_loglbol"] else None
    return None


LOG_LBOL_BASS = _bass_log_lbol()   # BASS DR2 log L_bol (erg/s); the value is not written in the code
KAPPA_2_10 = 20.0           # assumed 2-10 keV bolometric correction (cross-check only)
F_NEED_OK, F_NEED_TIGHT = 0.6, 1.0

HIGH, LOW, SAME, UNDECIDED, CHANGED = "높다", "낮다", "같다", "판정불가", "달라졌다"


def add_systematic(value: float, lo: float, hi: float, sys_dex: float) -> tuple[float, float]:
    """90% interval [lo, hi] around value widened by a systematic term in quadrature.

    The systematic is given as a 1-sigma size in dex; it is scaled to 90% (1.645).
    """
    s = 1.645 * sys_dex
    return (value - math.hypot(value - lo, s), value + math.hypot(hi - value, s))


def luminosity_label(delta: float, lo: float, hi: float) -> str:
    """Section 6.1. delta = log10 ratio, [lo, hi] its 90% interval (systematics included)."""
    # An interval may be open on one side (lo = -inf or hi = +inf) when the fit
    # statistic does not rise enough in that direction; the rules still apply.
    if math.isnan(delta) or math.isnan(lo) or math.isnan(hi):
        return UNDECIDED
    if lo > 0 and delta >= DELTA_L_MIN:
        return HIGH
    if hi < 0 and delta <= -DELTA_L_MIN:
        return LOW
    if abs(delta) < DELTA_L_MIN and lo > -DELTA_L_SAME_BOUND and hi < DELTA_L_SAME_BOUND:
        return SAME
    return UNDECIDED


def absorption_label(delta: float, lo: float, hi: float, *, measurable: bool = True) -> str:
    """Section 6.2. delta = log10 [N_H(2022) / N_H(2010)]."""
    if not measurable or math.isnan(delta) or math.isnan(lo) or math.isnan(hi):
        return UNDECIDED
    if (lo > 0 or hi < 0) and abs(delta) >= DELTA_NH_MIN:
        return CHANGED
    if abs(delta) < DELTA_NH_MIN and lo > -DELTA_NH_SAME_BOUND and hi < DELTA_NH_SAME_BOUND:
        return SAME
    return UNDECIDED


def overall_verdict(c1: str, c2: str, nh: str) -> dict:
    """Section 6.3. c1, c2: luminosity labels against the BAT average and 2010; nh: absorption label."""
    bright = (HIGH in (c1, c2)) and (LOW not in (c1, c2))
    if bright and nh == SAME:
        verdict, note = "H1 지지", "AGN 자체의 증광 후보"
    elif bright and nh == UNDECIDED:
        verdict, note = "H1 지지, 가림 변화는 배제 못 함", "AGN 자체의 증광 후보"
    elif bright and nh == CHANGED:
        verdict, note = "증광과 가림 변화 둘 다", ""
    elif c2 == SAME and c1 != HIGH and nh == CHANGED:
        verdict, note = "H2 지지", "가림 변화"
    else:
        verdict = UNDECIDED
        if c1 == SAME and c2 == SAME and nh == SAME:
            note = "2022.33의 X선에서 변화를 찾지 못함 (H3과 H4를 가르지 못함)"
        elif LOW in (c1, c2):
            note = "광도가 낮다는 비교가 있음"
        else:
            note = "정밀도 부족"
    return {"verdict": verdict, "note": note, "c1": c1, "c2": c2, "absorption": nh}


def energy_budget(delta1: float, log_lint_2022: float | None = None) -> dict:
    """Section 6.4: covering factor needed by the dust peak with the 2022 luminosity.

    The registered quantities start from the BASS bolometric luminosity. Without
    it (LOG_LBOL_BASS is None) they are returned as None and only the cross-check
    from the fitted 2-10 keV luminosity is computed.
    """
    if LOG_LBOL_BASS is None:
        out = {"log_lbol_2022": None, "f_need": None, "label": None, "f_need_with_bass_lbol": None,
               "covering_factor_from_excess_only": None}
        if log_lint_2022 is not None:
            log_lbol_k = log_lint_2022 + math.log10(KAPPA_2_10)
            out["log_lbol_2022_from_kappa20"] = log_lbol_k
            out["f_need_from_kappa20"] = 10 ** (LOG_L_DUST_PEAK - log_lbol_k)
        return out
    log_lbol = LOG_LBOL_BASS + delta1
    f_need = 10 ** (LOG_L_DUST_PEAK - log_lbol)
    label = "맞음" if f_need <= F_NEED_OK else ("빠듯함" if f_need <= F_NEED_TIGHT else "안 맞음")
    out = {"log_lbol_2022": log_lbol, "f_need": f_need, "label": label,
           "f_need_with_bass_lbol": 10 ** (LOG_L_DUST_PEAK - LOG_LBOL_BASS)}
    excess = 10 ** log_lbol - 10 ** LOG_LBOL_BASS
    out["covering_factor_from_excess_only"] = (10 ** LOG_L_DUST_PEAK / excess) if excess > 0 else None
    if log_lint_2022 is not None:
        log_lbol_k = log_lint_2022 + math.log10(KAPPA_2_10)
        out["log_lbol_2022_from_kappa20"] = log_lbol_k
        out["f_need_from_kappa20"] = 10 ** (LOG_L_DUST_PEAK - log_lbol_k)
    return out


def choose_model(aic: dict[str, float], order=("M0", "M1", "M2")) -> str:
    """Section 4.1: step to the next model only if AIC falls by at least AIC_STEP."""
    chosen = order[0]
    for nxt in order[1:]:
        if aic[chosen] - aic[nxt] >= AIC_STEP:
            chosen = nxt
        else:
            break
    return chosen


def powerlaw_band_flux(norm: float, gamma: float, e1: float, e2: float) -> float:
    """Energy flux (erg cm^-2 s^-1) of N(E) = norm E^-gamma between e1 and e2 keV."""
    kev = 1.602176634e-9
    if abs(gamma - 2.0) < 1e-9:
        return norm * math.log(e2 / e1) * kev
    return norm * (e2 ** (2 - gamma) - e1 ** (2 - gamma)) / (2 - gamma) * kev


def zpowerlw_rest_luminosity(norm: float, gamma: float, z: float, d_l_mpc: float,
                             e1: float = 2.0, e2: float = 10.0) -> float:
    """Rest-frame e1-e2 keV luminosity (erg/s) of an XSPEC zpowerlw component.

    zpowerlw: N(E) = norm [E (1+z)]^-gamma photons/keV/cm^2/s in the observed
    frame. The rest-frame band e1-e2 is observed at e/(1+z).
    """
    k_obs = norm * (1 + z) ** (-gamma)
    flux = powerlaw_band_flux(k_obs, gamma, e1 / (1 + z), e2 / (1 + z))
    d_cm = d_l_mpc * 3.0856775814913673e24
    return 4 * math.pi * d_cm ** 2 * flux


def net_rate(src_counts: float, bkg_counts: float, area_ratio: float, exposure_s: float) -> tuple[float, float]:
    """Net count rate and its Gaussian error; area_ratio = source area / background area."""
    net = src_counts - bkg_counts * area_ratio
    err = math.sqrt(src_counts + bkg_counts * area_ratio ** 2)
    return net / exposure_s, err / exposure_s


def hardness_ratio(hard: float, soft: float, hard_err: float, soft_err: float) -> tuple[float, float]:
    """(H - S) / (H + S) with first-order error propagation."""
    tot = hard + soft
    hr = (hard - soft) / tot
    err = 2.0 / tot ** 2 * math.sqrt((soft * hard_err) ** 2 + (hard * soft_err) ** 2)
    return hr, err


def interval_coverage(truth: float, lo, hi) -> float:
    lo, hi = np.asarray(lo, dtype=float), np.asarray(hi, dtype=float)
    ok = np.isfinite(lo) & np.isfinite(hi)
    return float(np.mean((lo[ok] <= truth) & (truth <= hi[ok]))) if ok.any() else float("nan")


# --- stage 7b (post hoc supplement, docs/XRAY_PLAN.md section 9) ----------------------------------

def bat_channels_in_band(e_min, e_max, lo: float = 14.0, hi: float = 50.0) -> list[int]:
    """1-based numbers of the channels that lie entirely inside [lo, hi] keV.

    Channels are chosen by number, not by an energy 'ignore': XSPEC also drops
    the channel that contains the limit itself (stage 7 lost channel 1 that way).
    """
    return [i + 1 for i, (a, b) in enumerate(zip(e_min, e_max)) if a >= lo - 1e-6 and b <= hi + 1e-6]


def binomial_se(p: float, n: int) -> float:
    return math.sqrt(p * (1.0 - p) / n)


def autocorr_time(chains, c: float = 5.0) -> float:
    """Integrated autocorrelation time (steps) of an ensemble chain, array [n_steps, n_walkers].

    Mean of the per-walker autocorrelation functions, summed with Sokal's
    automatic window (first M with M >= c * tau).
    """
    x = np.asarray(chains, dtype=float)
    n = x.shape[0]
    acf = np.zeros(n)
    for w in range(x.shape[1]):
        y = x[:, w] - x[:, w].mean()
        f = np.fft.rfft(y, n=2 * n)
        a = np.fft.irfft(f * np.conjugate(f))[:n]
        acf += a / a[0] if a[0] > 0 else 0.0
    acf /= x.shape[1]
    taus = 2.0 * np.cumsum(acf) - 1.0
    window = np.arange(n) >= c * taus
    m = int(np.argmax(window)) if window.any() else n - 1
    return float(taus[m])


def gelman_rubin(chains) -> float:
    """Potential scale reduction factor with every walker split in two halves."""
    x = np.asarray(chains, dtype=float)
    half = x.shape[0] // 2
    parts = np.concatenate([x[:half], x[half:2 * half]], axis=1)      # [half, 2 * n_walkers]
    n = parts.shape[0]
    within = parts.var(axis=0, ddof=1).mean()
    between = n * parts.mean(axis=0).var(ddof=1)
    return float(np.sqrt(((n - 1) / n * within + between / n) / within))


def poisson_net_samples(n_src: int, n_bkg: int, area_ratio: float, size: int, rng) -> np.ndarray:
    """Posterior samples of the net source counts (>= 0) in a source region.

    Background: N_bkg counts in a region 1/area_ratio times larger, flat prior
    on its mean. Source: flat prior on s >= 0 with N_src ~ Poisson(s + b).
    """
    from scipy.stats import gamma

    b = rng.gamma(n_bkg + 1.0, 1.0, size) * area_ratio
    u = rng.uniform(gamma.cdf(b, n_src + 1.0), 1.0, size)
    total = gamma.ppf(np.clip(u, 0.0, 1.0 - 1e-15), n_src + 1.0)
    return np.maximum(total - b, 0.0)


def gti_overlap_seconds(gti_a, gti_b) -> float:
    """Total overlap of two lists of (start, stop) intervals given on the same clock."""
    total = 0.0
    for a0, a1 in gti_a:
        for b0, b1 in gti_b:
            total += max(0.0, min(a1, b1) - max(a0, b0))
    return float(total)
