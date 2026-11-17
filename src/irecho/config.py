"""Constants and fixed analysis choices.

Every number here is either quoted from a cited document or is an analysis
choice that was fixed before the known-answer tests were run. The rationale
and citations are in docs/PIPELINE.md.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CACHE_DIR = PROJECT_ROOT / "data" / "cache"
MANIFEST_PATH = PROJECT_ROOT / "data" / "manifest.json"
RESULTS_DIR = PROJECT_ROOT / "results"
TARGETS_DIR = PROJECT_ROOT / "targets"

# --- IRSA ------------------------------------------------------------------
TAP_SYNC_URL = "https://irsa.ipac.caltech.edu/TAP/sync"
TABLE_ALLWISE_MEP = "allwise_p3as_mep"  # AllWISE Multiepoch Photometry Table
TABLE_NEOWISE = "neowiser_p1bs_psd"  # NEOWISE-R Single Exposure (L1b) Source Table
TABLE_ALLWISE_CAT = "allwise_p3as_psd"  # AllWISE Source Catalog

# Cone radius for matching single-exposure detections to a target.
# NEOWISE Explanatory Supplement II.3.b suggests "a small radius (e.g. 3 arcsec)".
CONE_RADIUS_ARCSEC = 3.0

# --- Photometric system ----------------------------------------------------
# Zero-magnitude flux densities for a F_nu = const spectrum, in Jy.
# WISE All-Sky Explanatory Supplement IV.4.h, Table 1 (Jarrett et al. 2011).
F0_JY = {"w1": 309.540, "w2": 171.787}

# Isophotal wavelengths in micron (Jarrett et al. 2011; Wright et al. 2010).
LAMBDA_ISO_UM = {"w1": 3.3526, "w2": 4.6028}

# Onset of saturation in single exposures, Vega mag. NEOWISE Explanatory
# Supplement II.2 quotes W1~8.0 and W2~7.5; the AllWISE Explanatory Supplement
# III.2.b quotes W1<8 and W2<7. We keep a margin and warn below these values.
SATURATION_WARN_MAG = {"w1": 8.5, "w2": 8.0}

# Angular resolution (FWHM) in arcsec, Wright et al. 2010.
FWHM_ARCSEC = {"w1": 6.1, "w2": 6.4}

BANDS = ("w1", "w2")

# --- Fixed analysis choices -------------------------------------------------
SEED = 20261112


@dataclass(frozen=True)
class QualityCuts:
    """Frame- and detection-level cuts. See docs/PIPELINE.md section 2."""

    min_qual_frame: float = 0.0  # NEOWISE only: keep qual_frame > 0
    min_qi_fact_neowise: float = 0.0  # keep qi_fact > 0
    require_qi_fact_allwise: float = 1.0  # keep qi_fact == 1
    min_saa_sep_deg: float = 5.0  # keep saa_sep > 5
    require_moon_unmasked: bool = True  # band character of moon_masked == '0'
    require_cc_clean: bool = True  # band character of cc_flags == '0'
    require_unsaturated: bool = True  # NEOWISE only: w?sat == 0
    reject_sso: bool = True  # NEOWISE only: sso_flg == 0


@dataclass(frozen=True)
class EpochBinning:
    """Visit grouping and per-visit statistics. See docs/PIPELINE.md section 3."""

    gap_days: float = 30.0  # a gap longer than this starts a new visit
    clip_sigma: float = 5.0  # one-pass outlier rejection inside a visit
    min_frames: int = 5  # a band needs at least this many frames per visit


@dataclass(frozen=True)
class FlareCriteria:
    """Brightening-event detector. See docs/PIPELINE.md section 5."""

    min_baseline_epochs: int = 4
    n_sigma: float = 5.0  # required in both W1 and W2
    min_amp_mag: float = 0.2  # required in both W1 and W2


@dataclass(frozen=True)
class SlowRiseCriteria:
    """Slow-rise detector (b). See docs/PIPELINE.md section 8. Fixed before the
    stage-4 injection tests were run; NGC 2273 had already been seen in stage 2."""

    min_visits: int = 8
    min_baseline_visits: int = 4  # flat visits required before a hinge onset
    min_ramp_visits: int = 3  # rising visits required after a hinge onset
    p_max: float = 1e-3  # F-test p-value, required in both W1 and W2
    min_cum_mag: float = 0.1  # cumulative brightening, required in both W1 and W2
    tau_grid_yr: tuple[float, ...] = (0.5, 1.0, 2.0, 4.0)  # exponential e-folding times


@dataclass(frozen=True)
class InjectionGrid:
    """Stage-4 injection grid. See docs/PIPELINE.md section 9."""

    amplitudes_mag: tuple[float, ...] = (0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0)
    rise_years: tuple[float, ...] = (0.0, 1.0, 3.0)  # 0 = instantaneous step
    baseline_visits: tuple[int, ...] = (4, 8, 12)  # visits before the onset
    post_visits: int = 7  # visits kept from the onset on (about 3.5 yr)


QUALITY = QualityCuts()
BINNING = EpochBinning()
FLARE = FlareCriteria()
SLOW_RISE = SlowRiseCriteria()
INJECTION = InjectionGrid()
