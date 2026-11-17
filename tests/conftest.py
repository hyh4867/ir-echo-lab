"""Shared synthetic data for the offline unit tests (no network, no cache)."""

from __future__ import annotations

import numpy as np
import pytest

RA0, DEC0 = 150.0, 20.0


def make_neowise_table(visit_mjds, n_frames=12, w1=10.3, w2=10.0, sigma=0.02,
                       seed=1, dmag=None):
    """Synthetic neowiser_p1bs_psd-like table: ``n_frames`` exposures per visit
    spread over ~1.5 days, Gaussian noise ``sigma`` mag. ``dmag`` (one value per
    visit, positive = fainter) shifts both bands."""
    rng = np.random.default_rng(seed)
    visit_mjds = np.asarray(visit_mjds, dtype=float)
    dmag = np.zeros(len(visit_mjds)) if dmag is None else np.asarray(dmag, dtype=float)
    mjd = np.concatenate([t + np.linspace(0, 1.5, n_frames) for t in visit_mjds])
    shift = np.repeat(dmag, n_frames)
    n = mjd.size
    return {
        "ra": RA0 + rng.normal(0, 3e-5, n),
        "dec": DEC0 + rng.normal(0, 3e-5, n),
        "mjd": mjd,
        "scan_id": np.array([f"{10000 + i:05d}a" for i in range(n)]),
        "frame_num": np.full(n, 100.0),
        "source_id": np.array([f"{10000 + i:05d}a100-000001" for i in range(n)]),
        "w1mpro": w1 + shift + rng.normal(0, sigma, n),
        "w1sigmpro": np.full(n, sigma),
        "w1snr": np.full(n, 50.0),
        "w1rchi2": np.full(n, 1.0),
        "w2mpro": w2 + shift + rng.normal(0, sigma, n),
        "w2sigmpro": np.full(n, sigma),
        "w2snr": np.full(n, 50.0),
        "w2rchi2": np.full(n, 1.0),
        "w1sat": np.zeros(n),
        "w2sat": np.zeros(n),
        "nb": np.ones(n),
        "na": np.zeros(n),
        "cc_flags": np.full(n, "0000"),
        "det_bit": np.full(n, 3.0),
        "ph_qual": np.full(n, "AA"),
        "sso_flg": np.zeros(n),
        "qual_frame": np.full(n, 10.0),
        "qi_fact": np.ones(n),
        "saa_sep": np.full(n, 60.0),
        "moon_masked": np.full(n, "00"),
        "w1frtr": np.full(n, "0000"),
        "w2frtr": np.full(n, "0000"),
        "allwise_cntr": np.full(n, "1"),
    }


def make_allwise_table(visit_mjds, n_frames=12, w1=10.3, w2=10.0, sigma=0.02, seed=2):
    """Synthetic allwise_p3as_mep-like table."""
    rng = np.random.default_rng(seed)
    mjd = np.concatenate([t + np.linspace(0, 1.5, n_frames) for t in visit_mjds])
    n = mjd.size
    return {
        "cntr_mf": np.full(n, "1"),
        "source_id_mf": np.full(n, "1500p200_ac51-000001"),
        "ra": np.full(n, RA0),
        "dec": np.full(n, DEC0),
        "nb": np.ones(n),
        "na": np.zeros(n),
        "cc_flags": np.full(n, "0000"),
        "frame_id": np.array([f"{i:05d}a100" for i in range(n)]),
        "mjd": mjd,
        "w1mpro_ep": w1 + rng.normal(0, sigma, n),
        "w1sigmpro_ep": np.full(n, sigma),
        "w1rchi2_ep": np.ones(n),
        "w2mpro_ep": w2 + rng.normal(0, sigma, n),
        "w2sigmpro_ep": np.full(n, sigma),
        "w2rchi2_ep": np.ones(n),
        "qi_fact": np.ones(n),
        "saa_sep": np.full(n, 60.0),
        "moon_masked": np.full(n, "0000"),
    }


@pytest.fixture
def neowise_visits():
    """20 visits, ~182 days apart, starting 2014."""
    return 56700.0 + 182.0 * np.arange(20)


@pytest.fixture
def allwise_visits():
    return np.array([55340.0, 55530.0])
