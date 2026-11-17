"""Stage 6A: flare measures, the consistency table and other-band readers
(synthetic inputs; the last tests use the committed results when present)."""

import json

import numpy as np
import pytest

from irecho import config, flares, otherbands, physics
from irecho.photometry import mag_to_flux_jy


def visits(dm_w1, dm_w2, err=0.01):
    n = len(dm_w1)
    mjd = 56650.0 + 182.6 * np.arange(n)
    v = {"mjd": mjd}
    for b, base, dm in (("w1", 11.0, dm_w1), ("w2", 10.5, dm_w2)):
        mag = base - np.asarray(dm, dtype=float)
        flux = mag_to_flux_jy(mag, b)
        v[f"{b}_mag"], v[f"{b}_err"] = mag, np.full(n, err)
        v[f"{b}_flux_jy"] = flux
        v[f"{b}_flux_err_jy"] = flux * np.log(10) / 2.5 * err
    return v


def bump(n, start, peak, tau):
    k = np.arange(n, dtype=float)
    out = np.zeros(n)
    rise = (k >= start) & (k <= peak)
    out[rise] = (k[rise] - start + 1) / (peak - start + 1)
    out[k > peak] = np.exp(-(k[k > peak] - peak) / tau)
    return out


def test_episodes_and_main_episode():
    rng = np.random.default_rng(1)
    shape = 0.3 * bump(21, 6, 7, 1.0) + 0.8 * bump(21, 13, 15, 3.0)
    shape[10:12] = 0.0
    v = visits(0.6 * shape + rng.normal(0, 0.003, 21), shape + rng.normal(0, 0.003, 21))
    base, base_err = flares.baseline(v, 5)
    ex = flares.excess(v, base, base_err)
    runs = flares.episodes(ex, 5)
    assert len(runs) == 2 and runs[0][0] == 6
    start, end, peak = flares.main_episode(ex, runs)
    assert peak == 15 and start >= 12
    assert flares.main_episode(ex, []) is None


def test_decay_fit_recovers_an_exponential():
    t = np.arange(0.0, 5.0, 0.5)
    y = 2.0 * np.exp(-t / 1.5)
    fits = flares.decay_fits(t, y, np.full(t.size, 0.01))
    assert fits["exponential"]["tau_yr"] == pytest.approx(1.5, rel=1e-3)
    assert fits["exponential"]["chi2"] < fits["power_law_5_3"]["chi2"]
    assert "exponential" not in flares.decay_fits(t[:3], y[:3], np.full(3, 0.01))


def test_time_scales_and_rebrightening():
    shape = 0.8 * bump(21, 9, 10, 4.0)
    v = visits(0.6 * shape, shape)
    base, base_err = flares.baseline(v, 5)
    ex = flares.excess(v, base, base_err)
    ep = flares.main_episode(ex, flares.episodes(ex, 5))
    ts = flares.time_scales(v, ex, 5, ep)
    assert ts["t_rise_yr"] == pytest.approx(1.0, abs=0.01)  # two visits
    assert ts["t_decay_efold_yr"] > ts["t_rise_yr"] and not ts["rebrightening_after_peak"]
    shape[16] += 0.3
    v = visits(0.6 * shape, shape)
    ex = flares.excess(v, *flares.baseline(v, 5))
    assert flares.time_scales(v, ex, 5, ep)["rebrightening_after_peak"]


def test_nu_l_nu_and_energy():
    d = 100.0
    val = flares.nu_l_nu(1e-3, "w2", d)
    nu = 2.99792458e8 / 4.6028e-6
    assert val == pytest.approx(4 * np.pi * (d * 3.0856775814913673e24) ** 2 * nu * 1e-26, rel=1e-9)
    v = {"mjd": np.array([0.0, 100.0, 200.0])}
    draws = np.full((3, 10), 1e42)
    e = flares.energy(v, draws, 0, 2, z=0.0)
    assert e["e_dust_erg"] == pytest.approx(1e42 * 200 * 86400.0)
    assert flares.energy(v, draws, 0, 2, z=1.0)["e_dust_erg"] == pytest.approx(1e42 * 100 * 86400.0)


def test_visit_physics_returns_blackbody_of_the_excess():
    shape = bump(12, 6, 7, 3.0)
    v = visits(0.5 * shape, shape, err=0.002)
    out = flares.visit_physics(v, 5, 0.02, 0)
    row = next(r for r in out["visits"] if r["index"] == 7)
    base = out["baseline_flux_jy"]
    want = physics.blackbody_quantities(v["w1_flux_jy"][7] - base["w1"], v["w2_flux_jy"][7] - base["w2"],
                                        0.02, out["d_l_mpc"])
    assert row["t_rest_k"] == pytest.approx(want["t_rest_k"])
    assert row["t_rest_k_p16"] < row["t_rest_k"] < row["t_rest_k_p84"]


BASE = {"o1_checks": {"a": True, "b": True, "c": None}, "nu_l_nu_w2_peak": 5e41, "e_dust": 5e49,
        "n_episodes": 1, "dm_w2_minus_dm_w1": 0.0, "dm_diff_err": 0.05,
        "position_shift_arcsec": 0.03, "position_shift_err_arcsec": 0.02,
        "rebrightening_after_peak": False, "t_rise_yr": 0.5, "t_decay_efold_yr": 2.0}


def verdicts(**kw):
    return {k: v["verdict"] for k, v in flares.classify({**BASE, **kw}).items()}


def test_table_rules_follow_the_plan():
    v = verdicts()
    assert v == {"artefact": "불일치", "supernova": "부합", "obscuration_change": "부합",
                 "tde_in_agn": "부합", "agn_flare": "부합"}
    assert verdicts(o1_checks={"a": True, "b": False})["artefact"] == "부합"
    assert verdicts(nu_l_nu_w2_peak=2e42)["supernova"] == "불일치"
    assert verdicts(e_dust=2e51)["supernova"] == "불일치"
    assert verdicts(n_episodes=2)["supernova"] == "불일치"
    assert verdicts(dm_w2_minus_dm_w1=0.3)["obscuration_change"] == "불일치"
    assert verdicts(dm_w2_minus_dm_w1=0.09)["obscuration_change"] == "부합"  # below 2 sigma
    for change in ({"n_episodes": 2}, {"rebrightening_after_peak": True}, {"t_rise_yr": 2.0},
                   {"e_dust": 1e53}):
        assert verdicts(**change)["tde_in_agn"] == "불일치"
    assert verdicts(t_decay_efold_yr=float("nan"))["tde_in_agn"] == "판정불가"
    off = verdicts(position_shift_arcsec=1.5, position_shift_err_arcsec=0.1)
    assert off["tde_in_agn"] == "불일치" and off["agn_flare"] == "불일치"
    assert verdicts(position_shift_arcsec=1.5, position_shift_err_arcsec=1.0)["agn_flare"] == "판정불가"


def _fits(path, time, rate):
    def card(k, v):
        return f"{k:<8}= {v:>20}".ljust(80)
    def block(cards):
        text = "".join(cards) + "END".ljust(80)
        return (text + " " * (-len(text) % 2880)).encode("ascii")
    primary = block([card("SIMPLE", "T"), card("BITPIX", 8), card("NAXIS", 0)])
    data = np.zeros(len(time), dtype=[("TIME", ">f8"), ("RATE", ">f8"), ("RATE_ERR", ">f8")])
    data["TIME"], data["RATE"], data["RATE_ERR"] = time, rate, 0.1
    cards = [card("XTENSION", "'BINTABLE'"), card("BITPIX", 8), card("NAXIS", 2),
             card("NAXIS1", 24), card("NAXIS2", len(time)), card("PCOUNT", 0), card("GCOUNT", 1),
             card("TFIELDS", 3)]
    for i, name in enumerate(("TIME", "RATE", "RATE_ERR"), 1):
        cards += [card(f"TTYPE{i}", f"'{name}'"), card(f"TFORM{i}", "'D'")]
    raw = data.tobytes()
    path.write_bytes(primary + block(cards) + raw + b"\0" * (-len(raw) % 2880))


def test_fits_reader_and_bat_yearly(tmp_path):
    time = (np.array([53200.0, 53800.0, 53830.0]) - otherbands.SWIFT_MET_MJD0) * 86400.0
    _fits(tmp_path / "x.lc", time, [1.0, 2.0, 4.0])
    t = otherbands.read_fits_table(tmp_path / "x.lc")
    assert np.allclose(t["RATE"], [1.0, 2.0, 4.0]) and t["TIME"].dtype.kind == "f"
    yearly = otherbands.bat_yearly(t)
    assert [y["year"] for y in yearly] == [2004, 2006]
    assert yearly[1]["rate"] == pytest.approx(3.0) and yearly[1]["n_months"] == 2
    with pytest.raises(ValueError):
        otherbands.read_fits_table(tmp_path / "x.lc", hdu=3)


def test_asassn_reader_and_optical_change(tmp_path):
    rng = np.random.default_rng(3)
    shape = bump(21, 9, 11, 4.0)
    v = visits(0.6 * shape, shape)
    rows = ["jd,mag,mag_err,limit,quality,phot_filter,camera"]
    for k, t in enumerate(v["mjd"]):
        for _ in range(8):
            rows.append(f"{t + rng.uniform(-30, 30) + 2400000.5:.4f},{15.0 - 0.3 * shape[k] + rng.normal(0, 0.02):.4f},"
                        "0.02,17.0,G,g,ba")
    rows.append("2456700.5,99.99,99.99,15.0,B,g,ba")  # bad row, dropped
    p = tmp_path / "a.csv"
    p.write_text("\n".join(rows) + "\n")
    lc = otherbands.read_asassn(p)
    assert set(lc) == {"g"} and lc["g"]["mjd"].size == 21 * 8
    res = otherbands.optical_change(v, lc["g"]["mjd"], lc["g"]["mag"])
    assert res["status"] == "usable" and res["dm_opt"] == pytest.approx(0.2, abs=0.06)
    assert res["upper_limit_3sigma"] > res["dm_opt"]
    few = otherbands.optical_change(v, lc["g"]["mjd"][:24], lc["g"]["mag"][:24])
    assert few["status"] == "fewer than 6 visits with an optical bin"


RESULT = config.RESULTS_DIR / "s6a_flares.json"


@pytest.mark.skipif(not RESULT.exists(), reason="stage 6A results not present")
def test_committed_table_follows_from_the_committed_measurements():
    d = json.loads(RESULT.read_text())
    assert set(d) == {"bat0463", "bat0582", "bat0074", "bat0686", "bat0560"}
    for name, o in d.items():
        m = {k: (float("nan") if v is None else v) for k, v in o["measured_for_table"].items()}
        again = flares.classify(m)
        assert {k: v["verdict"] for k, v in again.items()} == \
            {k: v["verdict"] for k, v in o["table"].items()}, name
    # the reference object (a published changing-look AGN) must not be called inconsistent
    assert d["bat0560"]["table"]["agn_flare"]["verdict"] != "불일치"
