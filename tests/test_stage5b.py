"""Stage 5b decision rules and per-object measures on synthetic inputs (no data)."""

import csv
import json

import numpy as np
import pytest

from irecho import stage5b as s5


def visits(n=21, w1=None, w2=None, err=0.02):
    mjd = 56650.0 + 182.6 * np.arange(n)
    w1 = np.full(n, 10.0) if w1 is None else np.asarray(w1, dtype=float)
    w2 = np.full(n, 9.5) if w2 is None else np.asarray(w2, dtype=float)
    e = np.full(n, err)
    return {"mjd": mjd, "w1_mag": w1, "w1_err": e, "w2_mag": w2, "w2_err": e.copy(),
            "color": w1 - w2, "color_err": np.hypot(e, e)}


def test_intervals_known_values():
    lo, hi = s5.wilson(16, 225)
    assert (lo, hi) == pytest.approx((0.0442, 0.1124), abs=5e-4)
    lo, hi = s5.clopper_pearson(0, 10)
    assert lo == 0.0 and hi == pytest.approx(1 - 0.025 ** 0.1, abs=1e-9)
    assert s5.clopper_pearson(10, 10)[1] == 1.0
    assert np.isnan(s5.wilson(0, 0)[0])


@pytest.mark.parametrize("lo, hi, expected", [
    (0.06, 0.15, "지지"), (0.05, 0.20, "지지"), (0.044, 0.112, "판정불가"),
    (0.01, 0.049, "기각"), (0.21, 0.30, "기각"), (0.15, 0.25, "판정불가"),
])
def test_verdict_q1(lo, hi, expected):
    assert s5.verdict_q1(lo, hi) == expected


def test_verdict_majority():
    assert s5.verdict_majority(18, 20)[0] == "지지"
    assert s5.verdict_majority(2, 20)[0] == "기각"
    assert s5.verdict_majority(11, 20)[0] == "판정불가"
    assert s5.verdict_majority(9, 9)[0] == "판정불가 (표본 부족)"  # fewer than 10 objects


def test_detectors_and_q1_class():
    flat = s5.run_detectors(visits())
    assert not flat["brightened"] and not flat["q1"] and flat["start_index"] is None
    rng = np.random.default_rng(1)
    ramp = np.concatenate([np.zeros(9), np.linspace(0.0, 0.5, 12)])
    noise = rng.normal(0.0, 0.01, (2, 21))
    v = visits(w1=10.0 - ramp + noise[0], w2=9.5 - ramp + noise[1])
    d = s5.run_detectors(v)
    assert d["b_flagged"] and d["q1"] and d["b_rise_years"] > 3
    small = np.concatenate([np.zeros(9), np.linspace(0.0, 0.14, 12)])
    v = visits(w1=10.0 - small + noise[0], w2=9.5 - small + noise[1], err=0.01)
    d = s5.run_detectors(v)
    assert d["b_flagged"] and not d["q1"]  # flagged at 0.1 mag but below the 0.2 mag class


def test_colour_change_sign_and_baseline_rule():
    ramp = np.concatenate([np.zeros(9), np.linspace(0.0, 0.5, 12)])
    v = visits(w1=10.0 - 0.6 * ramp, w2=9.5 - ramp)  # W2 brightens more: redder
    c = s5.colour_change(v, start_index=9)
    assert c["n_base"] == 9 and c["delta_color"] == pytest.approx(0.4 * np.mean(ramp[-3:]), abs=1e-9)
    assert c["delta_color_err"] > 0 and c["dm_w2"] > c["dm_w1"] > 0
    assert s5.colour_change(v, start_index=0)["n_base"] == 3
    assert np.isnan(s5.colour_change(v, start_index=19)["delta_color"])  # fewer than 3 visits left


def _write_ztf(path, mjd, mag):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["oid", "mjd", "mag", "magerr", "catflags"])
        for t, m in zip(mjd, mag):
            w.writerow([1, f"{t:.4f}", f"{m:.4f}", "0.02", 0])
    path.with_suffix(".csv.meta.json").write_text(json.dumps({"status": "ok"}))


def test_ztf_band_verdicts(tmp_path):
    rng = np.random.default_rng(2)
    ramp = np.concatenate([np.zeros(9), np.linspace(0.0, 0.5, 12)])
    v = visits(w1=10.0 - ramp, w2=9.5 - ramp)
    t = np.repeat(v["mjd"][9:], 8) + rng.uniform(-40, 40, 12 * 8)
    follow = 16.0 - 0.3 * np.repeat(ramp[9:], 8) + rng.normal(0, 0.02, t.size)
    flat = 16.0 + rng.normal(0, 0.02, t.size)
    for name, mag, expected in (("follow", follow, "accompanies"), ("flat", flat, "no change"),
                                ("anti", 32.0 - follow, "opposite")):
        path = tmp_path / f"{name}.csv"
        _write_ztf(path, t, mag)
        res = s5.ztf_band(v, path)
        assert res["status"] == "usable" and res["n_bins"] == 12 and res["verdict"] == expected
    quiet = visits()
    path = tmp_path / "quiet.csv"
    _write_ztf(path, t, flat)
    assert s5.ztf_band(quiet, path)["status"] == "IR change below 0.1 mag in the ZTF era"
    assert s5.ztf_band(v, tmp_path / "missing.csv")["status"] == "not fetched"
    few = tmp_path / "few.csv"
    _write_ztf(few, t[:24], follow[:24])
    assert s5.ztf_band(v, few)["status"] == "fewer than 6 visits with a ZTF bin"


def test_ztf_compare_labels_missing_coverage():
    v = visits()
    out = s5.ztf_compare(v, "nothing", {}, dec_deg=-45.0)
    assert out["ztf_g_status"] == "outside the ZTF declination range" and not out["ztf_usable"]
    out = s5.ztf_compare(v, "nothing", {}, dec_deg=20.0)
    assert out["ztf_r_status"] == "no ZTF object within 3 arcsec" and out["ztf_verdict"] == ""


def test_compare_groups_rules():
    rng = np.random.default_rng(3)
    above, below = rng.normal(0.0, 0.2, 12), rng.normal(1.2, 0.2, 12)
    res = s5.compare_groups(above, below, 1)
    assert res["verdict"] == "지지" and res["delta_r"] > 0.5 and res["p_one_sided"] < 0.05
    same = s5.compare_groups(rng.normal(0.0, 0.05, 30), rng.normal(0.0, 0.05, 30), 2)
    assert same["verdict"] == "기각" and same["delta_r_ci95"][1] < 0.5
    wide = s5.compare_groups(rng.normal(0.0, 1.5, 6), rng.normal(0.4, 1.5, 6), 3)
    assert wide["verdict"] == "판정불가"
    assert s5.compare_groups(above[:4], below, 4)["verdict"] == "판정불가 (표본 부족)"
    assert s5.compare_groups([], [], 5)["verdict"] == "판정불가 (표본 부족)"


def test_hysteresis_rules():
    rng = np.random.default_rng(4)
    assert s5.hysteresis_test(rng.normal(0.2, 0.1, 20), 0.0245, 1)["verdict"] == "지지"
    assert s5.hysteresis_test(rng.normal(-0.05, 0.02, 20), 0.0245, 2)["verdict"] == "기각"
    assert s5.hysteresis_test(rng.normal(0.02, 0.3, 10), 0.0245, 3)["verdict"] == "판정불가"
    assert s5.hysteresis_test(rng.normal(0.2, 0.1, 7), 0.0245, 4)["verdict"] == "판정불가 (표본 부족)"


def test_corrected_fraction_and_helpers():
    res = s5.corrected_fraction([1.0, 2.0, 5.0], 100, 1)
    assert res["fraction"] == pytest.approx(0.08)
    assert res["ci95"][0] < 0.08 < res["ci95"][1]
    assert s5.nearest([0.2, 0.3, 0.5, 0.75, 1.0], 0.38, log=True) == 0.3
    assert s5.nearest([0.2, 0.3, 0.5, 0.75, 1.0], 0.39, log=True) == 0.5  # log midpoint 0.387
    assert s5.nearest([0, 10, 20, 50], 36.3) == 50.0
    assert s5.fisher_2x2(5, 9, 6, 0) == pytest.approx(0.0141, abs=5e-4)
    assert s5.fisher_2x2(0, 0, 3, 4) is None
    m, e = s5.weighted_mean([1.0, 3.0, np.nan], [1.0, 1.0, 1.0])
    assert m == pytest.approx(2.0) and e == pytest.approx(2 ** -0.5)
