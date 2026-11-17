"""Stage 6A (A4): optical and hard X-ray behaviour of the flare candidates.

    python scripts/run_flares_other.py

Inputs: data/cache/s6/asassn/<name>.csv (ASAS-SN Sky Patrol v2),
data/cache/s6/bat157/<name>.lc (Swift-BAT 157-month crab-weighted monthly
light curves), the frozen WISE visits. Output: results/s6a_other_bands.json.
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from irecho import config, otherbands as ob  # noqa: E402
from irecho import stage5b as s5  # noqa: E402

NAMES = ("bat0463", "bat0582", "bat0074", "bat0686")


def main() -> None:
    flares = json.loads((config.RESULTS_DIR / "s6a_flares.json").read_text())
    out = {}
    for name in NAMES:
        v = s5.read_frozen(name)
        f = flares[name]
        res = {"counterpart": f["counterpart"], "asassn": {}, "bat157": {}}
        lc = ob.read_asassn(ob.S6_DIR / "asassn" / f"{name}.csv")
        for band, d in lc.items():
            cam = Counter(d["camera"].tolist()).most_common(1)[0][0]
            sel = d["camera"] == cam
            r = ob.optical_change(v, d["mjd"][sel], d["mag"][sel])
            r.update(camera=cam, n_good_points_all_cameras=int(d["mjd"].size),
                     n_points_used=int(sel.sum()),
                     mjd_range=[float(d["mjd"][sel].min()), float(d["mjd"][sel].max())])
            res["asassn"][band] = r
        t = ob.read_fits_table(ob.S6_DIR / "bat157" / f"{name}.lc")
        mjd = ob.SWIFT_MET_MJD0 + t["TIME"] / 86400.0
        yearly = ob.bat_yearly(t)
        allmean, allerr = ob.weighted_mean(t["RATE"], t["RATE_ERR"])
        chi2 = float(np.sum(((np.array([y["rate"] for y in yearly]) - allmean)
                             / np.array([y["rate_err"] for y in yearly])) ** 2))
        year = 2000.0 + (mjd - 51544.5) / 365.25
        early = (year >= 2012.0) & (year < 2016.0)
        late = (year >= 2016.0) & (year < 2018.0)
        m_e, e_e = ob.weighted_mean(t["RATE"][early], t["RATE_ERR"][early])
        m_l, e_l = ob.weighted_mean(t["RATE"][late], t["RATE_ERR"][late])
        compare = {"rate_2012_2015": m_e, "rate_2012_2015_err": e_e, "rate_2016_2017": m_l,
                   "rate_2016_2017_err": e_l, "difference_sigma": (m_l - m_e) / float(np.hypot(e_e, e_l)),
                   "note": "same windows for every object; descriptive, not in the plan"}
        res["bat157"] = {"compare_last_two_years": compare,
                         "mjd_range": [float(mjd.min()), float(mjd.max())],
                         "mean_rate": allmean, "mean_rate_err": allerr, "yearly": yearly,
                         "yearly_chi2_constant": chi2, "yearly_dof": len(yearly) - 1,
                         "covers_main_episode": bool(mjd.max() >= f["main_episode"]["start_mjd"]),
                         "main_episode_start_mjd": f["main_episode"]["start_mjd"]}
        out[name] = res
        a = res["asassn"]
        print(name, f["counterpart"], "| ASAS-SN",
              {b: (a[b]["status"], a[b]["n_bins"], round(a[b].get("dm_opt", float("nan")), 3),
                   round(a[b].get("dm_opt_err", float("nan")), 3),
                   round(a[b].get("ir_change_w2", float("nan")), 2)) for b in a},
              "| BAT chi2/dof", round(chi2, 1), len(yearly) - 1, "covers", res["bat157"]["covers_main_episode"])
        print("    BAT 2016-17 vs 2012-15:", {k: (round(1e5 * x, 1) if "rate" in k else round(x, 2))
                                           for k, x in compare.items() if k != "note"})
        print("    BAT yearly (rate/err):", [(y["year"], round(y["rate"] / y["rate_err"], 1)) for y in yearly])
        print("    BAT yearly rate x1e5:", [round(1e5 * y["rate"], 1) for y in yearly])
    (config.RESULTS_DIR / "s6a_other_bands.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
