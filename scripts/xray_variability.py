"""Stage 7 (check added after the fits, not in docs/XRAY_PLAN.md): did the source
vary during the 2022 observations? Source counts per NuSTAR orbit and per XRT
good-time interval.

    python scripts/xray_variability.py      (inside the HEASoft environment)

Output: results/s7_variability.json.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from irecho import config, xray  # noqa: E402

WORK = xray.S7_DIR / "work"


def main() -> None:
    reg = json.loads((config.RESULTS_DIR / "s7_regions.json").read_text())
    out = {"note": "computed after the registered results were seen; rates are counts per on-time, background included"}
    per_module = {}
    for mod, key in (("A", "fpma"), ("B", "fpmb")):
        with fits.open(WORK / "nustar" / "pipe" / f"nu60160365002{mod}01_cl.evt") as h:
            d, g, hd = h["EVENTS"].data, h["GTI"].data, h["EVENTS"].header
            mjdref = hd["MJDREFI"] + hd["MJDREFF"]
            cx, cy = reg[key]["centroid_xy"]
            t = d["TIME"][np.hypot(d["X"] - cx, d["Y"] - cy) * reg[key]["pixel_scale_arcsec"] < 50.0]
            orbits = [[g["START"][0], g["STOP"][0], g["STOP"][0] - g["START"][0]]]
            for a, b in zip(g["START"][1:], g["STOP"][1:]):
                if a - orbits[-1][1] < 300:
                    orbits[-1][1], orbits[-1][2] = b, orbits[-1][2] + (b - a)
                else:
                    orbits.append([a, b, b - a])
            per_module[key] = [{"mjd_start": float(mjdref + a / 86400), "mjd_stop": float(mjdref + b / 86400),
                                "ontime_s": float(ont), "counts": int(((t >= a) & (t <= b)).sum())}
                               for a, b, ont in orbits]
    out["nustar_orbits"] = per_module
    n = np.array([a["counts"] + b["counts"] for a, b in zip(per_module["fpma"], per_module["fpmb"])], dtype=float)
    ont = np.array([a["ontime_s"] + b["ontime_s"] for a, b in zip(per_module["fpma"], per_module["fpmb"])])
    rate, err = n / ont, np.sqrt(n) / ont
    mean = float(np.sum(rate / err ** 2) / np.sum(1 / err ** 2))
    out["nustar_summed"] = {"rate": [float(x) for x in rate], "rate_err": [float(x) for x in err], "weighted_mean": mean,
                            "chi2_constant": float(np.sum(((rate - mean) / err) ** 2)), "dof": int(rate.size - 1),
                            "max_over_min": float(rate.max() / rate.min())}
    with fits.open(WORK / "xrt22" / "pipe" / "sw00081051001xpcw3po_cl.evt") as h:
        d, g, hd = h["EVENTS"].data, h["GTI"].data, h["EVENTS"].header
        mjdref = hd["MJDREFI"] + hd["MJDREFF"]
        cx, cy = reg["xrt22"]["centroid_xy"]
        e = d["PI"] * 0.01
        t = d["TIME"][(np.hypot(d["X"] - cx, d["Y"] - cy) < 20.0) & (e >= 2) & (e <= 10)]
        rows = [{"mjd_start": float(mjdref + a / 86400), "mjd_stop": float(mjdref + b / 86400), "duration_s": float(b - a),
                 "counts_2_10": int(((t >= a) & (t <= b)).sum())} for a, b in zip(g["START"], g["STOP"])]
    n = np.array([r["counts_2_10"] for r in rows], dtype=float)
    dur = np.array([r["duration_s"] for r in rows])
    expect = n.sum() * dur / dur.sum()
    out["xrt22_gtis"] = rows
    out["xrt22_summary"] = {"chi2_constant_pearson": float(np.sum((n - expect) ** 2 / expect)), "dof": int(n.size - 1)}
    (config.RESULTS_DIR / "s7_variability.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({"nustar": out["nustar_summed"], "xrt22": out["xrt22_summary"]}, indent=1))


if __name__ == "__main__":
    main()
