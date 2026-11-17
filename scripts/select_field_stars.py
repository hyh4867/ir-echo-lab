"""Select comparison stars around each control galaxy (stage 3-4, step 0).

    python scripts/select_field_stars.py

The stars are used only to measure the visit-to-visit zero-point offset of
each control field (differential photometry). They obey the same criteria
R1-R13 as the NGC 3786 reference stars (scripts/select_targets.py), with two
differences fixed here before any of their light curves was downloaded:

  * W1 range 9.5-11.5 (the brightness range of the control nuclei, W1 9.0-10.3,
    is covered as far as the saturation margin allows; fainter stars have
    larger errors and add little to a zero-point estimate),
  * 5 stars from [9.5, 10.5) and 5 from [10.5, 11.5) per field.

Output: targets/field_stars.csv, targets/field_stars_funnel.json
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import select_targets as st  # noqa: E402
from irecho import config  # noqa: E402

W1_RANGE = (9.5, 11.5)
MAG_BINS = ((9.5, 10.5), (10.5, 11.5))
STARS_PER_BIN = 5


def main() -> None:
    with open(config.TARGETS_DIR / "control_agn.csv", newline="") as fh:
        controls = [r for r in csv.DictReader(fh) if r["selected"] == "yes"]
    rows, funnels = [], {}
    for field_index, c in enumerate(controls):
        ra, dec = float(c["ra"]), float(c["dec"])
        faint_limit = W1_RANGE[1] + st.ISOLATION_DMAG
        field = st.allwise_cone(
            f"field_{c['name']}", ra, dec,
            st.FIELD_RADIUS_DEG * 3600.0 + 2 * st.ISOLATION_ARCSEC,
            extra_where=f"w1mpro < {faint_limit:.1f}",
        )
        ok, funnel, dist, pm = st.star_candidates(field, ra, dec, W1_RANGE)
        w1, w2 = field["w1mpro"], field["w2mpro"]
        rng = np.random.default_rng([config.SEED, field_index])
        chosen: list[int] = []
        for lo, hi in MAG_BINS:
            idx = np.flatnonzero(ok & (w1 >= lo) & (w1 < hi))
            idx = idx[np.argsort(field["designation"][idx])]
            take = min(STARS_PER_BIN, idx.size)
            picked = rng.choice(idx, size=take, replace=False) if take else np.array([], int)
            funnel.append((f"bin {lo}-{hi}: candidates / drawn", f"{idx.size} / {take}"))
            chosen.extend(sorted(int(i) for i in picked))
        funnels[c["name"]] = funnel
        for rank, i in enumerate(chosen, start=1):
            rows.append({
                "name": f"{c['name']}_s{rank:02d}", "field": c["name"],
                "allwise_designation": str(field["designation"][i]),
                "ra": f"{field['ra'][i]:.7f}", "dec": f"{field['dec'][i]:.7f}",
                "w1mpro": f"{w1[i]:.3f}", "w2mpro": f"{w2[i]:.3f}",
                "w1rchi2": f"{field['w1rchi2'][i]:.2f}", "w2rchi2": f"{field['w2rchi2'][i]:.2f}",
                "var_flg": str(field["var_flg"][i]),
                "j_ks": f"{field['j_m_2mass'][i] - field['k_m_2mass'][i]:.3f}",
                "pm_mas_yr": f"{pm[i]:.0f}",
                "dist_to_target_deg": f"{dist[i] / 3600.0:.3f}",
            })
        print(f"{c['name']}: {funnel[-3][1]} candidates, drew {len(chosen)}", flush=True)
    with open(config.TARGETS_DIR / "field_stars.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    (config.TARGETS_DIR / "field_stars_funnel.json").write_text(
        json.dumps({"seed": [config.SEED, "field index"], "w1_range": W1_RANGE,
                    "funnel": funnels}, indent=1) + "\n")
    print(f"wrote {len(rows)} stars")


if __name__ == "__main__":
    main()
