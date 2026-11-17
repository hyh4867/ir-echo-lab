"""Stage 5b figures from results/s5b_* and the frozen light curves.

    python scripts/plot_stage5b.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from irecho import config, plots5b  # noqa: E402


def main() -> None:
    res = config.RESULTS_DIR
    obj, summ, val = res / "s5b_objects.csv", res / "s5b_summary.json", res / "s5b_validation.json"
    made = [
        plots5b.plot_validation(val, res / "fig_s5b_validation.png"),
        plots5b.plot_lightcurves(obj, res / "fig_s5b_lightcurves.png"),
        plots5b.plot_q1(summ, res / "fig_s5b_q1.png"),
        plots5b.plot_colour(obj, res / "fig_s5b_colour.png"),
        plots5b.plot_ztf(obj, res / "fig_s5b_ztf.png"),
        plots5b.plot_slope(obj, summ, res / "fig_s5b_slope.png"),
        plots5b.plot_hysteresis(obj, val, res / "fig_s5b_hysteresis.png"),
    ]
    for path in made:
        print(path.relative_to(config.PROJECT_ROOT))


if __name__ == "__main__":
    main()
