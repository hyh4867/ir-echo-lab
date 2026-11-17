"""Stage 6 figures from results/s6a_* and results/s6b_*.

    python scripts/plot_stage6.py [a|b]

Without an argument both figures are made. The 6B figure needs
results/s6b_objects.csv and results/s6b_summary.json (scripts/run_stage6b.py).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from irecho import config, plots5b  # noqa: E402


def main() -> None:
    which = sys.argv[1] if len(sys.argv) > 1 else "ab"
    res = config.RESULTS_DIR
    made = []
    if "a" in which:
        made.append(plots5b.plot_flares(res / "s6a_flares.json", res / "s6a_other_bands.json",
                                        res / "fig_s6a_flares.png"))
    if "b" in which:
        made.append(plots5b.plot_d_luminosity(res / "s6b_objects.csv", res / "s5b_objects.csv",
                                              res / "s6b_power_final.json", res / "s6b_summary.json",
                                              res / "fig_s6b_d_luminosity.png"))
    for path in made:
        print(path.relative_to(config.PROJECT_ROOT))


if __name__ == "__main__":
    main()
