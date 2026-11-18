"""Stage 7 figure from results/s7_summary.json and results/s7_spectra_plot.json.

    python scripts/plot_stage7.py

Needs scripts/run_stage7.py and scripts/xray_fit.py plotdata to have run.
The event-image figure (results/fig_s7_xray_images.png) is made by
scripts/xray_regions.py.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from irecho import config, plots7  # noqa: E402


def main() -> None:
    res = config.RESULTS_DIR
    path = plots7.plot_xray(res / "s7_summary.json", res / "s7_spectra_plot.json", res / "fig_s7_xray_spectra.png")
    print(path.relative_to(config.PROJECT_ROOT))


if __name__ == "__main__":
    main()
