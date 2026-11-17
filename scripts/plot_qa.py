"""Stage 5a: QA figures from the committed results/qa_* files (no target light curve).

    python scripts/plot_qa.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from irecho import config, plots5, sample_io  # noqa: E402


def main() -> None:
    res = config.RESULTS_DIR
    made = [
        plots5.plot_star_chi2(res / "qa_stars.csv", res / "fig_qa_star_chi2.png"),
        plots5.plot_offsets(res / "qa_offsets.csv", res / "fig_qa_offsets.png"),
        plots5.plot_pair_scatter(res / "qa_stars.csv", res / "qa_targets.json",
                                 sample_io.SAMPLE_PATH, res / "fig_qa_pair_scatter.png"),
    ]
    if (res / "qa_injection_stars.json").exists():
        made.append(plots5.plot_star_injection(res / "qa_injection_stars.json",
                                               res / "fig_qa_injection_stars.png"))
    for path in made:
        print(path.relative_to(config.PROJECT_ROOT))


if __name__ == "__main__":
    main()
