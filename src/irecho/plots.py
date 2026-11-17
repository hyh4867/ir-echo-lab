"""Static figures for results/. Labels are in English so that the figures do
not depend on a Korean font being installed."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from .analysis import PAPER_DATA_END_MJD, run_detector, visit_arrays  # noqa: E402
from .lightcurve import LightCurve, mjd_to_year  # noqa: E402

# Two-series categorical palette (blue / orange), checked for colour-vision
# deficiency separation against the light surface.
COLOR = {"w1": "#2a78d6", "w2": "#eb6834"}
MARKER = {"w1": "o", "w2": "s"}  # second encoding besides colour
LABEL = {"w1": "W1 (3.4 um)", "w2": "W2 (4.6 um)"}
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e4e3df"
NEUTRAL = "#8a8983"
SHADE = "#f3e3c3"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": NEUTRAL, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK_2, "ytick.color": INK_2, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False,
    "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold",
    "axes.titlelocation": "left", "legend.frameon": False, "lines.linewidth": 2.0,
    "figure.dpi": 110, "savefig.dpi": 160,
})


def _save(fig, path: Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_reference_stars(star_result: dict, path: Path) -> Path:
    """Left: distribution of reduced chi-square per star. Right: observed
    visit-to-visit RMS against the median adopted visit error."""
    rows = star_result["per_star"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.3))
    all_chi = np.array([[r[f"{b}_chi2_red_err"] for r in rows] for b in ("w1", "w2")])
    bins = np.logspace(np.log10(min(0.3, all_chi.min() * 0.9)),
                       np.log10(max(30.0, all_chi.max() * 1.1)), 19)
    for b in ("w1", "w2"):
        chi = np.array([r[f"{b}_chi2_red_err"] for r in rows])
        ax1.hist(chi, bins=bins, histtype="stepfilled", alpha=0.25, color=COLOR[b])
        ax1.hist(chi, bins=bins, histtype="step", color=COLOR[b], linewidth=2,
                 label=f"{LABEL[b]}: median {np.median(chi):.1f}")
    ax1.axvline(1.0, color=INK, linestyle="--", linewidth=1.2)
    ax1.annotate("errors explain\nthe scatter", xy=(1.0, 0.97), xycoords=("data", "axes fraction"),
                 xytext=(-6, 0), textcoords="offset points", ha="right", va="top",
                 color=INK_2, fontsize=9)
    ax1.set_xscale("log")
    ax1.set_xlabel("reduced chi-square of a constant model (NEOWISE visits, adopted errors)")
    ax1.set_ylabel("number of reference stars")
    ax1.set_title(f"Reference stars (n = {len(rows)}): is the scatter explained?")
    ax1.set_ylim(0, ax1.get_ylim()[1] * 1.3)  # headroom so the legend clears the bars
    ax1.legend(loc="upper right")

    lim = 0.0
    for b in ("w1", "w2"):
        err = 1e3 * np.array([r[f"{b}_median_err_mag"] for r in rows])
        rms = 1e3 * np.array([r[f"{b}_rms_mag"] for r in rows])
        ax2.scatter(err, rms, s=42, color=COLOR[b], marker=MARKER[b], edgecolor=SURFACE,
                    linewidth=1.0, label=LABEL[b], zorder=3)
        lim = max(lim, err.max(), rms.max())
    lim *= 1.08
    ax2.plot([0, lim], [0, lim], color=INK, linestyle="--", linewidth=1.2)
    ax2.annotate("RMS = error", xy=(lim * 0.93, lim * 0.93), xytext=(-4, -12),
                 textcoords="offset points", ha="right", color=INK_2, fontsize=9)
    ax2.set_xlim(0, lim)
    ax2.set_ylim(0, lim)
    ax2.set_xlabel("median adopted visit error (mmag)")
    ax2.set_ylabel("observed visit-to-visit RMS (mmag)")
    ax2.set_title("Observed scatter vs. error model")
    ax2.legend(loc="upper left")
    return _save(fig, path)


def plot_common_mode(common: dict, path: Path) -> Path:
    """Median residual of the reference stars at each NEOWISE visit."""
    fig, ax = plt.subplots(figsize=(9, 3.8))
    year = np.array(common["visit_year"])
    for b, dx in (("w1", -0.03), ("w2", 0.03)):
        off = 1e3 * np.array(common[b]["offset_mag"])
        err = 1e3 * np.array(common[b]["offset_err_mag"])
        ax.errorbar(year + dx, off, yerr=err, fmt=MARKER[b], color=COLOR[b], markersize=6,
                    markeredgecolor=SURFACE, elinewidth=1.2, capsize=0,
                    label=f"{LABEL[b]}: std {1e3 * common[b]['offset_std_mag']:.1f} mmag")
    ax.axhline(0.0, color=INK, linewidth=1.0)
    ax.set_xlabel("year")
    ax.set_ylabel("median star residual (mmag)\npositive = fainter")
    ax.set_title("Visit-to-visit offset shared by the reference stars (NGC 3786 field)")
    ax.legend(loc="lower right", ncol=1)
    return _save(fig, path)


def plot_ngc3786(lc: LightCurve, summary: dict, path: Path) -> Path:
    """Three stacked panels: W1, W2 and W1-W2 against time."""
    e = lc.epochs
    year = mjd_to_year(e["mjd"])
    fig, axes = plt.subplots(3, 1, figsize=(9.5, 8.4), sharex=True,
                             gridspec_kw={"height_ratios": [1, 1, 0.8], "hspace": 0.12})
    flagged_years = [row["year"] for row in summary["visit_table"] if row["flagged"]]
    paper_end = float(mjd_to_year(PAPER_DATA_END_MJD))
    for ax in axes:
        if flagged_years:
            ax.axvspan(min(flagged_years) - 0.2, max(flagged_years) + 0.2, color=SHADE,
                       zorder=0, linewidth=0)
        ax.axvline(paper_end, color=INK_2, linestyle=":", linewidth=1.2, zorder=1)
    for ax, b in zip(axes[:2], ("w1", "w2")):
        m = lc.masks[b]
        ax.scatter(mjd_to_year(lc.obs["mjd"][m]), lc.obs[f"{b}mag"][m], s=7, color=NEUTRAL,
                   alpha=0.45, linewidth=0, zorder=2, label="single exposures (after cuts)")
        ax.errorbar(year, e[f"{b}_mag"], yerr=e[f"{b}_err"], fmt=MARKER[b], color=COLOR[b],
                    markersize=7, markeredgecolor=SURFACE, elinewidth=1.4, zorder=3,
                    label="visit median")
        ax.invert_yaxis()
        ax.set_ylabel(f"{b.upper()} (Vega mag)")
        ax.legend(loc="upper left", ncol=2, fontsize=9)
    axes[2].errorbar(year, e["color"], yerr=e["color_err"], fmt="D", color=INK_2, markersize=6,
                     markeredgecolor=SURFACE, elinewidth=1.4, zorder=3)
    axes[2].set_ylabel("W1 - W2 (mag)")
    axes[2].set_xlabel("year")
    top = axes[0]
    top.set_title("NGC 3786: WISE/NEOWISE profile-fit photometry, visit medians")
    if flagged_years:
        axes[2].annotate("visits flagged as brightening",
                         xy=(min(flagged_years) - 0.1, 0.05),
                         xycoords=("data", "axes fraction"), fontsize=9, color=INK_2,
                         va="bottom")
    top.annotate("end of data in\nSon et al. (2022) Fig. 1", xy=(paper_end, 0.36),
                 xycoords=("data", "axes fraction"), xytext=(-6, 0), textcoords="offset points",
                 ha="right", va="center", fontsize=9, color=INK_2)
    return _save(fig, path)


def plot_controls(target: LightCurve, controls: list[LightCurve], path: Path) -> Path:
    """Small multiples: change relative to each object's own NEOWISE median, on
    one shared scale, with NGC 3786 first for reference."""
    panels = [target] + list(controls)
    ncols = 3
    nrows = int(np.ceil(len(panels) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(11, 2.7 * nrows), sharex=True, sharey=True)
    axes = np.atleast_2d(axes)
    for ax, lc in zip(axes.flat, panels):
        v = visit_arrays(lc, "neowise")
        res = run_detector(v)
        year = mjd_to_year(v["mjd"])
        for b in ("w1", "w2"):
            delta = np.median(v[f"{b}_mag"]) - v[f"{b}_mag"]  # positive = brighter
            ax.errorbar(year, delta, yerr=v[f"{b}_err"], fmt=MARKER[b] + "-", color=COLOR[b],
                        markersize=4.5, markeredgecolor=SURFACE, markeredgewidth=0.6,
                        linewidth=1.2, elinewidth=1.0, label=LABEL[b])
        ax.axhline(0.0, color=INK, linewidth=0.8)
        n_flag = int(res.flagged.sum())
        tag = f"{n_flag} flagged" if n_flag else "none flagged"
        ax.set_title(f"{lc.name}  ({tag})", fontsize=10)
    for ax in axes.flat[len(panels):]:
        ax.set_visible(False)
    for ax in axes[:, 0]:
        ax.set_ylabel("brightening (mag)")
    for ax in axes[-1, :]:
        ax.set_xlabel("year")
    axes.flat[0].legend(loc="upper left", fontsize=8)
    fig.suptitle("NGC 3786 and the Seyfert 2 control sample: change relative to each "
                 "object's NEOWISE median", x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return _save(fig, path)
