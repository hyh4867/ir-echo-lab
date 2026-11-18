"""Stage 7 figure: X-ray spectra of CGCG 312-012 and the registered comparisons.

Same visual system as plots.py. Colour follows the epoch (2010, 2022, BAT
average); the instrument is shown by the marker.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402,F401
import numpy as np  # noqa: E402

from . import xray  # noqa: E402
from .plots import GRID, INK, INK_2, NEUTRAL, _save  # noqa: E402

C2010, C2022, CBAT = "#2a78d6", "#eb6834", "#1baf7a"
SERIES = {  # key: (label, colour, marker, filled)
    "xrt10": ("Swift XRT, 2010-12-29", C2010, "o", True),
    "xrt22": ("Swift XRT, 2022-05-01", C2022, "o", True),
    "fpma": ("NuSTAR FPMA, 2022-05-01", C2022, "s", True),
    "fpmb": ("NuSTAR FPMB, 2022-05-01", C2022, "^", False),
    "bat": ("Swift-BAT, 2004-2017 average", CBAT, "D", True),
}
MODEL_TITLE = {"M0": "M0: absorbed power law", "M1": "M1: + scattered power law",
               "M2": "M2: + reflection and Fe line"}


LABEL_EN = {xray.HIGH: "higher", xray.LOW: "lower", xray.SAME: "same", xray.UNDECIDED: "undecided",
            xray.CHANGED: "changed"}


def _merge(x, y, yerr, xerr, n_bins: int = 14):
    """Merge points into logarithmic energy bins (inverse-variance mean). For display only."""
    x, y, yerr, xerr = (np.asarray(v, dtype=float) for v in (x, y, yerr, xerr))
    ok = np.isfinite(y) & np.isfinite(yerr) & (yerr > 0)
    x, y, yerr, xerr = x[ok], y[ok], yerr[ok], xerr[ok]
    if x.size <= n_bins:
        return x, y, yerr, xerr
    edges = np.geomspace((x - xerr).min(), (x + xerr).max() * 1.0001, n_bins + 1)
    idx = np.clip(np.digitize(x, edges) - 1, 0, n_bins - 1)
    out = []
    for i in range(n_bins):
        m = idx == i
        if not m.any():
            continue
        w = 1.0 / yerr[m] ** 2
        lo, hi = (x[m] - xerr[m]).min(), (x[m] + xerr[m]).max()
        centre = float(np.sqrt(lo * hi))
        out.append((centre, float(np.sum(w * y[m]) / w.sum()), float(1.0 / np.sqrt(w.sum())), 0.5 * (hi - lo)))
    return tuple(np.array(v) for v in zip(*out))


def _points(ax, x, y, yerr, xerr, key, label=True, floor=None):
    name, colour, marker, filled = SERIES[key]
    x, y, yerr, xerr = (np.asarray(v, dtype=float) for v in (x, y, yerr, xerr))
    ok = np.isfinite(y) & (y > (floor if floor is not None else -np.inf))
    ax.errorbar(x[ok], y[ok], yerr=yerr[ok], xerr=xerr[ok], fmt=marker, color=colour, markersize=4.5,
                markerfacecolor=colour if filled else "white", markeredgecolor=colour, markeredgewidth=1.0,
                elinewidth=0.8, linewidth=0, label=name if label else None, zorder=3)


def _spectrum_panels(ax, axr, arrays: dict, bat: dict, model: str, show_legend: bool) -> None:
    for key in ("xrt10", "xrt22", "fpma", "fpmb"):
        a = arrays[key]
        _points(ax, *_merge(a["energy"], a["eeuf"], a["eeuf_err"], a["energy_err"]), key, floor=0)
        order = np.argsort(a["energy"])
        colour = SERIES[key][1]
        if key in ("xrt10", "fpma", "xrt22"):
            ax.plot(np.array(a["energy"])[order], np.array(a["eeuf_model"])[order], color=colour, linewidth=1.4,
                    alpha=0.9, zorder=2)
        _points(axr, *_merge(a["energy"], a["ratio"], a["ratio_err"], a["energy_err"]), key, label=False)
    _points(ax, bat["energy"], bat["eeuf"], bat["eeuf_err"], bat["energy_err"], "bat", floor=0)
    order = np.argsort(bat["energy"])
    ax.plot(np.array(bat["energy"])[order], np.array(bat["eeuf_model"])[order], color=CBAT, linewidth=1.4,
            linestyle="--", zorder=2)
    _points(axr, bat["energy"], bat["ratio"], bat["ratio_err"], bat["energy_err"], "bat", label=False)
    ax.set_xscale("log")
    ax.set_yscale("log")
    axr.set_xscale("log")
    ax.set_xlim(0.3, 200)
    axr.set_xlim(0.3, 200)
    ax.set_ylim(2e-5, 4e-2)
    axr.set_ylim(0, 4.0)
    axr.axhline(1.0, color=INK, linewidth=0.9)
    ax.set_title(MODEL_TITLE[model], fontsize=10, loc="left", color=INK)
    axr.set_xlabel("observed energy (keV)")
    ax.tick_params(labelbottom=False)
    for a_ in (ax, axr):
        a_.set_xticks([0.5, 1, 2, 5, 10, 20, 50, 100])
        a_.set_xticklabels(["0.5", "1", "2", "5", "10", "20", "50", "100"])
    if show_legend:
        ax.legend(loc="upper left", fontsize=7.5, frameon=False, ncol=1, handletextpad=0.3)


def _interval_panel(ax, summary: dict, key: str, title: str, lim_min: float, lim_same: float, xlim) -> None:
    models = ["M0", "M1", "M2"]
    ax.axvspan(-lim_same, lim_same, color=GRID, alpha=0.6, zorder=0)
    ax.axvspan(-lim_min, lim_min, color=NEUTRAL, alpha=0.25, zorder=0)
    ax.axvline(0.0, color=INK, linewidth=0.9)
    rows = []
    for m in models:
        c = summary["by_model"][m]
        if key == "delta_nh":
            rows.append((f"{m}", c[key]["value"], c[key]["stat_90"], c[key]["label"], C2010))
        else:
            rows.append((f"{m}", c[key]["value"], c[key]["with_systematic_90"], c[key]["label"],
                         CBAT if key == "delta1" else C2010))
    for i, (name, v, (lo, hi), label, colour) in enumerate(rows):
        y = len(rows) - 1 - i
        lo, hi = (float("-inf") if lo is None else lo), (float("inf") if hi is None else hi)
        registered = " (registered)" if name == summary["chosen_model"] else ""
        if not np.isfinite(lo) and not np.isfinite(hi):
            ax.plot(xlim, [y, y], color=colour, linewidth=1.2, linestyle=":")
            ax.text(xlim[1], y + 0.28, f"{name}{registered}: not bounded  {LABEL_EN[label]}", ha="right",
                    va="bottom", fontsize=8, color=INK)
            continue
        ax.plot([max(lo, xlim[0]), min(hi, xlim[1])], [y, y], color=colour, linewidth=2.2, solid_capstyle="butt")
        ax.plot([v], [y], marker="o", color=colour, markersize=6)
        ax.text(xlim[1], y + 0.28, f"{name}{registered}: {v:+.2f} [{lo:+.2f}, {hi:+.2f}]  {LABEL_EN[label]}",
                ha="right", va="bottom", fontsize=8, color=INK)
    ax.set_ylim(-0.6, len(rows) - 0.1)
    ax.set_xlim(*xlim)
    ax.set_yticks([])
    ax.grid(axis="y", visible=False)
    ax.set_title(title, fontsize=9.5, loc="left", color=INK)


def plot_xray(summary_json: Path, plot_json: Path, path: Path) -> Path:
    summary = json.loads(Path(summary_json).read_text())
    plots = json.loads(Path(plot_json).read_text())
    chosen = summary["chosen_model"]
    other = "M2" if chosen != "M2" else "M0"
    fig = plt.figure(figsize=(15.5, 7.6))
    gs = fig.add_gridspec(3, 3, height_ratios=[2.4, 1.0, 0.02], width_ratios=[1.0, 1.0, 0.95], hspace=0.08,
                          wspace=0.2, left=0.055, right=0.985, top=0.87, bottom=0.09)
    for col, model in enumerate((chosen, other)):
        ax = fig.add_subplot(gs[0, col])
        axr = fig.add_subplot(gs[1, col])
        _spectrum_panels(ax, axr, plots[model]["c2"], plots[model]["bat"], model, show_legend=(col == 0))
        if col == 0:
            ax.set_ylabel(r"unfolded $E\,F_E$ (keV$^2$ photons cm$^{-2}$ s$^{-1}$ keV$^{-1}$)")
            axr.set_ylabel("data / model")
        tag = "registered model" if model == chosen else "lowest AIC" if model == summary["aic_lowest_model"] else ""
        ax.set_title(f"{MODEL_TITLE[model]} ({tag})" if tag else MODEL_TITLE[model], fontsize=10, loc="left", color=INK)
    sub = gs[0:2, 2].subgridspec(3, 1, hspace=0.75)
    ax1, ax2, ax3 = (fig.add_subplot(sub[i]) for i in range(3))
    _interval_panel(ax1, summary, "delta1", r"$\Delta_1$ = log L(2022) / L(BAT average)", xray.DELTA_L_MIN,
                    xray.DELTA_L_SAME_BOUND, (-0.6, 1.3))
    _interval_panel(ax2, summary, "delta2", r"$\Delta_2$ = log L$_{int}$(2022) / L$_{int}$(2010)", xray.DELTA_L_MIN,
                    xray.DELTA_L_SAME_BOUND, (-0.6, 1.3))
    _interval_panel(ax3, summary, "delta_nh", r"$\Delta_N$ = log N$_H$(2022) / N$_H$(2010)", xray.DELTA_NH_MIN,
                    xray.DELTA_NH_SAME_BOUND, (-0.6, 1.3))
    ax3.set_xlabel("dex (90% interval; dark band: 'same' size limit, light band: 'same' interval limit)", fontsize=8)
    fig.suptitle("CGCG 312-012: X-ray spectra in 2010, 2022 and the Swift-BAT average, with the registered comparisons",
                 x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.text(0.01, 0.925, "Unfolded spectra depend on the fitted model. Lines: model of each epoch (BAT: the 2022 model "
             "scaled by the fitted constant). Points are merged into wider energy bins for display only.", fontsize=8.5, color=INK_2)
    return _save(fig, path)


def plot_supplement(summary_json: Path, fits_json: Path, plot_json: Path, path: Path) -> Path:
    """Stage 7b (post hoc): hard-band flux ratio, Delta2 and DeltaN under the readings of the
    XRT constant, and the XRT residuals of M2 with the constant fixed and free."""
    s = json.loads(Path(summary_json).read_text())
    fits = json.loads(Path(fits_json).read_text())
    plots = json.loads(Path(plot_json).read_text())
    fig, axes = plt.subplots(1, 4, figsize=(16.5, 4.6), gridspec_kw={"width_ratios": [1.25, 1.0, 1.0, 1.0], "wspace": 0.2, "top": 0.8, "left": 0.02, "right": 0.985, "bottom": 0.13})

    def rows_panel(ax, rows, title, xlim, colour, bands=None):
        if bands:
            ax.axvspan(-bands[1], bands[1], color=GRID, alpha=0.6, zorder=0)
            ax.axvspan(-bands[0], bands[0], color=NEUTRAL, alpha=0.25, zorder=0)
        ax.axvline(0.0, color=INK, linewidth=0.9)
        for i, (name, v, lo, hi) in enumerate(rows):
            y = len(rows) - 1 - i
            ax.plot([lo, hi], [y, y], color=colour, linewidth=2.2, solid_capstyle="butt")
            ax.plot([v], [y], marker="o", color=colour, markersize=6)
            ax.text(xlim[1], y + 0.25, f"{name}: {v:+.2f} [{lo:+.2f}, {hi:+.2f}]", ha="right", va="bottom", fontsize=7.6, color=INK)
        ax.set_ylim(-0.6, len(rows) - 0.05)
        ax.set_xlim(*xlim)
        ax.set_yticks([])
        ax.grid(axis="y", visible=False)
        ax.set_title(title, fontsize=9.5, loc="left", color=INK)
        ax.set_xlabel("dex (90% interval)", fontsize=8)

    f = s["flux_ratio_14_50"]
    c1 = s["c1_shape_tied"]
    rows = [("BAT channels 1-4, free index", f["bat_channels_1_4_free_gamma"]["log_ratio"], f["bat_channels_1_4_free_gamma"]["lo"], f["bat_channels_1_4_free_gamma"]["hi"]),
            ("BAT channels 1-4, index fixed to NuSTAR", f["bat_channels_1_4_gamma_fixed_to_nustar"]["log_ratio"], f["bat_channels_1_4_gamma_fixed_to_nustar"]["lo"], f["bat_channels_1_4_gamma_fixed_to_nustar"]["hi"]),
            ("stage 7: BAT all channels, integrated 14-50 keV", f["stage7_value_renamed"]["log_ratio"], f["stage7_value_renamed"]["lo"], f["stage7_value_renamed"]["hi"]),
            ("C1, shape tied to 2022 model M0", c1["M0"]["delta1"], *c1["M0"]["stat_90"]),
            ("C1, shape tied to 2022 model M2", c1["M2"]["delta1"], *c1["M2"]["stat_90"])]
    rows_panel(axes[0], rows, "NuSTAR 2022 / BAT average,\nlog ratio of the 14-50 keV flux", (-0.1, 1.25), CBAT)
    v = s["m0_delta2_delta_nh"]
    names = (("stage 7 (XRT constant 1)", "stage7_xrt_constant_1"), ("(a) one constant on both XRT", "a_constant_on_both_xrt"),
             ("(b) constant on 2022 XRT only", "b_constant_on_2022_xrt_only"))
    rows_panel(axes[1], [(n, v[k]["delta2"], *v[k]["delta2_90"]) for n, k in names],
               r"M0: $\Delta_2$ =" "\n" r"log L$_{int}$(2022) / L$_{int}$(2010)", (-0.45, 1.35), C2010, (xray.DELTA_L_MIN, xray.DELTA_L_SAME_BOUND))
    rows_panel(axes[2], [(n, v[k]["delta_nh"], *v[k]["delta_nh_90"]) for n, k in names],
               r"M0: $\Delta_N$ =" "\n" r"log N$_H$(2022) / N$_H$(2010)", (-0.45, 1.35), C2010, (xray.DELTA_NH_MIN, xray.DELTA_NH_SAME_BOUND))
    ax = axes[3]
    fixed = plots["M2"]["c2"]["xrt22"]
    free = fits["m2_xrt_free"]["plot"]["xrt22"]
    for arr, marker, filled, label in ((fixed, "o", True, "XRT constant fixed at 1"), (free, "s", False, "XRT constant free (0.65)")):
        x, y, ye, xe = (np.asarray(arr[k], dtype=float) for k in ("energy", "ratio", "ratio_err", "energy_err"))
        ax.errorbar(x, y, yerr=ye, xerr=xe, fmt=marker, color=C2022, markersize=5, markerfacecolor=C2022 if filled else "white",
                    markeredgecolor=C2022, elinewidth=0.9, linewidth=0, label=label)
    ax.axhline(1.0, color=INK, linewidth=0.9)
    ax.set_xscale("log")
    ax.set_xlim(1.5, 10)
    ax.set_xticks([2, 3, 5, 7, 10])
    ax.set_xticklabels(["2", "3", "5", "7", "10"])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_ylim(0, 2.6)
    ax.set_xlabel("observed energy (keV)", fontsize=8)
    ax.set_title("2022 XRT, data / model M2\n(bins merged for display)", fontsize=9.5, loc="left", color=INK)
    ax.legend(loc="upper left", fontsize=7.5, frameon=False)
    fig.suptitle("CGCG 312-012, stage 7b (post hoc supplement): what the X-ray comparison depends on",
                 x=0.01, ha="left", fontweight="bold", fontsize=11)
    return _save(fig, path)
