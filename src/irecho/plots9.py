"""Stage 9 figures: identifiability on synthetic samples, example objects with the
fitted kernel, predicted against measured branch offset. Same visual system as plots.py.

Colour follows the band (W1, W2), the fitted model, or the sample; text stays in ink.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from .lightcurve import mjd_to_year  # noqa: E402
from .plots import GRID, INK, INK_2, NEUTRAL, _save  # noqa: E402

BAND_W1, BAND_W2, MODEL = "#2a78d6", "#eb6834", "#1baf7a"
SAMPLE_COLOUR = {"v1": "#2a78d6", "v2": "#eb6834"}
SAMPLE_LABEL = {"v1": "sample_v1 (z <= 0.04)", "v2": "sample_v2 (0.04 < z <= 0.30)"}
DRIVER_LABEL = {"fast": "driver 200 d", "slow": "driver 1000 d"}


def _title(fig, text: str) -> None:
    fig.suptitle(text, x=0.01, ha="left", fontweight="bold", fontsize=11)


def plot_identifiability(ident: dict, observed: dict | None, path: Path) -> Path:
    """Recovered lag and branch offset against the lag put in; residual after a fixed kernel by truth model."""
    fig, axes = plt.subplots(1, 3, figsize=(17.5, 5.6), gridspec_kw={"wspace": 0.3, "width_ratios": [1, 1, 1.0], "top": 0.8, "right": 0.985, "left": 0.05})
    pos = axes[2].get_position()
    axes[2].set_position([pos.x0 + 0.075, pos.y0, pos.width - 0.075, pos.height])
    delays = np.array([0.0, 20.0, 40.0, 80.0, 160.0])
    names = ("T0", "T1_20", "T1_40", "T1_80", "T1_160")
    for tag in ("v1", "v2"):
        for driver, marker, filled, dx in (("fast", "o", True, -2.0), ("slow", "s", False, 2.0)):
            truths = ident["configs"][f"{tag}_{driver}"]["truths"]
            for ax, key, scale in ((axes[0], "median_tau", 1.0), (axes[1], "median_d_colour_obs", 1e3)):
                med = np.array([truths[n][key]["median"] for n in names]) * scale
                rng = np.array([truths[n][key]["null_range"] for n in names]) * scale
                ax.errorbar(delays + dx + (1.0 if tag == "v2" else -1.0), med, yerr=[med - rng[:, 0], rng[:, 1] - med], fmt=marker,
                            color=SAMPLE_COLOUR[tag], markerfacecolor=SAMPLE_COLOUR[tag] if filled else "white", markersize=6,
                            elinewidth=1.0, linewidth=0, label=f"{SAMPLE_LABEL[tag]}, {DRIVER_LABEL[driver]}")
    ax = axes[0]
    ax.plot([0, 160], [0, 160], color=INK, linewidth=0.8)
    ax.set_xlabel("delay of W2 behind W1 put in (days, observed frame)")
    ax.set_ylabel("sample median of the fitted lag (days)")
    ax.set_title("B1: the lag is recovered at sample level\n(points: median of 50 synthetic samples; bars: null range)", fontsize=9.5, loc="left", color=INK)
    ax.legend(fontsize=7.5, frameon=False, loc="upper left")
    ax = axes[1]
    ax.axhline(0.0, color=INK, linewidth=0.8)
    if observed:
        for tag in ("v1", "v2"):
            y = 1e3 * observed[tag]
            ax.axhline(y, color=SAMPLE_COLOUR[tag], linewidth=1.6, linestyle="--")
            ax.text(162, y, f"observed\n{tag}: {y:+.1f}", fontsize=8, color=INK, va="center", ha="left")
    ax.set_xlim(-10, 195)
    ax.set_xlabel("delay of W2 behind W1 put in (days, observed frame)")
    ax.set_ylabel("sample median of the branch offset (mmag)\nnegative = bluer when rising")
    ax.set_title("B1: branch offset produced by a pure delay\n(visits used by the fit)", fontsize=9.5, loc="left", color=INK)
    ax = axes[2]
    show = ("T0", "T1_40", "T1_lyu", "T1_160", "T2_80", "T0b", "T3_1yr", "T3_3yr")
    y0 = np.arange(len(show))[::-1].astype(float)
    k = 0
    for tag in ("v1", "v2"):
        for driver, marker, filled in (("fast", "o", True), ("slow", "s", False)):
            truths = ident["configs"][f"{tag}_{driver}"]["truths"]
            med = np.array([truths[n]["median_d_colour_res_t2"]["median"] for n in show]) * 1e3
            rng = np.array([truths[n]["median_d_colour_res_t2"]["null_range"] for n in show]) * 1e3
            y = y0 + 0.27 - 0.18 * k
            ax.errorbar(med, y, xerr=[med - rng[:, 0], rng[:, 1] - med], fmt=marker, color=SAMPLE_COLOUR[tag],
                        markerfacecolor=SAMPLE_COLOUR[tag] if filled else "white", markersize=5, elinewidth=1.0, linewidth=0)
            k += 1
    ax.axvline(0.0, color=INK, linewidth=0.8)
    ax.set_yticks(y0)
    ax.set_yticklabels(["same response (T0)", "shift 40 d", "shift from Lyu+19", "shift 160 d", "tail, mean 80 d", "moving edge,\nno re-formation delay",
                        "re-formation 1 yr,\nno lag", "re-formation 3 yr,\nno lag"], fontsize=8.5)
    ax.set_xlabel("sample median of the residual offset after the fitted kernel T2 (mmag)")
    ax.set_title("B1: what a fixed kernel leaves over, by truth model\n(re-formation without a lag is absorbed by the fitted lag)", fontsize=9.5, loc="left", color=INK)
    for a in axes:
        a.grid(color=GRID, linewidth=0.6)
        a.set_axisbelow(True)
    fig.suptitle("Stage 9, B1: identifiability with the real visit times and errors (synthetic light curves)", x=0.01, y=0.97, ha="left", fontweight="bold", fontsize=11)
    return _save(fig, path)


def plot_examples(objects: list[dict], path: Path) -> Path:
    """W1 and W2 light curves with the fitted model of W2 and the kernel.

    objects: dicts with label, mjd, f1, s1, f2, s2 (Jy), idx (visits used), pred_t0, pred_t1, pred_t2 (Jy at idx),
    tau_t1, tau_lo, tau_hi, tau_t2, width_t2, ds_01, d_colour_obs, d_colour_t2.
    """
    n = len(objects)
    fig, axes = plt.subplots(n, 3, figsize=(17.5, 4.1 * n), gridspec_kw={"width_ratios": [1.35, 1.0, 0.75], "wspace": 0.3, "hspace": 0.55, "top": 0.9, "left": 0.05, "right": 0.985})
    axes = np.atleast_2d(axes)
    for i, o in enumerate(objects):
        yr = mjd_to_year(o["mjd"])
        idx = o["idx"]
        ax = axes[i, 0]
        z1 = o["f1"] / np.median(o["f1"])
        z2 = o["f2"] / np.median(o["f2"])
        ax.errorbar(yr, z1, yerr=o["s1"] / np.median(o["f1"]), fmt="o", color=BAND_W1, markerfacecolor="white", markersize=5,
                    elinewidth=0.8, linewidth=0.8, label="W1 (input)")
        ax.errorbar(yr, z2, yerr=o["s2"] / np.median(o["f2"]), fmt="s", color=BAND_W2, markersize=5, elinewidth=0.8, linewidth=0,
                    label="W2 (observed)")
        ax.plot(yr[idx], o["pred_t2"] / np.median(o["f2"]), color=MODEL, linewidth=1.8, label="W2 from W1, fitted kernel (T2)", zorder=4)
        ax.plot(yr[idx], o["pred_t0"] / np.median(o["f2"]), color=NEUTRAL, linewidth=1.2, linestyle="--", label="W2 from W1, no lag (T0)")
        ax.axvspan(yr[0] - 0.4, yr[idx[0]] - 0.2, color=GRID, alpha=0.6, zorder=0)
        ax.axvspan(yr[idx[-1]] + 0.2, yr[-1] + 0.5, color=GRID, alpha=0.6, zorder=0)
        ax.set_xlim(yr[0] - 0.4, yr[-1] + 0.5)
        ax.set_xlabel("year (shaded: visits not predicted, the kernel would need extrapolation)")
        ax.set_ylabel("flux / median")
        ax.set_title(o["label"], fontsize=9.5, loc="left", color=INK)
        if i == 0:
            ax.legend(fontsize=7.5, frameon=False, loc="best", ncol=2)
        ax = axes[i, 1]
        res0 = (o["f2"][idx] - o["pred_t0"]) / o["s2"][idx]
        res2 = (o["f2"][idx] - o["pred_t2"]) / o["s2"][idx]
        ax.axhline(0.0, color=INK, linewidth=0.8)
        ax.plot(yr[idx], res0, "o--", color=NEUTRAL, markersize=4, linewidth=1.0, label="no lag (T0)")
        ax.plot(yr[idx], res2, "o-", color=MODEL, markersize=4, linewidth=1.4, label="fitted kernel (T2)")
        ax.set_xlabel("year")
        ax.set_ylabel("(W2 observed - model) / W2 error")
        ax.set_title(f"residuals; T0 to T1 improves S by {o['ds_01']:.1f}\noffset: measured {1e3 * o['d_colour_obs']:+.1f}, "
                     f"kernel {1e3 * o['d_colour_t2']:+.1f} mmag", fontsize=9.5, loc="left", color=INK)
        if i == 0:
            ax.legend(fontsize=7.5, frameon=False, loc="best")
        ax = axes[i, 2]
        tau, w = o["tau_t2"], o["width_t2"]
        ax.axvline(0.0, color=INK, linewidth=0.8)
        if w > 0:
            ax.fill_between([tau - w / 2, tau - w / 2, tau + w / 2, tau + w / 2], [0, 1 / w, 1 / w, 0], color=MODEL, alpha=0.35, linewidth=0)
            ax.plot([tau - w / 2, tau - w / 2, tau + w / 2, tau + w / 2], [0, 1 / w, 1 / w, 0], color=MODEL, linewidth=1.6)
            top = 1 / w
        else:
            top = 0.012
            ax.plot([tau, tau], [0, top], color=MODEL, linewidth=2.2)
        ax.plot([o["tau_lo"], o["tau_hi"]], [1.15 * top, 1.15 * top], color=INK_2, linewidth=2.0)
        ax.plot([o["tau_t1"]], [1.15 * top], "o", color=INK_2, markersize=5)
        ax.set_xlim(-170, 370)
        ax.set_ylim(0, 1.45 * top)
        ax.set_yticks([])
        ax.set_xlabel("lag of W2 behind W1 (days, observed frame)")
        ax.set_title(f"T2: lag {tau:+.0f} d, width {w:.0f} d\nT1: {o['tau_t1']:+.0f} d [{o['tau_lo']:+.0f}, {o['tau_hi']:+.0f}] (grey)",
                     fontsize=9.5, loc="left", color=INK)
        for a in axes[i]:
            a.grid(color=GRID, linewidth=0.6)
            a.set_axisbelow(True)
    fig.suptitle("Stage 9, B2: W2 predicted from W1 with a fixed kernel (linear interpolation of W1, no extrapolation)", x=0.01, y=0.975, ha="left", fontweight="bold", fontsize=11)
    return _save(fig, path)


def plot_offsets(rows: dict, verdict: dict, lagl: dict, path: Path) -> Path:
    """Predicted against measured branch offset, the residual against its null, and lag against luminosity.

    rows[tag]: list of dicts with d_colour_obs, d_colour_t2, sigma_colour_obs (comparable objects).
    lagl: {"log_lbol": [...], "tau_rest": [...], "well": [...], "curve_x": [...], "lyu": [...], "mandal": [...]}.
    """
    fig, grid_axes = plt.subplots(2, 2, figsize=(13.5, 12.0), gridspec_kw={"wspace": 0.28, "hspace": 0.36, "top": 0.9, "bottom": 0.07,
                                                                          "left": 0.08, "right": 0.98})
    axes = grid_axes.ravel()
    for ax, tag in zip(axes[:2], ("v1", "v2")):
        obs = 1e3 * np.array([r["d_colour_obs"] for r in rows[tag]])
        pred = 1e3 * np.array([r["d_colour_t2"] for r in rows[tag]])
        lim = 160.0
        ax.plot([-lim, lim], [-lim, lim], color=INK, linewidth=0.8)
        ax.axhline(0.0, color=GRID, linewidth=0.8)
        ax.axvline(0.0, color=GRID, linewidth=0.8)
        ax.scatter(np.clip(obs, -lim, lim), np.clip(pred, -lim, lim), s=14, color=SAMPLE_COLOUR[tag], alpha=0.75, linewidth=0)
        e = verdict["samples"][tag]["explained"]["t2"]
        ax.plot([1e3 * e["m_obs_mag"]], [1e3 * e["m_pred_mag"]], marker="+", color=INK, markersize=16, markeredgewidth=2.0)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_xlabel("measured branch offset (mmag); negative = bluer when rising")
        ax.set_ylabel("offset given by the fitted kernel T2 (mmag)")
        ax.set_title(f"{SAMPLE_LABEL[tag]}: {len(obs)} objects\nmedians (cross): measured {1e3 * e['m_obs_mag']:+.1f}, kernel {1e3 * e['m_pred_mag']:+.1f},\n"
                     f"residual {1e3 * e['m_res_mag']:+.1f} mmag", fontsize=9.5, loc="left", color=INK)
    ax = axes[2]
    y = 0
    ticks, names = [], []
    for tag in ("v2", "v1"):
        for k, name in (("t2", "T2"), ("t1", "T1")):
            e = verdict["samples"][tag]["explained"][k]
            for j, (driver, b) in enumerate(e["null_fitted"].items()):
                ax.plot(1e3 * np.array(b["m_res"]), [y + 0.12 - 0.24 * j] * 2, color=NEUTRAL, linewidth=5, solid_capstyle="butt",
                        label="null range, truth = pure shift by the fitted lags" if y == 0 and j == 0 else None)
            ax.plot([1e3 * e["m_res_mag"]], [y], "o", color=SAMPLE_COLOUR[tag], markersize=8, label="observed" if y == 0 else None)
            ticks.append(y)
            names.append(f"{tag}, {name}")
            y += 1
    ax.axvline(0.0, color=INK, linewidth=0.8)
    ax.set_yticks(ticks)
    ax.set_yticklabels(names, fontsize=9)
    ax.set_ylim(-0.6, y - 0.4)
    ax.set_xlabel("sample median of the residual offset (mmag)")
    ax.set_title("residual after the kernel, against synthetic samples\n(upper bar: driver 200 d, lower bar: driver 1000 d)", fontsize=9.5, loc="left", color=INK)
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    ax.set_ylim(-0.6, y + 0.5)
    ax = axes[3]
    lum, tau, well = (np.array(lagl[k]) for k in ("log_lbol", "tau_rest", "well"))
    ax.axhline(0.0, color=INK, linewidth=0.8)
    ax.scatter(lum[~well], tau[~well], s=12, color=NEUTRAL, alpha=0.7, linewidth=0, label="fitted lag (T1), rest frame")
    ax.scatter(lum[well], tau[well], s=18, color=MODEL, alpha=0.9, linewidth=0, label="'well determined' by the plan's rule")
    ax.plot(lagl["curve_x"], lagl["lyu"], color=INK, linewidth=1.6, label="Lyu+19 Eqs. 15-16: W2 lag - W1 lag")
    ax.plot(lagl["curve_x"], lagl["mandal"], color=INK, linewidth=1.2, linestyle="--", label="Mandal+24 Table 5 (post hoc)")
    ax.set_ylim(-170, 370)
    ax.set_xlabel("log L_bol (erg/s), BASS DR2")
    ax.set_ylabel("lag of W2 behind W1 (days, rest frame)")
    ax.set_title(f"B3: lag against luminosity, both samples ({len(lum)} objects)", fontsize=9.5, loc="left", color=INK)
    ax.legend(fontsize=7.5, frameon=False, loc="upper left")
    for a in axes:
        a.grid(color=GRID, linewidth=0.6)
        a.set_axisbelow(True)
    fig.suptitle("Stage 9, B2-B3: does a fixed kernel W1 -> W2 give the measured colour asymmetry?", x=0.01, y=0.975, ha="left", fontweight="bold", fontsize=11)
    return _save(fig, path)
