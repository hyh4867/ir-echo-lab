"""Stage 8 figures: null test of the hysteresis statistic, the direct branch
comparison, and two example objects. Same visual system as plots.py.

Colour follows the branch (rising, declining) or the sample; text stays in ink.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from . import hyst8  # noqa: E402
from .lightcurve import mjd_to_year  # noqa: E402
from .plots import GRID, INK, INK_2, NEUTRAL, _save  # noqa: E402

RISE, FADE, OBS = "#eb6834", "#2a78d6", "#1baf7a"
SAMPLE_LABEL = {"v1": "sample_v1 (z <= 0.04)", "v2": "sample_v2 (0.04 < z <= 0.30)"}
NULL_LABEL = {"linear": "null (a): W1 linear in W2", "single_t": "null (b): one temperature, no lag"}


def _title(fig, text: str) -> None:
    fig.suptitle(text, x=0.01, ha="left", fontweight="bold", fontsize=11)


def plot_null(null_json: Path, cache_dir: Path, path: Path) -> Path:
    """Null distributions of the sample median of D against the observed medians."""
    null = json.loads(Path(null_json).read_text())
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 7.2), gridspec_kw={"width_ratios": [1, 1, 1.25], "wspace": 0.25, "hspace": 0.45})
    for i, tag in enumerate(("v1", "v2")):
        rec = null["samples"][tag]
        for j, kind in enumerate(("linear", "single_t")):
            ax = axes[i, j]
            rows = json.loads((Path(cache_dir) / f"null_{tag}_{kind}_primary.json").read_text())
            med = np.array([r["median_d"] for r in rows if r["median_d"] is not None])
            ax.hist(med, bins=np.linspace(-0.05, 0.05, 26), color=NEUTRAL, alpha=0.8, edgecolor="white", linewidth=0.5)
            s = rec["nulls"][kind]["primary"]
            obs = s["observed_median_d"]
            ax.axvline(obs, color=OBS, linewidth=2.0)
            ax.axvline(0.0, color=INK, linewidth=0.8)
            ax.set_xlim(-0.06, 0.32)
            ax.text(obs - 0.005, ax.get_ylim()[1] * 0.92, f"observed\nmedian D = +{obs:.3f}", ha="right", va="top", fontsize=8.5, color=INK)
            q = s["median_d_p2p5_p16_p50_p84_p97p5_p99"]
            ax.text(0.03, ax.get_ylim()[1] * 0.55, f"null: {q[2]:+.3f}\n[{q[0]:+.3f}, {q[4]:+.3f}]\n{s['n_realisations']} realisations",
                    ha="left", va="center", fontsize=8.5, color=INK_2)
            ax.set_title(f"{SAMPLE_LABEL[tag]}\n{NULL_LABEL[kind]}", fontsize=9.5, loc="left", color=INK)
            ax.set_xlabel("sample median of D (measurable objects)")
            if j == 0:
                ax.set_ylabel("realisations")
        ax = axes[i, 2]
        names = list(null["settings"]["variants"])
        for k, (kind, marker, filled) in enumerate((("linear", "o", True), ("single_t", "s", False))):
            y = np.arange(len(names))[::-1] + (0.14 if k == 0 else -0.14)
            vals = [rec["nulls"][kind][v]["median_d_p2p5_p16_p50_p84_p97p5_p99"] for v in names]
            ax.errorbar([v[2] for v in vals], y, xerr=[[v[2] - v[0] for v in vals], [v[4] - v[2] for v in vals]], fmt=marker,
                        color=INK_2, markerfacecolor=INK_2 if filled else "white", markersize=5, elinewidth=1.0, linewidth=0,
                        label=NULL_LABEL[kind])
        ax.axvline(0.0, color=INK, linewidth=0.8)
        ax.axvspan(-0.03, 0.03, color=GRID, alpha=0.7, zorder=0)
        ax.set_yticks(np.arange(len(names))[::-1])
        ax.set_yticklabels([n.replace("_", " ") for n in names], fontsize=8)
        ax.set_xlim(-0.06, 0.06)
        ax.grid(axis="y", visible=False)
        ax.set_xlabel("null median of D (2.5-97.5%); band: 'no bias' limit +-0.03")
        ax.set_title(f"{SAMPLE_LABEL[tag]}: estimator variants", fontsize=9.5, loc="left", color=INK)
        if i == 0:
            ax.legend(loc="lower left", bbox_to_anchor=(0.0, -0.36), fontsize=7.5, frameon=False, ncol=2)
    _title(fig, "Stage 8, B1: the hysteresis statistic D on light curves without a colour loop (real visit times, errors and W2 light curves)")
    return _save(fig, path)


def _rows(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def plot_direct(objects_csv: Path, direct_json: Path, null_json: Path, path: Path) -> Path:
    """Direct comparison: colour offset between the branches at the same W2 flux."""
    rows = _rows(objects_csv)
    direct = json.loads(Path(direct_json).read_text())
    null = json.loads(Path(null_json).read_text())
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 7.6), gridspec_kw={"wspace": 0.27, "hspace": 0.45})
    for i, tag in enumerate(("v1", "v2")):
        sub = [r for r in rows if r["sample"] == tag and r["comparable"] == "1"]
        col = np.array([_f(r["d_colour"]) for r in sub])
        sig = np.array([_f(r["sigma_colour"]) for r in sub])
        st = direct["samples"][tag]["direct"]
        ax = axes[i, 0]
        order = np.argsort(col)
        ax.errorbar(np.arange(col.size), 1e3 * col[order], yerr=1e3 * sig[order], fmt="o", color=OBS, markersize=3, elinewidth=0.6, linewidth=0)
        ax.axhline(0.0, color=INK, linewidth=0.9)
        ax.axhline(1e3 * st["median_d_colour_mag"], color=OBS, linewidth=1.2, linestyle="--")
        ax.set_xlabel("objects, sorted")
        ax.set_ylabel(r"$\Delta$(W1$-$W2), rise $-$ fade (mmag)")
        ax.set_title(f"{SAMPLE_LABEL[tag]}: {st['n']} comparable objects\n"
                     f"negative (rise bluer) in {st['n_delta_positive']} ({100 * st['frac_delta_positive']:.0f}%), "
                     f"median {1e3 * st['median_d_colour_mag']:+.1f} mmag", fontsize=9.5, loc="left", color=INK)
        lim = np.percentile(np.abs(1e3 * col), 97) * 1.5
        ax.set_ylim(-lim, lim)
        # sign fraction against the null
        ax = axes[i, 1]
        for kind, colour_alpha in (("linear", 0.9), ("single_t", 0.45)):
            q = null["samples"][tag]["nulls"][kind]["primary"]["direct"]["frac_delta_positive_p2p5_p50_p97p5_p99"]
            y = 1 if kind == "linear" else 0
            ax.plot([q[0], q[2]], [y, y], color=NEUTRAL, linewidth=6, solid_capstyle="butt", alpha=colour_alpha)
            ax.plot([q[1]], [y], marker="|", color=INK, markersize=12)
            short = "null (a)" if kind == "linear" else "null (b)"
            ax.text(q[0] - 0.015, y, f"{short}: {q[1]:.2f} [{q[0]:.2f}, {q[2]:.2f}]", ha="right", va="center", fontsize=8, color=INK_2)
        entries = [("pipeline photometry", direct["samples"][tag]["direct"], 3.2),
                   ("with a scan-parity term", direct["samples"][tag]["direct_with_scan_parity_term"], 2.6),
                   ("unTimely photometry", direct["samples"][tag]["untimely"], 2.0)]
        for name, e, y in entries:
            if e.get("n", 0) < 8:
                continue
            ax.plot([e["frac_delta_positive"]], [y], marker="o", color=OBS, markersize=7)
            ax.text(e["frac_delta_positive"] - 0.015, y, f"{name}: {e['frac_delta_positive']:.2f} (n = {e['n']})", ha="right", va="center", fontsize=8, color=INK)
        ax.axvline(0.5, color=INK, linewidth=0.8)
        ax.set_xlim(0.2, 1.02)
        ax.set_ylim(-0.6, 3.8)
        ax.set_yticks([])
        ax.grid(axis="y", visible=False)
        ax.set_xlabel("fraction of objects with the rising branch bluer")
        ax.set_title(f"{SAMPLE_LABEL[tag]}: sign fraction\n(grey: 2.5-97.5% of the null realisations)", fontsize=9.5, loc="left", color=INK)
        # D against the direct colour offset
        ax = axes[i, 2]
        both = [r for r in sub if r["hyst_measurable"] == "1"]
        ax.errorbar([_f(r["d_hyst"]) for r in both], [-1e3 * _f(r["d_colour"]) for r in both],
                    yerr=[1e3 * _f(r["sigma_colour"]) for r in both], xerr=[_f(r["sigma_d"]) for r in both],
                    fmt="o", color=OBS, markersize=3.5, elinewidth=0.5, linewidth=0)
        ax.axhline(0.0, color=INK, linewidth=0.8)
        ax.axvline(0.0, color=INK, linewidth=0.8)
        ax.set_xlabel("D (stage 5b / 6B statistic)")
        ax.set_ylabel(r"$-\Delta$(W1$-$W2) (mmag); positive: rise bluer")
        wd = direct["samples"][tag]["with_d"]
        ax.set_title(f"{SAMPLE_LABEL[tag]}: D against the direct offset\n{wd['n_same_sign']} of {wd['n_both']} objects with the same sign",
                     fontsize=9.5, loc="left", color=INK)
    _title(fig, "Stage 8, B2: colour of the rising and the declining branch at the same W2 flux (no reference visits, no shared peak, no blackbody)")
    return _save(fig, path)


def plot_examples(objects: list[dict], path: Path) -> Path:
    """Light curves and the colour - W2 flux plane of example objects.

    objects: dicts with name, label, mjd, f1, s1, f2, s2 (Jy), lo, hi (common W2 range), d, d_colour.
    """
    n = len(objects)
    fig, axes = plt.subplots(n, 2, figsize=(13.5, 4.3 * n), gridspec_kw={"width_ratios": [1.25, 1.0], "wspace": 0.22, "hspace": 0.42})
    axes = np.atleast_2d(axes)
    for i, o in enumerate(objects):
        yr = mjd_to_year(o["mjd"])
        rising = hyst8.branch_rising(o["f1"], o["f2"])
        ax = axes[i, 0]
        for b, marker, name in (("1", "o", "W1"), ("2", "s", "W2")):
            f, s = o[f"f{b}"] * 1e3, o[f"s{b}"] * 1e3
            ax.plot(yr, f, color=NEUTRAL, linewidth=0.8, zorder=1)
            for lab, colour, text in ((True, RISE, "rising"), (False, FADE, "declining")):
                m = rising == lab
                ax.errorbar(yr[m], f[m], yerr=s[m], fmt=marker, color=colour, markersize=5, elinewidth=0.8, linewidth=0,
                            markerfacecolor=colour if b == "2" else "white", label=f"{name}, {text}", zorder=3)
            ax.text(yr[-1] + 0.15, f[-1], name, va="center", fontsize=9, color=INK)
        ax.set_xlabel("year")
        ax.set_ylabel("flux (mJy)")
        ax.set_xlim(yr[0] - 0.4, yr[-1] + 0.9)
        ax.set_title(f"{o['label']}\nD = {o['d']:+.2f}; direct offset {1e3 * o['d_colour']:+.1f} mmag", fontsize=9.5, loc="left", color=INK)
        if i == 0:
            ax.legend(loc="lower center", fontsize=7.5, frameon=False, ncol=2)
        ax = axes[i, 1]
        colour = -2.5 * np.log10(o["f1"] / o["f2"])
        cerr = 2.5 / np.log(10) * np.hypot(o["s1"] / o["f1"], o["s2"] / o["f2"])
        colour = colour - np.median(colour)
        x = o["f2"] * 1e3
        ax.axvspan(o["lo"] * 1e3, o["hi"] * 1e3, color=GRID, alpha=0.7, zorder=0)
        ax.plot(x, colour, color=NEUTRAL, linewidth=0.8, zorder=1)
        for k in range(len(x) - 1):
            ax.annotate("", xy=(x[k + 1], colour[k + 1]), xytext=(x[k], colour[k]),
                        arrowprops={"arrowstyle": "->", "color": NEUTRAL, "lw": 0.7, "shrinkA": 4, "shrinkB": 4}, zorder=2)
        for lab, c, text in ((True, RISE, "rising"), (False, FADE, "declining")):
            m = rising == lab
            ax.errorbar(x[m], colour[m], yerr=cerr[m], xerr=o["s2"][m] * 1e3, fmt="o", color=c, markersize=5, elinewidth=0.7,
                        linewidth=0, label=text, zorder=3)
        ax.invert_yaxis()
        ax.set_xlabel("W2 flux (mJy); band: common range of the two branches")
        ax.set_ylabel("W1 - W2 relative to the median (mag); up = bluer")
        ax.set_title("colour against W2 flux, arrows follow time", fontsize=9.5, loc="left", color=INK)
        if i == 0:
            ax.legend(loc="best", fontsize=8, frameon=False)
    _title(fig, "Stage 8, B4: example objects (largest D, and D closest to the sample median; sample_v2)")
    return _save(fig, path)
