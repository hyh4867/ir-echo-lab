"""Stage 5a QA figures. Same visual system as plots.py. None of these figures
shows a target light curve: the inputs are comparison-star statistics, field
zero-point offsets, and the scatter of the difference between two photometries."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from .lightcurve import mjd_to_year  # noqa: E402
from .plots import COLOR, INK, INK_2, LABEL, MARKER, NEUTRAL, SURFACE, _save  # noqa: E402

BANDS = ("w1", "w2")
RISE_COLOR = {0.0: "#2a78d6", 1.0: "#eb6834", 3.0: "#1baf7a"}  # categorical slots 1-3
RISE_MARKER = {0.0: "o", 1.0: "s", 3.0: "D"}
RISE_LABEL = {0.0: "instant step", 1.0: "1-yr linear rise", 3.0: "3-yr linear rise"}


def _rows(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def plot_star_chi2(stars_csv: Path, path: Path) -> Path:
    """Reduced chi-square of every comparison star before and after the
    leave-one-out field correction."""
    rows = _rows(stars_csv)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    bins = np.logspace(np.log10(0.15), np.log10(40.0), 41)
    for ax, b in zip(axes, BANDS):
        c0 = np.array([float(r[f"{b}_chi2_red_before"]) for r in rows])
        c1 = np.array([float(r[f"{b}_chi2_red_after"]) for r in rows])
        ax.hist(np.clip(c0, bins[0], bins[-1]), bins=bins, histtype="stepfilled", color=NEUTRAL,
                alpha=0.35, linewidth=0)
        ax.hist(np.clip(c0, bins[0], bins[-1]), bins=bins, histtype="step", color=NEUTRAL,
                linewidth=2, label=f"before: median {np.median(c0):.2f}")
        ax.hist(np.clip(c1, bins[0], bins[-1]), bins=bins, histtype="stepfilled", color=COLOR[b],
                alpha=0.25, linewidth=0)
        ax.hist(np.clip(c1, bins[0], bins[-1]), bins=bins, histtype="step", color=COLOR[b],
                linewidth=2, label=f"after: median {np.median(c1):.2f}")
        ax.axvline(1.0, color=INK, linestyle="--", linewidth=1.2)
        ax.set_xscale("log")
        ax.set_xticks([0.2, 0.5, 1, 2, 5, 10, 30])
        ax.set_xticklabels(["0.2", "0.5", "1", "2", "5", "10", "30"])
        ax.minorticks_off()
        ax.set_ylim(0, ax.get_ylim()[1] * 1.25)
        ax.set_xlabel("reduced chi-square of a constant model")
        ax.set_title(f"{LABEL[b]}: {len(c0)} comparison stars")
        ax.legend(loc="upper right")
    axes[0].set_ylabel("number of stars")
    fig.suptitle("sample_v1 fields: field zero-point correction (leave-one-out)",
                 x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return _save(fig, path)


def plot_offsets(offsets_csv: Path, path: Path) -> Path:
    """Zero-point offsets of all fields against time: median and 16-84 % range
    per half-year, and the single most extreme field visit in each bin."""
    rows = [r for r in _rows(offsets_csv) if r["w1_offset_mmag"] and r["w2_offset_mmag"]]
    year = np.array([float(mjd_to_year(float(r["mjd"]))) for r in rows])
    edges = np.arange(2013.75, 2025.01, 0.5)
    centre = 0.5 * (edges[:-1] + edges[1:])
    which = np.digitize(year, edges) - 1
    fig, axes = plt.subplots(2, 1, figsize=(9.5, 6.2), sharex=True, gridspec_kw={"hspace": 0.12})
    for ax, b in zip(axes, BANDS):
        off = np.array([float(r[f"{b}_offset_mmag"]) for r in rows])
        med, lo, hi, ext = [], [], [], []
        for k in range(len(centre)):
            v = off[which == k]
            if v.size < 5:
                med.append(np.nan), lo.append(np.nan), hi.append(np.nan), ext.append(np.nan)
                continue
            p16, p50, p84 = np.percentile(v, [16, 50, 84])
            med.append(p50), lo.append(p16), hi.append(p84)
            ext.append(v[np.argmax(np.abs(v))])
        ax.fill_between(centre, lo, hi, color=COLOR[b], alpha=0.22, linewidth=0,
                        label="16-84 % of fields")
        ax.plot(centre, med, color=COLOR[b], marker=MARKER[b], markersize=6,
                markeredgecolor=SURFACE, markeredgewidth=0.8, label="median over fields")
        ax.plot(centre, ext, linestyle="none", marker="x", color=INK_2, markersize=5,
                label="most extreme field visit")
        ax.axhline(0.0, color=INK, linewidth=1.0, linestyle="--")
        ax.set_ylabel("offset (mmag, + = fainter)")
        ax.set_title(f"{LABEL[b]}: visit zero-point offset from comparison stars")
        ax.legend(loc="upper left", ncol=3, fontsize=8.5)
    axes[1].set_xlabel("year")
    return _save(fig, path)


def plot_pair_scatter(stars_csv: Path, targets_json: Path, sample_csv: Path, path: Path) -> Path:
    """Scatter of (pipeline - unTimely) over matched visits against AllWISE W1,
    for comparison stars and for sample targets. A difference of two
    measurements of the same visit; not a variability measure."""
    stars = _rows(stars_csv)
    targets = json.loads(Path(targets_json).read_text())
    w1_cat = {r["name"]: float(r["w1mpro"]) for r in _rows(sample_csv)}
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
    for ax, b in zip(axes, BANDS):
        sx = [float(r["w1_allwise"]) for r in stars if r[f"{b}_untimely_usable"] == "True"]
        sy = [float(r[f"{b}_untimely_pair_scatter_mmag"]) for r in stars
              if r[f"{b}_untimely_usable"] == "True"]
        tx, ty = [], []
        for name, row in targets.items():
            u = row.get("untimely", {})
            if u.get("status") == "ok" and u[b]["usable"]:
                tx.append(w1_cat[name])
                ty.append(u[b]["pair_scatter_mmag"])
        ax.plot(sx, sy, linestyle="none", marker=".", markersize=4, color=NEUTRAL, alpha=0.6,
                label=f"comparison stars ({len(sx)}): median {np.median(sy):.1f} mmag")
        ax.plot(tx, np.clip(ty, 0.5, 400), linestyle="none", marker=MARKER[b], markersize=5,
                color=COLOR[b], markeredgecolor=SURFACE, markeredgewidth=0.6,
                label=f"sample galaxies ({len(tx)}): median {np.median(ty):.1f} mmag")
        ax.set_yscale("log")
        ax.set_yticks([1, 3, 10, 30, 100, 300])
        ax.set_yticklabels(["1", "3", "10", "30", "100", "300"])
        ax.minorticks_off()
        ax.set_xlabel("AllWISE W1 (mag)")
        ax.set_title(f"{LABEL[b]}")
        ax.legend(loc="upper left", fontsize=8.5)
    axes[0].set_ylabel("scatter of pipeline - unTimely (mmag)")
    fig.suptitle("Two photometries of the same visits: visit-to-visit scatter of the difference",
                 x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return _save(fig, path)


def plot_star_injection(injection_json: Path, path: Path) -> Path:
    """Completeness of detector (a) or (b) on comparison stars, by star brightness."""
    doc = json.loads(Path(injection_json).read_text())
    rows = doc["completeness_pooled_over_baseline"]
    amps = doc["grid"]["amplitudes_mag"]
    groups = [g for g in ("9.5-10.5", "10.5-11.5") if g in doc["null"]]
    fig, axes = plt.subplots(1, len(groups), figsize=(11, 4.3), sharey=True)
    for ax, g in zip(np.atleast_1d(axes), groups):
        for rise in doc["grid"]["rise_years"]:
            sel = sorted((r for r in rows if r["star_w1_bin"] == g and r["rise_yr"] == rise),
                         key=lambda r: r["amp_mag"])
            ax.plot([r["amp_mag"] for r in sel], [r["completeness_a_or_b"] for r in sel],
                    color=RISE_COLOR[rise], marker=RISE_MARKER[rise], markersize=6,
                    markeredgecolor=SURFACE, markeredgewidth=0.8, linewidth=1.8,
                    label=RISE_LABEL[rise])
        ax.set_xscale("log")
        ax.set_xticks(amps)
        ax.set_xticklabels([f"{a:g}" for a in amps], fontsize=8.5)
        ax.minorticks_off()
        ax.set_ylim(-0.03, 1.03)
        ax.set_xlabel("injected brightening (mag, both bands)")
        ax.set_title(f"host stars with W1 {g} ({doc['null'][g]['n_stars']} stars)")
    np.atleast_1d(axes)[0].set_ylabel("fraction recovered by (a) or (b)")
    np.atleast_1d(axes)[0].legend(loc="lower right")
    fig.suptitle("Injection test on comparison stars of the sample fields "
                 "(all baseline lengths pooled)", x=0.01, ha="left", fontweight="bold",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return _save(fig, path)
