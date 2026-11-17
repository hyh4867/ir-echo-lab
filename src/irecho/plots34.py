"""Figures for stages 0, 3 and 4. Same visual system as plots.py."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from .lightcurve import mjd_to_year  # noqa: E402
from .plots import (COLOR, GRID, INK, INK_2, LABEL, MARKER, NEUTRAL, SHADE, SURFACE,  # noqa: E402
                    _save)
from .slowrise import model_curve  # noqa: E402
from .stats import constant_chi2  # noqa: E402

# Categorical slots 1-4 of the same palette (validated as a set on the light surface).
METHOD_COLOR = {"a": "#2a78d6", "b": "#eb6834", "a_or_b": "#1baf7a", "ml_matched": "#eda100"}
METHOD_MARKER = {"a": "o", "b": "s", "a_or_b": "D", "ml_matched": "^"}
METHOD_LABEL = {"a": "(a) fast flare", "b": "(b) slow rise", "a_or_b": "(a) or (b)",
                "ml_matched": "isolation forest"}


def plot_differential(before: list[dict], after: list[dict], path: Path) -> Path:
    """Per-star reduced chi-square before and after the field-offset correction."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    bins = np.logspace(np.log10(0.2), np.log10(8.0), 21)
    for ax, b in zip(axes, ("w1", "w2")):
        c0 = np.array([constant_chi2(v[f"{b}_mag"], v[f"{b}_err_scatter"]).chi2_red for v in before])
        c1 = np.array([constant_chi2(v[f"{b}_mag"], v[f"{b}_err"]).chi2_red for v in after])
        ax.hist(c0, bins=bins, histtype="stepfilled", color=NEUTRAL, alpha=0.35, linewidth=0)
        ax.hist(c0, bins=bins, histtype="step", color=NEUTRAL, linewidth=2,
                label=f"before: median {np.median(c0):.2f}")
        ax.hist(c1, bins=bins, histtype="stepfilled", color=COLOR[b], alpha=0.25, linewidth=0)
        ax.hist(c1, bins=bins, histtype="step", color=COLOR[b], linewidth=2,
                label=f"after: median {np.median(c1):.2f}")
        ax.axvline(1.0, color=INK, linestyle="--", linewidth=1.2)
        ax.set_xscale("log")
        ax.set_xticks([0.2, 0.5, 1, 2, 5])
        ax.set_xticklabels(["0.2", "0.5", "1", "2", "5"])
        ax.minorticks_off()
        ax.set_ylim(0, ax.get_ylim()[1] * 1.35)
        ax.set_xlabel("reduced chi-square of a constant model")
        ax.set_title(f"{LABEL[b]}: {len(c0)} reference stars")
        ax.legend(loc="upper right")
    axes[0].set_ylabel("number of stars")
    fig.suptitle("Field zero-point correction (leave-one-out), in-visit scatter errors",
                 x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return _save(fig, path)


def plot_untimely(v: dict, u: dict, path: Path) -> Path:
    """NGC 3786: profile-fit visit medians and unTimely coadd photometry."""
    fig, axes = plt.subplots(2, 1, figsize=(9.5, 6.4), sharex=True,
                             gridspec_kw={"hspace": 0.1})
    for ax, b in zip(axes, ("w1", "w2")):
        ax.errorbar(mjd_to_year(v["mjd"]), 1e3 * v[f"{b}_flux_jy"],
                    yerr=1e3 * v[f"{b}_flux_err_jy"], fmt=MARKER[b], color=COLOR[b],
                    markersize=7, markeredgecolor=SURFACE, elinewidth=1.2, zorder=3,
                    label="this pipeline: single-exposure profile fits, visit median")
        ax.errorbar(mjd_to_year(u["mjd"]), 1e3 * u[f"{b}_flux_jy"],
                    yerr=1e3 * u[f"{b}_flux_err_jy"], fmt="D", color=INK, markersize=6,
                    markerfacecolor=SURFACE, markeredgewidth=1.6, elinewidth=1.2, zorder=4,
                    label="unTimely: PSF fit on the per-visit unWISE coadd")
        ax.set_ylabel(f"{b.upper()} flux density (mJy)")
        ax.set_ylim(0, None)
    axes[0].legend(loc="upper left", fontsize=9)
    axes[1].set_xlabel("year")
    axes[0].set_title("NGC 3786: two photometries of the same visits (unTimely ends in 2020)")
    return _save(fig, path)


def plot_physics(phys: dict, path: Path) -> Path:
    """Temperature, luminosity and blackbody radius of the variable component."""
    rows = phys["flagged_visits"]
    year = np.array([r["year"] for r in rows])
    fig, axes = plt.subplots(3, 1, figsize=(9.5, 9.0), sharex=True,
                             gridspec_kw={"hspace": 0.12})

    def series(ax, key, scale, color, marker):
        y = np.array([r[key] for r in rows]) * scale
        lo = np.array([r[f"{key}_p16"] for r in rows]) * scale
        hi = np.array([r[f"{key}_p84"] for r in rows]) * scale
        ax.plot(year, y, color=color, linewidth=1.4, zorder=2)
        ax.errorbar(year, y, yerr=[np.clip(y - lo, 0, None), np.clip(hi - y, 0, None)],
                    fmt=marker, color=color, markersize=7, markeredgecolor=SURFACE,
                    elinewidth=1.4, zorder=3)

    series(axes[0], "t_rest_k", 1.0, INK_2, "D")
    axes[0].set_ylabel("T of the variable\ncomponent (K)")
    series(axes[1], "l_dust_erg_s", 1e-42, INK_2, "D")
    axes[1].set_ylabel("L_dust (1e42 erg/s)")
    axes[1].set_ylim(0, None)
    series(axes[2], "r_bb_light_days", 1.0, INK_2, "D")
    shown = {}
    for lit in phys["sublimation_radius"]:
        if lit["date"] is None:
            continue
        corrected = "extinction-corrected" in lit["label"]
        key = "corrected" if corrected else "observed"
        yr = float(mjd_to_year(_date_mjd(lit["date"])))
        h = axes[2].scatter([yr], [lit["r_sub_graphite_light_days"]], s=70,
                            marker="s" if corrected else "o",
                            color=COLOR["w2"] if corrected else COLOR["w1"],
                            edgecolor=SURFACE, linewidth=1.0, zorder=4)
        shown.setdefault(key, h)
    axes[2].set_ylabel("radius (light-days)")
    axes[2].set_ylim(0, None)
    axes[2].set_xlabel("year")
    axes[2].legend(
        [plt.Line2D([], [], color=INK_2, marker="D", markeredgecolor=SURFACE, markersize=7),
         shown["observed"], shown["corrected"]],
        ["blackbody radius R_bb of the variable component (lower limit on size)",
         "graphite R_sub from L_bol (broad H-alpha, Son et al. 2026)",
         "graphite R_sub from extinction-corrected L_bol (Son et al. 2026)"],
        loc="upper left", fontsize=8.5)
    axes[0].set_title("NGC 3786: blackbody fit to the excess flux over the pre-flare baseline")
    return _save(fig, path)


def _date_mjd(date: str) -> float:
    from datetime import datetime
    return (datetime.strptime(date, "%Y-%m-%d") - datetime(1858, 11, 17)).days


def plot_optical(label: str, v: dict, base_w1_flux_mjy: float, optical: dict,
                 path: Path) -> Path:
    """ZTF g and r (medians around each WISE visit) and the W1 excess flux."""
    fig, axes = plt.subplots(3, 1, figsize=(9.5, 8.0), sharex=True,
                             gridspec_kw={"hspace": 0.12})
    styles = {"g": ("#1baf7a", "o", "ZTF g"), "r": ("#e34948", "s", "ZTF r")}
    for ax, band in zip(axes[:2], ("g", "r")):
        info = optical["bands"].get(band, {})
        color, marker, band_label = styles[band]
        if info.get("status") == "ok":
            bins = info["bins"]
            yr = [float(mjd_to_year(_date_mjd(b["wise_visit"]) + b["offset_days"])) for b in bins]
            ax.errorbar(yr, [b["mag"] for b in bins], yerr=[b["err"] for b in bins],
                        fmt=marker, color=color, markersize=7, markeredgecolor=SURFACE,
                        elinewidth=1.3, zorder=3,
                        label=f"{band_label}: median of the points within 91 d of a WISE visit")
            pre = [b["mag"] for b in bins if not b["post_outburst"]]
            ax.axhline(float(np.median(pre)), color=INK_2, linewidth=1.0, linestyle="--",
                       zorder=1)
        ax.invert_yaxis()
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo, hi - 0.22 * (lo - hi))  # headroom for the legend (axis is inverted)
        ax.set_ylabel(f"{band_label} (mag)")
        ax.legend(loc="upper left", fontsize=8.5)
    excess = 1e3 * v["w1_flux_jy"] - base_w1_flux_mjy
    axes[2].errorbar(mjd_to_year(v["mjd"]), excess, yerr=1e3 * v["w1_flux_err_jy"],
                     fmt=MARKER["w1"], color=COLOR["w1"], markersize=7, markeredgecolor=SURFACE,
                     elinewidth=1.3, zorder=3,
                     label="WISE W1 excess over the level before 2020-06 (MJD 59000)")
    axes[2].axhline(0.0, color=INK, linewidth=0.8)
    axes[2].set_ylabel("W1 excess flux (mJy)")
    axes[2].set_xlabel("year")
    axes[2].legend(loc="upper left", fontsize=8.5)
    if v["mjd"].size:
        axes[2].set_xlim(2017.8, float(mjd_to_year(v["mjd"].max())) + 0.4)
    for ax in axes:
        ax.axvline(float(mjd_to_year(59000.0)), color=INK_2, linestyle=":", linewidth=1.2)
    axes[0].set_title(f"{label}: optical (ZTF) and mid-infrared; dashed = optical level before 2020-06")
    return _save(fig, path)


def plot_injection(inj: dict, path: Path) -> Path:
    """Recovery fraction against injected amplitude, by host type and rise time."""
    rises = inj["grid"]["rise_years"]
    amps = inj["grid"]["amplitudes_mag"]
    kinds = (("control", "hosts: control galaxies"), ("star", "hosts: comparison stars"))
    fig, axes = plt.subplots(2, len(rises), figsize=(11.5, 6.6), sharex=True, sharey=True)
    for i, (kind, kind_label) in enumerate(kinds):
        for j, rise in enumerate(rises):
            ax = axes[i, j]
            rows = [r for r in inj["recovery"] if r["kind"] == kind
                    and r["baseline_visits"] == "all" and r["rise_years"] == rise]
            rows.sort(key=lambda r: r["amp_mag"])
            for m in ("a", "b", "a_or_b", "ml_matched"):
                frac = [r[m] / r["n"] if r["n"] else np.nan for r in rows]
                ax.plot(amps, frac, color=METHOD_COLOR[m], marker=METHOD_MARKER[m],
                        markersize=6, markeredgecolor=SURFACE, markeredgewidth=0.8,
                        linewidth=1.8, label=METHOD_LABEL[m])
            ax.set_xscale("log")
            ax.set_xticks(amps)
            ax.set_xticklabels([f"{a:g}" for a in amps], fontsize=8)
            ax.minorticks_off()
            ax.set_ylim(-0.03, 1.03)
            rise_label = "instant step" if rise == 0 else f"{rise:g}-yr linear rise"
            ax.set_title(f"{rise_label}; {kind_label}", fontsize=9.5)
            if j == 0:
                ax.set_ylabel("fraction recovered")
            if i == 1:
                ax.set_xlabel("injected brightening (mag, both bands)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", ncol=4, fontsize=9,
               bbox_to_anchor=(0.995, 0.995))
    fig.suptitle("Injection test: recovery by detector (all baseline lengths pooled)",
                 x=0.01, y=0.985, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return _save(fig, path)


def plot_slow_rise(name: str, v: dict, result, path: Path) -> Path:
    """One light curve with the model chosen by detector (b)."""
    fig, ax = plt.subplots(figsize=(9.5, 4.4))
    year = mjd_to_year(v["mjd"])
    k = result.onset_index if result.onset_index is not None else 0
    n_ref = max(k, 3)
    for b in ("w1", "w2"):
        ref = np.median(v[f"{b}_mag"][:n_ref])
        ax.errorbar(year, ref - v[f"{b}_mag"], yerr=v[f"{b}_err"], fmt=MARKER[b],
                    color=COLOR[b], markersize=7, markeredgecolor=SURFACE, elinewidth=1.2,
                    zorder=3, label=LABEL[b])
        ax.plot(year, ref - model_curve(v["mjd"], v[f"{b}_mag"], v[f"{b}_err"], result),
                color=COLOR[b], linewidth=1.8, zorder=2)
    ax.axhline(0.0, color=INK, linewidth=0.8)
    if result.onset_index:
        ax.axvline(float(mjd_to_year(result.onset_mjd)), color=INK_2, linestyle=":",
                   linewidth=1.2)
        ax.annotate("fitted start of the rise", xy=(float(mjd_to_year(result.onset_mjd)), 0.97),
                    xycoords=("data", "axes fraction"), xytext=(-6, 0),
                    textcoords="offset points", ha="right", va="top", fontsize=9, color=INK_2)
    ax.set_xlabel("year")
    ax.set_ylabel("brightening relative to the\nvisits before the fitted start (mag)")
    verdict = "flagged" if result.flagged else "not flagged"
    label = name.upper().replace("NGC", "NGC ") if name.startswith("ngc") else name
    ax.set_title(f"{label}: detector (b), {result.model} model, {verdict}")
    ax.legend(loc="upper left")
    return _save(fig, path)


__all__ = ["plot_differential", "plot_untimely", "plot_physics", "plot_injection",
           "plot_slow_rise", "GRID", "SHADE"]
