"""Stage 10 figures: identifiability on synthetic samples, example objects with the optical driver
and the two fitted kernels, the lag difference against the stage 9 lag, lag against luminosity.
Same visual system as plots.py and plots9.py.

Colour follows the band (W1, W2) or the driver setting; text stays in ink. No axis carries two scales.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from .lightcurve import mjd_to_year  # noqa: E402
from .plots import GRID, INK, INK_2, NEUTRAL, _save  # noqa: E402
from .plots9 import BAND_W1, BAND_W2, MODEL, _title  # noqa: E402

DRIVER_COLOUR = {"fast": "#2a78d6", "slow": "#eb6834"}
DRIVER_LABEL = {"fast": "driver 200 d", "slow": "driver 1000 d"}
OPTICAL = "#5b6470"
# the labels of docs/PREDICTIONS_10.md in English (the figure font has no Hangul)
LABEL_EN = {"맞음": "consistent", "다름": "different", "판정불가": "undetermined", "이 자료로 시험할 수 없음": "not testable with these data",
            "판정불가(표본 부족)": "undetermined (sample too small)"}


def _bars(ax, ident: dict, truths, key: str, s: str = "d", centre: bool = True) -> None:
    """Median and null range of one statistic over the truth models, one marker per driver setting."""
    for driver, dx, filled in (("fast", -0.16, True), ("slow", 0.16, False)):
        tab = ident["drivers"][driver]["truths"]
        xs, med, lo, hi = [], [], [], []
        for i, name in enumerate(truths):
            b = tab[name][s].get(key, {})
            if "median" not in b:
                continue
            xs.append(i + dx)
            med.append(b["median"])
            lo.append(b["lo"])
            hi.append(b["hi"])
        med, lo, hi = np.array(med), np.array(lo), np.array(hi)
        ax.errorbar(xs, med, yerr=[med - lo, hi - med], fmt="o" if filled else "s", color=DRIVER_COLOUR[driver],
                    markerfacecolor=DRIVER_COLOUR[driver] if filled else "white", markersize=5.5, elinewidth=1.0, linewidth=0,
                    label=DRIVER_LABEL[driver])
    ax.set_xticks(range(len(truths)))
    ax.set_xticklabels(truths, rotation=45, ha="right", fontsize=8)
    ax.grid(axis="y", color=GRID, linewidth=0.6)
    if centre:
        ax.axhline(0.0, color=INK, linewidth=0.8)


def plot_identifiability(ident: dict, path: Path, observed: dict | None = None) -> Path:
    fig, axes = plt.subplots(2, 3, figsize=(17.5, 9.6), gridspec_kw={"wspace": 0.24, "hspace": 0.52, "top": 0.9, "right": 0.985,
                                                                      "left": 0.05, "bottom": 0.1})
    lag_truths = ["A_20", "A_50", "A_100", "A_200", "A_300", "A_lyu", "A_mandal", "A_s9", "G_lyu", "E_lyu", "W_thin", "W_wide"]
    for ax, key, label in ((axes[0, 0], "bias_tau1", "W1"), (axes[0, 1], "bias_tau2", "W2")):
        _bars(ax, ident, lag_truths, key)
        ax.axhspan(-10, 10, color=NEUTRAL, alpha=0.18, linewidth=0)
        ax.set_ylabel("median fitted lag - median lag put in (days)")
        ax.set_title(f"Lag of {label} behind the optical driver, detected set\n(points: median of 50 synthetic samples; bars: null range; band: +-10 d)",
                     fontsize=9.5, loc="left", color=INK)
    axes[0, 0].legend(fontsize=8, frameon=False, loc="best")
    ax = axes[0, 2]
    _bars(ax, ident, ["A0"] + lag_truths + ["B_20", "B_40", "B_both", "C_0", "C_1yr", "C_3yr"], "bias_dtau")
    ax.axhspan(-8, 8, color=NEUTRAL, alpha=0.18, linewidth=0)
    ax.set_ylabel("median fitted (tau2 - tau1) - median put in (days)")
    ax.set_title("Lag difference W2 - W1\n(band: +-8 d, the criterion of section 5.1)", fontsize=9.5, loc="left", color=INK)
    ax = axes[1, 0]
    _bars(ax, ident, ["X0", "Z0", "A0", "A_100", "A_lyu", "A_mandal", "B_40", "C_1yr"], "frac_detected", s="all", centre=False)
    ax.set_ylabel("share of objects detected in both bands")
    ax.set_ylim(0, 1)
    ax.set_title("Detection (r_w above the per-object null of X0)\nX0: infrared made from an independent driver", fontsize=9.5, loc="left", color=INK)
    if observed:
        ax.axhline(observed["frac_detected"], color=MODEL, linewidth=1.6, linestyle="--")
        ax.text(0.99, observed["frac_detected"] + 0.015, f"real light curves: {observed['frac_detected']:.2f} (null: drivers of other objects)",
                transform=ax.get_yaxis_transform(), fontsize=8, color=INK, ha="right")
    ax = axes[1, 1]
    _bars(ax, ident, ["A_lyu", "A_mandal", "G_lyu", "W_wide", "B_20", "B_40", "B_both", "C_0", "C_1yr", "C_3yr"], "med_alpha_w1")
    for y in (0.2, 0.4):
        ax.axhline(y, color=NEUTRAL, linewidth=0.8, linestyle=":")
    ax.set_ylabel("sample median of the undelayed share of W1 (alpha)")
    ax.set_title("Undelayed share in W1 (model MA)\nput in: 0.2 (B_20), 0.4 (B_40), 0.3 (B_both), 0 elsewhere", fontsize=9.5, loc="left", color=INK)
    if observed:
        ax.axhline(observed["med_alpha_w1"], color=MODEL, linewidth=1.6, linestyle="--")
        ax.text(0.01, observed["med_alpha_w1"] + 0.006, f"real light curves: {observed['med_alpha_w1']:.2f}", transform=ax.get_yaxis_transform(),
                fontsize=8, color=INK)
    ax = axes[1, 2]
    _bars(ax, ident, ["A0", "A_100", "A_lyu", "A_mandal", "A_s9", "G_lyu", "E_lyu", "B_40", "C_0", "C_1yr", "C_3yr"], "m_res")
    ax.set_ylabel("median residual branch offset (mag)")
    ax.set_title("Colour offset left after two fixed kernels (m_res)\nC_1yr, C_3yr: response that depends on the history", fontsize=9.5, loc="left", color=INK)
    over = ident["overall"]
    text = "  |  ".join(f"{k}: {'yes' if v else 'no'}" for k, v in over.items())
    fig.text(0.05, 0.02, f"Criteria of section 5.1 (both driver settings): {text}", fontsize=8.5, color=INK_2)
    _title(fig, f"Stage 10: what the optical-driver fit can recover ({ident['n_objects']} objects, real ZTF nights and WISE visits, synthetic light curves)")
    return _save(fig, path)


def plot_examples(objects: list[dict], path: Path) -> Path:
    """Per object: optical nights with the conditional mean, W1 and W2 with the fitted models, the two kernels."""
    n = len(objects)
    fig = plt.figure(figsize=(16.5, 4.6 * n + 0.6))
    gs = fig.add_gridspec(2 * n, 2, width_ratios=[3.2, 1.0], hspace=0.55, wspace=0.14, top=1 - 0.9 / (4.6 * n + 0.6),
                          left=0.055, right=0.985, bottom=0.06)
    for k, o in enumerate(objects):
        ax_o = fig.add_subplot(gs[2 * k, 0])
        ax_i = fig.add_subplot(gs[2 * k + 1, 0], sharex=ax_o)
        ax_k = fig.add_subplot(gs[2 * k:2 * k + 2, 1])
        yr = mjd_to_year(o["opt_t"])
        ax_o.plot(mjd_to_year(o["grid_t"]), o["grid_mean"], color=OPTICAL, linewidth=1.0, label="conditional mean (damped random walk)")
        ax_o.fill_between(mjd_to_year(o["grid_t"]), o["grid_mean"] - o["grid_sd"], o["grid_mean"] + o["grid_sd"], color=OPTICAL,
                          alpha=0.18, linewidth=0)
        ax_o.plot(yr, o["opt_y"], ".", color=INK, markersize=2.5, label=f"ZTF {o['band']} nights")
        ax_o.set_ylabel("optical flux\n(relative)")
        ax_o.set_title(o["title"], fontsize=9.5, loc="left", color=INK)
        ax_o.legend(fontsize=7.5, frameon=False, loc="upper left", ncol=2)
        ax_o.grid(color=GRID, linewidth=0.6)
        idx = np.asarray(o["idx"])
        used = np.zeros(o["mjd"].size, dtype=bool)
        used[idx] = True
        for f, s, pred, colour, label in ((o["f1"], o["s1"], o["pred_w1"], BAND_W1, "W1"), (o["f2"], o["s2"], o["pred_w2"], BAND_W2, "W2")):
            med = float(np.median(f[idx]))
            t = mjd_to_year(o["mjd"])
            ax_i.errorbar(t[used], f[used] / med, yerr=s[used] / med, fmt="o", color=colour, markersize=5, elinewidth=0.9, linewidth=0,
                          label=f"{label} (visits used)")
            ax_i.errorbar(t[~used], f[~used] / med, yerr=s[~used] / med, fmt="o", color=colour, markerfacecolor="white", markersize=4,
                          elinewidth=0.7, linewidth=0, alpha=0.6)
            ax_i.plot(t[used], np.asarray(pred) / med, "-", color=colour, linewidth=1.3, alpha=0.9)
        ax_i.set_ylabel("infrared flux /\nmedian of used visits")
        ax_i.set_xlabel("year")
        ax_i.set_xlim(mjd_to_year(o["opt_t"][0]) - 0.3, mjd_to_year(max(o["opt_t"][-1], o["mjd"][-1])) + 0.2)
        ax_i.legend(fontsize=7.5, frameon=False, loc="upper left", ncol=2)
        ax_i.grid(color=GRID, linewidth=0.6)
        ax_i.text(0.995, 0.04, "lines: fitted model M1 at the used visits; open symbols: visits outside the window",
                  transform=ax_i.transAxes, fontsize=7.5, color=INK_2, ha="right")
        for tau, colour, label, lo, hi in ((o["tau_w1"], BAND_W1, "W1", o["q_w1"][0], o["q_w1"][1]),
                                           (o["tau_w2"], BAND_W2, "W2", o["q_w2"][0], o["q_w2"][1])):
            a, b = sorted((0.5 * tau, 1.5 * tau))
            h = 1.0 / max(b - a, 5.0)
            ax_k.fill_between([a, b], [h, h], color=colour, alpha=0.3, linewidth=0)
            ax_k.plot([a, a, b, b], [0, h, h, 0], color=colour, linewidth=1.4, label=f"{label}: lag {tau:.0f} d (16-84%: {lo:.0f} to {hi:.0f})")
            ax_k.axvline(tau, color=colour, linewidth=0.9, linestyle=":")
        ax_k.axvline(0.0, color=INK, linewidth=0.8)
        ax_k.set_xlim(-110, 710)
        ax_k.set_xlabel("delay behind the optical driver (days, observed frame)")
        ax_k.set_ylabel("fitted kernel (1 / day)")
        ax_k.legend(fontsize=7.8, frameon=False, loc="upper right")
        ax_k.set_title(o["kernel_title"], fontsize=9, loc="left", color=INK)
        ax_k.grid(color=GRID, linewidth=0.6)
    _title(fig, "Stage 10: optical driver (ZTF), W1 and W2 with the fitted fixed kernels (top-hat of width = lag)")
    return _save(fig, path)


def plot_dtau(rows: list[dict], verdict: dict, path: Path) -> Path:
    """Lag difference from the optical driver against the stage 9 lag of the same objects."""
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 6.2), gridspec_kw={"wspace": 0.3, "top": 0.86, "right": 0.985, "left": 0.06,
                                                                      "width_ratios": [1.1, 1]})
    ax = axes[0]
    rest = [r for r in rows if not r["in_d"] and np.isfinite(r["tau9"])]
    det = [r for r in rows if r["in_d"] and np.isfinite(r["tau9"])]
    ax.plot([r["tau9"] for r in rest], [r["dtau"] for r in rest], "o", color=NEUTRAL, markerfacecolor="white", markersize=4.5,
            label=f"other fitted objects ({len(rest)})")
    if det:
        x, y = np.array([r["tau9"] for r in det]), np.array([r["dtau"] for r in det])
        lo, hi = np.array([r["dtau_q16"] for r in det]), np.array([r["dtau_q84"] for r in det])
        ax.errorbar(x, y, yerr=[np.clip(y - lo, 0, None), np.clip(hi - y, 0, None)], fmt="o", color=MODEL, markersize=6, elinewidth=0.9,
                    linewidth=0, label=f"detected in both bands ({len(det)}); bars: 16-84% of the joint fit")
        ax.plot(np.median(x), np.median(y), "D", color=INK, markersize=9, label=f"medians: ({np.median(x):.0f}, {np.median(y):.0f}) d")
    lim = (-160, 360)
    ax.plot(lim, lim, color=INK, linewidth=0.8)
    ax.axhline(0, color=GRID, linewidth=0.8)
    ax.axvline(0, color=GRID, linewidth=0.8)
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    n_out = sum(1 for r in rest + det if not lim[0] <= r["dtau"] <= lim[1])
    if n_out:
        ax.text(0.99, 0.02, f"{n_out} objects outside the plotted range", transform=ax.transAxes, fontsize=8, color=INK_2, ha="right")
    ax.set_xlabel("stage 9: lag of W2 behind W1, fitted without a driver (days, observed frame)")
    ax.set_ylabel("stage 10: tau(W2) - tau(W1) from the optical driver (days)")
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    ax.grid(color=GRID, linewidth=0.6)
    ax.set_title("Same objects, two ways to the delay of W2 behind W1 (line: equal)", fontsize=9.5, loc="left", color=INK)
    ax = axes[1]
    q1 = verdict["q1"]
    names, rng, med = [], [], []
    for d in ("fast", "slow"):
        b = q1["null_A_s9"][d]
        if "lo" in b:
            names.append(f"truth = stage 9 lag\n{DRIVER_LABEL[d]}")
            rng.append((b["lo"], b["hi"]))
            med.append(b["median"])
    for i, ((lo, hi), m) in enumerate(zip(rng, med)):
        ax.plot([lo, hi], [i, i], color=NEUTRAL, linewidth=7, solid_capstyle="butt")
        ax.plot(m, i, "|", color=INK, markersize=14)
    if q1["q"] is not None:
        ax.axvline(q1["q"], color=MODEL, linewidth=2.0)
        ci = q1["q_bootstrap_ci95"]
        if ci:
            ax.axvspan(ci[0], ci[1], color=MODEL, alpha=0.12, linewidth=0)
        ax.text(q1["q"], len(names) - 0.4, f" measured q = {q1['q']:.2f}", color=INK, fontsize=9, va="bottom")
    ax.plot([0.7, 1.6], [-0.45, -0.45], color=INK_2, linewidth=2.0)
    ax.text(0.7, -0.4, "registered prediction 0.7-1.6", fontsize=8, color=INK_2, va="bottom")
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=8.5)
    ax.set_ylim(-0.7, len(names) + 0.1)
    ax.set_xlabel("q = median(tau2 - tau1) / median(stage 9 lag), detected set")
    ax.grid(axis="x", color=GRID, linewidth=0.6)
    ax.set_title(f"Q1, label by the registered rule: '{LABEL_EN.get(q1['label'], q1['label'])}'\n(grey: null range of 50 synthetic samples; shaded: bootstrap 95% of q)", fontsize=9.5, loc="left", color=INK)
    _title(fig, "Stage 10, Q1: does the lag difference from the optical driver match the stage 9 lag?")
    return _save(fig, path)


def plot_lag_luminosity(rows: list[dict], relations: dict, verdict: dict, path: Path) -> Path:
    """Rest-frame lag of W1 and W2 against log L_bol with the two literature relations (points unlabelled)."""
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 6.4), sharey=True, gridspec_kw={"wspace": 0.06, "top": 0.8, "right": 0.985, "left": 0.06})
    x = np.linspace(43.0, 46.6, 50)
    for ax, b, colour in ((axes[0], "w1", BAND_W1), (axes[1], "w2", BAND_W2)):
        rest = [r for r in rows if not r["in_d"] and np.isfinite(r["logl"]) and r[f"rest_{b}"] > 0]
        det = [r for r in rows if r["in_d"] and np.isfinite(r["logl"]) and r[f"rest_{b}"] > 0]
        ax.plot([r["logl"] for r in rest], [r[f"rest_{b}"] for r in rest], "o", color=NEUTRAL, markerfacecolor="white", markersize=4.5,
                label=f"other fitted objects with a positive lag ({len(rest)})")
        ax.plot([r["logl"] for r in det], [r[f"rest_{b}"] for r in det], "o", color=colour, markersize=6.5,
                label=f"detected in both bands ({len(det)})")
        ax.plot(x, relations["lyu"][b](x), color=INK, linewidth=1.6, label="Lyu, Rieke & Smith 2019, Eqs. 15-16")
        ax.plot(x, relations["mandal"][b](x), color=INK, linewidth=1.6, linestyle="--", label="Mandal et al. 2024, Table 5 (disc-corrected)")
        ax.plot(x, relations["mandal_non_adc"][b](x), color=INK_2, linewidth=1.1, linestyle=":", label="Mandal et al. 2024, Table 5 (no disc correction)")
        ax.plot(x, relations["edge"](x), color=GRID, linewidth=1.0)
        ax.text(46.55, float(relations["edge"](np.array([46.55]))[0]) * 1.04, "end of the lag grid (465 d observed, z = 0.05)", fontsize=7.5,
                color=INK_2, ha="right", va="bottom")
        ax.set_yscale("log")
        ax.set_ylim(4, 1500)
        ax.set_xlim(43.0, 46.6)
        ax.set_xlabel("log L_bol (erg/s; BASS DR2)")
        ax.grid(color=GRID, linewidth=0.6)
        o2 = verdict["q2"]
        ax.set_title(f"{b.upper()}: median lag / relation = {o2['lyu']['observed'][f'k_lyu_{b}']:.2f} (Lyu), "
                     f"{o2['mandal']['observed'][f'k_mandal_{b}']:.2f} (Mandal);  Spearman rho = {o2['lyu']['observed'][f'rho_{b}']:+.2f}",
                     fontsize=9.5, loc="left", color=INK)
    axes[0].set_ylabel("lag behind the optical driver / (1 + z)  (days)")
    axes[0].legend(fontsize=7.8, frameon=False, loc="lower right")
    fig.text(0.06, 0.875, f"Labels by the registered rule: Lyu '{LABEL_EN.get(verdict['q2']['lyu']['label'])}', Mandal '{LABEL_EN.get(verdict['q2']['mandal']['label'])}'. "
             "Objects are not labelled; lags that are negative or at the end of the grid are not in the detected set.", fontsize=8.5, color=INK_2)
    _title(fig, "Stage 10, Q2: lag behind the optical driver against luminosity (same definition as the two relations)")
    return _save(fig, path)
