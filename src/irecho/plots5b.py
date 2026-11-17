"""Stage 5b figures. Same visual system as plots.py. These are made after the
predictions were committed; target light curves are shown here for the first
time."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from . import stage5b  # noqa: E402
from .lightcurve import mjd_to_year  # noqa: E402
from .plots import COLOR, INK, INK_2, LABEL, NEUTRAL, _save  # noqa: E402

TYPE_COLOR = {"1": "#2a78d6", "2": "#eb6834", "other": NEUTRAL}
TYPE_LABEL = {"1": "type 1 (Sy1-1.8)", "2": "type 2 (Sy1.9-2)", "other": "other"}
TYPE_MARKER = {"1": "o", "2": "s", "other": "D"}
GREEN = "#1baf7a"


def _rows(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def _f(x) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def _title(fig, text: str) -> None:
    fig.suptitle(text, x=0.01, ha="left", fontweight="bold", fontsize=11)


def plot_lightcurves(objects_csv: Path, path: Path) -> Path:
    """Visit light curves of every object flagged by detector (a) or (b)."""
    rows = [r for r in _rows(objects_csv) if r["brightened"] == "1"]
    ncol = 4
    nrow = int(np.ceil(len(rows) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(13, 2.6 * nrow), squeeze=False)
    for ax, r in zip(axes.ravel(), rows):
        v = stage5b.read_frozen(r["name"])
        year = mjd_to_year(v["mjd"])
        n_base = max(int(_f(r["start_index"])), 3)
        for b in ("w1", "w2"):
            dm = np.median(v[f"{b}_mag"][:n_base]) - v[f"{b}_mag"]
            ax.errorbar(year, dm, yerr=v[f"{b}_err"], color=COLOR[b], marker="o" if b == "w1" else "s",
                        markersize=3.5, linewidth=1.2, elinewidth=0.8, label=LABEL[b])
        ax.axhline(0.0, color=INK_2, linewidth=0.8, linestyle="--")
        tags = [t for t, k in (("a", "a_flagged"), ("b", "b_flagged"), ("slow rise", "q1"))
                if r[k] == "1"]
        note = "" if r["sample"] == "main" else " [separate]"
        ax.set_title(f"{r['counterpart'][:18]} ({r['type']}): {', '.join(tags)}{note}", fontsize=8.5)
        ax.tick_params(labelsize=8)
    for ax in axes.ravel()[len(rows):]:
        ax.set_visible(False)
    axes[0, 0].legend(loc="upper left", fontsize=8)
    fig.supylabel("brightening relative to the baseline visits (mag)", fontsize=10)
    _title(fig, "Objects flagged by detector (a) or (b): W1 and W2 visit light curves")
    fig.tight_layout(rect=(0.01, 0, 1, 0.97))
    return _save(fig, path)


def plot_q1(summary_json: Path, path: Path) -> Path:
    """Fractions with Wilson intervals against the predicted 5-20% band."""
    s = json.loads(summary_json.read_text())["floor"]
    fig, ax = plt.subplots(figsize=(9, 4.2))
    ax.axvspan(5, 20, color=GREEN, alpha=0.12, linewidth=0)
    ax.text(12.5, 1.02, "shaded: prediction, 5-20%", transform=ax.get_xaxis_transform(),
            ha="center", va="bottom", fontsize=9, color=INK_2)
    labels, y = [], 0
    series = (("main_all", COLOR["w1"], "o", "main sample"),
              ("main_clean", COLOR["w2"], "s", "main sample without QA flags"))
    items = (("slow rise: (b), fitted rise >= 0.2 mag in both bands", lambda q: (q["fraction"], q["ci95"])),
             ("completeness-corrected slow rise", lambda q: (q["completeness_corrected"]["fraction"],
                                                             q["completeness_corrected"]["ci95"])
              if "completeness_corrected" in q else None),
             ("detector (b) flag (>= 0.1 mag)", lambda q: (q["b_flag_0p1"]["fraction"], q["b_flag_0p1"]["ci95"])),
             ("detector (a) flag", lambda q: (q["a_flag"]["fraction"], q["a_flag"]["ci95"])),
             ("(a) or (b)", lambda q: (q["a_or_b"]["fraction"], q["a_or_b"]["ci95"])))
    for label, get in items:
        drawn = False
        for k, (key, color, marker, name) in enumerate(series):
            val = get(s[key]["q1"])
            if val is None:
                continue
            f, (lo, hi) = val
            yy = y + (0.16 if k == 0 else -0.16)
            ax.plot([100 * lo, 100 * hi], [yy, yy], color=color, linewidth=2)
            ax.plot(100 * f, yy, marker=marker, color=color, markersize=7, markeredgecolor="white",
                    label=name if y == 0 else None)
            ax.text(100 * hi + 0.4, yy, f"{100 * f:.1f}%", va="center", fontsize=8.5, color=INK)
            drawn = True
        if drawn:
            labels.append((y, label))
            y -= 1
    ax.set_yticks([p for p, _ in labels])
    ax.set_yticklabels([t for _, t in labels])
    ax.set_xlim(0, 25)
    ax.set_xlabel("fraction of the sample (%), 95% interval")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, fontsize=9)
    _title(fig, "Q1: how many AGN show a slow rise")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return _save(fig, path)


def plot_colour(objects_csv: Path, path: Path) -> Path:
    """Colour change of every brightened main-sample object, sorted."""
    rows = [r for r in _rows(objects_csv)
            if r["brightened"] == "1" and r["sample"] == "main" and r["delta_color"]]
    rows.sort(key=lambda r: _f(r["delta_color"]))
    fig, ax = plt.subplots(figsize=(9, 0.3 * len(rows) + 1.6))
    seen = set()
    for i, r in enumerate(rows):
        g = r["type_group"]
        ax.errorbar(_f(r["delta_color"]), i, xerr=_f(r["delta_color_err"]), color=TYPE_COLOR[g],
                    marker=TYPE_MARKER[g], markersize=6, markeredgecolor="white", elinewidth=1.5,
                    label=TYPE_LABEL[g] if g not in seen else None)
        seen.add(g)
    ax.axvline(0.0, color=INK, linewidth=1.0, linestyle="--")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r["counterpart"][:24] for r in rows], fontsize=8.5)
    ax.set_xlabel("change of W1 - W2 from the baseline to the three brightest visits (mag); "
                  "positive = redder")
    ax.legend(loc="lower right", fontsize=9)
    _title(fig, "Q2: colour change of the brightened objects")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return _save(fig, path)


def plot_ztf(objects_csv: Path, path: Path) -> Path:
    """Optical change against the W2 change between the IR-faint and IR-bright
    visits, for every main-sample object with a usable ZTF band."""
    rows = [r for r in _rows(objects_csv) if r["sample"] == "main"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
    for ax, band in zip(axes, ("g", "r")):
        seen = set()
        for r in rows:
            if r[f"ztf_{band}_status"] != "usable":
                continue
            g = r["type_group"]
            flagged = r["brightened"] == "1"
            ax.errorbar(_f(r[f"ztf_{band}_ir_change_w2"]), _f(r[f"ztf_{band}_dm_opt"]),
                        yerr=_f(r[f"ztf_{band}_dm_opt_err"]), color=TYPE_COLOR[g],
                        marker=TYPE_MARKER[g], markersize=8 if flagged else 5,
                        markeredgecolor=INK if flagged else "white", elinewidth=1.0, linestyle="",
                        label=TYPE_LABEL[g] if g not in seen else None)
            seen.add(g)
        ax.axhline(0.0, color=INK, linewidth=1.0, linestyle="--")
        ax.axhline(0.05, color=INK_2, linewidth=0.8, linestyle=":")
        ax.set_xlabel("W2 change between IR-faint and IR-bright visits (mag)")
        ax.set_title(f"ZTF {band}")
        ax.legend(loc="upper left", fontsize=9)
    axes[0].set_ylabel("optical change (mag); positive = brighter")
    _title(fig, "Q3: does the optical follow the infrared? (dark outline: flagged by (a) or (b))")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return _save(fig, path)


def plot_slope(objects_csv: Path, summary_json: Path, path: Path) -> Path:
    """Slope beta against the Eddington ratio, and the stacked means."""
    rows = [r for r in _rows(objects_csv) if r["sample"] == "main" and r["sl_status"] == "ok"
            and r["logedd"]]
    s = json.loads(summary_json.read_text())
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw={"width_ratios": [1.5, 1]})
    ax = axes[0]
    x = np.array([_f(r["logedd"]) for r in rows])
    beta = np.array([_f(r["sl_beta"]) for r in rows])
    sig = np.array([_f(r["sl_sigma_beta"]) for r in rows])
    meas = np.array([r["sl_measurable"] == "1" for r in rows])
    good = sig < 0.15
    ax.errorbar(x[good & ~meas], beta[good & ~meas], yerr=sig[good & ~meas], color=NEUTRAL,
                marker="o", markersize=3.5, linestyle="", elinewidth=0.7, alpha=0.7,
                label=f"slope error < 0.15 ({int((good & ~meas).sum())} objects)")
    ax.errorbar(x[meas], beta[meas], yerr=sig[meas], color=COLOR["w2"], marker="s", markersize=7,
                markeredgecolor=INK, linestyle="", elinewidth=1.5,
                label=f"'measurable' ({int(meas.sum())} objects)")
    ax.axvline(np.log10(stage5b.LAMBDA_EDD_THRESHOLD), color=INK, linewidth=1.2, linestyle="--")
    ax.axhline(0.0, color=INK_2, linewidth=0.8)
    ax.text(np.log10(stage5b.LAMBDA_EDD_THRESHOLD) + 0.03, 0.97, "threshold 0.02",
            transform=ax.get_xaxis_transform(), va="top", fontsize=8.5, color=INK)
    ax.set_ylim(-0.45, 0.45)
    ax.set_xlabel("log Eddington ratio (BASS DR2)")
    ax.set_ylabel("beta = d log T / d log (W2 excess)")
    ax.set_title("individual objects")
    ax.legend(loc="lower left", fontsize=8.5)
    ax = axes[1]
    pos = 0
    ticks = []
    for mode, sub, name in (("floor", "main_all", "main sample"),
                            ("floor", "main_clean", "without QA flags"),
                            ("nofloor", "main_all", "errors without floor")):
        prim = s[mode][sub]["h_breath"]["primary"]
        for k, (grp, color, marker) in enumerate((("below", COLOR["w1"], "o"),
                                                  ("above", COLOR["w2"], "s"))):
            g = prim[grp]
            yy = pos + (0.15 if k == 0 else -0.15)
            ax.errorbar(g["stack_beta"], yy, xerr=g["stack_beta_err"], color=color, marker=marker,
                        markersize=7, markeredgecolor="white", elinewidth=2,
                        label=f"{grp} the threshold" if pos == 0 else None)
            ax.plot(g["stack_beta_in_place_min"], yy, marker="|", color=color, markersize=14,
                    markeredgewidth=2, label="lower bound for heating in place" if pos == 0 and k == 0 else None)
        ticks.append((pos, name))
        pos -= 1
    ax.axvline(0.0, color=INK, linewidth=1.0, linestyle="--")
    ax.set_yticks([p for p, _ in ticks])
    ax.set_yticklabels([t for _, t in ticks])
    ax.set_ylim(pos + 0.4, 0.6)
    ax.set_xlabel("weighted mean beta of all objects with a slope")
    ax.set_title("stacked (reference, not used for the verdict)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), fontsize=8, ncol=1)
    _title(fig, "H-breath: slope of the variable component against the Eddington ratio")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return _save(fig, path)


def plot_hysteresis(objects_csv: Path, validation_json: Path, path: Path) -> Path:
    """Observed D against what the toy models give."""
    rows = [r for r in _rows(objects_csv) if r["sample"] == "main" and r["sl_hyst_measurable"] == "1"]
    d = np.array([_f(r["sl_d_hyst"]) for r in rows])
    val = json.loads(validation_json.read_text())["hysteresis"]
    pooled = {h["model"]: h["d_p16_p50_p84"] for h in val if h["driver"] == "all"}
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
    ax = axes[0]
    ax.hist(d, bins=np.arange(-0.3, 0.85, 0.05), color=COLOR["w1"], alpha=0.85, edgecolor="white")
    ax.axvline(0.0, color=INK, linewidth=1.0, linestyle="--")
    ax.axvline(np.median(d), color=INK, linewidth=1.5)
    ax.text(np.median(d) + 0.01, 0.95, f"median {np.median(d):+.2f}", transform=ax.get_xaxis_transform(),
            va="top", fontsize=9)
    ax.set_xlabel("D = beta(declining) - beta(rising)")
    ax.set_ylabel("number of objects")
    ax.set_title(f"sample: {d.size} objects with a measurable D")
    ax = axes[1]
    order = [("tau 0", "no hysteresis"), ("tau 0.25", "re-formation 0.25 yr"),
             ("tau 0.5", "0.5 yr"), ("tau 1", "1 yr"), ("tau 2", "2 yr"), ("tau 5", "5 yr"),
             ("tau 10", "10 yr"), ("tau inf", "never"), ("lag 10", "W2 lags W1 by 10 d"),
             ("lag 30", "by 30 d"), ("lag 90", "by 90 d")]
    for i, (key, label) in enumerate(order):
        lo, mid, hi = pooled[key]
        color = COLOR["w2"] if key.startswith("lag") else (NEUTRAL if key == "tau 0" else GREEN)
        ax.plot([lo, hi], [-i, -i], color=color, linewidth=2)
        ax.plot(mid, -i, marker="s" if key.startswith("lag") else "o", color=color, markersize=7,
                markeredgecolor="white")
    lo, mid, hi = np.percentile(d, [16, 50, 84])
    ax.axvspan(lo, hi, color=COLOR["w1"], alpha=0.12, linewidth=0)
    ax.axvline(mid, color=COLOR["w1"], linewidth=1.5)
    ax.text(mid + 0.01, 0.97, "sample median, 16-84%", transform=ax.get_xaxis_transform(),
            fontsize=8.5, color=INK, va="top")
    ax.axvline(0.0, color=INK, linewidth=1.0, linestyle="--")
    ax.set_yticks([-i for i in range(len(order))])
    ax.set_yticklabels([label for _, label in order], fontsize=8.5)
    ax.set_xlabel("D in the toy models (median, 16-84%)")
    ax.set_title("toy models with the real cadence and errors")
    _title(fig, "Hysteresis statistic: the sample against re-formation delay and inter-band lag")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return _save(fig, path)


def plot_validation(validation_json: Path, path: Path) -> Path:
    """How often an object is 'measurable', and what r it gives, by toy model."""
    cells = [c for c in json.loads(validation_json.read_text())["cells"] if c["group"] == "slope"]
    drivers = ["drw0.15", "drw0.3", "drw0.6", "ramp1.0", "ramp2.0"]
    names = ["random walk\n0.15 mag", "random walk\n0.3 mag", "random walk\n0.6 mag",
             "ramp\n1 mag", "ramp\n2 mag"]
    style = {"boundary": (COLOR["w1"], "o", "boundary shift"),
             "in place": (COLOR["w2"], "s", "heating in place"),
             "mix 50%": (GREEN, "D", "50% mixture")}
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
    for model, (color, marker, label) in style.items():
        for host, ls in ((0.0, "-"), (0.5, "--")):
            sel = {c["driver"]: c for c in cells if c["model"] == model and c["host"] == host}
            frac = [100.0 * sel[d]["n_measurable"] / sel[d]["n"] for d in drivers]
            axes[0].plot(range(len(drivers)), frac, color=color, marker=marker, linestyle=ls,
                         linewidth=2, markersize=6, markeredgecolor="white",
                         label=f"{label}, host {int(100 * host)}% of W1")
            r = [sel[d]["measurable_r_p16_p50_p84"] if sel[d]["n_measurable"] >= 10 else None
                 for d in drivers]
            xs = [i + (-0.12 if host == 0 else 0.12) for i, v in enumerate(r) if v]
            mid = [v[1] for v in r if v]
            err = [[v[1] - v[0] for v in r if v], [v[2] - v[1] for v in r if v]]
            axes[1].errorbar(xs, mid, yerr=err, color=color, marker=marker, linestyle="",
                             markersize=6, markeredgecolor="white", elinewidth=1.5,
                             alpha=1.0 if host == 0 else 0.55)
    for ax in axes:
        ax.set_xticks(range(len(drivers)))
        ax.set_xticklabels(names, fontsize=8.5)
        ax.set_xlabel("driver of the synthetic light curve")
    axes[0].set_ylabel("'measurable' objects (%)")
    axes[0].legend(loc="upper left", fontsize=8)
    axes[1].axhline(0.0, color=INK_2, linewidth=0.8)
    axes[1].axhline(1.0, color=INK_2, linewidth=0.8, linestyle="--")
    axes[1].set_ylabel("breathing index r of measurable objects (median, 16-84%)")
    axes[1].set_title("faded markers: host 50% of W1", fontsize=9)
    _title(fig, "Toy-model test of the slope estimator (400 synthetic objects per point)")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return _save(fig, path)


def plot_flares(flares_json: Path, other_json: Path, path: Path) -> Path:
    """Stage 6A: brightening, dust temperature and dust luminosity of the
    variable component for each flare candidate."""
    d = json.loads(flares_json.read_text())
    other = json.loads(other_json.read_text())
    fig, axes = plt.subplots(len(d), 3, figsize=(13, 2.7 * len(d)), squeeze=False)
    for row, (name, o) in zip(axes, d.items()):
        v = stage5b.read_frozen(name)
        year = mjd_to_year(v["mjd"])
        n_base = o["n_baseline_visits"]
        ax = row[0]
        for b in ("w1", "w2"):
            dm = np.median(v[f"{b}_mag"][:n_base]) - v[f"{b}_mag"]
            ax.errorbar(year, dm, yerr=v[f"{b}_err"], color=COLOR[b], marker="o" if b == "w1" else "s",
                        markersize=3.5, linewidth=1.2, elinewidth=0.8, label=LABEL[b])
        g = other.get(name, {}).get("asassn", {}).get("g", {})
        if g.get("status") == "usable":
            gy = mjd_to_year(np.array(g["bins_mjd"]))
            gm = np.array(g["bins_mag"])
            ax.plot(gy, np.median(gm) - gm, color=GREEN, marker="D", markersize=3.5, linewidth=1.0,
                    label="ASAS-SN g (binned)")
        for ep in o["episodes"]:
            ax.axvspan(mjd_to_year(ep["start_mjd"]) - 0.1, mjd_to_year(ep["end_mjd"]) + 0.1,
                       color=NEUTRAL, alpha=0.15, linewidth=0)
        ax.axhline(0.0, color=INK_2, linewidth=0.8, linestyle="--")
        ax.set_title(f"{o['counterpart']} ({o['type']}): brightening (mag)", fontsize=9)
        vis = [r for r in o["physics"]["visits"]
               if o["main_episode"]["start_mjd"] <= r["mjd"] <= o["main_episode"]["end_mjd"]
               and r["t_rest_k"] is not None]
        vy = mjd_to_year(np.array([r["mjd"] for r in vis]))
        for ax, key, title, log in ((row[1], "t_rest_k", "dust temperature of the excess (K)", False),
                                    (row[2], "l_dust_erg_s", "dust luminosity of the excess (erg/s)", True)):
            mid = np.array([r[key] for r in vis])
            lo = np.array([r[f"{key}_p16"] for r in vis])
            hi = np.array([r[f"{key}_p84"] for r in vis])
            ax.errorbar(vy, mid, yerr=[np.clip(mid - lo, 0, None), np.clip(hi - mid, 0, None)],
                        color=INK, marker="o", markersize=4, linewidth=1.2, elinewidth=0.9)
            if log:
                ax.set_yscale("log")
            ax.set_xlim(row[0].get_xlim())
            ax.set_title(title, fontsize=9)
        for ax in row:
            ax.tick_params(labelsize=8)
    axes[0, 0].legend(loc="upper left", fontsize=7.5)
    _title(fig, "Infrared flare candidates: variable component (shaded: episodes of the plan)")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    return _save(fig, path)


def plot_d_luminosity(objects_csv: Path, v1_csv: Path, power_json: Path, summary_json: Path,
                      path: Path) -> Path:
    """Stage 6B: hysteresis statistic D against log[(1+z) L^(1/2)] for the
    independent sample, with the toy-model expectations."""
    rows = [r for r in _rows(objects_csv) if r["hyst_measurable"] == "1" and r["x"]]
    v1 = [r for r in _rows(v1_csv) if r["sample"] == "main" and r["sl_hyst_measurable"] == "1"
          and r["loglbol"]]
    models = json.loads(power_json.read_text())["models"]
    summ = json.loads(summary_json.read_text())["d_vs_luminosity"]
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6), gridspec_kw={"width_ratios": [1.6, 1]})
    ax = axes[0]
    for key, color, label in (("model1", COLOR["w2"], "model 1: inter-band lag grows with (1+z) L^0.5"),
                              ("model2", GREEN, "model 2: re-formation delay, no lag")):
        b = models[key]["d_by_x_bin"]
        xm = [q["x_mid"] for q in b]
        ax.fill_between(xm, [q["d_p16"] for q in b], [q["d_p84"] for q in b], color=color,
                        alpha=0.15, linewidth=0)
        ax.plot(xm, [q["d_p50"] for q in b], color=color, linewidth=2, label=label)
    x1 = [np.log10(1 + _f(r["z"])) + 0.5 * _f(r["loglbol"]) for r in v1]
    ax.plot(x1, [_f(r["sl_d_hyst"]) for r in v1], linestyle="", marker="o", markersize=3.5,
            color=NEUTRAL, alpha=0.7, label=f"sample_v1 (stage 5b, {len(v1)} objects; already seen)")
    ax.errorbar([_f(r["x"]) for r in rows], [_f(r["d_hyst"]) for r in rows],
                yerr=[_f(r["sigma_d"]) for r in rows], linestyle="", marker="s", markersize=5,
                color=COLOR["w1"], markeredgecolor="white", elinewidth=0.9,
                label=f"sample_v2 ({len(rows)} objects)")
    ax.axhline(0.0, color=INK, linewidth=0.9, linestyle="--")
    ax.set_xlabel("x = log10(1+z) + 0.5 log10 L_bol (erg/s)")
    ax.set_ylabel("D = beta(declining) - beta(rising)")
    ax.legend(loc="upper left", fontsize=8)
    ax = axes[1]
    labels = []
    for i, (key, color, name) in enumerate((("model1", COLOR["w2"], "model 1"),
                                            ("model1_half", COLOR["w2"], "model 1, lag x 0.5"),
                                            ("model1_double", COLOR["w2"], "model 1, lag x 2"),
                                            ("model2", GREEN, "model 2"))):
        grid = sorted(int(k) for k in models[key]["by_n"])
        lo = np.interp(summ["n"], grid, [models[key]["by_n"][str(g)]["rho_p5"] for g in grid])
        mid = np.interp(summ["n"], grid, [models[key]["by_n"][str(g)]["rho_p50"] for g in grid])
        hi = np.interp(summ["n"], grid, [models[key]["by_n"][str(g)]["rho_p95"] for g in grid])
        ax.plot([lo, hi], [-i, -i], color=color, linewidth=2)
        ax.plot(mid, -i, marker="o", color=color, markersize=7, markeredgecolor="white")
        labels.append((-i, name))
    ax.axvline(summ["rho"], color=COLOR["w1"], linewidth=2)
    ax.axvspan(summ["rho_ci95"][0], summ["rho_ci95"][1], color=COLOR["w1"], alpha=0.12, linewidth=0)
    ax.axvline(summ["q2"], color=INK, linewidth=1.0, linestyle=":")
    ax.text(summ["q2"], 0.55, " q2", fontsize=8.5, color=INK)
    ax.text(summ["rho"], -3.75, f" observed {summ['rho']:+.2f}", fontsize=8.5, color=INK)
    ax.set_yticks([p for p, _ in labels])
    ax.set_yticklabels([t for _, t in labels])
    ax.set_ylim(-4.0, 0.8)
    ax.set_xlabel(f"Spearman rho of D with x (5-95% of toy samples, n = {summ['n']})")
    _title(fig, "Stage 6B: does the hysteresis statistic grow with luminosity in the independent sample?")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return _save(fig, path)
