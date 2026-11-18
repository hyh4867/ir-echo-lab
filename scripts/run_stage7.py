"""Stage 7: verdict table from the spectral fits (rules: docs/XRAY_PLAN.md section 6).

    python scripts/run_stage7.py

Needs only committed results (results/s7_fits.json, s7_regions.json,
s7_arf.json and, if present, s7_sim_validation.json and s7_mcmc.json); no
HEASoft. Output: results/s7_summary.json.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from irecho import config, xray  # noqa: E402

RES = config.RESULTS_DIR
INF = float("inf")


def _num(x, default):
    return default if x is None else float(x)


def comparison(fits: dict, model: str, sims: dict | None) -> dict:
    c1, c2 = fits["c1"][model], fits["c2"][model]
    out = {"model": model}
    # Delta1: 2022 against the BAT average
    d1 = c1["delta1"]
    lo, hi = _num(d1["lo"], -INF), _num(d1["hi"], INF)
    slo, shi = xray.add_systematic(d1["value"], lo, hi, xray.SYS_NUSTAR_BAT_DEX)
    out["delta1"] = {"value": d1["value"], "stat_90": [lo, hi], "with_systematic_90": [slo, shi],
                     "label": xray.luminosity_label(d1["value"], slo, shi)}
    # Delta2: 2022 against 2010
    d2 = c2["delta2"]
    lo, hi = _num(d2["lo"], -INF), _num(d2["hi"], INF)
    source = "profile of the fit statistic"
    if sims and sims.get("model") == model and not sims["pass"]["delta2"]:
        t = sims["truth"]["delta2"]
        p05, p95 = sims["p05_p95"]["delta2"]
        lo, hi = d2["value"] - (p95 - t), d2["value"] - (p05 - t)
        source = "distribution of the synthetic spectra (the profile interval failed the validation)"
    slo, shi = xray.add_systematic(d2["value"], lo, hi, xray.SYS_XRT_NUSTAR_DEX)
    out["delta2"] = {"value": d2["value"], "stat_90": [lo, hi], "with_systematic_90": [slo, shi],
                     "interval_source": source, "label": xray.luminosity_label(d2["value"], slo, shi)}
    # DeltaN
    dn = c2["delta_nh"]
    lo, hi = _num(dn["lo"], -INF), _num(dn["hi"], INF)
    source = "profile of the fit statistic"
    if sims and sims.get("model") == model and not sims["pass"]["delta_nh_coverage_only"]:
        t = sims["truth"]["delta_nh"]
        p05, p95 = sims["p05_p95"]["delta_nh"]
        lo, hi = dn["value"] - (p95 - t), dn["value"] - (p05 - t)
        source = "distribution of the synthetic spectra (the profile interval failed the validation)"
    iv = c2["intervals_90"]
    widths = {}
    for name in ("nh", "nh10"):
        a, b = iv[name][0], iv[name][1]
        widths[name] = math.log10(b / a) if a > 0 and b > 0 else INF
    out["nh_interval_width_dex"] = widths
    out["delta_nh"] = {"value": dn["value"], "stat_90": [lo, hi], "interval_source": source,
                       "delta_stat_if_tied": c2["delta_stat_nh_tied"]}
    out["values"] = {"nh_2022_1e22": c2["values"]["nh_1e22"], "nh_2022_90": iv["nh"][:2],
                     "nh_2010_1e22": c2["values"]["nh10_1e22"], "nh_2010_90": iv["nh10"][:2],
                     "gamma": c2["values"]["gamma"], "gamma_90": iv["gamma"][:2],
                     "log_lint_2022": c2["values"]["log_lint_2_10"], "log_lint_2022_90": iv["lum"][:2],
                     "log_lint_2010": c2["values"]["log_lint10_2_10"], "log_lint_2010_90": iv["lum10"][:2],
                     "cstat": c2["statistic"], "n_bins": c2["n_bins"], "n_free": c2["n_free"]}
    return out


def count_rates(regions: dict, arf: dict) -> dict:
    """Section 4.4 item 3b: XRT net count rates in the same region sizes."""
    out = {}
    for key in ("xrt10", "xrt22"):
        r, a = regions[key], arf[key]
        ratio = a["backscal_source"] / a["backscal_background"]
        rec = {"exposure_s": r["exposure_s"], "area_ratio": ratio}
        for band in ("0p3_2", "2_10"):
            rate, err = xray.net_rate(r[f"source_counts_{band}"], r[f"background_counts_{band}"], ratio, r["exposure_s"])
            rec[f"net_rate_{band}"], rec[f"net_rate_{band}_err"] = rate, err
            rec[f"net_counts_{band}"] = rate * r["exposure_s"]
        rec["net_counts_0p3_10"] = rec["net_counts_0p3_2"] + rec["net_counts_2_10"]
        # Stage 7 stored a hardness ratio here with a negative net soft rate cut to zero. Stage 7b
        # replaced it by a Poisson treatment with a lower limit (results/s7b_light.json, hardness_ratio).
        rec["mean_area_cm2"] = a["mean_area_cm2"]
        out[key] = rec
    a, b = out["xrt22"], out["xrt10"]
    ratio = a["net_rate_2_10"] / b["net_rate_2_10"]
    rel = math.hypot(a["net_rate_2_10_err"] / a["net_rate_2_10"], b["net_rate_2_10_err"] / b["net_rate_2_10"])
    area = a["mean_area_cm2"]["2-10"] / b["mean_area_cm2"]["2-10"]
    out["ratio_2_10_2022_over_2010"] = {"value": ratio, "err": ratio * rel, "log10": math.log10(ratio),
                                        "log10_err": rel / math.log(10), "effective_area_ratio": area,
                                        "significance_of_change_sigma": (a["net_rate_2_10"] - b["net_rate_2_10"])
                                        / math.hypot(a["net_rate_2_10_err"], b["net_rate_2_10_err"])}
    d = a["net_rate_0p3_2"] - b["net_rate_0p3_2"]
    out["difference_0p3_2_2022_minus_2010"] = {"value": d, "err": math.hypot(a["net_rate_0p3_2_err"], b["net_rate_0p3_2_err"])}
    return out


def main() -> None:
    if xray.LOG_LBOL_BASS is None:
        raise SystemExit("targets/sample_v1_props.csv has no BASS luminosity for the target (blank in the public copy). "
                         "Run scripts/restore_catalog_columns.py first; results/s7_summary.json was not changed.")
    fits = json.loads((RES / "s7_fits.json").read_text())
    regions = json.loads((RES / "s7_regions.json").read_text())
    arf = json.loads((RES / "s7_arf.json").read_text())
    sims = json.loads((RES / "s7_sim_validation.json").read_text()) if (RES / "s7_sim_validation.json").exists() else None
    mcmc = json.loads((RES / "s7_mcmc.json").read_text()) if (RES / "s7_mcmc.json").exists() else None
    chosen = xray.choose_model(fits["aic_2022"])
    assert chosen == fits["chosen_model"]
    rates = count_rates(regions, arf)
    out = {"object": xray.TARGET["counterpart"], "chosen_model": chosen, "aic_2022": fits["aic_2022"],
           "aic_lowest_model": min(fits["aic_2022"], key=fits["aic_2022"].get),
           "galactic_nh_1e22": fits["galactic_nh_1e22"], "xrt_count_rates": rates, "by_model": {}}
    enough = all(rates[k]["net_counts_0p3_10"] >= xray.MIN_XRT_COUNTS for k in ("xrt10", "xrt22"))
    for model in ("M0", "M1", "M2"):
        c = comparison(fits, model, sims)
        measurable = enough and all(w <= xray.MAX_NH_INTERVAL_DEX for w in c["nh_interval_width_dex"].values())
        dn = c["delta_nh"]
        dn["measurable"] = bool(measurable)
        dn["label"] = xray.absorption_label(dn["value"], dn["stat_90"][0], dn["stat_90"][1], measurable=measurable)
        c["verdict"] = xray.overall_verdict(c["delta1"]["label"], c["delta2"]["label"], dn["label"])
        c["energy_budget"] = xray.energy_budget(c["delta1"]["value"], fits["fits_2022"][model]["values"]["log_lint_2_10"])
        c["fit_2022"] = {"cstat": fits["fits_2022"][model]["statistic"], "n_bins": fits["fits_2022"][model]["n_bins"],
                         "n_free": fits["fits_2022"][model]["n_free"], "aic": fits["fits_2022"][model]["aic"],
                         "values": fits["fits_2022"][model]["values"],
                         "intervals_90": fits["fits_2022"][model]["intervals_90"],
                         "by_spectrum": fits["fits_2022"][model]["by_spectrum"]}
        out["by_model"][model] = c
    out["registered"] = out["by_model"][chosen]
    out["registered_verdict"] = out["registered"]["verdict"]
    out["verdict_same_for_all_models"] = len({out["by_model"][m]["verdict"]["verdict"] for m in out["by_model"]}) == 1
    sc = fits["simple_checks"]
    out["independent_checks"] = {
        "nustar_chi2_grouped_20": {m: {k: fits["sensitivity_nustar_chi2"][m][k] for k in (
            "statistic", "values", "intervals_90", "n_bins", "nustar_chi2", "nustar_bins", "n_free")}
            for m in ("M0", "M2")},
        "xrt_constant_free": fits["sensitivity_xrt_constant_free"],
        "flux_14_50_powerlaw": dict(sc["delta_14_50_powerlaw"], description_corrected_in_7b=(
            "BAT was fitted over all 8 channels (14-195 keV) and the flux integrated over 14-50 keV; it is not a fit "
            "of the 14-50 keV channels only. The fit of channels 1-4 is in results/s7b_light.json.")),
        "bat_powerlaw": sc["bat_powerlaw"],
        "nustar_14_50_powerlaw": sc["nustar_14_50_powerlaw"],
        "xrt_count_rate_ratio_2_10": rates["ratio_2_10_2022_over_2010"],
        "c3_reference_only": {m: {"values": fits["c3"][m]["values"], "intervals_90": fits["c3"][m]["intervals_90"],
                                  "statistic": fits["c3"][m]["statistic"]} for m in fits["c3"]},
        "mcmc": mcmc,
        "synthetic_validation": None if sims is None else {k: sims[k] for k in (
            "n_sims", "model", "truth", "recovered_median", "p05_p95", "coverage_90", "bias", "pass")}}
    (RES / "s7_summary.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    print("chosen", chosen, "| AIC", {k: round(v, 1) for k, v in fits["aic_2022"].items()})
    for model, c in out["by_model"].items():
        print(model, "D1 %+.3f [%+.3f, %+.3f] %s" % (c["delta1"]["value"], *c["delta1"]["with_systematic_90"], c["delta1"]["label"]),
              "| D2 %+.3f [%+.3f, %+.3f] %s" % (c["delta2"]["value"], *c["delta2"]["with_systematic_90"], c["delta2"]["label"]),
              "| DN %+.3f [%+.3f, %+.3f] %s" % (c["delta_nh"]["value"], *c["delta_nh"]["stat_90"], c["delta_nh"]["label"]),
              "|", c["verdict"]["verdict"], "| f_need %.3f %s" % (c["energy_budget"]["f_need"], c["energy_budget"]["label"]))
    r = rates["ratio_2_10_2022_over_2010"]
    print("XRT 2-10 keV rate ratio %.2f +- %.2f (%.1f sigma)" % (r["value"], r["err"], r["significance_of_change_sigma"]))


if __name__ == "__main__":
    main()
