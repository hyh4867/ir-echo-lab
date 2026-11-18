"""Stage 7b: table of old and revised values from the committed results (post hoc supplement).

    python scripts/run_stage7b.py

Reads results/s7_summary.json, s7_fits.json, s7b_light.json, s7b_fits.json;
writes results/s7b_summary.json and results/fig_s7b_supplement.png. No HEASoft.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from irecho import config, plots7, xray  # noqa: E402

RES = config.RESULTS_DIR


def main() -> None:
    s7 = json.loads((RES / "s7_summary.json").read_text())
    light = json.loads((RES / "s7b_light.json").read_text())
    fits = json.loads((RES / "s7b_fits.json").read_text())
    bat = light["bat_common_band"]
    out = {"note": "stage 7b, post hoc: computed after the stage 7 results were seen. The registered verdict of stage 7 is unchanged.",
           "registered_verdict_stage7": s7["registered_verdict"], "chosen_model_stage7": s7["chosen_model"]}
    out["flux_ratio_14_50"] = {
        "stage7_value_renamed": dict(bat["bat_all_channels_integrated_14_50"]["nustar_over_bat"],
                                     name="BAT fitted over all 8 channels (14-195 keV), flux integrated over 14-50 keV"),
        "bat_channels_1_4_free_gamma": bat["bat_ch1_4_free_gamma"]["nustar_over_bat"],
        "bat_channels_1_4_gamma_fixed_to_nustar": bat["bat_ch1_4_gamma_fixed_to_nustar"]["nustar_over_bat"],
        "joint_same_shape": {k: bat["joint_same_shape_14_50"][k] for k in ("log_ratio", "lo", "hi")}}
    out["c1_shape_tied"] = {m: {"delta1": s7["by_model"][m]["delta1"]["value"], "stat_90": s7["by_model"][m]["delta1"]["stat_90"],
                                "assumption": "BAT average has the spectral shape of the 2022 model; not a model-independent luminosity change"}
                            for m in ("M0", "M2")}
    variants = {"stage7_xrt_constant_1": {"delta2": s7["by_model"]["M0"]["delta2"]["value"], "delta2_90": s7["by_model"]["M0"]["delta2"]["stat_90"],
                                          "delta_nh": s7["by_model"]["M0"]["delta_nh"]["value"], "delta_nh_90": s7["by_model"]["M0"]["delta_nh"]["stat_90"]}}
    for key, name in (("c2_m0_const_both_xrt", "a_constant_on_both_xrt"), ("c2_m0_const_2022_only", "b_constant_on_2022_xrt_only")):
        r = fits[key]
        variants[name] = {"delta2": r["delta2"]["value"], "delta2_90": [r["delta2"]["lo"], r["delta2"]["hi"]],
                          "delta_nh": r["delta_nh"]["value"], "delta_nh_90": [r["delta_nh"]["lo"], r["delta_nh"]["hi"]],
                          "const_xrt": r["const_xrt22"], "const_xrt_90": r["const_xrt22_90"][:2], "statistic": r["statistic"]}
    for v in variants.values():
        v["label_delta2_stat_only"] = xray.luminosity_label(v["delta2"], *v["delta2_90"])
        v["label_delta_nh"] = xray.absorption_label(v["delta_nh"], *v["delta_nh_90"])
        v["factor"] = [10 ** v["delta2_90"][0], 10 ** v["delta2"], 10 ** v["delta2_90"][1]]
    out["m0_delta2_delta_nh"] = variants
    lo = min(v["delta2_90"][0] for v in variants.values())
    hi = max(v["delta2_90"][1] for v in variants.values())
    out["m0_delta2_envelope"] = {"lo": lo, "hi": hi, "factor": [10 ** lo, 10 ** hi],
                                 "note": "union of the 90% intervals of the readings of the XRT constant; M0 only"}
    m2 = fits["m2_xrt_free"]
    out["m2_xrt_constant_free"] = {k: m2[k] for k in ("statistic", "const_xrt22", "const_xrt22_90", "pegged", "values", "by_spectrum",
                                                      "reference_const_fixed", "nh_90", "gamma_90", "lum_90", "refl_norm_90", "line_norm_90")}
    out["observed"] = {"xrt_rate_ratio_2_10": s7["xrt_count_rates"]["ratio_2_10_2022_over_2010"],
                       "hardness_ratio": {k: light["hardness_ratio"][k] for k in ("xrt10", "xrt22")}}
    (RES / "s7b_summary.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    path = plots7.plot_supplement(RES / "s7b_summary.json", RES / "s7b_fits.json", RES / "s7_spectra_plot.json",
                                  RES / "fig_s7b_supplement.png")
    for k, v in out["flux_ratio_14_50"].items():
        print(k, "%+.3f [%+.3f, %+.3f]" % (v["log_ratio"], v["lo"], v["hi"]))
    for k, v in variants.items():
        print(k, "D2 %+.3f [%+.3f, %+.3f] %s | DN %+.3f [%+.3f, %+.3f] %s" % (v["delta2"], *v["delta2_90"], v["label_delta2_stat_only"],
                                                                              v["delta_nh"], *v["delta_nh_90"], v["label_delta_nh"]))
    print(path.relative_to(config.PROJECT_ROOT))


if __name__ == "__main__":
    main()
