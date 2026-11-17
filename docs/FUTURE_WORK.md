# Future work and open questions

This file collects the next steps and the open problems recorded at the end of stages 5b, 6A and 6B. Wording follows the rule used throughout the project: a result is a **candidate** until it is confirmed, and "not found" means **not found in this search**.

The first two items come before everything else: they test whether the hysteresis observation itself holds.

## 0. Checks in progress

1. **Null test of D.** Run light curves with no colour loop through the full pipeline: one case with ΔF_W1 = a ΔF_W2 exactly, and one with a single temperature T(t) shared by both bands. Apply the real cadence, the error covariance between the bands, the faintest-three baseline, the peak selection and the measurability cut. Report the distribution of D and the fraction with D > 0. W2 enters both the flux axis and the colour, so correlated errors could shift D if the two branches differ in signal-to-noise, flux range or distance from the peak.
2. **Direct comparison without the shared peak.** Where the two branches overlap in W2 excess flux, estimate F_W1,rise − F_W1,fade at matched W2 flux, with the W2 errors included and no extrapolation. A stable host flux cancels in the difference. Repeat with unTimely photometry. Report the common flux range, the number of visits per branch and the effect as Δ(W1−W2).

## (a) A model of the different W1 and W2 response widths

- Stage 6B could not separate the explanations. The baseline prediction of the pure-delay model (a strong positive correlation of D with luminosity) does not match, but a lag that depends weakly on luminosity, a difference in the W1 and W2 response that does not depend on luminosity, a re-formation delay, and combinations of these all remain.
- Build W1 and W2 transfer functions for a fixed, radially extended torus, starting from TORMAC (Almeyda et al. 2017, 2020). Compare in order: a pure shift; fixed kernels of different width and tail; kernels that change with the past illumination. Pass each through the real cadence and the same analysis, and compare D, the fraction with D > 0, the correlation with luminosity and the measurable fraction. Two symmetric kernels with the same centroid do not by themselves give a directed loop for a symmetric flare, so the rise and decay times of the driver matter.
- Rerun the lag model with a driver whose timescale and amplitude scale with black-hole mass and luminosity (observed variability timescales grow with mass; Burke et al. 2021), with scatter in the lag–luminosity relation, a luminosity exponent of 0.4 as well as 0.5, and the 210 → 94 → 70 selection applied to the simulated samples. Then check whether ρ = +0.13 still falls outside the model.
- Not yet varied in the toy models: driver amplitude and timescale, host fraction, shape of the dust response.

## (b) Dependence on density and Eddington ratio: a third sample

- The dust-breathing hypothesis could not be tested: the temperature-response index was measurable in only 3 objects, none below the Eddington-ratio threshold. With two bands, six-month sampling and an error floor near 20 mmag, only objects whose variable flux is comparable to their baseline dust flux can be measured.
- A third sample, selected for large variability or with more wavelengths, is needed. It has to be defined, frozen and pre-registered in the same way as the first two, and it must not reuse objects already inspected.
- Eddington ratio, type, absorbing column and luminosity are tied together in a hard-X-ray-selected sample; the design has to break that, for example by matching in luminosity. The threshold λ_Edd = 0.02 is for gas with N_H = 10²² cm⁻², and a line-of-sight column does not measure the gas density where grains grow.
- Open from 6B: the weak positive trends (median D higher by 0.07 in the brighter half and by 0.04 in sample_v2) could come from luminosity, from redshift (shorter rest wavelength, (1+z) time dilation) or from the mix of types; the rank correlation of luminosity with redshift in sample_v2 is 0.78. Only 17 type 2 objects had both a measurable D and a luminosity.
- Re-formation within one AGN is a cleaner test than a sample trend: find a sharp optical/UV drop followed by a sustained low state, and test whether the hot-dust lag shortens or the W1 response recovers, after checking how much a fixed kernel driven by the observed light curve already explains. NGC 4151 and Mrk 590 are the starting cases.

## (c) Denser infrared data

- Six-month visits cannot measure the lag of W2 behind W1, and both the hysteresis statistic and the slope move if W2 responds later. The stacked cross-correlation estimator failed validation and was not applied.
- Needed: infrared light curves with a cadence of weeks, in at least two bands, for a subset of the measurable objects; or archival or ground-based near-infrared monitoring that overlaps the WISE visits. No such data set has been identified yet for these objects.
- More wavelengths would also remove the single-temperature assumption and separate host light and, in type 1 objects, accretion-disc light in W1.

## (d) X-ray analysis of CGCG 312-012 (in progress)

- CGCG 312-012 is a mid-infrared nuclear flare **candidate**; a report of it was **not found in this search** (SIMBAD references, two flare catalogues, a web search; TNS, ADS and ATel were not queried directly).
- NuSTAR and Swift XRT observed it in 2022.33, close to the infrared peak, and a Swift observation from 2010 gives a comparison. The plan, with predictions and decision rules for each interpretation, was committed before any X-ray data were downloaded: [XRAY_PLAN.md](XRAY_PLAN.md). The analysis is running; no result is reported here.
- Infrared photometry alone cannot separate a tidal disruption inside an AGN from AGN variability; spectroscopy and X-rays are needed. The energy quoted for the event is a lower limit, because it had not ended when WISE stopped in 2024.
- Also open from 6A: whether the NGC 513 signal is real; the 2022 X-ray data of 2MASX J12005792+0648226; the 100 K difference in dust temperature between this work and the MIRONG list for that object, which comes from the emissivity assumption.

## Other open problems

- Q1 (slow rises): the interval ends at the 5% boundary; a sample about twice as large could decide it. Excluding quality-flagged objects moves the fraction from 7.1% to 4.7%, and whether the flag is independent of variability was not checked.
- Q2 and Q3 rest on 20 brightened objects.
- The baseline definition does not suit every object; for some the excess flux and energy are underestimated.
- A full rebuild on a second machine has not been tested.
- A systematic ADS search is needed before any statement about what has or has not been reported.

## Questions on which outside advice would help

- How does the same-flux colour statistic relate to the brightening/dimming amplitude asymmetry of Son, Kim & Ho (2026)? Could both come from the same process in the torus?
- Which transfer-function model would be the accepted benchmark for (a)?
