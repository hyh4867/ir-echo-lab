# Future work and open questions

This file collects the open problems recorded at the end of stages 5b, 6A, 6B, 7 (with its supplement 7b), 8 and 9, and the candidates for later work. **The project stops at stage 9.** Nothing below is planned or in progress; the items are listed so that the limits of the results are visible and a later study has a starting point. Wording follows the rule used throughout the project: a result is a **candidate** until it is confirmed, and "not found" means **not found in this search**.

## Where stage 9 left the question

Stage 8 showed that the colour asymmetry is in the data. Stage 9 asked whether a lag of W2 behind W1 and the width of the response give its size object by object ([PREDICTIONS_9.md](PREDICTIONS_9.md), Korean; the result is in the README, "Conclusions in three layers"). The answer has two parts.

- A fixed linear response from W1 to W2, with a median lag of 40 d and 50 d in the observed frame, predicts 95% and 99% of the median asymmetry. No evidence was found for a component beyond it. The label by the committed rule is "depends on the model" (two cells "explained", two "undetermined", none "a part remains").
- With six-month visits this does not identify the cause. On synthetic samples a model with delayed re-formation only, and no lag between the bands, is fitted as an apparent lag of 15–25 d and is absorbed by the fixed kernel; a pure shift is not separated from a kernel of different width or tail; undelayed accretion-disc light in W1 leaves a fixed linear relation between the bands as well. The statement is limited to "consistent with an inter-band-lag model; not distinguished from re-formation".

The fitted lags are about 1.4 times the W2 − W1 lag difference of Lyu, Rieke & Smith (2019) and depend weakly on luminosity (ρ = +0.38). This hints that the fitted lag is not only a lag between the dust seen in the two bands, but the published relations differ by a factor of 2.5 at log L_bol = 44.5 (Lyu et al. against Mandal et al. 2024), and no published direct measurement of the W1-to-W2 lag was found in this search.

## Candidates for later work

Two candidates follow directly from what stage 9 could not separate. Neither has been started.

**1. Separating undelayed disc light in W1 with ZTF.** The design was written in stage 9 and not run. For objects with ZTF g and r light curves, average the optical flux in a window around each WISE visit, subtract c times that flux from W1, and repeat the transfer fit; c runs over a grid from zero to the largest value the spectral energy distribution allows, and the question is how the fitted lag and the residual change with c. To be counted first: how many objects keep at least 10 usable W2 visits after 2018, when ZTF starts, and for what fraction the seasonal gaps of the optical light curve miss the WISE visits. The method does not work where host light dilutes the optical variability, which is the case for most type 2 objects. An exploratory split in stage 9 found no difference in the median lag between types 1 and 2 (37.5 and 40 d in `sample_v1`, p = 0.68; 50 and 50 d in `sample_v2`, p = 0.18) and a larger asymmetry in type 1; this split alone does not decide anything, because the size of the asymmetry also depends on the amplitude and speed of the variability.

**2. Denser infrared data.** The separation of a fixed kernel from re-formation needs one of: (i) light curves at 3–5 μm with a cadence of weeks, in at least two bands, for a subset of the objects; (ii) an optical driver for the same objects, so that the response of each band can be fitted instead of the relation between the bands; (iii) one AGN with a sharp optical/UV drop followed by a long low state, in which a shortening of the hot-dust lag can be looked for after checking how much a fixed kernel driven by the observed light curve already explains (NGC 4151 and Mrk 590 are the starting cases). No such data set has been identified for the objects of this project. More wavelengths would also remove the single-temperature assumption and separate host light.

Smaller items left open by stage 9:

- The synthetic samples are rougher than the real light curves (residuals beyond 3σ in 3–4 and 9–10 objects of the two samples, against medians of 16–43 in the synthetic samples). The driver's strength on time scales shorter than the visit spacing cannot be set from WISE; optical structure functions could set it.
- The fit used a linear interpolation of W1. Another interpolation (a Gaussian process) was left out of the plan and was not tried. The fit was not repeated with the unTimely photometry.
- Why the lags exceed the Lyu et al. difference is not determined. The PDF of Mandal et al. (2024) and the ratio of the W2 and W1 lags in Tomar et al. (2025) were not checked in the papers; ADS was not queried.
- The mock for the lag–luminosity test of stage 8 comes closer to the data with a driver matched in amplitude (79–100 measurable objects against 70; 30–34 before) but still does not reproduce the sample, so its ρ test was not used. Host fraction and the shape of the response were not varied.

## (a) The response-width question

- Stage 6B could not separate the explanations, and stage 8 showed that its verdict against the pure-delay model depends on assuming no scatter in the lag–luminosity relation. A lag that depends weakly on luminosity, a difference in the W1 and W2 response that does not depend on luminosity, a re-formation delay, undelayed disc light in W1, and combinations of these all remain.
- Response functions for a fixed, radially extended torus exist in the literature (TORMAC; Almeyda et al. 2017, 2020). Only the abstracts were read; the per-band response-weighted delays and widths have not been taken from the papers, and whether the response functions are available as data was not checked. The kernels of stage 9 are parametric (a top-hat with a lag and a width) and are not TORMAC results.
- Not yet varied in the toy models: host fraction, shape of the dust response, a luminosity exponent of 0.4 as well as 0.5 in the lag relation.

## (b) Dependence on density and Eddington ratio: a third sample

- The dust-breathing hypothesis could not be tested: the temperature-response index was measurable in only 3 objects, none below the Eddington-ratio threshold. With two bands, six-month sampling and an error floor near 20 mmag, only objects whose variable flux is comparable to their baseline dust flux can be measured.
- A third sample, selected for large variability or with more wavelengths, is needed. It has to be defined, frozen and pre-registered in the same way as the first two, and it must not reuse objects already inspected.
- Eddington ratio, type, absorbing column and luminosity are tied together in a hard-X-ray-selected sample; the design has to break that, for example by matching in luminosity. The threshold λ_Edd = 0.02 is for gas with N_H = 10²² cm⁻², and a line-of-sight column does not measure the gas density where grains grow.
- Open from 6B: the weak positive trends (median D higher by 0.07 in the brighter half and by 0.04 in sample_v2) could come from luminosity, from redshift (shorter rest wavelength, (1+z) time dilation) or from the mix of types; the rank correlation of luminosity with redshift in sample_v2 is 0.78. Only 17 type 2 objects had both a measurable D and a luminosity.
- Re-formation within one AGN is a cleaner test than a sample trend: find a sharp optical/UV drop followed by a sustained low state, and test whether the hot-dust lag shortens or the W1 response recovers, after checking how much a fixed kernel driven by the observed light curve already explains. NGC 4151 and Mrk 590 are the starting cases.

## (c) What six-month visits allow

- The stacked cross-correlation estimator failed validation and was not applied. The lags of stage 9 come from a fit of W2 to linearly interpolated W1: the sample median is recovered within 10 d in the synthetic tests, the lag of a single object only to 21–68 d (rms).
- Both the hysteresis statistic and the slope move if W2 responds later. Archival or ground-based near-infrared monitoring that overlaps the WISE visits would be the cheapest addition; see candidate 2 above.

## (d) X-ray analysis of CGCG 312-012 (done; what the data cannot decide)

- CGCG 312-012 is a mid-infrared nuclear flare **candidate**; a report of it was **not found in this search** (SIMBAD references, two flare catalogues, web searches; TNS, ADS and ATel were not queried directly).
- Stage 7 and its post hoc supplement 7b are summarised in the README. The observed rise (a factor of 3.3 in the XRT 2–10 keV rate against 2010, 3.4 in the 14–50 keV flux against the BAT average) stands; the factor by which the intrinsic luminosity rose (4.1 or 6.0 with an absorbed power law, depending on how the XRT normalisation is read) and the state of the source in 2010 (unconstrained once a reflection component is allowed) are not determined. No further model was added: more parameters do not add counts to the 58-count spectrum of 2010.
- What observations would decide it:
  - A broad-band X-ray spectrum now (a soft X-ray telescope together with NuSTAR) would show whether the flux has returned to the BAT average and whether the reflection component and the 6.4 keV line followed the direct continuum. The second point tests the assumption of the reflection model that these stayed constant, and with it how the 2010 data should be read.
  - Hard X-ray monitoring after 2017 (this work used the BAT average up to the 157-month release) would date the start and the length of the brightening.
  - The state in 2010 cannot be re-observed. No other X-ray archive data from around 2010 were found in the HEASARC lists of Swift, XMM-Newton, Chandra and NuSTAR; other missions were not searched.
  - Optical spectroscopy (broad lines, coronal lines) and 3–5 μm photometry after 2024 remain open from stage 6A.
- Unexplained: the XRT/NuSTAR normalisation of 0.65–0.68 on the same day. The values of this object in the INTEGRAL 17-year survey and the Palermo BAT catalogue were not looked up; they would be a hard X-ray average independent of BAT.
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
- Which transfer-function model would be the accepted benchmark for a test like stage 9, and are per-band response functions for the WISE bands available as data?
- Is there a published direct measurement of the lag of W2 behind W1 (not the difference of two lags behind the optical) for Seyferts of this luminosity?
