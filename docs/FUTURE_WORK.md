# Future work and open questions

This file collects the next steps and the open problems recorded at the end of stages 5b, 6A, 6B, 7 (with its supplement 7b) and 8. Wording follows the rule used throughout the project: a result is a **candidate** until it is confirmed, and "not found" means **not found in this search**.

The two checks that came before everything else (a null test of D and a direct comparison without the shared peak) were done in stage 8; the result is in the README ("Conclusions in three layers") and the criteria are in [VALIDATION_8.md](VALIDATION_8.md). The colour asymmetry is in the data; its cause is open. The next stage follows from that.

## Stage 9 (planned): is the asymmetry explained by an inter-band lag alone?

**Question.** The direct comparison gives a median Δ(W1−W2) of −17.5 mmag (`sample_v1`) and −29.8 mmag (`sample_v2`) between the rising and the fading branch. Is this, object by object, what a lag of W2 behind W1 and the width of the response give, or is a part left over? The design note is [TRANSFER_DESIGN.md](TRANSFER_DESIGN.md) (Korean).

The stage has not been run. Its models, pass/fail criteria and interpretation rule are to be committed before any model is fitted to the real light curves, and the document will state that the results up to stage 8 had already been seen.

- **B0. Pre-registration.** Four fixed models in order of complexity: T0, the same response in both bands (no loop); T1, a pure shift of the W2 response; T2, fixed responses of different width and tail; T3, a response that changes with the past illumination (the inner edge recedes and returns late). Every model is compared with T0 directly, not only with its predecessor: the stage 7 rule that stepped from one model to the next missed the better-fitting one.
- **B1. Identifiability first.** With the real visit times and errors, simulate T0–T3 and measure how well a six-month cadence separates them and how accurately a lag is recovered. If they cannot be separated, that is the result, and the scope of B2 is reduced accordingly.
- **B2. Band-to-band transfer without a driver.** For each object fit F_W2(t) = a + b·(K ∗ F_W1)(t), using W1 as the input of W2, for a lag and a kernel width; the interpolation and the error propagation are to be stated and no extrapolation is allowed. Predict the direct offset δ of stage 8 from the fitted kernel and compare it with the measured δ per object and per sample: the ratio of predicted to measured, and the distribution and sign fraction of the residual. A residual consistent with zero reads "the asymmetry is explained by an inter-band lag"; otherwise its size and the objects that carry it are reported.
- **B3. Lag against luminosity.** Compare the fitted rest-frame lags with the relation of Lyu, Rieke & Smith (2019) and other published relations, with equations and sources given from the papers themselves; report how many objects have a well-determined lag and what selection that implies.
- **B4. Undelayed light in W1 (optional).** For objects with ZTF light curves, compare with a model that adds an immediate (accretion-disc) component to W1. If this is too heavy for the available machine, only the design is written.
- **B5. A better mock for the lag–luminosity test (optional).** The stage 8 mock samples had half the measurable objects of the real sample and a larger D. Find a driver that first reproduces the number of measurable objects and the distribution of D, then repeat the ρ test.

Expectation recorded in the design note: with six-month visits a pure shift (T1) and a history-dependent response (T3) were hard to tell apart in the stage 6B simulations, and T1 against T2 is expected to be harder.

## (a) What remains of the response-width question after stage 9

- Stage 6B could not separate the explanations, and stage 8 showed that its verdict against the pure-delay model depends on assuming no scatter in the lag–luminosity relation. A lag that depends weakly on luminosity, a difference in the W1 and W2 response that does not depend on luminosity, a re-formation delay, undelayed disc light in W1, and combinations of these all remain.
- Response functions for a fixed, radially extended torus exist in the literature (TORMAC; Almeyda et al. 2017, 2020). Only the abstracts were read; the per-band response-weighted delays and widths have not been taken from the papers, and whether the response functions are available as data was not checked. Until then the kernels of stage 9 are parametric and are to be described as such, not as TORMAC results.
- Not yet varied in the toy models: host fraction, shape of the dust response, a luminosity exponent of 0.4 as well as 0.5 in the lag relation.

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
- Which transfer-function model would be the accepted benchmark for stage 9, and are per-band response functions for the WISE bands available as data?
