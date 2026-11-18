# Dust Colour Hysteresis in AGN

Repository name: `ir-echo-lab`.

A pre-registered test, on public WISE/NEOWISE photometry, of whether the mid-infrared colour of variable dust emission in active galactic nuclei (AGN) differs between brightening and fading at the same W2 flux. The pipeline is first checked on a known case (the 2020 flare of NGC 3786), then applied to two disjoint samples of hard-X-ray-selected AGN from the Swift-BAT 105-month catalogue. **Current status: the colour asymmetry is in the data. It survives a null test of the statistic, a direct comparison at matched W2 flux and a second photometry. Its cause is not determined: every process in which W1 responds before W2 gives the same sign, and the inter-band-lag model is no longer excluded. The next step, a test of whether an inter-band lag alone accounts for the size of the asymmetry, is described in [docs/FUTURE_WORK.md](docs/FUTURE_WORK.md).**

- Theory (English): [docs/THEORY.md](docs/THEORY.md)
- Pre-registration files: [docs/PREDICTIONS.md](docs/PREDICTIONS.md), [docs/PREDICTIONS_5b.md](docs/PREDICTIONS_5b.md), [docs/PREDICTIONS_6.md](docs/PREDICTIONS_6.md), [docs/FLARES_PLAN.md](docs/FLARES_PLAN.md), [docs/XRAY_PLAN.md](docs/XRAY_PLAN.md) (sections 1–8 registered, section 9 a post hoc supplement), [docs/VALIDATION_8.md](docs/VALIDATION_8.md); the order of the commits is in [docs/PREREGISTRATION_LOG.md](docs/PREREGISTRATION_LOG.md)
- Next steps and open questions: [docs/FUTURE_WORK.md](docs/FUTURE_WORK.md); design note for the transfer-function comparison (Korean, no computation yet): [docs/TRANSFER_DESIGN.md](docs/TRANSFER_DESIGN.md)
- Method, sample, data and literature notes (Korean): [docs/PIPELINE.md](docs/PIPELINE.md), [docs/SAMPLE.md](docs/SAMPLE.md), [docs/DATA.md](docs/DATA.md), [docs/LITERATURE.md](docs/LITERATURE.md). They refer to internal stage reports that are not part of this repository.
- A browser toy model of a dust echo, independent of the pipeline: [sim/](sim/)

The questions, predictions and analysis decisions are mine. The code was written and run with AI coding agents under my direction; I reviewed each stage and had the key numbers recomputed independently of the pipeline.

## Conclusions in three layers

**1. The colour asymmetry is in the data.** (Stage 8, an additional analysis whose criteria were committed before it was run; [docs/VALIDATION_8.md](docs/VALIDATION_8.md).)

- *Null test.* Synthetic light curves without a colour loop (W1 a linear function of W2, or one blackbody whose temperature alone changes) were given the real visit times, the real errors and the W1–W2 error correlation estimated per object (median r = 0.42 and 0.33 in the two samples), and passed through the same estimator, 200 times per case. The median D of the null samples is 0.000 (2.5–97.5%: −0.012 to +0.013 for `sample_v1`, −0.010 to +0.010 for `sample_v2`) and the fraction with D > 0 is 0.50. All 36 combinations (two samples, two null cases, nine variants of the estimator settings) meet the "unbiased" criterion. The observed medians, +0.228 and +0.268, are far above the 99th percentile of the null (+0.014, +0.012).
- *Direct comparison.* Without reference visits, excess fluxes, a shared peak temperature or a blackbody: within the range of W2 flux covered by both branches, W1 is brighter on the rising branch than on the fading branch in 180 of 211 objects (85%) in `sample_v1` and 194 of 200 (97%) in `sample_v2`. The median colour difference Δ(W1−W2), rising minus fading, is −17.5 mmag (95% interval −20.9 to −13.2) and −29.8 mmag (−33.4 to −25.8). On the loop-free null data the same estimator gives a positive fraction of 0.50 (0.42–0.57).
- *Second photometry.* With the unTimely catalogue the fractions are 170 of 218 (78%) and 173 of 192 (90%).
- A term for the scan direction changes the median by 0.0% and 1.7%. A linear time term, added after the plan was committed, leaves 179 of 211 and 189 of 200.

**2. The cause is not determined.** Any case in which W1 responds before W2 gives this sign: a lag between the bands from a fixed dust distribution, a W2 response that is broader, undelayed accretion-disc light in W1, or delayed re-formation of dust. The evidence for re-formation has not increased. The stage 6B verdict against the inter-band-lag model rested on an assumption: when the simulated driver follows the mass and luminosity scaling of MacLeod et al. (2010), the observed ρ = +0.129 lies outside the model (1.5–2.0% of mock samples at or below it) only if the lag–luminosity relation has no scatter; with 0.15 dex of scatter the fraction is 7.5–8.5% and with 0.30 dex it is 32%. The size of the real scatter was not found in the literature checked, so the lag model is no longer excluded. The mock samples do not reproduce the sample well either: they have 30–34 measurable objects where the data have 70, and a median D of 0.32–0.48 where the data have 0.27. The verdicts in the table below are left as they were recorded.

**3. CGCG 312-012 in X-rays.** (Stage 7, registered in [docs/XRAY_PLAN.md](docs/XRAY_PLAN.md) before the data were downloaded; "7b" is a post hoc supplement written after the results and an outside review, section 9 of the same file.)

- *Observed.* On 2022-05-01 the Swift XRT 2–10 keV count rate was 3.3 ± 0.6 times that of 2010-12-29 (6.5σ). The NuSTAR 14–50 keV flux of the same day was 3.4 times the Swift-BAT 2004–2017 average in the same band (+0.53 dex, 90% interval +0.44 to +0.63, BAT channels 1–4 only).
- *Model-dependent.* With an absorbed power law (M0, the model picked by the registered rule) the intrinsic 2–10 keV luminosity rose by a factor of 4.1 (2.7–6.2) if the fitted XRT/NuSTAR normalisation of 0.68 applies to both XRT observations, or 6.0 (4.1–8.7) if it applies only to 2022; the two fits have the same statistic and cannot be told apart. With a reflection component and an iron line (M2, AIC lower by 36.7; the registered sequential rule stopped before testing it) the direct continuum of 2010 is not constrained. A change of the absorbing column is undecided (+0.13 dex, −0.06 to +0.35).
- *Interpretation.* The X-ray source was in a high state 33 days after the brightest WISE visit; when the infrared flare actually peaked between the six-month visits is not known. This is what a brightening of the AGN itself would give, and the registered verdict was "brightening of the AGN supported; a change of obscuration not excluded". It remains a **candidate** dust echo. The factor by which the intrinsic luminosity rose and the state of the source in 2010 are not determined by these data; that is the final result of stage 7. The dust luminosity at the brightest visit is about 4% of a bolometric luminosity estimated from the fitted 2–10 keV luminosity with an assumed correction of 20: an energy requirement, not a measurement of the covering factor.

## What was done

1. **Pipeline and known-answer test.** WISE and NEOWISE single-exposure photometry (2010–2024, W1 3.4 μm and W2 4.6 μm) is binned into visits about six months apart. A per-visit zero point from 10 comparison stars per field (differential photometry) and an error floor for extended galaxies, measured against the unTimely catalogue, are applied. The pipeline was run on NGC 3786, whose 2020 flare was reported by Son et al. (2022).
2. **Detection limits.** Injection–recovery on 2,310 comparison stars (about 1.1 million trials).
3. **Two frozen samples.** `sample_v1`: 231 Swift-BAT AGN at z ≤ 0.040 (225 enter the statistics). `sample_v2`: 210 AGN at 0.040 < z ≤ 0.30, with no object in common. For each sample the light curves were downloaded, quality-checked without any target variability statistic, frozen and hashed (`data/freeze_v1.json`, `data/freeze_v2.json`); predictions and decision rules were committed; estimators were validated on toy models; only then were the frozen files opened.
4. **Analyses.** Descriptive predictions Q1–Q3, a "dust breathing" hypothesis against Eddington ratio, a hysteresis statistic D, a test of an inter-band-lag model against luminosity in the second sample, and a classification of four infrared flare candidates against criteria committed beforehand.
5. **X-ray follow-up of one flare candidate.** Archival NuSTAR and Swift XRT observations of CGCG 312-012 from 2022, a Swift XRT observation from 2010 and the Swift-BAT 157-month average spectrum were reprocessed with HEASoft and fitted with XSPEC, against predictions and decision rules committed before the download; a post hoc supplement followed an outside review.
6. **Validation of the hysteresis statistic.** A null test on loop-free light curves, a direct comparison of the two branches that shares none of the steps of D, and a re-test of the lag model with a driver that depends on mass and luminosity.

## Results

**Known-answer test (NGC 3786).** W1 brightened by 0.560 mag and W2 by 0.855 mag (Son et al. 2022: about 0.5 and 0.8 mag), with the onset between the 2020-05 and 2020-12 visits (published: mid-2020) and W1−W2 going from 0.21 to 0.50. A second photometry (unTimely) and an independent script that does not use the pipeline code agree; the flux excess agrees to 1.2% in W1 and 6.9% in W2.

**Injection–recovery.** Rises of 0.3 mag or more were recovered in essentially all trials; a 0.2 mag rise in 76–100% depending on rise time.

**Pre-registered verdicts.**

| Test | Result | Verdict by the committed rule |
|---|---|---|
| Q1: 5–20% of AGN show a slow rise ≥ 0.2 mag in both bands | 16/225 = 7.1% (95% interval 4.4–11.2%) | undetermined (interval crosses 5%) |
| Q2: brightening objects mostly redden in W1−W2 | 11/20 | undetermined |
| Q3: optical (ZTF) mostly brightens with the infrared | 9/15 | undetermined |
| Dust breathing vs. Eddington ratio | index measurable in 3 of 210 fitted objects, none below the threshold | undetermined (sample too small) |
| Hysteresis, sample_v1: D > 0 | 71 measurable objects, median D = +0.228 (95% interval 0.200–0.248); D above +2σ in 52, below −2σ in none | supported |
| Hysteresis, sample_v2 | 94 of 210 measurable, D > 0 in 94 of 94, median +0.268 (0.241–0.295) | sign reproduced |
| Lag model 1: D rises with x = log(1+z) + 0.5 log L_bol | Spearman ρ = +0.129 (95% bootstrap interval −0.120 to +0.368, n = 70); model 1 predicted +0.42 to +0.69 | "no model-1-strength correlation" supported (stage 8: this depends on assuming no scatter in the lag–luminosity relation; see layer 2 above) |

D is the difference between the temperature-response indices fitted to the brightening and the fading visits of one object, with a shared temperature at maximum W2 excess flux; D > 0 means a bluer (hotter) variable component while brightening than while fading at the same W2 excess flux.

**Additional analysis (stage 8), verdicts by criteria committed before it was run.**

| Test | Result | Label by the committed criterion |
|---|---|---|
| Null test of D (loop-free light curves, real cadence, errors and error correlation) | median D 0.000 in both samples, fraction with D > 0 0.50, all 36 combinations | unbiased; the observed medians are not explained by the null |
| Direct comparison at matched W2 flux | rising branch bluer in 180/211 (85%) and 194/200 (97%); median Δ(W1−W2) −17.5 and −29.8 mmag | asymmetry supported in both samples |
| Same with unTimely photometry | 170/218 (78%) and 173/192 (90%) | same direction |
| Scan-direction term | median changes by 0.0% and 1.7% | unrelated to scan direction |
| Lag model with a mass- and luminosity-dependent driver | observed ρ outside the model for no scatter in the lag–luminosity relation, inside for 0.15 and 0.30 dex | depends on the assumption |

D and the direct offset have the same sign in 67 of 71 and 94 of 94 objects where both are measured. One planned variant of the null test (meant to locate the source of a bias) was not run because no bias was found; the linear time term was not in the plan.

**What this does and does not show.** The colour asymmetry appears in two disjoint samples, is not produced by the estimator on loop-free data, and is seen without the reference visits, the shared peak temperature and the blackbody that D uses. Its cause is open: slow re-formation of dust was the motivating idea, but a W2 response that is delayed or broader than the W1 response from a fixed dust distribution, or undelayed disc light in W1, gives the same sign. Whether a lag and a kernel width fitted per object account for the size of the asymmetry is the question of the next stage.

**Infrared flare candidates (side result).** Four objects flagged by the sudden-flare detector were classified against pre-committed criteria. Two are consistent with an AGN flare and are known objects (2MASX J12005792+0648226, in the MIRONG list of Jiang et al. 2021; NGC 5273, a known changing-look AGN). For NGC 513 an instrumental artefact is not excluded. CGCG 312-012 (Seyfert 2) is a mid-infrared nuclear flare **candidate**: W1 +0.50 mag and W2 +0.82 mag at a peak near 2022.2, no optical change in ZTF or ASAS-SN at about 0.05–0.10 mag (3σ). **No report of it was found in this search**, which covered a SIMBAD reference list, two flare catalogues and a web search, but not TNS or ADS. Infrared photometry alone cannot separate a tidal disruption inside an AGN from AGN variability. The X-ray analysis of this object is summarised in layer 3 above; no paper analysing the 2022 NuSTAR and Swift observations was found in this search (titles of the SIMBAD reference list and web searches, not ADS).

Figures are in `results/` (for example `fig_ngc3786_lightcurve.png`, `fig_s5b_hysteresis.png`, `fig_s6b_d_luminosity.png`, `fig_s6a_flares.png`, `fig_s8_null.png`, `fig_s8_direct.png`, `fig_s8_examples.png`, `fig_s7_xray_spectra.png`, `fig_s7b_supplement.png`).

## Related work and how this statistic differs

The closest sample-level study is Son, Kim & Ho (2026, A&A 706, A122), from the group whose NGC 3786 paper this project started from. They compare ensemble structure functions of the brightening and dimming phases in ZTF and NEOWISE light curves of optically selected AGN at z = 0.15–0.4, find that bluer AGN vary more when brightening and redder AGN when dimming, and read this as hot dust being lost or gained. That is an asymmetry of variability *amplitude* across an ensemble. The statistic here is different: for each object, the W1−W2 colour (expressed as a colour temperature of the variable component) on the rising branch is compared with the fading branch *at the same W2 excess flux*, in a hard-X-ray-selected sample at lower redshift, and the result is then tested against a specified inter-band-lag model (a test that stage 8 found to depend on the assumed scatter of the lag–luminosity relation). Whether the two asymmetries come from the same process in the torus is an open question.

Other related work: Wu et al. (2026, arXiv:2611.04236) report a loop in the mid-infrared colour–magnitude diagram of one recurrent changing-look Seyfert, SDSS J234728.73+242745.7, read as the delayed response of an extended torus. Izviekova, Kompaniiets & Vavilova (2026, arXiv:2607.17815) find phase-resolved W1−W2 hysteresis in a mid-infrared-selected TDE candidate, UGC 11487, with W2 delayed by 224 ± 45 d and broader; this supports the delayed, broader W2 reading for a single transient. Wu et al. (2026, arXiv:2609.12390) study W1−W2 against W1 for 6,090 Seyfert 1s (bluer- and redder-when-brighter sources, host dilution) without separating branches. Within the literature checked, no study combines a hard-X-ray-selected sample, a same-flux comparison of rising and fading branches and a test of a specified lag model; the check used web searches and abstracts, not a systematic ADS query, so this is not a claim of priority.

## Reproducing

Requirements: Python ≥ 3.11 and the packages in `requirements.txt` (`pip install -r requirements.txt`, or `pip install -e ".[dev]"`).

```bash
python -m pytest                               # unit tests; tests that need the data cache or frozen files are skipped
python scripts/restore_catalog_columns.py      # download BAT 105-month and BASS tables, fill the blank catalogue columns
python scripts/fetch_data.py                   # stage 1-4 fields (NGC 3786 and controls) from IRSA
python scripts/run_analysis.py                 # known-answer test
python scripts/run_stage34.py                  # detectors, injection, physics, unTimely comparison
python scripts/fetch_sample.py                 # sample_v1 fields from IRSA (about 2 hours); IRECHO_SAMPLE=v2 for sample_v2
python scripts/fetch_untimely.py --sample
python scripts/run_qa.py                       # blind quality checks, error floors, frozen light curves
python scripts/run_stage5b.py                  # Q1-Q3, dust breathing, hysteresis (sample_v1)
python scripts/run_stage6b.py                  # hysteresis and luminosity (sample_v2); run with IRECHO_SAMPLE=v2
python scripts/run_stage8_noise.py             # stage 8: W1-W2 error correlation (needs the cache); again with IRECHO_SAMPLE=v2
python scripts/run_stage8.py null 2            # stage 8: null test (about 50 minutes with 2 workers)
python scripts/run_stage8.py direct            # stage 8: direct comparison
python scripts/run_stage8.py lag 2             # stage 8: lag model re-test (about 25 minutes)
python scripts/run_stage8.py select; python scripts/plot_stage8.py
python scripts/run_stage7.py                   # X-ray verdict table from the committed fit results (after restore_catalog_columns.py)
python scripts/run_stage7b.py                  # X-ray supplement: table and figure from the committed results
```

The X-ray fits themselves need HEASoft 6.37.1 (installed with conda into a git-ignored `.heasoft/` folder), the calibration files and the archive data (`scripts/fetch_xray.py`, 322 MB; `scripts/fetch_caldb_subset.py`); the steps are `scripts/xray_process.sh`, `scripts/xray_fit.py` and `scripts/xray_7b.py`, listed in docs/PIPELINE.md sections 15–17. The final fits are saved as XSPEC command files and logs in `results/s7b_xcm/`. None of the X-ray event files, spectra, response matrices, calibration files or the HEASoft environment is in this repository; `results/s7_spectra_plot.json` holds only the binned points drawn in the spectrum figure.

The exact order, options and run times are in [docs/PIPELINE.md](docs/PIPELINE.md). Raw IRSA downloads go to `data/cache/` and frozen light curves to `data/frozen/`; both are git-ignored and are **not** in this repository. `data/freeze_v1.json` and `data/freeze_v2.json` hold the SHA-256 of every frozen light curve, and the analysis scripts stop if a rebuilt file does not match. IRSA can reprocess or re-serve tables, so a fresh download is not guaranteed to reproduce the hashes; I have not tested a full rebuild on a second machine.

Tests: 287 pass in this repository as published and 58 are skipped (regression tests that need the local IRSA cache and the frozen light curves, and one that needs a catalogue value); with those present the suite runs all 345. The toy model in `sim/` has its own 9 tests (`cd sim && node --test model.test.mjs`).

## Data

No raw data are redistributed. The table lists every source, the terms I found (checked 2026-11-17 unless a date is given in docs/DATA.md) and the citation or acknowledgement it asks for.

| Data | Where from | Terms found | Citation / acknowledgement |
|---|---|---|---|
| AllWISE Multiepoch Photometry, NEOWISE-R Single Exposure Source Table, AllWISE Source Catalog | NASA/IPAC Infrared Science Archive (IRSA) | IRSA terms of use: "Most data served by IRSA is public, with no usage restrictions." | Wright et al. 2010, AJ 140, 1868; Mainzer et al. 2011, ApJ 731, 53; Mainzer et al. 2014, ApJ 792, 30; doi:10.26131/IRSA134, doi:10.26131/IRSA144. Acknowledgements below |
| unWISE Time-Domain Catalog (unTimely) | IRSA | same IRSA terms; no separate restriction found | Meisner et al. 2023, AJ 165, 36; doi:10.26131/IRSA580 |
| ZTF public light curves | IRSA | served without login; no restriction found on the IRSA ZTF page (ZTF's own policy document was not checked) | Bellm et al. 2019, PASP 131, 018002; Masci et al. 2019, PASP 131, 018003. Acknowledgement below |
| ASAS-SN Sky Patrol light curves | asas-sn.osu.edu | no licence text found; the site asks for citations | Shappee et al. 2014, ApJ 788, 48; Kochanek et al. 2017, PASP 129, 104502 |
| Swift-BAT 105-month catalogue | NASA GSFC, swift.gsfc.nasa.gov/results/bs105mon/ | no terms of use stated on the page | Oh et al. 2018, ApJS 235, 4 |
| Swift-BAT 157-month light curves, HEASARC archive records | NASA GSFC / HEASARC | no terms of use stated on the pages checked | HEASARC acknowledgement below |
| NuSTAR observation 60160365002, Swift XRT observations 00040946001 and 00081051001 | HEASARC archive (no login) | the HEASARC acknowledgement page states no restriction; the mission pages were not checked | Harrison et al. 2013, ApJ 770, 103; Burrows et al. 2005, Space Sci. Rev. 120, 165; HEASARC acknowledgement below |
| Swift-BAT 157-month average spectrum, response matrix and catalogue file | NASA GSFC, swift.gsfc.nasa.gov/bat_survey/bs157mon/ | no terms of use stated on the page | no paper for the 157-month release was checked; the page is cited as given |
| 2SXPS catalogue (one row looked up) | VizieR (IX/58) | not confirmed | Evans et al. 2020, ApJS 247, 54 |
| HEASoft 6.37.1 and calibration files (CALDB) | HEASARC | software and calibration, used locally and not redistributed | HEASARC acknowledgement below |
| BASS DR2 (table 9) and BASS DR1 X-ray catalogue | VizieR (J/ApJS/261/2, J/ApJS/233/17) | the BASS site says "© BASS collaboration. All rights reserved."; the VizieR licence page could not be opened | Koss et al. 2022, ApJS 261, 2; Ricci et al. 2017, ApJS 233, 17 |
| MIRONG catalogue (one row looked up) | VizieR (J/ApJS/252/32) | not confirmed | Jiang et al. 2021, ApJS 252, 32 |
| NED (names, positions and redshifts of the nine stage 1–4 objects) | ned.ipac.caltech.edu | no usage restriction found; NED's acknowledgement page did not load | "The NASA/IPAC Extragalactic Database (NED) is funded by the National Aeronautics and Space Administration and operated by the California Institute of Technology." |

**Catalogue values left out.** Because redistribution terms for the BAT and BASS tables were not confirmed (and the BASS site reserves all rights), the per-object values copied from them are blank in the tables and result files of this repository: BAT class, type and redshift; BASS type, redshift, distance, black-hole mass and its method, bolometric luminosity, Eddington ratio and column density; and the quantities computed directly from them per object (type group, threshold flags, the lag proxy x). Object names and catalogue IDs are kept. The same holds for the stage 7 and 8 files: type and redshift in `results/s8_direct_objects.csv`, the two objects outside types 1 and 2 in `results/s8_selection.json`, and the part of the X-ray energy budget that starts from the BASS luminosity in `results/s7_summary.json`. `scripts/restore_catalog_columns.py` downloads the tables from the Swift site and VizieR and fills all of these in place; with the files cached it reproduces my private copies byte for byte. Rows of three other tables looked up for CGCG 312-012 (the BASS X-ray catalogue, 2SXPS and the Swift-BAT 157-month catalogue file) are replaced by a note in `results/s7_literature.json` and `results/s7b_light.json`; `scripts/xray_literature.py` queries them again. Aggregate statistics (counts, medians, correlation coefficients) and figures are kept: a figure does not let the table be recovered, since points are not labelled with object names and positions can be read only to plotting accuracy.

Acknowledgements requested by the data providers:

> This research has made use of the NASA/IPAC Infrared Science Archive, which is funded by the National Aeronautics and Space Administration and operated by the California Institute of Technology.

> This publication makes use of data products from the Wide-field Infrared Survey Explorer, which is a joint project of the University of California, Los Angeles, and the Jet Propulsion Laboratory/California Institute of Technology, and NEOWISE, which is a project of the Jet Propulsion Laboratory/California Institute of Technology. WISE and NEOWISE are funded by the National Aeronautics and Space Administration.

> This publication makes use of data products from the Near-Earth Object Wide-field Infrared Survey Explorer (NEOWISE), which is a joint project of the Jet Propulsion Laboratory/California Institute of Technology and the University of California, Los Angeles. NEOWISE is funded by the National Aeronautics and Space Administration.

> ZTF: Supported by the National Science Foundation and a collaboration including current partners Caltech, IPAC, the Oskar Klein Center at Stockholm University, the University of Maryland, the University of California at Berkeley, the University of Wisconsin at Milwaukee, University of Warwick, Ruhr University, Cornell University, Northwestern University and Drexel University. Operations are conducted by COO, IPAC, and UW.

> This research has made use of data and/or software provided by the High Energy Astrophysics Science Archive Research Center (HEASARC), which is a service of the Astrophysics Science Division at NASA/GSFC.

## Limits

- Six-month sampling cannot measure an inter-band lag directly. A stacked W1–W2 cross-correlation estimator failed its validation on simulated data and was not applied. The direct comparison shows that the colours of the two branches differ; it does not say why.
- The variable component is treated as a single-temperature blackbody; host light is assumed constant; in type 1 AGN the W1 flux can include accretion-disc light.
- The baseline is the mean of the three faintest visits, which is not a true floor for objects that never faded in ten years.
- The "measurable" cut selects objects with large variability, and the selection 210 → 94 → 70 (measurable, with a bolometric luminosity) was not reproduced in the simulations.
- Eddington ratio, type, absorbing column and luminosity are correlated in the sample; stacked trends cannot be assigned to one of them. Catalogue values were used as published, without their uncertainties.
- Differences computed after the results were seen (for example a higher median D in the brighter half) are labelled post hoc in the result files and are not part of any verdict.
- Literature checks used web searches and abstracts; ADS was not queried systematically.
- The pre-registration protects against tuning on the target data. It does not make the estimator unbiased. The null test of stage 8 covers two kinds of loop-free light curves and starts from the visit level, not from single exposures; the error correlation it uses mixes the errors of two photometries.
- The stage 8 criteria were written after the results of stages 5b and 6B had been seen; the document lists what was already known.
- Lag-model re-test: the simulated samples have about half the measurable objects of the real one and a larger D, and the variability relations of MacLeod et al. (2010), fitted to SDSS Stripe 82 quasars, are applied to less luminous Seyferts.
- X-rays: the XRT/NuSTAR normalisation of 0.65–0.68 is unexplained (the good-time intervals of the two instruments overlap for 62 s of 6,963 s, so variability between them is not excluded). The synthetic-spectrum validation only recovered model M0 from M0, and its 90% interval for the photon index covered the truth in 85.5% of 200 trials (2.1σ below 90%). The cross-calibration terms of 0.06 and 0.04 dex are assumptions whose basis was not confirmed in the literature. Everything labelled 7b was done after the registered results were seen.

## References

Almeyda et al. 2017, ApJ 843, 3 · Almeyda et al. 2020, ApJ 891, 26 · Burke et al. 2021, Science 373, 789 · Burrows et al. 2005, Space Sci. Rev. 120, 165 · Evans et al. 2020, ApJS 247, 54 · Fabian, Vasudevan & Gandhi 2008, MNRAS 385, L43 · Harrison et al. 2013, ApJ 770, 103 · Izviekova, Kompaniiets & Vavilova 2026, arXiv:2607.17815 · Jiang et al. 2021, ApJS 252, 32 · Kokubo & Minezaki 2020, MNRAS 491, 4615 · Koshida et al. 2009, ApJ 700, L109 · Koss et al. 2022, ApJS 261, 1 and 2 · Lyu, Rieke & Smith 2019, ApJ 886, 33 · MacLeod et al. 2010, ApJ 721, 1014 · Madsen et al. 2017, AJ 153, 2 · Meisner et al. 2023, AJ 165, 36 · Oh et al. 2018, ApJS 235, 4 · Ricci et al. 2017, Nature 549, 488 · Ricci et al. 2017, ApJS 233, 17 · Son et al. 2022, ApJ 937, 3 · Son et al. 2026, arXiv:2604.17173 · Son, Kim & Ho 2026, A&A 706, A122 · Wu et al. 2026, arXiv:2609.12390 · Wu et al. 2026, arXiv:2611.04236. Volume and page numbers are as cited in the project notes and have not all been checked against ADS.

## Licence

Code: MIT (see [LICENSE](LICENSE)). The data belong to their providers under the terms above.
