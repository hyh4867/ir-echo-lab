# ir-echo-lab

A pre-registered test, on public WISE/NEOWISE photometry, of whether the mid-infrared colour of variable dust emission in active galactic nuclei (AGN) differs between brightening and fading at the same W2 flux. The pipeline is first checked on a known case (the 2020 flare of NGC 3786), then applied to two disjoint samples of hard-X-ray-selected AGN from the Swift-BAT 105-month catalogue. **Current status: a null test of the hysteresis statistic and a direct comparison at matched W2 flux are in progress; until they are done, the hysteresis result below is an observation whose cause, and whose freedom from estimator bias, are open.**

- Theory (English): [docs/THEORY.md](docs/THEORY.md)
- Pre-registration files: [docs/PREDICTIONS.md](docs/PREDICTIONS.md), [docs/PREDICTIONS_5b.md](docs/PREDICTIONS_5b.md), [docs/PREDICTIONS_6.md](docs/PREDICTIONS_6.md), [docs/FLARES_PLAN.md](docs/FLARES_PLAN.md), [docs/XRAY_PLAN.md](docs/XRAY_PLAN.md); the order of the commits is in [docs/PREREGISTRATION_LOG.md](docs/PREREGISTRATION_LOG.md)
- Next steps and open questions: [docs/FUTURE_WORK.md](docs/FUTURE_WORK.md)
- Method, sample, data and literature notes (Korean): [docs/PIPELINE.md](docs/PIPELINE.md), [docs/SAMPLE.md](docs/SAMPLE.md), [docs/DATA.md](docs/DATA.md), [docs/LITERATURE.md](docs/LITERATURE.md). They refer to internal stage reports that are not part of this repository.
- A browser toy model of a dust echo, independent of the pipeline: [sim/](sim/)

The questions, predictions and analysis decisions are mine. The code was written and run with AI coding agents under my direction; I reviewed each stage and had the key numbers recomputed independently of the pipeline.

## What was done

1. **Pipeline and known-answer test.** WISE and NEOWISE single-exposure photometry (2010–2024, W1 3.4 μm and W2 4.6 μm) is binned into visits about six months apart. A per-visit zero point from 10 comparison stars per field (differential photometry) and an error floor for extended galaxies, measured against the unTimely catalogue, are applied. The pipeline was run on NGC 3786, whose 2020 flare was reported by Son et al. (2022).
2. **Detection limits.** Injection–recovery on 2,310 comparison stars (about 1.1 million trials).
3. **Two frozen samples.** `sample_v1`: 231 Swift-BAT AGN at z ≤ 0.040 (225 enter the statistics). `sample_v2`: 210 AGN at 0.040 < z ≤ 0.30, with no object in common. For each sample the light curves were downloaded, quality-checked without any target variability statistic, frozen and hashed (`data/freeze_v1.json`, `data/freeze_v2.json`); predictions and decision rules were committed; estimators were validated on toy models; only then were the frozen files opened.
4. **Analyses.** Descriptive predictions Q1–Q3, a "dust breathing" hypothesis against Eddington ratio, a hysteresis statistic D, a test of an inter-band-lag model against luminosity in the second sample, and a classification of four infrared flare candidates against criteria committed beforehand.

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
| Lag model 1: D rises with x = log(1+z) + 0.5 log L_bol | Spearman ρ = +0.129 (95% bootstrap interval −0.120 to +0.368, n = 70); model 1 predicted +0.42 to +0.69 | "no model-1-strength correlation" supported |

D is the difference between the temperature-response indices fitted to the brightening and the fading visits of one object, with a shared temperature at maximum W2 excess flux; D > 0 means a bluer (hotter) variable component while brightening than while fading at the same W2 excess flux.

**What this does and does not show.** The colour asymmetry appears in two disjoint samples. Its cause is open. Slow re-formation of dust was the motivating idea, but a W2 response that is delayed or broader than the W1 response from a fixed dust distribution gives the same sign, and that case has not been modelled. The luminosity test disfavours only the specified pure-delay model, whose simulated driver had one timescale for every object. W2 enters both axes of the statistic, so correlated errors could bias D; the null test that would show this has not been run yet.

**Infrared flare candidates (side result).** Four objects flagged by the sudden-flare detector were classified against pre-committed criteria. Two are consistent with an AGN flare and are known objects (2MASX J12005792+0648226, in the MIRONG list of Jiang et al. 2021; NGC 5273, a known changing-look AGN). For NGC 513 an instrumental artefact is not excluded. CGCG 312-012 (Seyfert 2) is a mid-infrared nuclear flare **candidate**: W1 +0.50 mag and W2 +0.82 mag at a peak near 2022.2, no optical change in ZTF or ASAS-SN at about 0.05–0.10 mag (3σ). **No report of it was found in this search**, which covered a SIMBAD reference list, two flare catalogues and a web search, but not TNS or ADS. Infrared photometry alone cannot separate a tidal disruption inside an AGN from AGN variability. An analysis of archival NuSTAR and Swift data is planned in [docs/XRAY_PLAN.md](docs/XRAY_PLAN.md) and in progress.

Figures are in `results/` (for example `fig_ngc3786_lightcurve.png`, `fig_s5b_hysteresis.png`, `fig_s6b_d_luminosity.png`, `fig_s6a_flares.png`).

## Related work and how this statistic differs

The closest sample-level study is Son, Kim & Ho (2026, A&A 706, A122), from the group whose NGC 3786 paper this project started from. They compare ensemble structure functions of the brightening and dimming phases in ZTF and NEOWISE light curves of optically selected AGN at z = 0.15–0.4, find that bluer AGN vary more when brightening and redder AGN when dimming, and read this as hot dust being lost or gained. That is an asymmetry of variability *amplitude* across an ensemble. The statistic here is different: for each object, the W1−W2 colour (expressed as a colour temperature of the variable component) on the rising branch is compared with the fading branch *at the same W2 excess flux*, in a hard-X-ray-selected sample at lower redshift, and the result is then tested against a specified inter-band-lag model. Whether the two asymmetries come from the same process in the torus is an open question.

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
```

The exact order, options and run times are in [docs/PIPELINE.md](docs/PIPELINE.md). Raw IRSA downloads go to `data/cache/` and frozen light curves to `data/frozen/`; both are git-ignored and are **not** in this repository. `data/freeze_v1.json` and `data/freeze_v2.json` hold the SHA-256 of every frozen light curve, and the analysis scripts stop if a rebuilt file does not match. IRSA can reprocess or re-serve tables, so a fresh download is not guaranteed to reproduce the hashes; I have not tested a full rebuild on a second machine.

Tests: 238 pass in this repository as published and 56 are skipped (regression tests that need the local IRSA cache and the frozen light curves); with those present the suite runs all 294. The toy model in `sim/` has its own 9 tests (`cd sim && node --test model.test.mjs`).

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
| BASS DR2 (table 9) and BASS DR1 X-ray catalogue | VizieR (J/ApJS/261/2, J/ApJS/233/17) | the BASS site says "© BASS collaboration. All rights reserved."; the VizieR licence page could not be opened | Koss et al. 2022, ApJS 261, 2; Ricci et al. 2017, ApJS 233, 17 |
| MIRONG catalogue (one row looked up) | VizieR (J/ApJS/252/32) | not confirmed | Jiang et al. 2021, ApJS 252, 32 |
| NED (names, positions and redshifts of the nine stage 1–4 objects) | ned.ipac.caltech.edu | no usage restriction found; NED's acknowledgement page did not load | "The NASA/IPAC Extragalactic Database (NED) is funded by the National Aeronautics and Space Administration and operated by the California Institute of Technology." |

**Catalogue values left out.** Because redistribution terms for the BAT and BASS tables were not confirmed (and the BASS site reserves all rights), the per-object values copied from them are blank in the tables and result files of this repository: BAT class, type and redshift; BASS type, redshift, distance, black-hole mass and its method, bolometric luminosity, Eddington ratio and column density; and the quantities computed directly from them per object (type group, threshold flags, the lag proxy x). Object names and catalogue IDs are kept. `scripts/restore_catalog_columns.py` downloads the tables from the Swift site and VizieR and fills the columns in place; with the files cached it reproduces my private copies byte for byte. Aggregate statistics (counts, medians, correlation coefficients) and figures are kept: a figure does not let the table be recovered, since points are not labelled with object names and positions can be read only to plotting accuracy.

Acknowledgements requested by the data providers:

> This research has made use of the NASA/IPAC Infrared Science Archive, which is funded by the National Aeronautics and Space Administration and operated by the California Institute of Technology.

> This publication makes use of data products from the Wide-field Infrared Survey Explorer, which is a joint project of the University of California, Los Angeles, and the Jet Propulsion Laboratory/California Institute of Technology, and NEOWISE, which is a project of the Jet Propulsion Laboratory/California Institute of Technology. WISE and NEOWISE are funded by the National Aeronautics and Space Administration.

> This publication makes use of data products from the Near-Earth Object Wide-field Infrared Survey Explorer (NEOWISE), which is a joint project of the Jet Propulsion Laboratory/California Institute of Technology and the University of California, Los Angeles. NEOWISE is funded by the National Aeronautics and Space Administration.

> ZTF: Supported by the National Science Foundation and a collaboration including current partners Caltech, IPAC, the Oskar Klein Center at Stockholm University, the University of Maryland, the University of California at Berkeley, the University of Wisconsin at Milwaukee, University of Warwick, Ruhr University, Cornell University, Northwestern University and Drexel University. Operations are conducted by COO, IPAC, and UW.

> This research has made use of data and/or software provided by the High Energy Astrophysics Science Archive Research Center (HEASARC), which is a service of the Astrophysics Science Division at NASA/GSFC.

## Limits

- Six-month sampling cannot measure an inter-band lag directly. A stacked W1–W2 cross-correlation estimator failed its validation on simulated data and was not applied.
- The variable component is treated as a single-temperature blackbody; host light is assumed constant; in type 1 AGN the W1 flux can include accretion-disc light.
- The baseline is the mean of the three faintest visits, which is not a true floor for objects that never faded in ten years.
- The "measurable" cut selects objects with large variability, and the selection 210 → 94 → 70 (measurable, with a bolometric luminosity) was not reproduced in the simulations.
- Eddington ratio, type, absorbing column and luminosity are correlated in the sample; stacked trends cannot be assigned to one of them. Catalogue values were used as published, without their uncertainties.
- Differences computed after the results were seen (for example a higher median D in the brighter half) are labelled post hoc in the result files and are not part of any verdict.
- Literature checks used web searches and abstracts; ADS was not queried systematically.
- The pre-registration protects against tuning on the target data. It does not make the estimator unbiased; see the status line at the top.

## References

Burke et al. 2021, Science 373, 789 · Fabian, Vasudevan & Gandhi 2008, MNRAS 385, L43 · Izviekova, Kompaniiets & Vavilova 2026, arXiv:2607.17815 · Jiang et al. 2021, ApJS 252, 32 · Kokubo & Minezaki 2020, MNRAS 491, 4615 · Koshida et al. 2009, ApJ 700, L109 · Koss et al. 2022, ApJS 261, 1 and 2 · Lyu, Rieke & Smith 2019, ApJ 886, 33 · Meisner et al. 2023, AJ 165, 36 · Oh et al. 2018, ApJS 235, 4 · Ricci et al. 2017, Nature 549, 488 · Ricci et al. 2017, ApJS 233, 17 · Son et al. 2022, ApJ 937, 3 · Son et al. 2026, arXiv:2604.17173 · Son, Kim & Ho 2026, A&A 706, A122 · Wu et al. 2026, arXiv:2609.12390 · Wu et al. 2026, arXiv:2611.04236. Volume and page numbers are as cited in the project notes and have not all been checked against ADS.

## Licence

Code: MIT (see [LICENSE](LICENSE)). The data belong to their providers under the terms above.
