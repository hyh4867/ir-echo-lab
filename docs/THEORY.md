# THEORY: dust reverberation, sublimation radius, and what two WISE bands can measure

This document collects the physics and the statistics used in IR Echo Lab at the level of a fourth-year undergraduate. Every equation carries its source, and the notes say whether the source was read directly while writing this document or is quoted from memory of a standard text. All example numbers were produced by running the code in this repository (`src/irecho/theory.py`, `src/irecho/photometry.py`, `src/irecho/physics.py`, `scripts/run_stage34.py`) on 2026-11-12.

Contents: 1 the picture · 2 reverberation · 3 sublimation radius · 4 NGC 3786 example · 5 two-band blackbody temperature · 6 why the temperature is measured on the variable component · 7 luminosity and radius of the variable component · 8 differential photometry · 9 tests for a slow rise · 10 an error floor for galaxy nuclei

## 1. The picture

An active galactic nucleus (AGN) contains a supermassive black hole and an accretion disk that radiates mostly in the UV and optical. Farther out there is gas mixed with dust, usually called the dusty torus. Dust absorbs UV/optical light, heats up, and re-emits the energy as thermal infrared radiation. Cackett, Bentz & Kara (2021, Sec. 5) attribute the 2–10 μm emission of AGN to thermal emission from hot dust.

When the central source brightens, the dust gets hotter and the infrared brightens too. Light needs time to travel from the centre to the dust, so the infrared change arrives late. This delayed infrared response is a dust echo, and using the delay to measure the distance to the dust is dust reverberation mapping.

The WISE bands W1 (3.4 μm) and W2 (4.6 μm) are sensitive to dust near 1000 K (Section 5). That is why an echo from the nucleus can show up in WISE/NEOWISE light curves.

## 2. Reverberation: delay and smearing

### 2.1 Transfer function

Let ΔF_i(t) be the variable part of the central source and ΔF_r(t) the variable part of the reprocessed (here infrared) light. If reprocessing is linear,

```
ΔF_r(t) = ∫_0^τmax Ψ(τ) ΔF_i(t − τ) dτ                                (1)
```

Ψ(τ) is the transfer function (impulse response). It encodes the geometry of the reprocessing region: the infrared light curve is the central light curve delayed and smeared by Ψ.

- Source: Cackett, Bentz & Kara (2021) Eq. 1 (read directly). The framework goes back to Blandford & McKee (1982), as cited there.

### 2.2 Iso-delay surfaces and a thin shell

Consider dust at distance R from the centre, at angle θ from the line of sight. Light that goes from the centre to the dust and then to the observer arrives later than direct light by

```
τ = (R / c) (1 − cos θ)                                               (2)
```

Dust on the near side (θ = 0) has τ = 0 and dust on the far side (θ = π) has τ = 2R/c. For an optically thin spherical shell the delays are spread uniformly between 0 and 2R/c, so the mean delay is

```
⟨τ⟩ = R / c                                                           (3)
```

A lag multiplied by the speed of light is therefore a characteristic radius of the responding region. Equation (2) also says how fast a shell can respond: the whole shell has not answered until 2R/c after the central change.

- Source: Lyu, Rieke & Smith (2019) Appendix B, Eq. B8 (delay range 2R/c) and Eq. B9 (mean delay R/c), read directly. The general form (2) is the standard result of the reverberation-mapping review by Peterson (1993), which was not re-read for this document.

### 2.3 What this means for WISE

WISE/NEOWISE revisits a source about every six months. If Ψ is much narrower than six months the lag is unresolved, and the data only say between which two visits the brightening happened. Section 4 shows that NGC 3786 is such a case.

## 3. Sublimation radius: R_sub ∝ L^(1/2)

### 3.1 Derivation

A dust grain that gets too hot sublimates. Call the limiting temperature T_sub. Netzer (2015, Sec. 1.2.4) quotes T_sub ≃ 1800 K for ISM-like graphite grains and ≃ 1400 K for silicate grains.

For a spherical grain of radius a at distance R from a point source of luminosity L, radiative equilibrium (absorbed power = emitted power) reads

```
(L / 4πR²) · πa² · Q_abs = 4πa² · Q_em(T) · σT⁴                       (4)
```

with Q_abs and Q_em the wavelength-averaged absorption and emission efficiencies. Setting T = T_sub and solving for R,

```
R_sub = [ L Q_abs / (16π σ T_sub⁴ Q_em) ]^(1/2)  ∝  L^(1/2)           (5)
```

The key point is the square root: a source 100 times more luminous has its inner dust edge 10 times farther out.

- Equations (4)–(5) follow directly from the equilibrium condition. Netzer (2015) Eq. 1 is the numerical version and refers to Barvainis (1987).
- With Q_abs = Q_em = 1 (a perfect blackbody grain), L = 10^46 erg/s gives R = 0.27 pc at T = 1500 K and 0.19 pc at 1800 K (computed here from Eq. 5).

### 3.2 Real grains: Netzer (2015) Eq. 1

Small grains absorb UV well but emit poorly at infrared wavelengths longer than their own size (Q_em < Q_abs). They reach a given temperature farther out than a blackbody grain would, and the temperature dependence is steeper than T⁻². For a mean grain size a = 0.05 μm, Netzer (2015) Eq. 1 gives

```
⟨R_sub,C⟩  ≃ 0.5 L_46^(1/2) (1800 K / T_sub)^2.6 f(θ) pc    (graphite)   (6)
⟨R_sub,Si⟩ ≃ 1.3 L_46^(1/2) (1500 K / T_sub)^2.6 f(θ) pc    (silicate)   (7)
```

L_46 = L_AGN / 10^46 erg/s, and f(θ) describes the anisotropy of the central source (f = 1 if isotropic). The same section notes that the distance scales roughly as a⁻²: larger grains survive closer in.

- Source: Netzer (2015) Eq. 1 (read directly).
- For comparison, the review of dust echoes in tidal disruption events by van Velzen et al. (2021, Eq. 1) quotes R_sub = 0.15 pc × [L_45 / (a_0.1² T_1850^5.8)]^(1/2) for graphite-like grains, with a_0.1 the grain size in units of 0.1 μm (read directly). At L = 10^45 erg/s Eq. (6) gives 0.16 pc for 0.05 μm grains, and the van Velzen et al. expression gives 0.15 pc for 0.1 μm grains or 0.30 pc when scaled to 0.05 μm. The two prescriptions differ by a factor of about 2 at equal grain size, and the grain size itself is uncertain, so any R_sub is uncertain by a factor of a few before the luminosity is even considered.

### 3.3 Observations: the radius–luminosity relation

Lags measured by dust reverberation do follow L^(1/2). Cackett et al. (2021, Sec. 5) summarise:

- Dust lags follow a radius–luminosity relation close to R ∝ L^(1/2) over four orders of magnitude in luminosity.
- The inner wall is 2–3 times smaller than first expected for ISM-like grains; one interpretation is that large graphite grains, which resist sublimation best, form the inner wall.
- Dust lags are on average about four times the Hβ lag, so most of the broad-line region lies inside the dusty torus.

In the WISE bands, Lyu et al. (2019) found from 87 PG quasars

```
Δt_W1 / day = 10^(2.10±0.06) (L_AGN / 10^11 L_⊙)^(0.47±0.06)          (8)
Δt_W2 / day = 10^(2.20±0.06) (L_AGN / 10^11 L_⊙)^(0.45±0.05)          (9)
```

and size ratios R_K : R_W1 : R_W2 = 0.6 : 1 : 1.2. Longer wavelengths see cooler dust farther out.

- Source: Lyu et al. (2019) Eqs. 15, 16 and abstract (read directly).

## 4. Example: NGC 3786

Son et al. (2022, Sec. 3.3) give a bolometric luminosity L_bol = 10^42.89 erg/s before the flare and 10^43.35 erg/s after it, both estimated from the W1 magnitude. Their follow-up (Son et al. 2026, Sec. 4) derives L_bol from the broad Hα line at three epochs, 10^43.31, 10^43.24 and 10^43.34 erg/s without an extinction correction and 10^44.05, 10^44.39 and 10^44.64 erg/s with one (quoted uncertainty about 0.1 dex). Putting these into Eqs. (6)–(9):

| L_bol (erg/s) | source | R_sub,C / c | R_sub,Si / c | Eq. (8) W1 lag | Eq. (9) W2 lag |
|---|---|---|---|---|---|
| 10^42.89 | Son+2022, before the flare | 16.6 d | 43.1 d | 20.2 d | 27.4 d |
| 10^43.35 | Son+2022, after the flare | 28.2 d | 73.3 d | 33.2 d | 44.2 d |
| 10^43.31 | Son+2026, 2022-02-02, observed | 26.9 d | 70.0 d | | |
| 10^43.24 | Son+2026, 2023-06-10, observed | 24.8 d | 64.6 d | | |
| 10^43.34 | Son+2026, 2023-12-23, observed | 27.9 d | 72.4 d | | |
| 10^44.05 | Son+2026, 2022-02-02, extinction-corrected | 63.1 d | 164.0 d | | |
| 10^44.39 | Son+2026, 2023-06-10, extinction-corrected | 93.3 d | 242.6 d | | |
| 10^44.64 | Son+2026, 2023-12-23, extinction-corrected | 124.4 d | 323.6 d | | |

Computed with `theory.sublimation_radius_pc`, `theory.light_travel_days`, `theory.lyu2019_lag_days`; the rows from Son et al. (2026) are in `results/stage34_physics_ngc3786.json` (`sublimation_radius`).

How to read this:

- The expected lag is a few weeks to a few months, shorter than the six-month cadence for all but the largest luminosities. WISE alone cannot measure the lag of this object. That is why the onset of the flare is only bracketed between two visits.
- Equations (8)–(9) were calibrated on quasars with L_AGN ≈ 10^11–10^13 L_⊙. Using them for NGC 3786 (10^9.3–10^9.8 L_⊙) is an extrapolation.
- These are sizes of the inner dust edge. If dust extends farther out the echo lasts much longer; in Eq. (1) this is a wide Ψ.
- The literature luminosity spans 10^42.66–10^44.64 erg/s depending on the method and the extinction correction. Because R_sub ∝ L^(1/2), that is a factor of about 10 in radius. Any comparison with R_sub inherits this.

## 5. Two-band blackbody temperature

### 5.1 Equations

The specific intensity of a blackbody at temperature T is the Planck function

```
B_ν(T) = (2hν³ / c²) / [exp(hν / kT) − 1]                             (10)
```

(Rybicki & Lightman 1979, Ch. 1; a textbook equation, not re-checked for this document). The ratio of the flux densities that one region emits at two wavelengths λ1 < λ2 (frequencies ν1 > ν2) does not depend on the solid angle and is a function of temperature only:

```
F_ν(λ1) / F_ν(λ2) = (ν1 / ν2)³ · [exp(hν2 / kT) − 1] / [exp(hν1 / kT) − 1]      (11)
```

The ratio increases monotonically with T, so an observed ratio fixes T uniquely (`photometry.blackbody_temperature`, Brent's method). A temperature obtained this way is a colour temperature.

In the Wien regime (hν ≫ kT) the −1 in Eq. (11) can be dropped:

```
T ≈ (hc / k) (1/λ1 − 1/λ2) / ln[ (F_ν(λ2) / F_ν(λ1)) (λ2 / λ1)³ ]               (12)
```

with hc/k = 14387.8 μm K. At the isophotal wavelengths of W1 and W2, 3.3526 and 4.6028 μm (Jarrett et al. 2011), and T ≈ 1000 K, hν/kT ≈ 3–4, so Eq. (12) is only good to a few per cent. All numbers in this project solve Eq. (11) exactly.

### 5.2 Temperature and the W1−W2 colour

Converting the flux ratio into a Vega colour (zero points F0 in docs/PIPELINE.md Sec. 4),

```
W1 − W2 = −2.5 log10[ (F_ν(W1) / F_ν(W2)) · (F0_W2 / F0_W1) ]                    (13)
```

Values from `photometry.blackbody_flux_ratio`:

| T (K) | F_ν(W1)/F_ν(W2) | W1−W2 (mag) |
|---|---|---|
| 500 | 0.251 | 2.14 |
| 800 | 0.593 | 1.21 |
| 1000 | 0.782 | 0.91 |
| 1200 | 0.933 | 0.71 |
| 1500 | 1.105 | 0.53 |
| 1800 | 1.229 | 0.42 |

Stellar photospheres are on the Rayleigh–Jeans side at these wavelengths and have W1−W2 ≈ 0. Adding ~1000 K dust emission (W1−W2 ≈ 0.9) to a galaxy centre dominated by starlight therefore makes the total colour redder, which is what Son et al. (2022, Sec. 4) report: W1−W2 became redder as W1 brightened.

### 5.3 Limits of the estimate

1. Monochromatic approximation. W1 and W2 are broad bands; a steep 1000 K spectrum needs a colour correction (Wright et al. 2010, Table 1). Section 7.3 quantifies it: about +1% in temperature.
2. Dust is not a blackbody. With an emissivity Q_ν ∝ ν^β (β > 0) the true temperature for a given colour is lower.
3. One temperature. Real dust has a range of temperatures that two bands cannot separate.
4. Observed frame. The rest-frame temperature is (1+z) times higher; for NGC 3786 (z = 0.0091) that is 1%.

## 6. Why the temperature is measured on the variable component

### 6.1 The decomposition

What a WISE profile fit measures at the position of a nucleus is a sum,

```
F_obs(t) = F_host + F_AGN,0 + ΔF(t)                                    (14)
```

F_host is starlight of the galaxy inside the PSF (about 6″), F_AGN,0 is whatever the nucleus emitted in the reference state, and ΔF(t) is the change. The first two terms are unknown individually. Three things follow.

1. The colour of F_obs is not the colour of the dust. Starlight has W1−W2 ≈ 0 and dilutes the red dust emission. For NGC 3786 the total colour at the first flare visit is W1−W2 = 0.50, which Eq. (13) would turn into about 1570 K, while the excess alone has W1−W2 = 0.83 (1074 K). The first number mixes stars and dust and has no physical meaning as a dust temperature.
2. The amplitude in magnitudes depends on how much host light the photometry includes. A method that includes more host light reports a smaller Δm for the same event.
3. The difference ΔF(t) = F_obs(t) − F_obs(reference) cancels every constant term: host light, a constant aperture or PSF mismatch for an extended source, and the quiescent nucleus. It is the only quantity in Eq. (14) that the data determine without a model of the galaxy.

So every physical quantity in stage 3 is computed from ΔF in W1 and W2, never from the total flux. The same idea is used for dust echoes of tidal disruption events, where the pre-flare flux of the galaxy is subtracted before a blackbody is fitted (reviewed by van Velzen et al. 2021, Sec. 2; the specific subtraction recipe used here is ours).

### 6.2 An empirical check with a second photometry

The claim in point 3 can be tested, because NGC 3786 has two independent photometries of the same visits: this pipeline (median of single-exposure profile fits) and unTimely (a PSF fit on the coadded image of each visit; Meisner et al. 2023). The two disagree on the total flux and agree on the excess (`results/stage34_untimely.json`):

| | this pipeline | unTimely |
|---|---|---|
| W1 baseline (2014–2020-05, 13 visits) | 10.392 mag | 10.545 mag |
| W2 baseline | 10.189 mag | 10.287 mag |
| W1 brightening at 2020-12 | 0.560 mag | 0.619 mag |
| W2 brightening at 2020-12 | 0.855 mag | 0.866 mag |
| W1 excess flux at 2020-12 | 14.56 mJy | 14.39 mJy |
| W2 excess flux at 2020-12 | 17.29 mJy | 16.09 mJy |

unTimely assigns 13% (W1) and 9% (W2) less total flux to the source, so its amplitude in magnitudes is larger, as point 2 predicts. The excess fluxes differ by 1% in W1 and 7% in W2. The W2 difference is not explained and matters for the temperature: the unTimely excess gives 1144 K where the pipeline gives 1074 K. That 70 K is a fair estimate of the systematic uncertainty of the temperature, larger than the statistical interval.

### 6.3 What the subtraction assumes

- The host term is constant in the measurement, not only on the sky. For an extended source the profile-fit flux can depend on the scan direction, which alternates between visits; the subtraction does not remove that.
- ΔF is the change, not the dust emission. If the pre-existing dust (inside F_AGN,0) also changes temperature, ΔF is a difference of two spectra and need not look like a blackbody. Negative ΔF has no blackbody solution at all.
- The reference state is a choice. Section 7.3 varies it.

## 7. Luminosity and radius of the variable component

### 7.1 Equations

A sphere of radius R with uniform surface brightness B_ν(T) at distance D produces the flux density

```
F_ν = π B_ν(T) (R / D)²                                                (15)
```

(Rybicki & Lightman 1979, Sec. 1.3; textbook result, not re-checked). Integrating over frequency with ∫B_ν dν = σT⁴/π gives the Stefan–Boltzmann luminosity

```
L = 4π R² σ T⁴                                                         (16)
```

The code uses these in three steps (`physics.blackbody_quantities`). With the observed-frame colour temperature T_obs from Eq. (11),

```
N       = ΔF_ν / B_ν(T_obs)                    (an effective solid angle)
L_dust  = 4π D_L² · N σ T_obs⁴ / π                                     (17)
R_bb    = [ L_dust / (4π σ T_rest⁴) ]^(1/2),   T_rest = (1 + z) T_obs  (18)
```

D_L is the luminosity distance, 40.68 Mpc for z = 0.009083 (NED) in the flat ΛCDM cosmology of Son et al. (2022) (H0 = 67.4 km/s/Mpc, Ω_m = 0.315). The redshift terms change the results by 1–4%.

### 7.2 What R_bb means, and two other radii to compare it with

R_bb is the radius of a sphere that would emit L_dust as a perfect blackbody at T. Real dust covers only part of the sky seen from the centre and emits less efficiently than a blackbody, so a real emitting region with the same L_dust and T must be larger:

```
L_dust = 4π R² f σ T⁴  with  f ≤ 1   ⟹   R = R_bb / f^(1/2) ≥ R_bb     (19)
```

R_bb is therefore a lower limit on the size of the emitting region, not a measurement of where the dust is.

Two independent radii bracket it from the other side.

- Light travel time. By Eq. (2) a shell of radius R responds over 2R/c. If most of the brightening happened between two visits Δt apart, the responding dust is within roughly R ≲ c Δt. This is an order-of-magnitude bound, not a fit.
- Sublimation radius. Dust cannot survive inside R_sub of the current luminosity (Eqs. 6–7), so hot dust that responds to the nucleus should sit at or outside R_sub.

The ratio L_dust / L_bol is the fraction of the central luminosity that the variable dust component re-emits. For an echo it estimates the covering factor of the dust; van Velzen et al. (2021, Sec. 2.1) define the covering factor through the ratio of the energy radiated by the dust to that of the transient (read directly). Here it is a ratio of luminosities at one time, so it equals a covering factor only if the dust has had time to respond.

### 7.3 Results for NGC 3786

Input: differentially corrected visit fluxes (Section 8), baseline = median of the 13 visits from 2014-05-29 to 2020-05-10 (21.63 and 14.37 mJy in W1 and W2). Intervals are 16–84% of 4000 Monte Carlo draws (seed 20261112) over the visit and baseline flux errors; the baseline draws are shared between visits because that error is common to all of them. Source: `results/stage34_physics_ngc3786.json`.

| visit | ΔF W1, W2 (mJy) | T_rest (K) | L_dust (10^42 erg/s) | R_bb (light-days) |
|---|---|---|---|---|
| 2020-12-02 | 14.5, 17.4 | 1074 (1026–1125) | 3.50 (3.40–3.62) | 23.4 (21.6–25.4) |
| 2021-05-10 | 9.1, 12.3 | 961 (908–1019) | 2.28 (2.22–2.37) | 23.6 (21.2–26.5) |
| 2021-12-04 | 19.3, 24.7 | 1008 (974–1045) | 4.71 (4.62–4.83) | 30.9 (28.9–33.0) |
| 2022-05-10 | 26.4, 28.2 | 1214 (1177–1252) | 6.45 (6.28–6.63) | 24.9 (23.7–26.2) |
| 2022-12-03 | 21.2, 28.2 | 972 (945–1000) | 5.26 (5.18–5.35) | 35.1 (33.3–37.0) |
| 2023-05-11 | 11.4, 16.0 | 930 (888–975) | 2.92 (2.85–3.01) | 28.5 (26.1–31.3) |
| 2023-12-05 | 31.5, 34.4 | 1184 (1153–1215) | 7.63 (7.46–7.81) | 28.5 (27.3–29.8) |
| 2024-05-10 | 34.7, 41.6 | 1073 (1047–1099) | 8.35 (8.22–8.50) | 36.3 (34.8–37.9) |

Reading the table:

- The temperature stays between 930 and 1214 K (mean 1052 K, standard deviation 104 K) while L_dust changes by a factor of 3.7 between its lowest and highest visit. There is no monotonic heating or cooling trend. The two dips (2021-05, 2023-05) are the coolest visits, but the rank correlation between temperature and W1 excess flux over the eight visits is not significant (Spearman ρ = 0.67, p = 0.07; a post hoc check, `results/stage34_hypotheses.json`). Three hypotheses about this temperature history (temperature mixing, silicate, dust destruction) and the alternative of a brightening central source were tested with the decision rules of `docs/PREDICTIONS.md`; the verdicts are in README Section 6.2. The temperatures in this table had already been computed when those rules were written (PREDICTIONS.md Section 0).
- R_bb grows from 23 to 36 light-days, a factor 1.55 while L_dust grows by 2.39. A fit of log R_bb against log L_dust gives a slope of 0.21, shallower than the 0.5 of Eq. (5). This is a statement about a colour temperature and an equivalent radius, not a test of R_sub ∝ L^(1/2): R_bb is not R_sub, and by Eq. (18) its scatter is driven by the temperature.
- Compared with R_sub: for the observed broad-Hα luminosities of Son et al. (2026), graphite R_sub is 25–28 light-days and R_bb at the nearest visits is 1.02–1.15 times that; silicate R_sub is 65–72 light-days (ratio 0.39–0.44). For their extinction-corrected luminosities graphite R_sub is 63–124 light-days and the ratio is 0.23–0.49. By Eq. (19) R_bb ≤ R, so all of these are consistent with dust at or beyond the sublimation radius; none of them measures the radius. The nearest visits are 18–59 days from the spectroscopic epochs.
- L_dust / L_bol is 0.17–0.35 for the observed luminosities and 0.012–0.042 for the extinction-corrected ones. The first range would need a large covering factor for the varying dust alone; the second is a few per cent. The data here cannot choose between them, because the choice is the choice of L_bol.
- Light-travel bound: the first flagged visit is 206 days after the last unflagged one, so R ≲ 206 light-days for the dust that produced the initial jump. Together with Eq. (19): 23 ≲ R ≲ 206 light-days at onset.

Systematic checks, applied to every visit (same file):

| change | effect on T_rest | effect on R_bb |
|---|---|---|
| colour correction (Wright et al. 2010, Table 1, blackbody rows interpolated in log T) | +0.1% to +0.9% | −0.6% to +0.5% |
| baseline = faint half of the pre-flare visits (7 visits) | +0.1% to +0.7% | −0.4% to +0.2% |
| baseline = the last pre-flare visit only (2020-05-10) | +0.5% to +3.1% | −4.8% to −13.3% (L_dust −5% to −20%) |
| unTimely instead of the pipeline (onset visit only) | +6.5% | not computed |
| distance ±10% (a peculiar velocity of ~270 km/s at cz = 2720 km/s) | none | ±10% (L_dust ±20%) |

The statistical intervals in the table are a few per cent; the systematic terms are comparable or larger. The blackbody assumption itself (Section 5.3, items 2–3) is not in either list.

## 8. Differential photometry

### 8.1 Model

Let m_ik be the magnitude of star i in visit k. Assume

```
m_ik = M_i + o_k + ε_ik                                                (20)
```

M_i is the constant magnitude of the star, o_k is an error of the photometric zero point that is common to everything observed in visit k, and ε_ik is independent noise with variance σ_ik². Stage 2 showed that o_k exists: the residuals of 30 stars in the NGC 3786 field move together, with a standard deviation of 2.5 mmag in W1 and 7.8 mmag in W2 and a drift of about +2 mmag/yr in W2.

If o_k is ignored, the scatter of a constant star is σ² + var(o), larger than its error bars say, and reduced χ² exceeds 1. Worse for trend tests, a drift in o_k is a trend in every light curve of the field.

### 8.2 Estimator

For each star subtract its own median, r_ik = m_ik − median_k(m_ik), and estimate the offset of visit k from N stars as

```
ô_k = median_i(r_ik),    σ(ô_k) = sqrt(π/2) · s_k / sqrt(N)            (21)
```

with s_k the robust scatter (1.4826 × MAD) of the r_ik in that visit. The factor sqrt(π/2) is the asymptotic ratio of the standard error of a median to that of a mean for Gaussian noise. The corrected magnitude of a target observed in the same visits is

```
m'_k = m_k − ô_k,    σ'_k² = σ_k² + σ(ô_k)²                            (22)
```

The median, rather than a weighted mean, is used so that one variable star among the comparison stars cannot bias the offset (`tests/test_differential.py` checks this with a star that drifts by 0.5 mag).

Differential (ensemble) photometry against comparison stars in the same field is the standard way to remove frame-to-frame zero-point errors in time-series photometry; the classical formulation is Honeycutt (1992), and Lyu et al. (2019, Appendix A) use calibration stars to measure the stability of NEOWISE itself, finding σ ≈ 0.016 mag for NEOWISE and 0.029 mag for WISE single visits (the Lyu et al. numbers were read directly; Honeycutt 1992 was not re-read). Equations (20)–(22) are the simplest version: one additive term per visit, no colour or position terms.

### 8.3 Validation must be leave-one-out

If star j is included when ô_k is computed, part of its own noise is subtracted from it and its scatter after correction is too small. Each star is therefore corrected with offsets from the other N−1 stars (`differential.leave_one_out`). Result for the 30 reference stars, reduced χ² of a constant model, median over stars (`results/stage34_differential.json`):

| visit error model | W1 before → after | W2 before → after |
|---|---|---|
| catalogue errors only (`err_formal`) | 0.82 → 0.62 | 1.84 → 0.90 |
| in-visit scatter only (`err_scatter`) | 0.92 → 0.74 | 1.72 → 0.85 |
| larger of the two (`err`, stage 2) | 0.66 → 0.54 | 1.39 → 0.72 |

In W2 the correction removes the excess: the number of stars whose constant model is rejected at p < 0.01 drops from 13 to 1 (catalogue errors) and the visit-to-visit RMS from 11.5 to 8.1 mmag. In W1 reduced χ² was already below 1 and moves further below it. The correction is doing its job there too (RMS 5.7 → 5.3 mmag); what the numbers show is that the W1 visit errors are overestimated, by roughly 10–25% in σ depending on the error model. The stage-2 rule "take the larger of the two error estimates" overestimates in both bands.

The same leave-one-out test in the eight control fields (ten comparison stars each, in-visit scatter errors, same file, key `control_fields`) gives a different W1 picture:

| field | W1 before → after | W2 before → after |
|---|---|---|
| Mrk 573 | 1.54 → 0.94 | 4.11 → 0.84 |
| NGC 1386 | 1.30 → 0.85 | 2.77 → 1.02 |
| NGC 2273 | 1.03 → 1.09 | 1.80 → 0.78 |
| NGC 3081 | 1.30 → 1.03 | 1.77 → 0.87 |
| NGC 3393 | 0.96 → 0.71 | 2.08 → 0.99 |
| NGC 4941 | 1.21 → 0.88 | 2.88 → 1.15 |
| NGC 5347 | 0.81 → 0.66 | 3.03 → 0.74 |
| NGC 5728 | 1.33 → 1.05 | 3.93 → 1.15 |

After correction the field medians are 0.66–1.09 in W1 (median of the eight: 0.91) and 0.74–1.15 in W2 (0.93). Six of the eight fields had W1 above 1 before the correction, so the "W1 errors are overestimated" statement above holds for the NGC 3786 field and two others, not everywhere. The zero-point scatter also differs from field to field: the standard deviation of the W2 offsets is 7–17 mmag and single visits reach 48 and 57 mmag (Mrk 573 and NGC 5728 fields), against 7.8 and 21 mmag in the NGC 3786 field.

From stage 3 on the visit error is the in-visit scatter estimate plus σ(ô_k) in quadrature. It does not rely on catalogue errors, which are not trustworthy for extended sources (the profile-fit reduced χ² of NGC 3786 is 44 in W1), and after correction it is the model closest to 1 in the median of both bands. In the NGC 3786 field it still overestimates the noise of stars by roughly 10–15% in σ; in the control fields it ranges from overestimating by about 23% (reduced χ² 0.66) to underestimating by about 7% (1.15).

### 8.4 Assumptions that were checked, and one that was not

- The offset does not depend on brightness within W1 = 9.6–12.4: offsets from the brightest ten and the faintest ten stars differ by an RMS of 2.7 mmag in W1 and 6.8 mmag in W2, against 2.8 and 5.5 mmag expected from their errors.
- Ten stars are enough: offsets from ten randomly drawn stars, applied to the other twenty, give reduced χ² 0.95 → 0.74 (W1) and 1.59 → 0.80 (W2). The control fields use ten stars each.
- Not checked: that a galaxy shares the offset of the stars. An extended source is fitted with a point-source profile, and its error budget has terms (scan direction, seeing-like PSF changes) that a zero-point offset does not describe. The correction is necessary for galaxies but not sufficient.

## 9. Tests for a slow rise

### 9.1 Why a second detector

Detector (a) compares one visit with the earlier ones and flags a jump of at least 5σ and 0.2 mag in both bands. Every visit that is not flagged joins the baseline of the next one. A rise that never jumps by 5σ in a single step therefore raises the baseline scatter as it goes and is never flagged. The control galaxy NGC 2273 showed this in stage 2. A slow rise has to be tested on the light curve as a whole.

### 9.2 Linear and exponential trend models, and the F-test

Detector (b) fits each band with a constant and with "flat, then rising" models (`src/irecho/slowrise.py`):

```
constant        m(t) = a
hinge           m(t) = a + s · max(0, t − t_on)          t_on scanned over the visits      (23)
exponential     f(t) = a + c · exp[(t − t_end)/τ]         relative flux; τ = 0.5, 1, 2, 4 yr
```

The hinge with t_on at the first visit is a plain linear trend. The exponential describes an accelerating rise. Each rising model has one more linear parameter than the constant, so the two can be compared with the F-test for an additional term,

```
F = (χ²_const − χ²_model) / [ χ²_model / (n − 2) ],     p = P(F_{1, n−2} > F)             (24)
```

(Bevington & Robinson 2003, Sec. 11.4; textbook result, not re-read). The reason for F rather than Δχ² is that F divides by the scatter around the model. If all error bars are too large or too small by a common factor, F does not change. That matters here because the error bars of galaxies are not validated (Section 8.4).

### 9.3 Why p is not a false-alarm probability

Two conditions behind Eq. (24) fail.

1. The model was selected. t_on and τ are scanned and the best is kept, so the smallest p of many correlated trials is reported. The true chance of p ≤ 10⁻³ somewhere is larger than 10⁻³.
2. Residuals are not independent. AGN light curves are red noise: neighbouring visits are correlated, and a random walk produces stretches that look like trends. Optical quasar variability is commonly described as a damped random walk (Kelly et al. 2009; MacLeod et al. 2010; not re-read). For such a process a test that assumes white residuals overstates the significance of any trend.

For these reasons the threshold p ≤ 10⁻³ is used only as a ranking statistic, and the false-alarm rate is measured: on constant stars, where it tests point 1, and on the control galaxies, where it also includes point 2 but where a flag may be a real rise.

### 9.4 Cumulative change over several visits

Statistical significance is not enough. Ordinary Seyfert variability can produce a significant trend of a few hundredths of a magnitude. Detector (b) therefore also requires, in both bands,

```
fitted rise  Δm_fit ≥ 0.1 mag                                                               (25)
median(visits before t_on) − median(last 3 visits)  ≥ 0.1 mag
```

The second line is a model-free cumulative change. Using the median of the last three visits means that at least two of them must be bright: a single high point at the end is a flare, which detector (a) handles, not a slow rise.

### 9.5 Rank-based checks

Two distribution-free statistics are reported next to the fits, as a cross-check that does not depend on error bars at all.

- Mann–Kendall test (Mann 1945; Kendall 1975). S = Σ_{i<j} sign(m_j − m_i) counts concordant minus discordant pairs. For independent data without ties var(S) = n(n−1)(2n+5)/18 and Z = (S ∓ 1)/sqrt(var S) is approximately standard normal. It tests for any monotonic trend.
- Theil–Sen slope (Theil 1950; Sen 1968): the median of all pairwise slopes, robust to outliers.

Both assume independent data, so the red-noise caveat of Section 9.3 applies to them as well. (These four papers are the standard references and were not re-read.)

### 9.6 What the injection tests measure

Because none of the p-values is a calibrated probability, performance is measured by injection (`src/irecho/injection.py`): a brightening of known amplitude and rise time is added to real light curves of control galaxies and stars, and the fraction recovered is counted. The same segments without injection give the false-alarm rate. Results are in README Section 7.

## 10. An error floor for galaxy nuclei from two photometries

Added in stage 5a (2026-11-14). Code: `src/irecho/errorfloor.py`; rules: PIPELINE.md Section 12.4. This section gives the reasoning only. The floor was computed for the sample objects, but its effect on their light curves was not inspected in stage 5a.

### 10.1 The problem

The visit error used so far is the scatter of the single exposures inside a visit, divided by √n. It measures only noise that changes from exposure to exposure. A galaxy nucleus has error terms that are the same for all exposures of one visit and different in the next visit: how much host light the point-source profile fit picks up at that visit's scan angle, whether the source is split by deblending, and what is left of the zero-point offset after Section 8. These terms do not appear in the in-visit scatter, so the visit error of a galaxy is too small. Stage 3–4 found the visit-to-visit disagreement of two photometries of the same galaxy to be 15–40 mmag in W1 and 9–22 mmag in W2, against a quoted error of about 10 mmag.

### 10.2 Why the difference of two photometries isolates measurement error

Write the two measurements of visit k as

  m_A,k = μ_k + c_A + ε_A,k,  m_B,k = μ_k + c_B + ε_B,k,  (26)

where μ_k is the true (possibly variable) magnitude, c_A and c_B are constant offsets of the two methods (they pick up different fractions of the host), and ε are the measurement errors. The difference

  d_k = m_A,k − m_B,k = (c_A − c_B) + (ε_A,k − ε_B,k)  (27)

does not contain μ_k. However strongly the source varies, the visit-to-visit variance of d is

  s_d² = σ_A² + σ_B² − 2 Cov(ε_A, ε_B).  (28)

Here A is the pipeline (profile fit on single exposures, visit median, differentially corrected) and B is unTimely (PSF fit on the coadd of the same visit; Meisner et al. 2023). Both use the same exposures, so the photon noise of the two is correlated and the covariance term is positive: s_d underestimates σ_A² + σ_B² for the part of the noise the two share. That is acceptable here, because the shared part (photon noise) is exactly the part the in-visit scatter already measures. What s_d captures is the part that differs between two reasonable ways of measuring the same nucleus, which is the definition of the systematic term we are missing.

Equation (27) assumes that the offsets c_A, c_B do not depend on the brightness of the nucleus. If the two methods respond differently to a change of the nucleus-to-host ratio, a variable nucleus adds variance to d that is not measurement noise, and the floor is overestimated for strongly variable objects. This was not tested (it cannot be tested without looking at the variability).

### 10.3 The reference from stars

For a point source both methods are appropriate, and s_d is set by ordinary noise and by the zero-point errors of the two calibrations. This value, s_ref, is measured on the comparison stars, in two magnitude bins. The floor of a galaxy is its excess over the stars,

  σ_floor² = max(0, s_d² − s_ref²),  (29)

and the visit error becomes

  σ_visit² = σ_stat² + σ(ô_k)² + σ_floor².  (30)

Equation (29) assigns the whole excess to the pipeline measurement. Part of it belongs to unTimely, which also treats an extended nucleus as a point source, so σ_floor is biased high. Too large an error makes the detectors less sensitive but does not create false detections; the opposite bias would.

### 10.4 Precision of the estimate

s_d is estimated as 1.4826 × MAD, so that one bad coadd does not set the floor. The price is precision. A simulation of 200,000 Gaussian samples (seed 20261112, run while writing this section) gives:

| pairs n | mean of the estimate / true σ | relative standard deviation |
|---|---|---|
| 8 | 0.885 | 0.41 |
| 14 | 0.941 | 0.31 |
| 16 | 0.949 | 0.29 |

The relative scatter follows 1.16/√n, which is what the 37% efficiency of the MAD relative to the sample standard deviation implies. With about 14 pairs the scatter of a single object is therefore known to about 30%, and its floor to worse than that when s_d is close to s_ref, because the subtraction in Eq. (29) amplifies the error. The estimator is also biased low by 5–11% at these sample sizes. The bias affects the galaxy and the reference stars alike (they have similar numbers of pairs), so it partly cancels in Eq. (29); it was not corrected. Objects with fewer than eight pairs receive the sample median instead of their own value.

The estimate uses visits up to 2020 only (the neo7 release of unTimely ends there) and is applied to all visits.

## 11. Slope of the variable component: d log T / d log L

Added in stage 5b (2026-11-17), before the frozen light curves of the sample were opened. Code: `src/irecho/breath.py`; rules: PREDICTIONS_5b.md Sections 6–7 and 10.

### 11.1 Two limiting geometries

Write the hot-dust emission as a blackbody of temperature T that covers an effective area A:

```
L_dust = A σ T⁴.     (31)
```

Two things can change when the central luminosity L changes.

- Fixed area. The dust stays where it is and its temperature follows the heating. For grains in radiative equilibrium at a fixed radius, absorbed power ∝ L and emitted power ∝ T⁴ (blackbody grains), so T ∝ L^(1/4) and, with A constant, L_dust ∝ T⁴:

```
d log T / d log L_dust = 1/4.     (32)
```

- Receding boundary. The inner edge sits at the sublimation radius, R_sub ∝ L^(1/2) (Eq. 6). Dust at the edge is at the sublimation temperature whatever L is, so T is constant and the area changes, A ∝ R_sub² ∝ L:

```
d log T / d log L_dust = 0.     (33)
```

A one-parameter family joins the two: T ∝ L_dust^s and A ∝ L_dust^(1−4s), with s = 0 (boundary) and s = 1/4 (fixed area). Real grains do not radiate as blackbodies; with an emissivity ∝ ν^q the fixed-area exponent is 1/(4+q), about 0.18 for q = 1.6. That correction is not applied; the two reference values are 0 and 1/4 as in the prediction.

### 11.2 Why T and L must not be regressed on each other

With two bands, the temperature comes from the flux ratio and the luminosity from one flux and that temperature (Eqs. 15–18). Let γ_b = d ln B_ν / d ln T in band b (γ = x eˣ/(eˣ − 1), x = hν/kT). Small errors propagate as

```
δ ln T = (δ ln F₁ − δ ln F₂) / (γ₁ − γ₂),
δ ln L = δ ln F₂ + (4 − γ₂) δ ln T.     (34)
```

A W1 error alone moves a point in the (log L, log T) plane along a line of slope 1/(4 − γ₂); a W2 error alone moves it along a line of slope −1/(γ₁ − 4). At 1000 K (γ₁ = 4.35, γ₂ = 3.27) these slopes are +1.37 and −2.85. Both are far outside the range 0–0.25 that the hypothesis is about, so a regression of log T on log L of noisy visits measures the noise ellipse, not the dust. The same is true for any independent variable built from the fluxes that also give T.

### 11.3 The definition used: a curve in the flux–flux plane

The two fluxes of a visit are measured by different detectors and their errors are treated as independent. The fit is therefore done on the fluxes themselves.

1. Reference state A = (A₁, A₂): the weighted mean of the three visits with the smallest normalised sum F₁/median(F₁) + F₂/median(F₂). Using the sum, not one band, keeps the selection from favouring a noise dip in one band.
2. Variable component: e_b,k = F_b,k − A_b. Any constant contribution, host starlight or dust that does not vary, cancels in this difference.
3. Curve. The colour temperature of the variable component is a power law of its W2 flux,

```
T_c(e₂) = T_p (e₂ / e₂,p)^β,     e₁ = e₂ ρ(T_c),     ρ(T) = B_ν₁(T) / B_ν₂(T),     (35)
```

   and (T_p, β) minimise

```
χ² = Σ_k [e₁,k − e₂,k ρ(T_c(e₂,k))]² / [σ₁,k² + (de₁/de₂)² σ₂,k²] + (a₁/σ_A1)² + (a₂/σ_A2)²,     (36)
```

   with de₁/de₂ = ρ [1 + β(γ₁ − γ₂)]. The last two terms let the reference state move by (a₁, a₂) within its own uncertainty, so the error that all visits share through A is carried into σ_β.
4. No visit has a temperature measured from it. T enters only as a parameter of the curve, so the correlation of Section 11.2 does not arise. The independent variable is the W2 excess; W2 is chosen because the hot dust is brighter there relative to the host and to the accretion disc.
5. Conversion. With L_dust ∝ e₂ T⁴ / B_ν₂(T),

```
s = d log T / d log L_dust = β / [1 + β (4 − γ₂(T_p))].     (37)
```

   β = 0.306 at 1000 K corresponds to s = 0.25. Decisions are made on β; the conversion is monotonic.

An optical driver (ZTF) was not used as the independent variable: it exists for part of the sample and of the time span only, the nucleus is hidden in the optical in type 2 objects, and the catalogue photometry of an extended galaxy depends on seeing.

### 11.4 What the two geometries predict for β

- Receding boundary: e_b = (ΔA) B_ν_b(T_sub), so e₁/e₂ = ρ(T_sub) for every visit and **β = 0 exactly**. The result does not depend on the host, on constant dust, or on which visits define the reference state.
- Fixed area: e_b = A [B_ν_b(T) − B_ν_b(T_ref)]. This difference is not a blackbody at T. For a small change its colour is that of dB/dT, which is hotter than T_ref and does not depend on the size of the change; β → 0 as the amplitude → 0. Example for T_ref = 1000 K: when the W2 dust flux rises by 1%, 10%, 100% and 300% the dust is at 1003, 1030, 1263 and 1678 K and the colour temperature of the excess is 1379, 1404, 1637 and 2082 K. Over the first decade of excess T_c changes by 2%; only when the excess is comparable to the dust flux of the reference state does T_c climb.

So β is second order in the amplitude for heating in place, and a slope of 0.25 is not what this model gives for modest variability. What it gives depends on the ratio of the excess to the dust flux of the reference state, which is not known because the host fraction is not known. The smallest β is obtained when all of the reference flux is hot dust. That lower bound, β_A,min, is computed for each object from its own excess fluxes and fitted T_p, and the breathing index

```
r = β / β_A,min     (38)
```

is 0 for a receding boundary and ≥ 1 for heating in place (larger with more host). An object is "measurable" when σ_β ≤ β_A,min / 3.

### 11.5 Toy models

All toy models share one description of the inner edge. Let ℓ(t) be the luminosity seen by the dust in units of a reference value, R_in the inner radius in units of the sublimation radius at ℓ = 1.

```
R_sub = ℓ^(1/2),
R_in = R_sub                      when R_sub ≥ R_in (sublimation is immediate),
dR_in/dt = −(R_in − R_sub)/τ_re   when R_sub < R_in (re-formation),
T = T_sub (R_sub / R_in)^(1/2),     F_ν ∝ R_in² B_ν(T).     (39)
```

τ_re = 0 is the receding boundary. With τ_re = ∞ the edge stays at its largest past position and the dust heats and cools in place. The pure fixed-area model has R_in = 1 and T = T₀ ℓ^(1/4) with no upper limit. A two-layer model adds a fixed-area layer further out. A constant host of colour W1 − W2 = 0 can be added, and W2 can be made to respond later than W1. The driver is a damped random walk (time scale 1000 d) with an optional one-sided ramp; the dust sees it delayed and smoothed by 60 d.

### 11.6 What the toy tests showed

46,800 synthetic objects with the visit times and visit errors of the real sample (`results/s5b_validation.json`; PREDICTIONS_5b.md Section 10).

- Boundary model: median β between −0.015 and +0.002, β/σ_β has standard deviation 0.95–1.02, and 90–97% of measurable objects are classified on the "0 side".
- Fixed-area model: median r = 0.97–1.04 without host, 1.35–1.54 with half of W1 from the host.
- Measurability needs large changes. With a driver rms of 0.3 mag, 7% of boundary objects and 2% of fixed-area objects pass; with 0.15 mag none do.
- A 50% mixture is classified on the "0 side" in 59–82% of measurable cases. "0 side" means that a component of changing area is present, not that nothing is heated in place.
- **Inter-band lag.** If W2 responds 10 days later than W1, β is pushed to −0.03 to −0.04 in ramp-like light curves and 23–40% of measurable boundary objects are misclassified. The hysteresis statistic of Section 13 becomes +0.05 to +0.10, the size expected from slow re-formation. Lyu, Rieke & Smith (2019, abstract) give size ratios R_W1 : R_W2 = 1 : 1.2 from reverberation of PG quasars, so a W1 lag of 50 days or more implies a W2 − W1 lag of at least 10 days. With visits six months apart this lag cannot be measured here.

## 12. The effective Eddington limit for dusty gas

Added in stage 5b. This is the reasoning behind the threshold fixed in PREDICTIONS_5b.md Section 6.1. The derivation below is a textbook-level argument written for this document; the published curves come from photoionisation calculations that were not reproduced.

### 12.1 Classical limit

A proton–electron pair at distance r from a source of luminosity L feels the radiation force σ_T L / (4π r² c) on the electron and gravity G M m_p / r². They balance at

```
L_Edd = 4π G M m_p c / σ_T ≃ 1.26 × 10³⁸ (M / M_⊙) erg s⁻¹,     λ_Edd = L_bol / L_Edd.     (40)
```

### 12.2 Dusty gas, optically thin

If the gas contains dust that is coupled to it, the cross-section per hydrogen nucleus is no longer σ_T but a flux-weighted dust cross-section σ_d that is far larger in the ultraviolet. The same balance gives

```
L_Edd,eff = L_Edd / A,     A = σ_d / σ_T,     λ_eff = 1 / A.     (41)
```

Fabian, Vasudevan & Gandhi (2008, p. 1 abstract and p. 4) give A between unity and about 500, with the maximum "σd/σT ∼ 500" for gas that the radiation passes through unattenuated. With σ_T = 6.65 × 10⁻²⁵ cm² that maximum is σ_d ≈ 3.3 × 10⁻²² cm² per hydrogen nucleus, so the gas becomes optically thick to the driving radiation at N_H ≈ 1/σ_d ≈ 3 × 10²¹ cm⁻². Below that column the limit is λ_eff ≈ 1/500 = 0.002.

### 12.3 Dusty gas, optically thick

A column with N_H ≫ 1/σ_d absorbs the momentum of the radiation once. If a fraction f_abs of L is absorbed and the re-emitted infrared escapes, the force per unit area is f_abs L / (4π r² c) and the weight per unit area is G M m_p N_H / r² (taking a mass m_p per hydrogen nucleus). Balance gives

```
λ_eff = N_H σ_T / f_abs = (N_H / 1.5 × 10²⁴ cm⁻²) / f_abs.     (42)
```

The boost factor is A = f_abs / (N_H σ_T): it falls as 1/N_H, because a thicker column is heavier but cannot absorb more momentum. A reaches unity for f_abs = 1 at N_H = 1.5 × 10²⁴ cm⁻², the Compton-thick column; this is why Compton-thick gas is never above its limit at λ_Edd < 1 (Fabian et al. 2008, p. 2).

The published dividing line is N_H = 5 × 10²³ λ cm⁻² for N_H > 5 × 10²¹ cm⁻² (Fabian et al. 2008, abstract and p. 3), i.e. λ_eff = N_H / 5 × 10²³ cm⁻². Compared with Eq. (42) it corresponds to f_abs ≈ 1/3. The paper computes the absorbed part of the spectrum with CLOUDY for two AGN spectral energy distributions, "from absorption only (the input spectrum minus the transmitted one)", with the gas optically thin to its own infrared emission and no trapping (p. 2). At N_H = 10²² cm⁻² the line gives λ_eff = 0.02, the value Ricci et al. (2017, p. 6) quote from a later calculation for interstellar grain abundances, together with 0.15 at 10²³ cm⁻² (the line gives 0.2). Trapping of the reprocessed infrared lowers the limit further at high column (Ishibashi et al. 2018, abstract only).

The threshold adopted for the hypothesis is the single value λ_Edd = 0.02, the limit for N_H = 10²² cm⁻².

### 12.4 What the limit does and does not say about the hot-dust boundary

Two different things can move the inner edge of the dust.

- Sublimation. The edge is where grains reach their sublimation temperature; it moves out as L^(1/2) whatever λ_Edd is, on the sublimation time of Section 13.1.
- Radiation pressure. Above λ_eff, dusty gas is accelerated outward. The time to move it a distance comparable to its radius is at least the dynamical time R / v_K. For M = 10^7.7 M_⊙ (the sample median) and R = 0.03 pc, v_K ≈ 2700 km s⁻¹ and R / v_K ≈ 11 yr; the net acceleration is (λ_Edd/λ_eff − 1) times gravity, so the time is of this order unless λ_Edd ≫ λ_eff.

The light curves cover 10 years with visits every 6 months. Motion of the gas by radiation pressure is therefore too slow to follow the luminosity changes within one light curve. What the slope of Section 11 can respond to within a light curve is sublimation and re-formation. The hypothesis H-breath can still hold in the sense that objects above the limit have had their dusty gas cleared to the point where the remaining hot dust sits at the sublimation edge, while objects below it keep dust at larger radii that is only heated; but that is a statement about the long-term state, and this section gives no derivation of it.

## 13. Sublimation and re-formation times, and when hysteresis appears

### 13.1 Sublimation

Waxman & Draine (2000, Eqs. 10 and 12; read directly) write the sublimation rate of a grain of radius a as

```
da/dt = −(m/ρ)^(1/3) ν₀ exp(−B/kT),
t_surv(T) = a / |da/dt| = 7.7 a₋₅ exp[7 × 10⁴ K (1/T − 1/2300 K)] s,     (43)
```

with ν₀ = 10¹⁵ s⁻¹, B/k = 7 × 10⁴ K, ρ/m = 10²³ cm⁻³ as representative values for refractory grains and a₋₅ = a / 0.1 μm. For a 0.1 μm grain this gives 2.7 yr at 1500 K, 54 d at 1600 K, 4 d at 1700 K and 10 h at 1800 K. The temperature at which a grain lasts one visit interval (0.5 yr) is 1557 K. The dependence on T is so steep that "the sublimation temperature" changes by about 100 K when the time scale changes by a factor of 100, and a rise of the luminosity by a few tens of per cent at the edge removes grains within days to weeks. On the 6-month sampling of these light curves sublimation is immediate.

### 13.2 Re-formation

There is no equally simple expression for making grains. A lower bound on the time is growth of existing grain cores by sticking of gas-phase atoms, da/dt = S n_X m_X v̄ / (4 ρ_gr) with v̄ = (8kT/π m_X)^(1/2). For carbon with C/H = 2.7 × 10⁻⁴ all in the gas (an assumed solar value, not checked against a source here), sticking probability S = 1, gas at 1000 K and ρ_gr = 2.26 g cm⁻³,

```
da/dt ≃ 7.9 × 10⁻²³ n_H cm s⁻¹,     t_grow(0.1 μm) ≃ 4 × 10⁹ yr / (n_H / cm⁻³).     (44)
```

That is 0.4 yr at n_H = 10¹⁰ cm⁻³ and 4 yr at 10⁹ cm⁻³. This is an estimate made for this document. It ignores nucleation of new cores, which is slower, and the destruction of molecules by the ultraviolet field.

Observed values quoted in abstracts: about one year for dust formation in NGC 4151 (Koshida et al. 2009), less than four years for replenishment of the innermost dust in Mrk 590 (Kokubo & Minezaki 2020). Both are single objects.

### 13.3 Conditions for hysteresis

With Eq. (39), a slope that differs between brightening and fading needs all of the following.

1. Sublimation faster than the sampling: satisfied (Section 13.1).
2. Re-formation slower than the sampling, τ_re ≳ 0.5 yr. With Eq. (44) this means n_H ≲ 8 × 10⁹ cm⁻³ at the edge.
3. The luminosity must exceed its previous maximum within a time τ_re. Below the previous maximum the dust is heated and cooled in place on both branches and the track is reversible. If τ_re is much longer than the light curve the edge sits at its historical maximum and almost all variability is reversible; hysteresis is then weak again. The toy tests give the largest effect for τ_re between 1 and 5 years.
4. The two bands must respond to the same dust at the same time. If W2 comes from more extended or cooler dust and responds later, the track in the flux–flux plane opens into a loop of the same sense as hysteresis (bluer excess while rising). Section 11.6 shows that a 10-day difference is enough to imitate the effect, and the size ratios of Lyu et al. (2019) make such a difference likely.

The statistic is D = β_decline − β_rise, with one shared temperature at the largest excess and a separate exponent on each branch (PREDICTIONS_5b.md Section 7). Slow re-formation predicts D > 0. So does an inter-band lag. The two are not separable with this cadence.

## References

- Barvainis, R. 1987, ApJ, 320, 537 (cited by Netzer 2015 for Eqs. 6–7; not read directly)
- Bevington, P. R., & Robinson, D. K. 2003, Data Reduction and Error Analysis for the Physical Sciences, 3rd ed. (McGraw-Hill) (textbook; not re-read)
- Blandford, R. D., & McKee, C. F. 1982, ApJ, 255, 419 (cited by Cackett et al. 2021; not read directly)
- Cackett, E. M., Bentz, M. C., & Kara, E. 2021, iScience, 24, 102557, "Reverberation mapping of active galactic nuclei: from X-ray corona to dusty torus" (arXiv:2105.06926; review; Eq. 1 and Sec. 5 read)
- Fabian, A. C., Vasudevan, R. V., & Gandhi, P. 2008, MNRAS, 385, L43, "The effect of radiation pressure on dusty absorbing gas around AGN" (arXiv:0712.0277v1 read in full; page numbers are those of the arXiv version)
- Honeycutt, R. K. 1992, PASP, 104, 435, "CCD ensemble photometry on an inhomogeneous set of exposures" (not re-read)
- Ishibashi, W., Fabian, A. C., Ricci, C., & Celotti, A. 2018, MNRAS, 479, 3335, "Revisiting the 'forbidden' region: AGN radiative feedback with radiation trapping" (arXiv:1806.08387; abstract only; journal reference not checked)
- Jarrett, T. H., et al. 2011, ApJ, 735, 112 (WISE zero points and isophotal wavelengths; values checked in the All-Sky Explanatory Supplement IV.4.h)
- Kelly, B. C., Bechtold, J., & Siemiginowska, A. 2009, ApJ, 698, 895 (not re-read)
- Kendall, M. G. 1975, Rank Correlation Methods, 4th ed. (Griffin) (not re-read)
- Kokubo, M., & Minezaki, T. 2020, MNRAS, 491, 4615 (arXiv:1904.08946; abstract only)
- Koshida, S., et al. 2009, ApJ, 700, L109 (arXiv:0907.0573; abstract only)
- Lyu, J., Rieke, G. H., & Smith, P. S. 2019, ApJ, 886, 33, "Mid-IR Variability and Dust Reverberation Mapping of Low-z Quasars. I." (arXiv:1909.11101; Eqs. 15, 16, B8, B9, abstract and Appendix A read)
- MacLeod, C. L., et al. 2010, ApJ, 721, 1014 (not re-read)
- Mann, H. B. 1945, Econometrica, 13, 245 (not re-read)
- Meisner, A. M., Caselden, D., Schlafly, E. F., & Kiwy, F. 2023, AJ, 165, 36, "unTimely: a Full-sky, Time-domain unWISE Catalog" (arXiv:2209.14327; column definitions read)
- Netzer, H. 2015, ARA&A, 53, 365, "Revisiting the Unified Model of Active Galactic Nuclei" (arXiv:1505.00811; review; Eq. 1 and Sec. 1.2.4 read)
- Peterson, B. M. 1993, PASP, 105, 247, "Reverberation mapping of active galactic nuclei" (review; not read directly)
- Rybicki, G. B., & Lightman, A. P. 1979, Radiative Processes in Astrophysics (Wiley) (textbook; not re-checked)
- Sen, P. K. 1968, J. Am. Stat. Assoc., 63, 1379 (not re-read)
- Ricci, C., et al. 2017, Nature, 549, 488, "The close environments of accreting massive black holes are shaped by radiative feedback" (arXiv:1709.09651v1; main text, Methods and Extended Data captions read; page numbers are those of the arXiv version)
- Son, S., Kim, M., Ho, L. C., Kim, D., & Kim, T. 2022, ApJ, 937, 3, "A Mid-infrared Flare in the Seyfert Galaxy NGC 3786: A Changing-look Event Triggered by an Obscured Tidal Disruption Event?" (arXiv:2208.08062; read in full)
- Son, S., Kim, M., Ho, L. C., Kim, D., & Kim, T. 2026, "The Changing-look Phenomenon Accompanied by an Accretion Mode Transition in NGC 3786" (arXiv:2604.17173, accepted by ApJ according to the arXiv listing; Secs. 2 and 4 read)
- Theil, H. 1950, Proc. K. Ned. Akad. Wet., 53, 386 (not re-read)
- van Velzen, S., Pasham, D. R., Komossa, S., Yan, L., & Kara, E. A. 2021, Space Sci. Rev., 217, 63, "Reverberation in tidal disruption events: dust echoes, coronal emission lines, multi-wavelength cross-correlations, and QPOs" (arXiv:2107.12268; review; Sec. 2.1 read)
- Waxman, E., & Draine, B. T. 2000, ApJ, 537, 796, "Dust sublimation by GRBs and its implications" (arXiv:astro-ph/9909020; Eqs. 10–12 read)
- Wright, E. L., et al. 2010, AJ, 140, 1868 (WISE mission paper; Table 1 read from arXiv:1008.0031)
