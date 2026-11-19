# PIPELINE: 정제 규칙과 그 근거

파이프라인이 내리는 모든 선택과 그 근거 문서를 적는다. 수치 상수는 `src/irecho/config.py`에 있고, 이 문서의 절 번호가 코드 주석에 인용되어 있다. 2단계 known-answer test의 판정 기준(5절)과 대상 선정 기준(6절)은 해당 광도곡선을 보기 전에 커밋했다(commit `4f397cb`, `65da8eb`).

인용한 문서의 약칭:

- NEOWISE ES: *NEOWISE Data Release Explanatory Supplement* (IRSA, https://irsa.ipac.caltech.edu/data/WISE/docs/release/NEOWISE/expsup/). 2026-11-12에 II.1.a, II.2, II.3절을 읽었다. 해당 페이지의 개정 기록에는 2024-11-14 "NEOWISE Final Data Release" 갱신이 적혀 있다.
- AllWISE ES: *AllWISE Explanatory Supplement* (https://irsa.ipac.caltech.edu/data/WISE/docs/release/AllWISE/expsup/). III.2절.
- All-Sky ES: *WISE All-Sky Release Explanatory Supplement* (https://irsa.ipac.caltech.edu/data/WISE/docs/release/All-Sky/expsup/). IV.4.h절.

## 1. 조회

- 서비스: IRSA TAP sync (`https://irsa.ipac.caltech.edu/TAP/sync`), `FORMAT=csv`.
- 표: `neowiser_p1bs_psd` (NEOWISE-R Single Exposure Source Table), `allwise_p3as_mep` (AllWISE Multiepoch Photometry Table).
- 조건: 표적 좌표 중심 반지름 3″ cone, `ORDER BY mjd`. NEOWISE ES II.3.b는 정지한 천체의 모든 검출을 찾는 방법으로 "a simple positional cone search with a small radius (e.g. 3 arcsec)"를 든다.
- 표적 좌표는 NED 이름 조회로 얻은 위치에서 가장 가까운 AllWISE Source Catalog 항목의 좌표다(`targets/` 파일에 고정).
- ADQL 전문, 받은 시각, 행 수, SHA-256은 `data/manifest.json`에 기록한다. 캐시가 있으면 다시 받지 않고, 캐시에 기록된 ADQL이 요청과 다르면 오류를 낸다(`irsa.fetch_csv`).
- 조회는 한 번에 하나씩, 최소 1초 간격으로 보낸다.
- 한 frame 안에 cone에 드는 검출이 둘 이상이면 표적 좌표에 가장 가까운 것 하나만 쓴다(`quality.nearest_per_frame`).

## 2. 품질 플래그

밴드별로 적용한다. cc_flags와 moon_masked는 밴드마다 한 글자씩이므로, W1에는 첫째 글자, W2에는 둘째 글자만 본다. 그래서 한 frame에서 W1은 버리고 W2는 쓰는 경우가 있다.

| 기준 | NEOWISE-R | AllWISE MEP | 근거 |
|---|---|---|---|
| 측정값 존재 | `w?mpro`, `w?sigmpro`가 null이 아님 | `w?mpro_ep`, `w?sigmpro_ep`가 null이 아님 | AllWISE ES III.2.b: SNR<2이면 등급은 95% 상한이고 "the magnitude uncertainty, w?sigmpro_ep, is null". 상한은 측정값이 아니므로 뺀다. |
| frame 전체 품질 | `qual_frame > 0` | 컬럼 없음 | NEOWISE ES II.3.a: "Select database entries that have qual_frame>0 to obtain detections from framesets with good quality." AllWISE는 품질이 나쁜 frameset을 처리 단계에서 이미 제외했다(III.2.a). |
| image quality | `qi_fact > 0` | `qi_fact = 1` | NEOWISE ES II.3.a: "Select database entries with qi_fact>0 to exclude extractions from framesets that may have smeared images". AllWISE ES III.2.b: "limit entries to those with qi_fact=1 to avoid cases were low level smearing may lead a slight phtometric bias." 각 문서의 권고를 그대로 따랐다. |
| South Atlantic Anomaly | `saa_sep > 5` | `saa_sep > 5` | NEOWISE ES II.3.a와 AllWISE ES III.2.b는 `saa_sep>0`을 권한다. NEOWISE ES II.1.a의 컬럼 설명은 여기에 더해 "negative or small positive (< 5.0) values of saa_sep should be used with caution"이라고 적는다. 더 엄격한 5°를 택했다. |
| 달 산란광 | `moon_masked`의 해당 밴드 글자 = '0' | 같음 | NEOWISE ES II.3.a: "Select database entries where moon_masked='00'". AllWISE ES III.2.b: "Select database entries where moon_masked=0". |
| artifact 오염 | `cc_flags`의 해당 밴드 글자 = '0' | 같음 (단, 표적 단위 플래그) | NEOWISE ES II.3.b: "have a confusion and contamination flag value of "0" in each band. (cc_flags like '00%')". AllWISE MEP의 cc_flags는 frame별이 아니라 AllWISE 카탈로그 항목의 값이다. |
| 포화 | `w?sat = 0` | 컬럼 없음 | NEOWISE ES II.2: 포화된 천체의 profile-fit 밝기는 체계적으로 과대평가된다. 4절 참고. |
| 태양계 천체 | `sso_flg = 0` | 컬럼 없음 | NEOWISE ES II.2: sso_flg=1은 알려진 태양계 천체의 예측 위치와 겹친다는 뜻이다. 소행성 빛이 섞였을 수 있으므로 뺀다. |

일부러 쓰지 않은 기준:

- `w?rchi2` 상한: NEOWISE ES II.3.b는 w?rchi2가 작은 검출이 더 믿을 만하다고 하면서도 은하처럼 PSF보다 넓은 천체는 실제 검출인데도 값이 크다고 적는다. 표적이 은하핵이므로 자르지 않았다. 대신 표적별 중앙값을 결과에 적는다(NGC 3786은 W1 43.8, W2 15.6).
- `w?frtr` (transient 플래그): NEOWISE ES II.3.b의 권고 조건이 `w1mpro>11`, `w2mpro>10`인 어두운 천체용이다. 표적들은 그보다 밝고, cosmic ray는 방문 안의 outlier 제거(3절)로 다룬다.
- `allwise_cntr`/`tmass_key` 연관 요구: 표적 좌표 자체가 AllWISE 항목이므로 생략했다.

각 표적에서 기준별로 몇 행이 걸렸는지는 `LightCurve.cut_report`에 남고, NGC 3786의 값은 `results/ngc3786_summary.json`의 `cut_report`에 있다.

## 3. 방문(epoch) 묶기와 오차

- 시간순으로 정렬해 이웃한 측정의 간격이 30일을 넘으면 새 방문으로 나눈다. 이 표적들의 실제 방문 간격은 약 6개월이고 한 방문은 며칠 안에 끝난다. 황도 극 근처처럼 관측이 이어지는 곳에는 이 규칙을 쓸 수 없다.
- 통계는 flux 공간에서 낸다. 방문의 대푯값은 flux의 중앙값이다.
- outlier 제거는 한 번만 한다. 중앙값에서 5 × max(1.4826·MAD, 카탈로그 오차의 중앙값)보다 멀리 떨어진 측정을 버린다.
- 품질 기준을 통과한 측정이 한 밴드에서 5개 미만이면 그 방문의 그 밴드는 비운다(NaN).
- 방문 시각은 사용한 측정들의 평균 MJD다.

오차 모형은 두 가지를 모두 계산해 둔다.

```
err_formal  = sqrt(π/2) · sqrt(mean(σ_i²)) / sqrt(n)      카탈로그 오차가 예측하는 중앙값의 표준오차
err_scatter = sqrt(π/2) · s / sqrt(n)                     방문 안에서 실제로 관측된 표준편차 s를 쓴 표준오차
err         = max(err_formal, err_scatter)                 채택값
```

sqrt(π/2) ≈ 1.2533은 정규분포에서 중앙값의 표준오차가 평균의 표준오차보다 큰 비율(점근값)이다. n ≈ 5–15에서는 실제 비율이 조금 작아 약간 보수적이다(단위 테스트 `test_error_model_calibration_on_gaussian_noise`).

이 오차 모형이 실제 흩어짐을 설명하는지는 2단계 (a)에서 기준 별로 시험한다. 결과와 한계는 README에 적는다.

Son et al. (2022)은 방문별 오차를 "the root mean square sum of the measurement errors and the standard deviation of the measured brightnesses in the multiple observations in each visit"로 정의한다. sqrt(n)으로 나누는지 등 세부는 논문에서 확인하지 못했고, 이 파이프라인의 정의와 같지 않다.

## 4. 등급 → flux, 색, 포화, 혼입

### 변환

```
F_ν [Jy] = F_ν0 · 10^(−m / 2.5),     F_ν0(W1) = 309.540 Jy,  F_ν0(W2) = 171.787 Jy
σ_F = F · (ln 10 / 2.5) · σ_m
W1 − W2 = m_W1 − m_W2,   σ = sqrt(σ_W1² + σ_W2²)
```

- 출처: All-Sky ES IV.4.h Eq. 1과 Table 1(Jarrett et al. 2011). 표의 값을 2026-11-12에 직접 확인했다.
- F_ν0은 F_ν = 일정인 스펙트럼 기준이다. All-Sky ES IV.4.h는 스펙트럼이 가파르면 color correction이 필요하다고 적는다. 이 파이프라인은 color correction을 적용하지 않는다. 그래서 Jy 값은 "F_ν = 일정 가정의 flux density"이고, 붉은 스펙트럼에서는 몇 퍼센트 수준의 계통 차이가 있을 수 있다(정확한 크기는 이 프로젝트에서 계산하지 않았다).
- 등급은 Vega 기준 profile-fit 등급(`w?mpro`)이다. NEOWISE ES II.2는 aperture 측광보다 profile-fit 측광을 쓰라고 권한다("Defer to the profile-fit photometric measurements ... whenever possible").

### 포화 한계

- NEOWISE ES II.2: "Profile-fit measurements of sources brighter than the W1 and W2 saturation limits (W1~8.0 and W2~7.5 mags) are systematically overestimated".
- AllWISE ES III.2.b: post-cryo 측정은 "W1<8 and W2<7 mag"부터 cryo 측정보다 체계적으로 밝아지며, 이 차이가 "can give a false indication of flux variability".
- 그래서 대조군은 AllWISE 등급 W1 > 8.5, W2 > 8.0인 천체만 뽑았고(`SATURATION_WARN_MAG`), `w?sat = 0`인 측정만 쓴다.
- NGC 3786의 가장 밝은 방문(2024-05-10)은 W1 9.35, W2 8.73이다. 한계까지 W1 1.3 mag, W2 1.2 mag 여유가 있다.
- 대조군 중 NGC 1386은 AllWISE W2 = 8.09로 기준(8.0)을 0.09 mag 차이로 통과했다. 포화 한계에 가장 가까운 대상이다.

### 혼입 (blending)과 host 은하

- WISE의 angular resolution(FWHM)은 W1 6.1″, W2 6.4″다(Wright et al. 2010). 이보다 가까운 천체는 분리되지 않는다.
- NEOWISE ES II.2: profile-fit 측광은 점광원에 맞춰져 있어 "systematically underestimate the brightness of sources that are extended"이고, 그런 천체는 `w?rchi2>2`일 수 있다.
- NGC 3786의 single-exposure `w1rchi2` 중앙값은 43.8, `w2rchi2`는 15.6이다. 즉 측정값은 "핵 + PSF 안에 든 host 은하 빛"이고 점광원 모형과 잘 맞지 않는다.
- 결과 해석에 미치는 영향:
  1. 등급으로 잰 섬광 증가폭은 host 빛에 희석된 값이다. 핵만의 증가폭은 더 크다. host에 무관한 양은 flux 초과분(mJy)이다.
  2. 방문마다 scan 방향이 달라 넓은 천체의 profile-fit 값이 흔들릴 수 있다. 이 효과는 점광원인 기준 별로는 검증되지 않는다(README 한계 절).
- 주변 천체: AllWISE Source Catalog에서 NGC 3786 중심 36″ 안의 다른 항목은 모두 W1 ≥ 15.8로, 중심 항목(W1 10.28)보다 5.5 mag 이상 어둡다(2026-11-12에 반지름 36″ cone으로 직접 조회해 확인. 이 항목들은 은하 원반의 일부일 수 있다). 짝은하 NGC 3788은 NED 좌표 기준 1.38′ 떨어져 있어 PSF 밖이다.

## 5. 밝아짐(섬광 후보) 판정식

방문 단위 광도곡선(NEOWISE-R 방문만, 두 밴드 모두 측정된 방문만)에 적용한다. AllWISE 방문은 측정 방식이 달라(카탈로그 위치에서의 forced 측정) 판정에 쓰지 않고 그림에만 넣는다.

방문 k에 대해, 그보다 앞선 방문 중 이미 판정된 방문을 뺀 것을 baseline으로 한다. baseline이 4개 미만이면 판정하지 않는다.

```
m_base = baseline 등급의 중앙값
s_base = max(1.4826·MAD(baseline 등급), baseline 방문 오차의 중앙값)
Δm_k   = m_base − m_k                    (양수 = 밝아짐)
S_k    = Δm_k / sqrt(s_base² + err_k²)
판정:   W1과 W2 모두에서  S_k ≥ 5  이고  Δm_k ≥ 0.2 mag
```

- 시작 시기(onset)는 처음 판정된 방문이다. 6개월 간격 자료이므로 실제 시작은 "직전 방문과 그 방문 사이"라고만 말할 수 있다.
- 두 밴드 동시 조건은 한 밴드에만 나타나는 artifact를 막기 위한 것이다.
- 0.2 mag 하한은 통계적으로는 유의하지만 작은 변화(보정 오차, 보통의 AGN 변광)를 섬광으로 부르지 않기 위한 것이다.
- s_base에 baseline의 실제 흩어짐이 들어가므로, 3절의 오차가 과소평가되어 있어도 판정 통계량이 그만큼 부풀지는 않는다.
- 한계: baseline이 "그 이전 전체"이므로 몇 년에 걸친 느린 밝아짐도 판정될 수 있다. 이 판정식은 "과거 대비 유의한 밝아짐"을 찾는 것이지 섬광의 물리적 원인을 가리지 않는다. 0.2 mag보다 작은 섬광은 설계상 찾지 못한다.
- 실행 뒤에 덧붙임(2026-11-12, 기준값은 바꾸지 않았다): 실제로는 반대 방향의 한계도 드러났다. 판정되지 않은 방문은 다음 방문의 baseline에 들어가므로, 한 방문에 5σ를 넘지 않는 느린 밝아짐은 baseline의 흩어짐 s_base를 키워 끝까지 판정되지 않는다. 대조군 NGC 2273이 그 예다(README 3.3절). 즉 이 판정식은 6개월 안에 일어나는 급한 밝아짐에만 민감하다.

기준값(5σ, 0.2 mag, 4개)은 NGC 3786을 포함한 어떤 표적의 결과도 보기 전에 정했다.

## 6. 대상 선정 기준

`scripts/select_targets.py`에 상수로 들어 있다. 결과는 `targets/`에 고정했다.

### 기준 별 (R1–R13)

NGC 3786과 같은 방문에 관측되도록 같은 시야에서 고른다. 이렇게 하면 방문별 보정 오차를 별들의 공통 잔차로 따로 볼 수 있다.

| | 기준 | 값 | 이유 |
|---|---|---|---|
| R1 | 위치 | NGC 3786에서 1.5° 이내, 3′ 이내는 제외 | 같은 방문 공유. 은하 빛 회피 |
| R2 | 밝기 | 9.5 ≤ W1 < 12.5 | 포화와 멀고 single-exposure SNR이 수십 이상. NGC 3786(W1 ≈ 10.3)을 포함하는 범위 |
| R3 | 색 | −0.15 < W1−W2 < 0.15 | 별의 photosphere 색. AGN과 먼지 많은 천체 제외 |
| R4 | cc_flags | W1, W2 글자 모두 '0' | artifact 없음 |
| R5 | ext_flg | 0 | 점광원 |
| R6 | nb, na | nb = 1, na = 0 | deblending 없음 |
| R7 | 포화 | w1sat = w2sat = 0 | |
| R8 | var_flg | W1, W2 자리 모두 0–5 | AllWISE(2010–2011) 기준 변광 가능성 낮음. 'n'(자료 부족)은 제외 |
| R9 | ph_qual | W1, W2 모두 'A' | |
| R10 | rchi2 | w1rchi2, w2rchi2 < 3 | PSF와 잘 맞음 |
| R11 | 2MASS | n_2mass = 1, J−Ks ≤ 1.0 | 근적외선에서도 별 색 |
| R12 | 고유운동 | AllWISE pm < 150 mas/yr | 14년 동안 3″ cone을 벗어나지 않음 |
| R13 | 고립 | 15″ 안에 W1 < W1_별 + 4인 이웃 없음 | 혼입 2.5% 미만 |

통과 수(2026-11-12 실행, `targets/reference_stars_funnel.json`): 시야 조회 29,228 → R1 28,880 → R2 804 → R3 753 → R4 721 → R5 692 → R6 679 → R7 648 → R8 647 → R9 647 → R10 647 → R11 646 → R12 616 → R13 589.

589개 중에서 W1 구간 [9.5, 10.5), [10.5, 11.5), [11.5, 12.5)마다 10개씩, 모두 30개를 `numpy.random.default_rng(20261112)`로 뽑았다(구간별 후보 79, 178, 332개). 뽑은 뒤에는 결과를 보고 별을 빼거나 바꾸지 않는다.

주의: R8은 WISE 자신의 2010–2011년 자료로 만든 플래그다. 시험에 쓰는 NEOWISE-R 자료(2014–2024)와는 시기가 겹치지 않지만 같은 장비다. Gaia 같은 독립 자료로 변광 여부를 확인하지는 않았다.

### 대조군 AGN (C1–C4)

| | 기준 |
|---|---|
| C1 | 가까운 은하의 optical Seyfert 2 핵. 광도곡선을 보기 전에 문헌에서 이름을 골랐다 |
| C2 | 선정 시점에 changing-look 사건이나 중적외선 섬광 보고를 알지 못함 |
| C3 | AllWISE W1 > 8.5, W2 > 8.0 (포화 회피) |
| C4 | \|b\| > 20° (혼잡하지 않은 시야) |

후보 8개(Mrk 573, NGC 1386, NGC 2273, NGC 3081, NGC 3393, NGC 4941, NGC 5347, NGC 5728)가 모두 C3, C4를 통과했다(`targets/control_agn.csv`).

C2의 한계를 분명히 적는다. 2026-11-12에 이 여덟 이름과 "changing-look", "mid-infrared flare", "mid-infrared variability"를 묶어 웹 검색을 두 번 했고, 해당 천체를 다룬 결과는 나오지 않았다. 이것은 문헌 전수 조사가 아니다. "변광이 보고되지 않았다"가 아니라 "간단한 검색에서 보고를 찾지 못했다"가 정확한 표현이다. Seyfert 2도 중적외선에서 어느 정도 변할 수 있으므로, 대조군 시험이 묻는 것은 "전혀 변하지 않는가"가 아니라 "5절의 판정식이 이들을 섬광으로 판정하는가"다.

## 7. 차등 측광(differential photometry)과 비교 별

근거와 식은 THEORY.md 8절에 있다. 여기에는 정한 값만 적는다.

- 한 시야의 비교 별마다 자기 중앙값을 뺀 잔차를 만들고, 방문마다 별들 잔차의 중앙값을 그 방문의 영점 오프셋으로 삼는다. 오프셋의 오차는 `sqrt(π/2) × 1.4826·MAD / sqrt(N)`이다. 별이 3개 미만인 방문은 오프셋을 내지 않는다(`differential.MIN_STARS`).
- 표적의 방문은 15일 안에서 가장 가까운 시야 방문의 오프셋을 뺀다(`differential.MATCH_DAYS`). 해당 방문이 없으면 그 표적 방문은 버린다.
- 보정 뒤 방문 오차는 `sqrt(err_scatter² + offset_err²)`이다. `err_scatter`는 방문 안 흩어짐에서 온 오차(3절)다. 2단계의 "둘 중 큰 쪽" 대신 이것을 쓰는 이유는 README 5.1절에 있다.
- 별의 검증은 leave-one-out으로 한다. 별 하나를 보정할 때 그 별을 뺀 나머지로 오프셋을 계산한다.
- AllWISE 방문에는 적용하지 않는다. 비교 별의 AllWISE 자료를 대조군 시야에서는 받지 않았고, AllWISE 방문은 판정에도 쓰지 않는다.

### 대조군 시야의 비교 별

`scripts/select_field_stars.py`. 대조군 은하 8개마다 그 은하에서 1.5° 이내(3′ 이내 제외)의 AllWISE 항목에 6절의 R1–R13을 그대로 적용하되, 두 가지만 다르게 정했다. 별의 광도곡선을 받기 전에 정해 커밋했다(`11b21e4`).

| | 값 | 이유 |
|---|---|---|
| W1 범위 | 9.5 ≤ W1 < 11.5 | 대조군 핵의 밝기(W1 9.0–10.3)에 가깝게. 9.5보다 밝은 쪽은 포화 여유 때문에 쓰지 않는다. 더 어두운 별은 오차가 커서 영점 추정에 보탬이 적다 |
| 개수 | [9.5, 10.5)에서 5개, [10.5, 11.5)에서 5개 | 시야당 10개. 10개로 충분한지는 NGC 3786 시야에서 따로 시험했다(THEORY.md 8.4절) |

추첨은 `numpy.random.default_rng([20261112, 시야 번호])`다. 시야별 통과 수는 `targets/field_stars_funnel.json`에 있다(R13까지 통과한 후보: Mrk 573 338, NGC 1386 329, NGC 2273 610, NGC 3081 392, NGC 3393 397, NGC 4941 295, NGC 5347 221, NGC 5728 593). 뽑은 80개는 `targets/field_stars.csv`에 있다.

## 8. 느린 상승 판정식 (b)

5절의 판정식을 (a)라고 부른다. (a)는 한 방문의 급한 밝아짐만 찾는다. (b)는 광도곡선 전체에 "평탄하다가 오르는" 모형을 맞춘다. 근거는 THEORY.md 9절, 코드는 `src/irecho/slowrise.py`, 기준값은 `config.SlowRiseCriteria`다.

입력은 차등 측광으로 보정한 NEOWISE 방문(두 밴드 모두 측정된 방문, 8개 이상)이다.

```
모형 (밴드별, t는 년):
  상수          m = a
  hinge         m = a + s·max(0, t − t_on)     t_on = 방문 k−1의 시각. k = 0(전체 직선) 또는 4 ≤ k ≤ n−3
  지수          f = a + c·exp((t − t_end)/τ)   상대 flux. τ = 0.5, 1, 2, 4년

검정:  F = (χ²_상수 − χ²_모형) / (χ²_모형 / (n − 2)),  p = P(F_{1, n−2} > F)
       밝아지는 방향이 아니면 p = 1

선택:  두 밴드 중 나쁜 쪽의 p가 가장 작은 모형 하나를 고른다(두 밴드에 같은 모형, 같은 t_on).

판정:  W1과 W2 모두에서
         p ≤ 0.001
         모형이 말하는 누적 밝아짐 ≥ 0.1 mag
         중앙값(t_on 이전 방문) − 중앙값(마지막 3개 방문) ≥ 0.1 mag
       (지수 모형이나 k = 0이면 "이전 방문"은 처음 3개 방문)
```

- hinge는 "언제부터 올랐는가"를 방문 단위로 준다. 6개월 간격이므로 그 이상의 정밀도는 없다.
- p ≤ 0.001은 보정된 확률이 아니다. t_on과 τ를 훑어 가장 좋은 것을 고르고, AGN의 변광은 red noise이기 때문이다(THEORY.md 9.3절). 거짓 경보율은 9절의 방법으로 직접 잰다.
- 0.1 mag 하한은 보통의 Seyfert 변광(대조군의 방문 간 RMS 0.02–0.1 mag)을 느린 상승으로 부르지 않기 위한 것이다. 별의 방문 간 흩어짐(W1 약 5 mmag, W2 약 8 mmag)의 10배 이상이다.
- 마지막 방문 하나만 밝은 경우는 (b)가 아니라 (a)의 몫이다. 마지막 3개 방문의 중앙값을 쓰므로 그중 2개 이상이 밝아야 한다.
- Mann–Kendall Z와 Theil–Sen 기울기를 함께 기록하지만 판정에는 쓰지 않는다.

**이 판정식이 눈먼 시험이 아닌 부분.** (b)를 설계할 때 NGC 2273의 2단계 광도곡선(2020년 이후 느린 상승)을 이미 보았다. 기준값(0.001, 0.1 mag)은 NGC 2273의 값을 계산해 보기 전에 정했고 그 뒤 바꾸지 않았지만, "느린 상승을 찾는 판정식이 필요하다"는 동기 자체가 그 천체에서 왔다. 따라서 NGC 2273에 대한 (b)의 결과는 독립적인 확인이 아니다. 기준 별, 나머지 대조군, 주입 시험에 대해서는 실행 전에 정한 것이다(`11b21e4`).

## 9. 주입 시험과 ML 비교

코드는 `src/irecho/injection.py`, `src/irecho/mlcompare.py`, 격자는 `config.InjectionGrid`다. 실제 자료에 돌리기 전에 커밋했다(`11b21e4`).

**주입 대상(host).** 차등 측광으로 보정한 광도곡선: 대조군 은하 8개, NGC 3786 시야의 기준 별 30개와 대조군 시야의 비교 별 80개(별은 leave-one-out 보정). NGC 3786 자체에는 주입하지 않는다.

**구간.** host 광도곡선에서 "baseline 방문 L개 + 그 뒤 방문 7개"를 잘라 쓴다. 시작 방문 k를 옮겨 가며 가능한 구간을 모두 쓴다(방문 21개인 host라면 L = 4에서 11개, L = 8에서 7개, L = 12에서 3개).

**격자.**

| 축 | 값 |
|---|---|
| 증가폭 (두 밴드에 같은 값) | 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0 mag |
| 상승 시간 | 즉시(계단), 1년, 3년 (마지막 baseline 방문에서 시작하는 직선 상승, 그 뒤 유지) |
| baseline 길이 L | 4, 8, 12 방문 (약 2, 4, 6년) |

- 회수(recovery): 주입한 구간이 판정되면 회수로 센다. (a)는 주입 시작 이후의 방문이 판정된 경우만 센다. (b)는 구간 전체에 대한 판정이므로 판정 여부만 본다(찾은 시작 시점이 맞는지는 요구하지 않는다).
- 거짓 경보(false alarm): 주입하지 않은 같은 구간이 판정되는 비율. 별에서는 순수한 거짓 경보다. 대조군 은하에서는 실제 변광일 수 있으므로 "판정 비율"로 읽어야 한다.
- 두 밴드에 같은 증가폭을 넣는 것은 단순화다. 실제 dust echo는 W2가 더 크게 변한다(NGC 3786: 0.56, 0.86 mag). 여기의 "증가폭"은 작은 쪽 밴드의 증가폭으로 읽으면 된다.

**ML 비교: isolation forest 하나.** (Liu, Ting & Zhou 2008. scikit-learn 1.9.1의 `IsolationForest`, 나무 200개.)

- 입력은 구간마다 계산한 특징 18개다(`injection.FEATURE_NAMES`): 밴드별로 처음 4개 방문 중앙값 대비 최대 밝아짐(mag, σ 단위), Theil–Sen 기울기 × 기간, Mann–Kendall Z, 누적 변화(처음 4개 대 마지막 3개), 이웃한 방문 사이 최대 밝아짐, RMS/오차, von Neumann 비. 그리고 W1–W2 잔차의 상관계수와 색 변화.
- host를 5개 묶음으로 나눈다. 한 host의 구간이 학습과 평가에 함께 들어가지 않는다.
- 학습에는 학습 묶음 host의 **주입하지 않은 구간만** 쓴다. 주입한 구간은 학습에 들어가지 않는다. 즉 통계 판정식과 같은 주입 시험으로만 평가한다.
- 문턱값은 학습 묶음의 점수에서 정한다. "matched": (a) 또는 (b)가 같은 학습 구간에서 낸 거짓 경보율과 같은 분위수. 그 비율이 0이면 학습 구간 가운데 가장 이상한 것보다 더 이상해야 판정한다. 참고로 1%, 5% 분위수도 함께 낸다.
- isolation forest는 "다른 것과 다르다"를 찾는 것이지 "밝아졌다"를 찾는 것이 아니다. 그래서 어두워짐(−0.3 mag 계단)을 넣었을 때의 판정 비율도 따로 잰다.

## 10. unTimely와의 비교

`src/irecho/untimely.py`, `stage34.untimely_comparison`. 자료와 받는 방법은 DATA.md에 있다.

- 등급: `mag = 22.5 − 2.5 log10(FLUX)` (Vega nanomaggy. Meisner et al. 2023 Table 1). Jy로 바꿀 때는 4절과 같은 F0를 쓴다.
- 표적에서 5″ 안의 `primary = 1` 검출 가운데 가장 가까운 것을 그 epoch의 값으로 한다. 5″ 안의 검출 수를 함께 기록한다(은하가 여러 개로 쪼개졌는지 보려고).
- NGC 3786과 대조군은 모두 `flags_info = 66`(bit 1: HyperLeda의 큰 은하 안, bit 6: deblending을 하지 않는 영역. Schlafly et al. 2019 Table의 정의)이고 5″ 안의 검출은 epoch마다 하나였다.
- unTimely의 `DFLUX`는 통계 오차만이다(같은 표). 비교에서 오차로 쓰지 않고 그림에만 표시한다.
- 파이프라인의 NEOWISE 방문과 unTimely epoch은 15일 안에서 짝짓는다. 비교하는 양: 등급 차이의 중앙값, 그 차이의 흩어짐, 변화량끼리의 기울기. NGC 3786은 같은 baseline 방문 13개에 대해 증가폭과 초과 flux를 두 측광에서 따로 계산한다.
- neo7 판은 2020년 12월에서 끝난다. 2021년 이후의 상승은 비교할 수 없다.

## 11. 물리량 (3단계)

식과 해석은 THEORY.md 6–7절, 코드는 `src/irecho/physics.py`, `stage34.physics_analysis`다.

- 입력: 차등 측광으로 보정한 NGC 3786의 방문 flux. baseline은 (a)가 처음 판정한 방문 이전의 판정되지 않은 NEOWISE 방문 전부(13개)의 flux 중앙값, 그 오차는 `sqrt(π/2) × 표준편차 / sqrt(13)`이다.
- 초과 flux ΔF = F(방문) − F(baseline)을 (a)가 판정한 방문마다 계산한다. 판정되지 않은 방문은 ΔF가 오차에 비해 작아 계산하지 않는다.
- T: W1, W2의 ΔF 비에서 식 (11). L_dust, R_bb: 식 (17), (18). 거리는 NED z = 0.009083과 Son et al. (2022)의 우주론(H0 = 67.4, Ω_m = 0.315)으로 40.68 Mpc. 특이 속도 보정은 하지 않았다.
- 오차: 방문 flux와 baseline flux를 정규분포로 4000번 뽑는 Monte Carlo(seed 20261112)의 16–84% 구간. baseline 추출은 모든 방문이 공유한다. 거리, blackbody 가정, baseline 선택은 이 구간에 들어 있지 않고 따로 계통 오차로 보고한다.
- R_sub: THEORY.md 식 (6), (7)에 문헌의 L_bol을 넣는다. 출처는 `stage34.LITERATURE_LBOL`에 적었다(Son et al. 2022 3.3절, Son et al. 2026 4절).

### 가설 판정과 광학 비교

- 가설 H1–H3와 대안, 그 판정 규칙은 `docs/PREDICTIONS.md`에 있다. 판정 통계량을 계산하기 전, ZTF 자료를 받기 전에 커밋했다(`0aefbe4`). 온도 값 자체는 그 전에 이미 계산돼 있었다(같은 문서 0절).
- 규칙의 구현은 `stage34.hypothesis_tests`, `stage34.optical_comparison`이다. 문턱값은 그 문서의 값 그대로다(σ_sys 70 K, 100 K, 1000–1200 K, 1500 K, 3배, p < 0.05).
- 광학 자료(`src/irecho/ztf.py`): ZTF의 각 점을 가장 가까운 WISE 방문에 배정하고(91.3일 이내), 방문마다 점이 5개 이상이면 등급의 중앙값을 낸다. PREDICTIONS.md는 "관측 시즌(약 6개월) 중앙값"이라고만 적었고, 묶음을 WISE 방문에 맞추는 것은 ZTF 값을 보기 전에 정한 구현 세부다. 상관 계산에는 묶음의 중앙 시각이 WISE 방문에서 60일 이내인 것만 쓴다.
- `stage34.exploratory_temperature_checks`(온도와 밝기의 순위 상관)는 PREDICTIONS.md에 없는 사후 분석이다. 결과에 그렇게 표시한다.
- NGC 2273에 같은 광학 비교를 적용한 것도 사후에 추가한 것이다.

## 12. 5a단계: 표본, 시야 단위 조회, blind 품질 확인

표본의 선택 기준은 `docs/SAMPLE.md`에 있다. 이 절은 자료를 받고 확인하는 방법을 적는다.

### 12.1 한 번에 여러 위치를 조회하기

1–4단계는 천체 하나에 cone 조회 하나를 보냈고 조회 하나에 1–2분이 걸렸다. IRSA TAP은 cone 조건을 OR로 묶은 조회를 받지 않는다(1절). 5a단계는 **위치 목록을 표로 올려(TAP upload) 한 번에 조회**한다(`irsa.fetch_upload_csv`).

```sql
SELECT p.pid, t.ra, t.dec, t.mjd, ... FROM neowiser_p1bs_psd AS t, TAP_UPLOAD.pos AS p
WHERE CONTAINS(POINT('ICRS', t.ra, t.dec), CIRCLE('ICRS', p.pra, p.pdec, 0.0008333)) = 1
```

- 올리는 표는 `pid`(위치 번호), `pra`, `pdec` 세 열의 VOTable이다. 결과의 첫 열 `pid`가 각 행이 어느 위치의 것인지 알려 준다.
- 시야 하나(표적 + 비교 별)를 조회 하나로 받는다. `pid` 0이 표적, 1부터가 비교 별이다(`sample_io.field_positions`). NEOWISE-R과 AllWISE MEP에 각각 한 번씩, 시야당 조회 두 번이다.
- 조건(반지름 3″, 열 목록, 품질 컷)은 1–4단계와 같다. 같은 위치를 두 방식으로 받은 결과가 같은지는 `tests/test_regression_stage5.py`가 NGC 3786으로 확인한다.
- 서버가 행 순서를 보장하지 않으므로 받은 CSV를 `pid`, `mjd`, 식별자 순으로 정렬해 저장한다. 해시가 실행마다 달라지지 않게 하려는 것이다.
- 캐시: `data/cache/s5/neowise/<이름>.csv`, `data/cache/s5/allwise_mep/<이름>.csv`. 옆의 `.meta.json`과 `data/manifest.json`에 ADQL, 올린 위치 목록, 받은 시각, 행 수, SHA-256을 적는다. 받은 시야는 다시 받지 않으므로 중단한 뒤 이어 받을 수 있다.
- AllWISE Source Catalog 조회(표적 후보 주변 5′, 비교 별 시야 1.5°)도 같은 방식으로 40개, 5개씩 묶어 보낸다.

### 12.2 unTimely

`scripts/fetch_untimely.py --sample`. 표적이 속한 HEALPix partition 하나를 열어, 표적과 **같은 partition에 있는 비교 별**의 행을 함께 꺼낸다. 다른 partition에 있는 별은 받지 않는다(파일을 더 열지 않으려고). 별의 unTimely 측광은 오차 바닥의 기준으로 쓴다(12.4절).

### 12.3 blind 규칙과 품질 확인

5a단계에서는 표적 은하의 변광을 보지 않는다. 코드가 지키는 규칙은 다음과 같다(`src/irecho/qa.py`).

| | 표적 은하 | 비교 별 |
|---|---|---|
| 계산하고 저장하는 것 | 노출 수, 품질 컷별 탈락 수, 방문 수, 관측 기간과 간격, 포화 화소가 있는 노출의 비율, 포화 여유보다 밝은 방문의 **수**, 방문 오차의 크기, 두 측광 차이의 흩어짐 | 제한 없음: χ²_red, RMS, 방문별 오프셋, 주입 시험 |
| 계산하지 않는 것 | 판정식 (a), (b), 증가폭, 색 변화, 추세, ZTF 비교 | |
| 그리지 않는 것 | 광도곡선 | |

- 표적의 광도곡선은 방문 수를 세고 최종 파일을 만들기 위해 메모리에서 계산되지만, 화면과 `results/`에는 위 표의 값만 나간다. `tests/test_qa.py::test_target_counts_are_blind`는 1 mag 밝아진 합성 표적과 일정한 합성 표적이 같은 QA 값을 내는지 확인한다.
- 두 측광 차이의 흩어짐은 같은 방문을 두 방법으로 잰 값의 차이에서 나온다. 표적의 변광은 두 값에 똑같이 들어 있어 차이에서 빠진다(`tests/test_errorfloor.py::test_pair_scatter_ignores_source_variability`).
- 최종 광도곡선(차등 보정 + 오차 바닥)은 `data/frozen/sample_v1/<이름>.csv`에 쓰고 git에는 넣지 않는다. 파일별 SHA-256만 `data/freeze_v1.json`에 커밋한다. 다음 단계에서 판정식을 돌릴 때 이 해시와 같은 파일인지 확인할 수 있다.
- QA 플래그(`run_qa.py`의 상수)는 주의할 천체를 표시할 뿐이고 표본에서 빼지 않는다: 보정된 방문 15개 미만, 포화 화소 노출 10% 초과 또는 포화 여유보다 밝은 방문 있음, 시야의 별 χ²_red 중앙값 2 초과, 10일 넘는 방문, unTimely 없음·짝 부족·여러 개로 쪼개짐.

### 12.4 은하의 오차 바닥

3–4단계 보고서의 남은 문제 2번(은하의 방문 오차가 과소평가됨)에 대한 설계다. 코드는 `src/irecho/errorfloor.py`, 근거는 THEORY.md 10절.

1. 방문 k를 두 측광으로 잰 차이 d_k = m_pipe,k − m_unTimely,k를 만든다. 파이프라인 쪽은 차등 보정한 값이다.
2. 차이의 흩어짐 s = 1.4826 × MAD(d_k). 짝이 8개 이상일 때만 쓴다.
3. 같은 값을 비교 별에 대해 구해, W1 9.5–10.5와 10.5–11.5 두 구간의 중앙값을 기준 s_ref로 둔다. 별보다 밝거나 어두운 은하는 가까운 구간의 값을 쓴다.
4. 바닥 = sqrt(max(0, s² − s_ref²)). 밴드마다 천체마다 하나다.
5. 짝이 8개 미만이거나 unTimely 검출이 없는 천체는 표본의 중앙값을 쓴다.
6. 방문 오차 = sqrt(방문 안 흩어짐에서 구한 오차² + 오프셋 오차² + 바닥²) (`errorfloor.apply_floor`). 등급과 flux는 바꾸지 않는다.

- 바닥을 과대평가하는 쪽으로 치우친 설계다. 은하에서 생기는 초과 흩어짐의 일부는 unTimely의 몫일 텐데 전부 파이프라인에 돌린다.
- unTimely neo7이 2020년에서 끝나므로 바닥은 2014–2020년 방문에서 구해 모든 방문에 적용한다.
- 5a단계에서는 바닥 값을 구해 최종 파일에 넣는 것까지만 한다. 바닥을 넣은 뒤 표적의 χ²나 판정 결과가 어떻게 달라지는지는 보지 않았다.

### 12.5 별에 대한 주입 시험

`scripts/run_injection_stars.py`. 9절의 격자(`config.InjectionGrid`)를 표본 시야의 비교 별(leave-one-out 보정)에만 적용해 판정식 (a), (b)의 회수율과 거짓 경보율을 미리 표로 만든다. 표적의 광도곡선은 읽지 않는다. ML 모형은 돌리지 않는다.

- 결과: 별 2,310개, 시도 1,085,128번. 주입 없는 구간 49,324개에서 (a) 0개, (b) 8개(모두 별 하나). 표는 SAMPLE.md 9절과 `results/qa_injection_completeness.csv`.
- 변형 `--floor-mmag X`: 별을 오차 바닥이 X mmag인 은하처럼 만든다. 방문마다 두 밴드에 표준편차 X의 독립 Gaussian 잡음을 더하고, 같은 X를 방문 오차에 제곱합으로 더한다(`errorfloor.apply_floor`가 표적에 하는 것과 같다). 시야마다 별 2개(밝기 구간마다 첫 번째)만 쓴다. X = 10, 20, 50 mmag로 돌렸다. 이 값은 표적의 바닥 분포(12.4절)를 본 뒤에 골랐다. 실제 초과 흩어짐이 Gaussian이고 독립이라는 보장은 없다.
- 12.4절의 바닥은 0.01 mmag로 반올림한 값을 채택한다. `results/qa_error_floor.json`에 적힌 수가 동결 파일에 들어간 수와 정확히 같아서, 캐시와 이 파일만으로 동결 파일을 다시 만들 수 있다(`tests/test_regression_stage5.py::test_field_qa_and_frozen_file_rebuild_from_cache`).

## 13. 5b단계: 사전 등록 분석

규칙의 전문은 [PREDICTIONS_5b.md](PREDICTIONS_5b.md)에 있고 이 절은 절차만 적는다. 순서가 중요하다: 예측과 판정 규칙 커밋 → 동결 해시 확인 → 합성 자료로 추정량 시험과 기준 확정 커밋 → 그 뒤에 동결 파일을 처음 열었다.

| 단계 | 코드 | 산출물 |
|---|---|---|
| 동결 확인 | `scripts/run_stage5b.py`의 `check_freeze` (실행할 때마다) | `results/s5b_summary.json`의 `freeze_check` |
| 합성 자료 시험 | `src/irecho/breath.py`, `scripts/run_breath_validation.py` | `results/s5b_validation.json` |
| 물리량 표 | `scripts/build_props.py` (BASS DR2의 `logEdd`, BASS DR1의 `logNH`를 그대로 옮김) | `targets/sample_v1_props.csv` |
| Q1 완전도 | `scripts/run_injection_q1.py` (별에만 주입) | `results/s5b_q1_completeness.csv` |
| ZTF | `scripts/fetch_ztf.py --sample` | `data/cache/s5/ztf/` (git 밖), manifest |
| 분석 | `src/irecho/stage5b.py`, `scripts/run_stage5b.py` | `results/s5b_objects.csv`, `s5b_objects_nofloor.csv`, `s5b_summary.json` |
| 그림 | `src/irecho/plots5b.py`, `scripts/plot_stage5b.py` | `results/fig_s5b_*.png` |

- **동결 파일 읽기** (`stage5b.read_frozen`): 주 결과는 `w?_err`(바닥 포함). 참고 결과는 `w?_err_stat`과 `w?_offset_err`의 제곱합으로 오차를 다시 만들고 flux 오차도 그에 맞춘다. 등급과 flux는 같다.
- **탐지기**: `detect_flares`, `detect_slow_rise`를 기준값 그대로 부른다. "slow rise" = (b) 표시이고 `dm_fit`이 두 밴드 모두 0.2 mag 이상(`injection.Q1_MIN_MAG`).
- **기울기** (`breath.measure`): 기준 상태 3개 방문 → 초과 flux → flux 평면의 곡선 맞춤(매개변수 T_p, β와 기준 상태의 이동 두 개) → 제자리 가열의 최솟값 β_A,min → '측정 가능' 여부와 분류. 근거는 THEORY.md 11절.
- **hysteresis**: 같은 방문들을 밝아지는 구간과 어두워지는 구간으로 나눠 지수 두 개를 맞추고 D = β_D − β_R.
- **ZTF**: 서비스가 요청 하나에 1–4분 걸린다. 먼저 `ztf_objects_dr24` 표를 위치 목록으로 한 번 조회해(17초) ZTF 천체가 있는 표적만 광도곡선을 요청한다. `--order`로 밝아진 천체를 먼저 받고 `--workers`로 동시에 여러 개를 요청한다. 파일 단위로 이어 받는다.
- **표시와 오차의 네 가지 조합**: (주 표본 전체, 표시 없는 주 표본) × (바닥 포함, 바닥 제외)를 모두 계산한다. 판정은 바닥 포함의 두 표본으로 하고, 둘이 다르면 "갈림"이다.

## 14. 6단계: 섬광 후보 조사와 독립 표본 시험

### 14.1 A: 섬광 후보 (`docs/FLARES_PLAN.md`)

| 단계 | 코드 | 산출물 |
|---|---|---|
| 신호 확인, 시간 척도, 물리량, 분류 | `src/irecho/flares.py`, `scripts/run_flares.py` | `results/s6a_flares.json` |
| 다른 파장 (ZTF, ASAS-SN, Swift-BAT 157개월) | `src/irecho/otherbands.py`, `scripts/run_flares_other.py` | `results/s6a_other_bands.json` |
| 문헌과 목록 조회 기록 | (조회는 손으로 실행해 기록) | `results/s6a_literature.json` |
| 그림 | `plots5b.plot_flares`, `scripts/plot_stage6.py a` | `results/fig_s6a_flares.png` |

분류 기준(O1–O11)과 해석별 예측은 실제 조사 전에 커밋했다(`8793b03`).

### 14.2 B: hysteresis 통계량과 광도 (`docs/PREDICTIONS_6.md`)

순서가 중요하다: 예측·판정 규칙·표본 정의 커밋 → 표적과 비교 별 커밋 → 받기 → blind 품질 확인과 동결 → 실제 cadence와 오차로 문턱값 확정 커밋 → 그 뒤에 D를 처음 계산했다.

| 단계 | 코드 | 산출물 |
|---|---|---|
| 독립 표본 | `IRECHO_SAMPLE=v2 scripts/build_sample.py targets`, `stars` | `targets/sample_v2*.csv`, `sample_v2_funnel.json` |
| 받기 | `IRECHO_SAMPLE=v2 scripts/fetch_sample.py`, `fetch_untimely.py --sample` | `data/cache/s6b/` (git 밖), manifest |
| blind 품질 확인, 오차 바닥, 동결 | `IRECHO_SAMPLE=v2 scripts/run_qa.py` | `results/qa2_*`, `data/frozen/sample_v2/` (git 밖), `data/freeze_v2.json` |
| 두 모형의 예상 분포와 검정력 | `scripts/run_d_lum_power.py prelim` (받기 전), `final` (동결 뒤) | `results/s6b_power_prelim.json`, `s6b_power_final.json` |
| 밴드 간 지연 추정기의 검증 | `src/irecho/bandlag.py`, `scripts/run_bandlag_validation.py` | `results/s6b_bandlag_validation_prelim.json`, `_final.json` |
| 분석 | `scripts/run_stage6b.py` | `results/s6b_objects.csv`, `s6b_summary.json` |
| 그림 | `plots5b.plot_d_luminosity`, `scripts/plot_stage6.py b` | `results/fig_s6b_d_luminosity.png` |

- **표본 꼬리표** (`sample_io.SAMPLE_TAG`): 환경 변수 `IRECHO_SAMPLE`이 `v1`(기본)이면 5a단계의 경로(`targets/sample_v1*`, `data/cache/s5/`, `results/qa_*`, `data/frozen/sample_v1/`), `v2`이면 6B의 경로(`targets/sample_v2*`, `data/cache/s6b/`, `results/qa2_*`, `data/frozen/sample_v2/`)를 쓴다. 기준과 절차의 코드는 같다. v2에서는 z 규칙(표본 크기를 맞추는 z 상한 찾기)과 calibration 표시를 쓰지 않는다.
- **장난감 모형** (`run_d_lum_power.py`): 천체마다 실제 방문 시각과 오차, 적색편이, 광도를 쓴다. 먼지 반응의 지연과 뭉개짐은 (1+z) Δt_W1(L), 모형 1의 밴드 간 지연은 (1+z)[Δt_W2(L) − Δt_W1(L)](Lyu et al. 2019의 식 15–16), 모형 2는 밴드 간 지연 없이 재형성 시간 1, 2, 5년 중 하나다. 합성 표본마다 실제 자료와 같은 코드(`breath.measure`)로 D를 구해 ρ(D, x)의 분포를 만든다.
- **판정**: q2는 실제 D를 천체 사이에서 섞어 구한 ρ의 95번째 백분위, q1은 모형 1의 합성 표본에서 실제 '측정 가능' 수 n에 맞춘 ρ의 5번째 백분위다(n이 격자 사이면 선형 보간). 규칙은 PREDICTIONS_6.md 7절.
- **밴드 간 지연 추정기** (`bandlag.stacked_lag`): F2 = a + b[F1 − δ dF1/dt]를 천체마다 맞추고 표본 전체의 χ²(δ)를 더해 δ를 구한다. 합성 자료에서 통과 조건 셋을 모두 만족할 때만 실제 자료에 적용한다.

## 15. 7단계: CGCG 312-012의 X선 분석 (`docs/XRAY_PLAN.md`)

순서: 해석별 예측과 판정 기준 커밋(`e34ce17`) → 도구 설치 → 자료 받기 → 재처리 → 영역 → 스펙트럼 → 맞춤 → 판정. 문헌 조사는 계획을 커밋한 뒤에 했다.

| 단계 | 코드 | 산출물 |
|---|---|---|
| 판정 규칙과 수치 도우미 (HEASoft 없이 동작) | `src/irecho/xray.py` | 단위 테스트 `tests/test_xray.py` |
| 자료 받기 (HEASARC, Swift-BAT survey 페이지) | `scripts/fetch_xray.py` | `data/cache/s7/{nustar,swift,bat157}/` (git 밖), manifest |
| NuSTAR 보정 파일 일부만 받기 | `scripts/fetch_caldb_subset.py 2022-05-01` | `data/cache/s7/caldb/` (git 밖), `results/s7_caldb.json` |
| 환경 | `scripts/heasoft_env.sh` | (conda 환경 `.heasoft/`, git 밖) |
| 재처리 | `scripts/xray_process.sh xrt10`, `xrt22`, `nustar` | `data/cache/s7/work/` (git 밖) |
| 영역과 배경 광도곡선 확인 | `scripts/xray_regions.py` (`xray_process.sh regions`) | `results/s7_regions.json`, `results/fig_s7_xray_images.png` |
| 스펙트럼, 응답 행렬, 묶기, 우리 은하 N_H | `scripts/xray_process.sh products` | `data/cache/s7/work/*/prod*` (git 밖) |
| 맞춤 (PyXspec) | `scripts/xray_fit.py fits` | `results/s7_fits.json` |
| 광도 정의의 독립 확인, XRT 유효 면적, 그림용 배열 | `scripts/xray_fit.py lumcheck`, `arf`, `plotdata` | `results/s7_lum_check.json`, `s7_arf.json`, `s7_spectra_plot.json` |
| 합성 자료 검증, MCMC | `scripts/xray_fit.py sims 200`, `mcmc` | `results/s7_sim_validation.json`, `s7_mcmc.json` |
| 사후 확인 (계획에 없음) | `scripts/xray_fit.py posthoc` | `results/s7_posthoc.json` |
| 판정 표 (HEASoft 없이 동작) | `scripts/run_stage7.py` | `results/s7_summary.json` |
| 문헌·목록 조회 | `scripts/xray_literature.py` | `results/s7_literature.json` |
| 그림 | `src/irecho/plots7.py`, `scripts/plot_stage7.py` | `results/fig_s7_xray_spectra.png` |

- **도구**: HEASoft 6.37.1(NuSTARDAS 2.1.6, swxrtdas 3.7.2, XSPEC 13.0.1, PyXspec 2.2.0)을 conda로 프로젝트 폴더 안의 `.heasoft/`에 설치했다. 맞춤 스크립트는 그 환경의 Python에서 돌고, 판정 표·그림·테스트는 프로젝트의 `.venv`에서 돈다(커밋된 JSON만 읽는다).
- **경로 길이**: CALDB 조회(`quzcif`)와 `nuexpomap`은 파일 경로가 길면 실패한다("Error constructing path", "contains too many characters"). 그래서 (가) `CALDB` 환경 변수는 사용자 임시 폴더의 짧은 symlink(`$TMPDIR/irecho-s7-caldb` → `data/cache/s7/caldb`)를 가리키고, (나) `nuproducts`는 작업 폴더로 들어가 상대 경로로 부른다. 프로젝트 폴더 밖에 생기는 것은 이 symlink 하나와 XSPEC이 처음 실행될 때 만드는 설정 폴더(`~/.xspec`)다.
- **NuSTAR**: `nupipeline`을 기본 설정으로 한 번 돌린 뒤 배경 광도곡선(천체에서 120″ 밖, 200초 구간)을 봤다. 두 모듈 모두 290초짜리 짧은 구간 하나의 배경이 나머지의 1.5배였다. 계획대로 `saamode=optimized tentacle=yes`로 다시 처리했고(노출 23.05 → 22.64 ks, 22.82 → 22.44 ks) 이쪽을 썼다. 천체 영역은 영상 중심의 반지름 50″ 원, 배경은 같은 칩(DET0)에서 천체와 190″ 이상 떨어지고 시야 안에 완전히 들어가는 반지름 100″ 원 가운데 가장 먼 것(`xray_regions.nustar_background`)이다.
- **Swift XRT**: PC mode, grade 0–12, 천체 영역 반지름 20 화소, 배경 고리 60–110 화소. 고리 안에서 반지름 10 화소 칸을 옮겨 가며 Poisson 기대치를 10⁻⁵ 수준으로 넘는 곳을 찾았고 두 관측 모두 없었다. `xrtmkarf`에 exposure map을 넣었다.
- **묶기와 통계량**: XRT와 NuSTAR는 bin마다 1 count 이상으로 묶고(`ftgrouppha grouptype=min groupscale=1`) 배경 스펙트럼과 함께 `cstat`(W 통계량)으로, BAT의 8채널은 χ²로 맞춘다. 스펙트럼마다 data group을 따로 둔다.
- **모형**: `constant*TBabs*(zTBabs*clumin*zpowerlw + clumin*zpowerlw + pexrav + zgauss)`에서 M0는 첫 항만, M1은 둘째 항까지, M2는 전부다. `clumin`은 흡수를 걷어낸 정지계 2–10 keV 광도를 맞춤 매개변수(lg10Lum)로 만들어 준다. 그래서 광도의 구간이 N_H, Γ와 같은 방법(통계량 차이 2.706)으로 나온다. `pexrav`는 `rel_refl = −1`(반사만), cutoff 없음(`foldE`를 상한 10⁶ keV에 고정), cos i = 0.45, Γ는 주 성분과 같게 두고 norm을 자유로 둔다. 반사 세기 R은 norm의 비로 나중에 계산한다.
- **계산 격자**: 모형은 0.1–300 keV의 로그 격자 1500칸에서 계산해 응답 행렬의 격자로 보간한다(`AllModels.setEnergies`). 응답 행렬의 격자를 그대로 쓰면 `pexrav`가 든 맞춤 한 번이 2분 넘게 걸렸다. M0에서 두 방법의 차이는 통계량 0.08, N_H 0.1%다.
- **Δ2와 ΔN의 구간**: 두 매개변수의 차와 비이므로 `error`로 바로 구할 수 없다. 2010년의 값을 "2022년 값 − Δ"(또는 "× 10^−ΔN")로 묶어 Δ를 0.03 dex씩 옮기며 다시 맞추고, 통계량이 2.706 오르는 곳을 보간했다(`xray_fit.profile`).
- **PyXspec에서 조심할 점**: 둘째 이후 data group의 매개변수는 `Parameter.index`가 그 group 안의 번호를 준다. `error`나 chain의 열 이름에는 전체 번호(`startParIndex + index − 1`)를 써야 한다. 처음에 이를 놓쳐 2010년 매개변수의 구간이 2022년 것으로 계산됐고, 결과를 쓰기 전에 고쳐 다시 돌렸다.
- **판정**: `scripts/run_stage7.py`가 커밋된 맞춤 결과에 XRAY_PLAN.md 6절의 규칙(`xray.luminosity_label`, `absorption_label`, `overall_verdict`, `energy_budget`)을 적용한다.

## 16. 7b와 8단계

### 16.1 7b: X선 분석의 사후 보완 (`docs/XRAY_PLAN.md` 9절)

| 단계 | 코드 | 산출물 |
|---|---|---|
| 가벼운 확인 (BAT 채널 선택, XRT 보정 계수, GTI 겹침, chain 진단, hardness ratio) | `scripts/xray_7b.py light` | `results/s7b_light.json` |
| 맞춤 (`.xcm` 저장과 새 세션 재현, 3000칸 격자, XRT 상수의 두 해석, 기본 SAA 스펙트럼) | `scripts/xray_7b.py fits`, `scripts/xray_process.sh products_saa_none` | `results/s7b_fits.json`, `results/s7b_xcm/` |
| 표와 그림 | `scripts/run_stage7b.py` | `results/s7b_summary.json`, `results/fig_s7b_supplement.png` |

- BAT의 채널은 번호로 고른다(`xray.bat_channels_in_band`). XSPEC의 `ignore`에 실수를 주면 그 에너지를 포함한 채널까지 빠진다.
- XSPEC은 이름이 `none`으로 시작하는 응답 파일을 "응답 없음"으로 읽는다.

### 16.2 8단계: hysteresis 통계량의 검증 (`docs/VALIDATION_8.md`)

순서: 계획과 기준 커밋(`1e28fff`) → 오차 상관 추정 → null test → (직접 비교가 null에서 치우치지 않음을 확인한 뒤) 실제 자료의 직접 비교 → lag 모형 재시험.

| 단계 | 코드 | 산출물 |
|---|---|---|
| null 자료, 직접 비교의 추정량, driver 식 | `src/irecho/hyst8.py` | 단위 테스트 `tests/test_hyst8.py` |
| W1·W2 오차의 상관, unTimely 광도곡선 | `scripts/run_stage8_noise.py` (`IRECHO_SAMPLE=v2`로 한 번 더) | `results/s8_error_correlation_v1.json`, `_v2.json`, `data/cache/s8/untimely_lc_*.json` (git 밖) |
| B1 null test와 직접 비교의 null 시험 | `scripts/run_stage8.py null 2` | `results/s8_null.json` (회차별 값은 `data/cache/s8/null_*.json`) |
| B2 직접 비교, unTimely, scan 방향 | `scripts/run_stage8.py direct` | `results/s8_direct.json`, `results/s8_direct_objects.csv` |
| B3 lag 모형 재시험 | `scripts/run_stage8.py lag 2` | `results/s8_lag_model.json` |
| B4 선택 효과와 보고 보완 | `scripts/run_stage8.py select` | `results/s8_selection.json` |
| 그림 | `src/irecho/plots8.py`, `scripts/plot_stage8.py` | `results/fig_s8_null.png`, `fig_s8_direct.png`, `fig_s8_examples.png`, `s8_examples.json` |

- **오차 상관**: 같은 방문을 pipeline과 unTimely가 각각 잰 값의 차이를 W1과 W2에서 구하고, 방문들에 걸친 Pearson 상관을 천체마다 낸다(쌍 6개 이상). 두 측광의 오차 상관이 섞인 값이다.
- **null 자료**: 참값 W2 = 관측된 방문 flux. W1은 (a) 직선, (b) 넓이가 고정된 blackbody의 온도만 변하는 경우(`hyst8.truth_single_temperature`). 잡음은 동결 파일의 오차와 위의 상관으로 넣는다. 추정은 `breath.measure` 그대로다.
- **직접 비교** (`hyst8.direct_offset`): 두 구간의 W2 flux가 겹치는 구간 안의 방문만으로 F1 = c + m·F2 + δ·[밝아지는 구간]을 두 변수의 오차를 넣어 맞춘다(`hyst8.eiv_fit`). 기준 방문, 초과 flux, 공유 peak, blackbody를 쓰지 않는다. δ > 0가 D > 0에 대응한다.
- **scan 방향**: 방문 시각을 반년 단위로 반올림한 홀짝. unTimely의 `FORWARD` 표시와 맞는지 확인했다.
- 병렬 실행은 `multiprocessing`의 fork 방식을 쓴다(준비한 광도곡선을 worker가 물려받는다).

## 17. 9단계: 색 비대칭이 inter-band lag만으로 설명되는가 (`docs/PREDICTIONS_9.md`)

순서: 모형·기준·해석 규칙 커밋(`e099db9`, 실제 광도곡선에 맞추기 전) → B1 합성 실험과 그 결과 커밋 → 실제 자료 맞춤 → 맞춘 지연을 넣은 null 표본 → 판정.

| 단계 | 코드 | 산출물 |
|---|---|---|
| 보간 가중치, 일반화 최소제곱, 모형 사다리, 예측 δ, 합성 광도곡선 | `src/irecho/transfer9.py` | 단위 테스트 `tests/test_transfer9.py` |
| B1 식별 가능성 (진실 모형 11가지 × driver 2가지 × 표본 2개 × 50회) | `scripts/run_stage9.py ident 2` | `results/s9_identifiability.json` (회차별 값은 `data/cache/s9/ident_*.json`, git 밖) |
| B2 실제 자료 맞춤 | `scripts/run_stage9.py fit` | `results/s9_objects.csv`, `results/s9_transfer.json` |
| B2 null 표본 (맞춘 지연을 진실로) | `scripts/run_stage9.py null 2` | `data/cache/s9/fitted_*.json` (git 밖) |
| B2–B4 판정, 지연–광도 | `scripts/run_stage9.py verdict` (`scripts/stage9_verdict.py`) | `results/s9_verdict.json` |
| B5 진폭을 맞춘 driver의 mock | `scripts/run_stage9.py mock5 2` | `results/s9_mock5.json` |
| 그림 | `src/irecho/plots9.py`, `scripts/plot_stage9.py` | `results/fig_s9_identifiability.png`, `fig_s9_examples.png`, `fig_s9_offsets.png`, `s9_examples.json` |

- **모형**: F_W2(t_i) = a + b·(K ∗ F_W1)(t_i) [+ η·R_i]. K는 지연 τ를 중심으로 한 폭 w의 top-hat. 방문 사이의 W1은 이웃 방문을 잇는 직선이고, top-hat 구간의 평균은 그 직선의 적분으로 정확히 계산한다(`transfer9.kernel_weights`). 그래서 예측은 W1 방문 flux의 가중합이다.
- **외삽 금지**: kernel은 W2 방문의 과거 350일, 미래 150일까지만 닿는다. 그 구간이 W1 방문으로 덮이고 250일보다 긴 간격이 없는 W2 방문만 쓴다(`transfer9.usable_visits`). 보간 함수는 범위 밖의 시각을 받으면 예외를 낸다.
- **오차 전파**: 잔차의 공분산은 diag(σ2²) + b²·c·diag(σ1²)·cᵀ − b·(X + Xᵀ)다(`transfer9.residual_covariance`. 단위 테스트에서 Monte Carlo와 대조). 통계량은 S = χ² + ln det C. 선형 보간의 오차는 넣지 않았고 합성 자료로 쟀다.
- **모형 사다리**: T0(지연 없음), T1(τ), T2(τ, w), T3(η를 더함). 모든 모형을 T0와 직접 비교한다. 격자는 τ 5일 간격, w 0·100·200·300일의 284점.
- **예측 δ**: 맞춤에 쓴 방문에서 8단계의 `hyst8.direct_offset`을 관측된 W2와 모형의 W2로 각각 계산한다(`transfer9.offsets`). 가지 표시와 공통 범위 안의 방문은 관측 쪽 것을 쓴다(`direct_offset`에 `inside` 인자를 더했다. 기본 동작은 그대로다).
- **합성 표본**: driver는 damped random walk(시간척도 200일, 1000일), 진폭은 천체마다 관측된 W2의 상대 rms에 맞춘다. W1의 응답은 Lyu et al. (2019)의 W1 지연을 폭으로 하는 boxcar. 방문 시각·오차·오차 상관은 실제 값이다(`transfer9.mock_truth`).
- **null 범위**: 합성 표본 50개의 값에서 가장 작은 것과 가장 큰 것을 하나씩 뺀 범위(`run_stage9.band`).
- `results/s9_objects.csv`에는 목록 값(유형, z, L_bol)을 넣지 않았다. 이름과 맞춤 결과만 있다.

## 18. 10단계: 광학(ZTF) driver로 W1, W2의 지연을 따로 재기 (`docs/PREDICTIONS_10.md`)

코드: `src/irecho/driver10.py`(맞춤과 합성 광도곡선), `scripts/fetch_ztf10.py`, `scripts/run_stage10.py`, `scripts/stage10_analysis.py`, `scripts/stage10_verdict.py`, `scripts/plot_stage10.py`, `src/irecho/plots10.py`.

1. **광학 광도곡선**(`run_stage10.py build`): `catflags == 0`인 점이 가장 많은 `oid`의 측광점(`magerr` < 0.2) → 0.5일 안에 이어진 점을 역분산 가중 평균한 밤 단위 flux(등급 중앙값의 flux를 1로) → ±15일 안의 다른 밤의 중앙값에서 5 robust σ 넘게 벗어난 밤 제거(한 번).
2. **쓰는 방문**: WISE 방문 t는 [t − 700일, t + 100일]이 광학 광도곡선 안에 들고 그 안에 250일보다 긴 광학 공백이 없을 때만 쓴다. 모든 모형이 같은 방문을 쓴다. 8개 이상이어야 맞춘다.
3. **계절 공백**: 밤 단위 flux에 damped random walk + 백색 잡음(jitter) + 상수 평균(평평한 사전분포)의 Gaussian process를 최대가능도로 맞추고(천체마다 한 번), 2일 격자 위의 조건부 평균(밤 flux의 선형 결합)과 조건부 공분산 P를 구한다.
4. **맞춤**: F_W(t_i) = a + b (K ∗ D)(t_i). K는 평균 지연 τ의 kernel(기준: 폭이 |τ|인 top-hat). 예측의 불확도 b²·cPcᵀ를 잔차 공분산에 더하고 S = χ² + ln det로 τ를 5일 격자에서 고른다(일반화 최소제곱, b는 3번 갱신). 꼴을 바꾼 것(M2, MG), 지연 없는 항을 더한 것(MA), 가지 표시를 더한 것(ME)은 모두 기준 모형 M1과 직접 비교한다.
5. **Δτ의 구간**: 두 밴드를 함께 맞춘 10일 격자의 분포(driver 불확도의 공유 항과 같은 방문의 오차 상관 r 포함).
6. **검출**: 예측과 관측의 가중 상관 r_w가 다른 천체 20개의 driver로 바꿔 맞춘 값 모두보다 클 때. 검출 표본 D는 두 밴드 모두 검출이고 τ가 격자 끝이 아닌 천체.
7. **합성 표본**(`run_stage10.py ident`): 실제 ZTF 밤의 시각·오차·jitter와 실제 WISE 방문 시각·오차·오차 상관에, 진실 모형 21가지 × driver 2가지 × 50회차. 천체별 준비 계산(가중치 W와 공분산 V)은 관측 시각에만 달려 있어 회차마다 다시 쓰인다.
8. **판정**(`stage10_verdict.py`): 식별 가능성의 기준 → 실제 통계량을 합성 null 범위와 비교(`docs/PREDICTIONS_10.md` 5.2절, 6절).

## 19. 재현

```
uv venv --python 3.12 .venv && uv pip install --python .venv/bin/python -r requirements.txt
.venv/bin/python scripts/select_targets.py     # targets/ 생성 (NED, IRSA 조회)
.venv/bin/python scripts/fetch_data.py         # data/cache/ 채우기 (대상당 약 2분, 순차)
.venv/bin/python scripts/run_analysis.py       # results/ 생성 (네트워크 불필요)
.venv/bin/python scripts/independent_check.py  # NGC 3786 증가폭을 src/irecho 없이 다시 계산
.venv/bin/python scripts/select_field_stars.py # 대조군 주변 비교 별 목록 (3–4단계)
.venv/bin/python scripts/fetch_untimely.py     # unTimely 행 받기 (3–4단계)
.venv/bin/python scripts/fetch_ztf.py          # ZTF g, r 광도곡선 받기 (3단계 가설의 대안)
.venv/bin/python scripts/run_stage34.py        # 보정, 물리량, 판정식, 주입 시험 (네트워크 불필요, 약 3분)
.venv/bin/python scripts/build_sample.py targets   # 5a: 표본 (S1–S10, z 규칙)
.venv/bin/python scripts/build_sample.py stars     # 5a: 비교 별 (S11)
.venv/bin/python scripts/fetch_sample.py           # 5a: 시야별 NEOWISE-R, AllWISE MEP (이어 받기 가능)
.venv/bin/python scripts/fetch_untimely.py --sample            # 5a: 시야별 unTimely
.venv/bin/python scripts/fetch_untimely.py --register-sample   # 5a: unTimely 파일을 manifest에 등록
.venv/bin/python scripts/run_qa.py                 # 5a: blind 품질 확인, 오차 바닥, 광도곡선 동결
.venv/bin/python scripts/run_injection_stars.py    # 5a: 비교 별 주입 시험
.venv/bin/python scripts/run_injection_stars.py --floor-mmag 20   # 5a: 오차 바닥을 흉내 낸 변형 (10, 20, 50)
.venv/bin/python scripts/plot_qa.py                # 5a: QA 그림 (별과 집계값만)
.venv/bin/python scripts/run_breath_validation.py  # 5b: 합성 자료로 기울기 추정량 시험 (약 10분)
.venv/bin/python scripts/build_props.py            # 5b: Eddington 비와 N_H 표
.venv/bin/python scripts/run_injection_q1.py       # 5b: slow rise 정의의 완전도 (별, 약 25분)
.venv/bin/python scripts/fetch_ztf.py --sample --workers 4   # 5b: ZTF g, r (수 시간, 이어받기 가능)
.venv/bin/python scripts/fetch_ztf.py --register-sample      # 5b: ZTF 파일을 manifest에 등록
.venv/bin/python scripts/run_stage5b.py            # 5b: 사전 등록 분석 (약 10초)
.venv/bin/python scripts/plot_stage5b.py           # 5b: 그림
.venv/bin/python scripts/run_flares.py             # 6A: 섬광 후보의 신호 확인, 물리량, 분류
.venv/bin/python scripts/run_flares_other.py       # 6A: ZTF, ASAS-SN, Swift-BAT (캐시 필요)
IRECHO_SAMPLE=v2 .venv/bin/python scripts/build_sample.py targets   # 6B: 독립 표본 (0.04 < z <= 0.30)
IRECHO_SAMPLE=v2 .venv/bin/python scripts/build_sample.py stars     # 6B: 비교 별
.venv/bin/python scripts/run_d_lum_power.py prelim # 6B: 두 모형의 어림 검정력 (받기 전, 약 10분)
IRECHO_SAMPLE=v2 .venv/bin/python scripts/fetch_sample.py                     # 6B: NEOWISE-R, AllWISE MEP (약 2.6시간)
IRECHO_SAMPLE=v2 .venv/bin/python scripts/fetch_untimely.py --sample          # 6B: unTimely
IRECHO_SAMPLE=v2 .venv/bin/python scripts/fetch_untimely.py --register-sample # 6B: manifest 등록
IRECHO_SAMPLE=v2 .venv/bin/python scripts/run_qa.py                           # 6B: blind 품질 확인, 오차 바닥, 동결
.venv/bin/python scripts/run_d_lum_power.py final  # 6B: 실제 cadence와 오차로 문턱값 확정 (D 계산 전)
.venv/bin/python scripts/run_bandlag_validation.py final   # 6B: 밴드 간 지연 추정기의 검증
.venv/bin/python scripts/run_stage6b.py            # 6B: D와 광도의 관계, 판정
.venv/bin/python scripts/plot_stage6.py            # 6단계 그림
.venv/bin/python scripts/update_data_doc.py    # docs/DATA.md 해시 표 갱신
.venv/bin/python -m pytest -v                  # 단위 테스트 + 캐시 기반 회귀 테스트
bash scripts/xray_process.sh xrt10 && bash scripts/xray_process.sh xrt22   # 7: Swift XRT 재처리 (HEASoft 환경 필요, 아래 참고)
bash scripts/xray_process.sh nustar            # 7: NuSTAR 재처리 (NU_EXTRA="saamode=optimized tentacle=yes")
bash scripts/xray_process.sh regions && bash scripts/xray_process.sh products   # 7: 영역, 스펙트럼
python scripts/xray_fit.py fits                # 7: 맞춤 (HEASoft 환경의 python, 약 1시간)
python scripts/xray_fit.py plotdata; python scripts/xray_fit.py arf; python scripts/xray_fit.py lumcheck
python scripts/xray_fit.py mcmc; python scripts/xray_fit.py sims 200; python scripts/xray_fit.py posthoc
.venv/bin/python scripts/run_stage7.py         # 7: 판정 표 (커밋된 JSON만 읽음)
.venv/bin/python scripts/plot_stage7.py        # 7: 그림
python scripts/xray_7b.py light; python scripts/xray_7b.py fits   # 7b: 사후 보완 (HEASoft 환경)
.venv/bin/python scripts/run_stage7b.py        # 7b: 표와 그림
.venv/bin/python scripts/run_stage8_noise.py; IRECHO_SAMPLE=v2 .venv/bin/python scripts/run_stage8_noise.py   # 8: 오차 상관 (캐시 필요)
.venv/bin/python scripts/run_stage8.py null 2  # 8: null test (약 50분)
.venv/bin/python scripts/run_stage8.py direct  # 8: 직접 비교
.venv/bin/python scripts/run_stage8.py lag 2   # 8: lag 모형 재시험 (약 25분)
.venv/bin/python scripts/run_stage8.py select; .venv/bin/python scripts/plot_stage8.py   # 8: 보고 보완, 그림
.venv/bin/python scripts/run_stage9.py ident 2 # 9: B1 식별 가능성 (worker 2개로 약 70분)
.venv/bin/python scripts/run_stage9.py fit     # 9: 실제 자료 맞춤
.venv/bin/python scripts/run_stage9.py null 2  # 9: 맞춘 지연을 넣은 null 표본
.venv/bin/python scripts/run_stage9.py verdict; .venv/bin/python scripts/run_stage9.py mock5 2; .venv/bin/python scripts/plot_stage9.py   # 9: 판정, B5, 그림
.venv/bin/python scripts/fetch_ztf10.py objects; .venv/bin/python scripts/fetch_ztf10.py curves --workers 3; .venv/bin/python scripts/fetch_ztf10.py meta   # 10: ZTF (약 2시간)
.venv/bin/python scripts/run_stage10.py build            # 10: 밤 단위 광도곡선, 기준, 동결 (약 4분)
.venv/bin/python scripts/run_stage10.py ident --workers 3; .venv/bin/python scripts/run_stage10.py ident-summary   # 10: 식별 가능성 (약 37분)
.venv/bin/python scripts/run_stage10.py fit              # 10: 실제 자료 맞춤, driver 바꾸기, 시간 반전 (약 4분)
.venv/bin/python scripts/run_stage10.py fitted --workers 3; .venv/bin/python scripts/run_stage10.py fitted-summary   # 10: 맞춘 지연을 넣은 null 표본
.venv/bin/python scripts/stage10_verdict.py; .venv/bin/python scripts/run_stage10.py checks; .venv/bin/python scripts/plot_stage10.py   # 10: 판정, 민감도 (약 13분), 그림
```

난수는 기준 별·비교 별 추첨과 blackbody 온도의 Monte Carlo에만 쓴다. 모두 seed 20261112에서 나온다(`config.SEED`. 5a단계의 비교 별은 `[20261112, 5, BAT 목록 번호]`).

7단계의 도구 설치(프로젝트 폴더 안, git 밖):

```
conda create -y -p .heasoft --override-channels -c https://heasarc.gsfc.nasa.gov/FTP/software/conda/ -c conda-forge heasoft=6.37.1
.venv/bin/python scripts/fetch_xray.py         # 관측 자료 (322 MB)
# CALDB: data/cache/s7/caldb/ 에 caldb_setup_files.tar.Z, goodfiles_swift_mis.tar.Z, goodfiles_swift_xrt.tar.Z 를 풀고
source scripts/heasoft_env.sh && python scripts/fetch_caldb_subset.py 2022-05-01   # NuSTAR 보정 파일 94개 (1.95 GB)
```
