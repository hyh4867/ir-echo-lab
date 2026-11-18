# DATA: 자료 출처, 이용 조건, 받는 방법

원자료는 git에 넣지 않는다. 이 문서와 `data/manifest.json`에 받는 방법과 해시를 적는다.

## 1. 사용한 자료

| 자료 | IRSA 표 이름 | DOI | 이 프로젝트에서의 쓰임 |
|---|---|---|---|
| AllWISE Multiepoch Photometry Table | `allwise_p3as_mep` | [10.26131/IRSA134](https://doi.org/10.26131/IRSA134) | 2010–2011년 single-exposure 측광 (W1, W2) |
| NEOWISE-R Single Exposure (L1b) Source Table | `neowiser_p1bs_psd` | [10.26131/IRSA144](https://doi.org/10.26131/IRSA144) | 2014–2024년 single-exposure 측광 (W1, W2) |
| AllWISE Source Catalog | `allwise_p3as_psd` | 이 문서에서 DOI를 확인하지 않음 | 표적 좌표 고정, 기준 별 선정 |
| NED Object Lookup | (서비스) | | 은하 이름 → 좌표, redshift |
| ZTF 공개 광도곡선 (g, r) | (IRSA ZTF lightcurve 서비스) | 이 문서에서 DOI를 확인하지 않음 | 3단계: NGC 3786의 광학 변화를 중적외선과 비교(PREDICTIONS.md의 대안). NGC 2273에도 추가로 적용 |
| unWISE Time-Domain Catalog (unTimely, neo7) | (TAP 표 없음. S3 parquet) | [10.26131/IRSA580](https://doi.org/10.26131/IRSA580) | 3–4단계: NGC 3786과 대조군의 방문별 coadd 측광(2010–2020), profile-fit 결과와 비교 |

- 두 측광 표의 DOI와 설명은 IPAC DOI 페이지(https://www.ipac.caltech.edu/doi/irsa/10.26131/IRSA134, .../IRSA144)에서 2026-11-12에 확인했다.
- AllWISE MEP는 AllWISE Source Catalog와 Reject Table의 각 항목 위치에서 single-exposure 영상마다 잰 flux 모음이다(DOI 페이지 설명: "time-tagged fluxes measured on the individual Single-exposure image sets forced at the position of each deep source extraction").
- NEOWISE-R 표는 7.7초 W1, W2 single-exposure 영상에서 독립적으로 검출·측정한 항목 모음이다.
- 이 프로젝트가 받은 NGC 3786 자료의 실제 시간 범위는 AllWISE MEP 2010-05-28 ~ 2010-12-05, NEOWISE-R 2014-05-29 ~ 2024-05-17이다(캐시 파일의 mjd 최솟값·최댓값). 2011–2013년은 WISE가 동면 중이라 자료가 없다.

- unTimely(Meisner et al. 2023)는 WISE/NEOWISE의 방문별 unWISE coadd 영상에서 `crowdsource`로 PSF 측광한 목록이다. neo7 판은 2010년부터 2020년 12월까지를 담는다. IRSA DOI 페이지와 버킷의 `README.txt`(Version Date 2023-05-20)에서 2026-11-12에 확인했다. README는 전체 크기를 4.4 TB, 43,459,534,445행, 32열로 적는다.
- 3–4단계에서는 대조군 은하 8개 각각의 주변 비교 별 10개(모두 80개)의 NEOWISE-R 자료도 같은 방법(3″ cone)으로 받았다. 목록은 `targets/field_stars.csv`, 선정 기준은 PIPELINE.md 7절에 있다.

- ZTF 광도곡선은 IRSA의 ZTF lightcurve 서비스가 공개 Data Release에서 돌려주는 PSF 측광이다. IRSA ZTF 페이지(https://irsa.ipac.caltech.edu/Missions/ztf.html, 2026-11-12 확인)는 최신 공개판을 Data Release 24(2026-01-22)로 적는다. 받은 자료의 실제 범위는 NGC 3786이 2018-03-25 – 2025-10-09, NGC 2273이 2018-03-27 – 2025-10-20이다(캐시 파일의 mjd).

- 5a단계(2026-11-14)에서 더한 자료:

| 자료 | 받은 곳 | 이 프로젝트에서의 쓰임 |
|---|---|---|
| Swift-BAT 105-month 목록 (Oh et al. 2018) | NASA GSFC `https://swift.gsfc.nasa.gov/bat_survey/bs105mon/data/BAT_105m_catalog_07jul2019.txt` | 표본의 모집단: 이름, counterpart 좌표, 적색편이, 분류 |
| BASS DR2 표 (Koss et al. 2022, VizieR J/ApJS/261/2 table9) | VizieR CfA 미러 `https://vizier.cfa.harvard.edu/viz-bin/asu-tsv` | DR2 분류(Sy1, Sy1.9, Sy2, beamed), 적색편이, 거리, black hole 질량 |
| AllWISE Source Catalog | IRSA TAP | 표본 후보 주변 5′(기준 S6–S10), 비교 별 시야 1.5° |
| NEOWISE-R, AllWISE MEP | IRSA TAP (위치 목록을 올려 시야 단위로 조회) | 표본 231개 시야의 표적과 비교 별의 single-exposure 측광 |
| unTimely neo7 | 공개 S3 parquet | 표본 표적과, 같은 partition에 있는 비교 별 |

  표본의 정의와 두 목록의 이용 조건은 `docs/SAMPLE.md` 2절과 4절에 있다.

- 5b단계(2026-11-17)에서 더한 자료:

| 자료 | 받은 곳 | 이 프로젝트에서의 쓰임 |
|---|---|---|
| BASS DR2 표의 `logLbol`, `logEdd` 열 (Koss et al. 2022) | 5a단계에 받은 파일에 이미 있음 | Eddington 비(문턱 0.02 위/아래), 광도 |
| BASS DR1 X선 목록 (Ricci et al. 2017, ApJS 233, 17. VizieR J/ApJS/233/17의 `agns`, `Xraydata` 표) | VizieR CfA 미러 `https://vizier.cfa.harvard.edu/viz-bin/asu-tsv` | column density `logNH`. 위치 5″로 맞춤 |
| ZTF DR24 objects 표 (`ztf_objects_dr24`) | IRSA TAP (위치 목록을 올려 조회 하나) | 표본 231개 위치에 ZTF 천체가 있는지 |
| ZTF 공개 광도곡선 g, r | IRSA `nph_light_curves` | Q3: 적외선이 변할 때 광학이 같이 변하는지 |

  BASS DR1 목록의 이용 조건도 VizieR를 통해 받은 것이어서 5a단계와 같은 상태다(조건 전문 미확인, `docs/SAMPLE.md` 4절). 원본 표는 저장소에 넣지 않았고, `targets/sample_v1_props.csv`에는 천체별 `logEdd`, `logLbol`, `logMBH`, `logNH` 값만 출처와 함께 옮겼다.

- 6단계 A(2026-11-17)에서 더한 자료 (섬광 후보 4개에 대해서만):

| 자료 | 받은 곳 | 이 프로젝트에서의 쓰임 |
|---|---|---|
| ASAS-SN Sky Patrol v2 광도곡선 (V, g) | 공개 client `skypatrol`(pyasassn)의 `cone_search(..., download=True)`. 로그인 없음 | 섬광 시기의 광학 변화와 상한 |
| Swift-BAT 157개월 survey의 crab-weighted 월별 광도곡선 (14–195 keV, 2004-12 ~ 2017-12) | NASA GSFC `https://swift.gsfc.nasa.gov/bat_survey/bs157mon/crab_lc/` | 섬광 전후의 hard X-ray 변화 (대부분의 섬광은 이 기간 뒤라 덮지 못함) |
| HEASARC 관측 목록 (Swift, XMM-Newton, Chandra, NuSTAR의 master 표) | HEASARC TAP `https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync` | 공개 X선 관측이 있는지와 날짜. 자료 자체는 받지 않음 |
| SIMBAD 참고문헌 목록 | SIMBAD TAP 미러 `simbad.harvard.edu` | 천체별 문헌 제목 확인 |
| MIRONG 목록 (Jiang et al. 2021, VizieR J/ApJS/252/32), Flaires 목록 (Necker et al. 2025, VizieR J/A+A/695/A228) | VizieR CfA 미러 | 섬광 후보가 기존 중적외선 섬광 목록에 있는지 |

  - ASAS-SN client는 프로젝트의 `requirements.txt`에 넣지 않았다. 저장소 밖의 임시 환경에서 한 번 받아 `data/cache/s6/asassn/`에 CSV로 두었다(git 밖). 다시 받으려면 `pip install skypatrol` 뒤 `SkyPatrolClient().cone_search(ra_deg, dec_deg, radius=10, units='arcsec', catalog='master_list', download=True)`.
  - ASAS-SN과 BAT 157개월 자료의 이용 조건 전문은 확인하지 않았다. 원자료는 저장소에 넣지 않았고 요약 수치만 `results/s6a_other_bands.json`에 있다. 인용: ASAS-SN은 Shappee et al. (2014), Kochanek et al. (2017), Sky Patrol v2는 Hart et al. (2023). BAT 157개월은 survey 페이지가 가리키는 arXiv:2506.04109. 이 네 문헌은 서지 정보만 옮겼고 읽지 않았다.
  - TNS는 공개 검색 주소가 403을 돌려줘 조회하지 못했다. ATLAS forced photometry는 로그인이 필요해 쓰지 않았다.

- 6단계 B(2026-11-17)에서 더한 자료:

| 자료 | 받은 곳 | 이 프로젝트에서의 쓰임 |
|---|---|---|
| AllWISE Source Catalog | IRSA TAP | 독립 표본 후보 주변 5′(기준 S6–S10), 비교 별 시야 1.5° |
| NEOWISE-R, AllWISE MEP | IRSA TAP (위치 목록을 올려 시야 단위로 조회) | sample_v2 210개 시야의 표적과 비교 별 2,100개의 single-exposure 측광 |
| unTimely neo7 | 공개 S3 parquet | sample_v2 표적과, 같은 partition에 있는 비교 별 |

  모집단은 5a단계와 같은 두 목록(Swift-BAT 105개월, BASS DR2)이고 새로 받은 목록은 없다. 받는 방법은 5a단계와 같고 명령 앞에 `IRECHO_SAMPLE=v2`를 붙인다. 캐시 위치는 `data/cache/s6b/neowise/`, `s6b/allwise_mep/`, `s6b/untimely/`, `data/cache/allwise_cat/s6b_targets_<번호>.csv`, `s6b_fields_<번호>.csv`다. 표본의 정의는 `docs/PREDICTIONS_6.md` 4절, 실행 기록은 `docs/SAMPLE.md` 10절에 있다.

- 7단계(2026-11-17)에서 더한 자료: CGCG 312-012의 X선 관측. 모두 로그인 없이 받는 공개 자료다.

| 자료 | 관측 ID | 날짜 (UTC) | 노출 | 보관소의 처리 | 이 프로젝트에서 다시 한 처리 |
|---|---|---|---|---|---|
| Swift XRT (PC mode) | 00040946001 | 2010-12-29 00:21 – 16:25 | 10.31 ks | 처리 버전 3.17.06 | `xrtpipeline` (swxrtdas 3.7.2, HEASoft 6.37.1), XRT CALDB 20260327 |
| Swift XRT (PC mode) | 00081051001 | 2022-05-01 04:01 – 12:17 | 6.93 ks | 처리 버전 3.18.19 | 위와 같음 |
| NuSTAR FPMA, FPMB | 60160365002 | 2022-05-01 (MJD 59700.11–59700.59) | FPMA 22.64 ks, FPMB 22.44 ks (SAA 필터 뒤. 기본 설정으로는 23.05, 22.82 ks) | NuSTARDAS 2.1.1, CALDB 20220215 | `nupipeline` (NuSTARDAS 2.1.6, HEASoft 6.37.1, `saamode=optimized tentacle=yes`), CALDB 색인 20261102 |
| Swift-BAT 157개월 평균 스펙트럼 (8채널, 14–195 keV)과 응답 행렬, 목록 | SWIFT J0929.7+6232 | 2004-12 ~ 2017-12 | survey 평균 | survey 페이지의 제품 | 그대로 사용 |

  - 받은 곳: NuSTAR와 Swift는 HEASARC 보관소 `https://heasarc.gsfc.nasa.gov/FTP/nustar/data/obs/01/6/60160365002/`, `https://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/2010_12/00040946001/`, `.../2022_05/00081051001/`. BAT는 `https://swift.gsfc.nasa.gov/bat_survey/bs157mon/`의 `spectra/SWIFT_J0929.7+6232.pha`, `data/BAT_157m_swiftbat_survey_full_157m.rsp`, `data/BAT_157m_catalog_20250525.txt`.
  - 받는 방법: `python scripts/fetch_xray.py` (2026-11-17 04:35–04:49 UTC, 파일 123개, 321.8 MB). NuSTAR는 `auxil/`, `hk/`, `event_uf/`, `event_cl/`, Swift는 `auxil/`, `xrt/`만 받았다(UVOT와 BAT event 자료는 받지 않음). 캐시 위치는 `data/cache/s7/`이고 git에 넣지 않는다. 파일별 URL, 크기, 받은 시각, SHA-256은 `data/manifest.json`의 `s7/`로 시작하는 항목에 있다.
  - 관측 목록의 값(관측 ID, 날짜, 노출, 보관소 처리 버전)은 HEASARC TAP의 `swiftmastr`, `numaster` 표에서 읽었다.
  - 주요 입력 파일의 SHA-256:

| 캐시 파일 (data/cache/ 아래) | 크기 (byte) | SHA-256 |
|---|---|---|
| `s7/nustar/60160365002/event_uf/nu60160365002A_uf.evt.gz` | 60,789,253 | `eaba64787fcc37e09b8d7d099cf6b266311ef148972856114ca63193aee25cd4` |
| `s7/nustar/60160365002/event_uf/nu60160365002B_uf.evt.gz` | 91,143,832 | `99b90adb2787f967682acf930364df75906dc112323c3057d267aa269958a936` |
| `s7/nustar/60160365002/auxil/nu60160365002_att.fits.gz` | 12,575,727 | `27d422ba1e4ab8bf959886a611a51c5b3a277de8cb364b7f2886e2381ed93d73` |
| `s7/nustar/60160365002/auxil/nu60160365002_orb.fits.gz` | 2,005,012 | `fe24c7cc0f2f60d7aeca75eb1e71ed96313015facd6a783ea05fb5ba2ab8c6e1` |
| `s7/swift/00040946001/xrt/event/sw00040946001xpcw3po_uf.evt.gz` | 22,539,792 | `2cb4ff0d4e7cc2ce0e08ba087301e9644bc0ab5e7451b211a571740c7aed3bf4` |
| `s7/swift/00040946001/auxil/sw00040946001sat.fits.gz` | 198,046 | `00930518d6b4c43a1927f4ce2f6467dd0533e639cb5d58f4373115dd34b64cb5` |
| `s7/swift/00081051001/xrt/event/sw00081051001xpcw3po_uf.evt.gz` | 39,196,934 | `ca39034f4a22d6c3bf879ff24a5836b9581310fd2e768b7ff5fecf3731285c4d` |
| `s7/swift/00081051001/auxil/sw00081051001sat.fits.gz` | 158,803 | `e2fba0c2ddd7bcdfde0601f26c8c305ff88ecfa0e3848fd2ae612d1e61988159` |
| `s7/bat157/SWIFT_J0929.7+6232.pha` | 17,280 | `94fdc98c2086a4be7b8107ca2af0930db800b84da4344a79ba3390be1ee67794` |
| `s7/bat157/BAT_157m_swiftbat_survey_full_157m.rsp` | 31,680 | `f9113bc147d84508f01ebdf5e6d9b0314c757d68182b15ac4730d20dd918c427` |
| `s7/bat157/BAT_157m_catalog_20250525.txt` | 470,669 | `9901ec0dc86cd67781383242b4808f4d87e0ec396fb23eb134c30b4b177b0249` |

  - **도구와 보정 파일** (git에 넣지 않음):
    - HEASoft 6.37.1을 conda로 프로젝트 폴더 안의 `.heasoft/`에 설치했다(1.6 GB. 채널 `https://heasarc.gsfc.nasa.gov/FTP/software/conda/`와 conda-forge). XSPEC 13.0.1, PyXspec 2.2.0이 들어 있다.
    - CALDB는 `data/cache/s7/caldb/`에 두었다. Swift XRT는 `goodfiles_swift_xrt.tar.Z`(660,443,737 byte, 색인 20260327)와 `goodfiles_swift_mis.tar.Z` 전체. NuSTAR는 전체 묶음(2.3 GB)을 받지 않고, 색인(`caldb.indx`, 버전 20261102)에서 2022-05-01에 유효한 파일 94개(1.95 GB)만 받았다(`scripts/fetch_caldb_subset.py`. 'nocut' 응답 행렬 8개는 이 처리에서 쓰이지 않아 뺐다). 묶음과 파일별 SHA-256은 `results/s7_caldb.json`에 있다.
  - 문헌·목록 조회(`scripts/xray_literature.py`, `results/s7_literature.json`): SIMBAD TAP 미러(참고문헌 제목), VizieR CfA 미러의 J/ApJS/233/17(Ricci et al. 2017)과 IX/58(2SXPS, Evans et al. 2020).
  - 이용 조건: HEASARC의 안내(`https://heasarc.gsfc.nasa.gov/docs/heasarc/acknow.html`, 2026-11-17 확인)는 다음 문구를 쓰라고 한다. "This research has made use of data and/or software provided by the High Energy Astrophysics Science Archive Research Center (HEASARC), which is a service of the Astrophysics Science Division at NASA/GSFC." 그 쪽에는 이용 제한에 관한 문장이 없다. NuSTAR와 Swift 각 임무의 이용 조건·인용 안내 전문은 확인하지 않았다. 원자료와 처리된 event 파일, 스펙트럼, 응답 행렬은 저장소에 넣지 않았고 맞춤 결과의 수치와 그림만 `results/s7_*`에 있다.

- 7b와 8단계(2026-11-17)에서는 새 자료를 받지 않았다. 7b에서 2SXPS(VizieR IX/58, CfA 미러)의 이 천체 행을 노출·관측 수·밴드별 rate 열까지 다시 조회했고(공개 저장소의 결과 파일과 문서에는 목록의 값을 옮기지 않았고, 자체 측정값과 차이만 `results/s7b_light.json`의 `comparison_with_2sxps`에 있다. 다만 비교에 쓴 행의 값은 `scripts/xray_7b.py`에 적혀 있다), 8단계는 이미 받은 캐시(5a·6B의 NEOWISE-R, unTimely)와 동결 광도곡선만 썼다. 8단계가 캐시에서 다시 뽑은 unTimely 광도곡선은 `data/cache/s8/`에 있고 git에 넣지 않는다. 문헌의 식을 확인한 기록은 `docs/LITERATURE.md` 11절에 있다.

- 9단계(2026-11-18)에서도 새 관측 자료를 받지 않았다. 동결 광도곡선과 8단계의 오차 상관만 썼다. 합성 표본의 회차별 값은 `data/cache/s9/`에 있고 git에 넣지 않는다. 문헌의 식을 확인하려고 arXiv의 논문 PDF 한 편(1909.11101v2)을 받아 본문을 읽었고(저장소에 넣지 않음), 확인한 내용은 `docs/LITERATURE.md` 12절에 있다. `results/s9_objects.csv`에는 목록 값(유형, z, L_bol)이 없다.

## 2. 이용 조건

IRSA "Terms of Use of Data Hosted by IRSA" (https://irsa.ipac.caltech.edu/data_use_terms.html, 최종 갱신 2023-02-08, 2026-11-12 확인):

> 1.1 Public Data — Most data served by IRSA is public, with no usage restrictions.

WISE/NEOWISE 표는 로그인 없이 받을 수 있는 공개 자료다. 같은 문서는 저작권이 있는 자료로 Digital Sky Survey만 따로 든다. 이 프로젝트는 DSS를 쓰지 않는다.

같은 문서 2절은 출판물에 acknowledgment와 citation을 넣는 것을 "best practice"로 안내한다. 이 저장소는 논문이 아니지만 아래 문구를 그대로 옮겨 둔다.

unTimely도 IRSA가 호스팅하는 공개 자료이며 같은 약관이 적용된다. 공개 버킷 `nasa-irsa-wise`는 자격 증명 없이(anonymous) 읽을 수 있다고 IRSA "Cloud Data Access" 페이지(https://irsa.ipac.caltech.edu/cloud_access/)가 안내한다. unTimely에만 적용되는 별도의 이용 제한 문구는 DOI 페이지, 개요 페이지, 버킷 README에서 찾지 못했다.

ZTF 공개 Data Release도 IRSA가 로그인 없이 제공한다. IRSA ZTF 페이지에는 acknowledgement 문구와 대표 논문 두 편이 있고, 그 밖의 이용 제한 문구는 그 페이지에서 찾지 못했다. ZTF 프로젝트 자체의 자료 이용 정책 문서는 따로 확인하지 않았다.

## 3. 인용 문구

IRSA (https://irsa.ipac.caltech.edu/ack.html):

> This research has made use of the NASA/IPAC Infrared Science Archive, which is funded by the National Aeronautics and Space Administration and operated by the California Institute of Technology.

AllWISE (AllWISE Explanatory Supplement I.6.b):

> This publication makes use of data products from the Wide-field Infrared Survey Explorer, which is a joint project of the University of California, Los Angeles, and the Jet Propulsion Laboratory/California Institute of Technology, and NEOWISE, which is a project of the Jet Propulsion Laboratory/California Institute of Technology. WISE and NEOWISE are funded by the National Aeronautics and Space Administration.

NEOWISE (NEOWISE Explanatory Supplement I.6.b, 최종 갱신 2024-08-18):

> This publication makes use of data products from the Near-Earth Object Wide-field Infrared Survey Explorer (NEOWISE), which is a joint project of the Jet Propulsion Laboratory/California Institute of Technology and the University of California Los Angeles. NEOWISE is funded by the National Aeronautics and Space Administration.

함께 인용할 논문과 DOI:

- Wright, E. L., et al. 2010, AJ, 140, 1868 (WISE)
- Mainzer, A., et al. 2011, ApJ, 731, 53 (NEOWISE. AllWISE DOI 페이지가 인용)
- Mainzer, A., et al. 2014, ApJ, 792, 30 (NEOWISE Reactivation)
- AllWISE Multiepoch Photometry Table, WISE Team, IPAC, doi:10.26131/IRSA134
- NEOWISE-R Single Exposure (L1b) Source Table, NEOWISE Team, IPAC, doi:10.26131/IRSA144

unTimely (버킷 `README.txt`와 https://irsa.ipac.caltech.edu/data/WISE/unWISE/overview.html, 2026-11-12 확인):

> If you use the unWISE Time-Domain Catalog, please cite the dataset Digital Object Identifier (DOI): 10.26131/IRSA580 (https://www.ipac.caltech.edu/doi/irsa/10.26131/IRSA580). In addition, please follow the project's acknowledgement guidelines (https://irsa.ipac.caltech.edu/data/WISE/unWISE/overview.html) and the IRSA acknowledgement guidelines (https://irsa.ipac.caltech.edu/ack.html).

- Meisner, A. M., Caselden, D., Schlafly, E. F., & Kiwy, F. 2023, AJ, 165, 36, doi:10.3847/1538-3881/aca2ab (unTimely. 개요 페이지가 인용을 요청)
- unWISE Time-Domain Catalog, Meisner et al., IPAC, doi:10.26131/IRSA580
- Schlafly, E. F., Meisner, A. M., & Green, G. M. 2019, ApJS, 240, 30 (unWISE Catalog. `flags_info` 비트 정의의 출처)

ZTF (IRSA ZTF 페이지의 "Acknowledgement" 항목, 2026-11-12 확인):

> Supported by the National Science Foundation and a collaboration including current partners Caltech, IPAC, the Oskar Klein Center at Stockholm University, the University of Maryland, the University of California at Berkeley, the University of Wisconsin at Milwaukee, University of Warwick, Ruhr University, Cornell University, Northwestern University and Drexel University. Operations are conducted by COO, IPAC, and UW.

- Bellm, E. C., et al. 2019, PASP, 131, 018002 (ZTF 사양과 survey 설계. 같은 페이지의 "Canonical Papers")
- Masci, F. J., et al. 2019, PASP, 131, 018003 (ZTF 자료 처리. 같은 페이지의 "Canonical Papers")

NED를 쓴 부분(이름 조회)에 대한 NED의 공식 acknowledgment 문구는 이 문서를 쓰면서 확인하지 않았다.

## 4. 받는 방법

```
.venv/bin/python scripts/select_targets.py   # 표적 목록 (이미 targets/에 고정되어 있으면 생략)
.venv/bin/python scripts/fetch_data.py       # data/cache/에 CSV 저장
```

- 접근 경로: IRSA TAP sync, `https://irsa.ipac.caltech.edu/TAP/sync?QUERY=<ADQL>&FORMAT=csv`. 인증이 필요 없다.
- 표적 하나에 조회 두 번(NEOWISE-R, AllWISE MEP)을 보낸다. 조건은 표적 좌표 중심 반지름 3″ cone이다. 예:

```sql
SELECT ra, dec, mjd, scan_id, frame_num, source_id, w1mpro, w1sigmpro, w1snr, w1rchi2,
       w2mpro, w2sigmpro, w2snr, w2rchi2, w1sat, w2sat, nb, na, cc_flags, det_bit, ph_qual,
       sso_flg, qual_frame, qi_fact, saa_sep, moon_masked, w1frtr, w2frtr, allwise_cntr
FROM neowiser_p1bs_psd
WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', 174.9271290, 31.9094276, 0.0008333)) = 1
ORDER BY mjd
```

- 각 파일의 ADQL 전문은 `data/manifest.json`의 `adql` 항목과 캐시 옆의 `.meta.json`에 있다.
- 2026-11-12 실행에서는 조회 하나에 1–2분이 걸렸다. 조회는 순차로 보내며, 받은 파일은 다시 받지 않는다.
- 캐시 위치: `data/cache/neowise/<이름>.csv`, `data/cache/allwise_mep/<이름>.csv`, `data/cache/allwise_cat/<이름>.csv`, `data/cache/untimely/<이름>.csv`, `data/cache/ztf/<이름>_<밴드>.csv`. `.gitignore`로 제외한다.

### 5a단계: 시야 단위로 받는 방법

```
.venv/bin/python scripts/build_sample.py targets              # 두 목록과 AllWISE 목록 → targets/sample_v1.csv
.venv/bin/python scripts/build_sample.py stars                # 비교 별 → targets/sample_v1_stars.csv
.venv/bin/python scripts/fetch_sample.py                      # data/cache/s5/neowise/, s5/allwise_mep/
.venv/bin/python scripts/fetch_untimely.py --sample           # data/cache/s5/untimely/
.venv/bin/python scripts/fetch_untimely.py --register-sample  # unTimely 파일을 manifest에 등록
```

- 시야 하나(표적 + 비교 별 최대 10개)를 조회 하나로 받는다. 위치 목록을 VOTable로 올리고(`UPLOAD=pos,param:pos`) 다음 꼴의 ADQL을 POST로 보낸다. 조건은 위와 같은 반지름 3″다.

```sql
SELECT p.pid, t.ra, t.dec, t.mjd, ... FROM neowiser_p1bs_psd AS t, TAP_UPLOAD.pos AS p
WHERE CONTAINS(POINT('ICRS', t.ra, t.dec), CIRCLE('ICRS', p.pra, p.pdec, 0.0008333)) = 1
```

- `pid` 0이 표적, 1부터가 `targets/sample_v1_stars.csv`에 적힌 순서의 비교 별이다. manifest의 `upload` 항목에 올린 위치 목록 전체가 있다.
- 받은 CSV는 `pid`, `mjd`, 식별자 순으로 정렬해 저장한다(서버가 행 순서를 보장하지 않는다).
- 소요 시간과 방식별 비교는 `docs/SAMPLE.md` 5절에 있다.
- 캐시 위치: `data/cache/s5/neowise/<이름>.csv`, `data/cache/s5/allwise_mep/<이름>.csv`, `data/cache/s5/untimely/<이름>.csv`, `data/cache/allwise_cat/s5_targets_<번호>.csv`, `data/cache/allwise_cat/s5_fields_<번호>.csv`, `data/cache/catalogs/`.

### unTimely 받는 방법

```
.venv/bin/python scripts/select_field_stars.py   # 대조군 주변 비교 별 목록 (이미 있으면 생략)
.venv/bin/python scripts/fetch_data.py           # 비교 별 80개의 NEOWISE-R 표 추가 (별 하나에 약 2분)
.venv/bin/python scripts/fetch_untimely.py       # data/cache/untimely/<이름>.csv
```

- unTimely에는 cone search 서비스가 없다(2026-11-12에 IRSA TAP의 `TAP_SCHEMA.tables`를 조회했을 때 이 목록의 표가 없었다). IRSA가 만든 parquet 사본이 HEALPix order 5(NESTED)로 나뉘어 `s3://nasa-irsa-wise/unwise/neo7/catalogs/time_domain/healpix_k5/unwise-neo7-time_domain-healpix_k5.parquet/healpix_k0=<k0>/healpix_k5=<k5>/part0.snappy.parquet`에 있다(region us-west-2).
- `scripts/fetch_untimely.py`는 표적이 속한 partition 하나를 열어 좌표 열만 먼저 읽고, 표적에서 12″ 안의 행만 골라 나머지 열을 한 열씩 읽는다. partition 하나는 약 170–260만 행(NGC 3786의 partition 파일은 203 MB)이지만 저장하는 것은 표적당 32–34행이다. 표적 하나에 14–24초가 걸렸다.
- partition 번호는 `src/irecho/healpix.py`가 계산한다. NGC 3786은 `healpix_k5=7148`이고, 그 파일의 좌표 범위(RA 172.97–175.78°, Dec 30.00–32.79°)가 표적을 포함함을 parquet 통계에서 확인했다.
- 각 캐시 파일 옆의 `.meta.json`과 `data/manifest.json`에 partition 경로, 선택 조건, 행 수, SHA-256을 적는다.

### ZTF 받는 방법

```
.venv/bin/python scripts/fetch_ztf.py            # data/cache/ztf/<이름>_<밴드>.csv
```

- 접근 경로: `https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?POS=CIRCLE <ra> <dec> 0.0008333&BANDNAME=<g|r>&FORMAT=csv`. 반지름 3″, 인증 불필요. 밴드 하나에 10–50초가 걸렸다.
- 받은 행: NGC 3786 g 747, r 933. NGC 2273 g 630, r 633.
- 분석에 쓰는 행: `catflags = 0`이고, 그런 행이 가장 많은 object ID 하나(`src/irecho/ztf.py`). object ID가 다르면 field·CCD 조합이 달라 영점이 다르기 때문이다. NGC 3786은 g 588행, r 746행을 썼다.

### 5b단계: 표본의 ZTF 받는 방법

```
.venv/bin/python scripts/fetch_ztf.py --sample --workers 4 --order data/cache/s5/ztf_priority.txt
.venv/bin/python scripts/fetch_ztf.py --register-sample
```

- 먼저 `ztf_objects_dr24` 표를 표본 231개 위치(반지름 3″)로 한 번 조회한다(17초). ZTF 천체가 있는 표적만 광도곡선을 요청한다.
- 광도곡선 서비스는 요청 하나에 1–4분이 걸렸다. 천체 ID 여러 개를 한 요청에 넣어도 빨라지지 않았다(12개에 10분 넘게 걸려 중단). 그래서 요청을 동시에 여러 개 보내고(`--workers`), 판정에 필요한 천체를 먼저 받았다(`--order`).
- 파일은 `data/cache/s5/ztf/<이름>_<밴드>.csv`, 파일마다 조회 조건·시각·행 수·SHA-256을 담은 sidecar가 있고 `--register-sample`이 manifest에 옮긴다.

## 5. 해시

다시 받은 파일이 분석에 쓴 파일과 같은지는 `tests/test_regression_cached.py::test_cached_files_match_manifest_hashes`가 확인한다. IRSA가 표를 재처리하거나 행 순서·숫자 표기가 바뀌면 해시가 달라질 수 있다. 그때는 행 수와 값부터 비교해야 한다.

아래 표는 `scripts/update_data_doc.py`가 `data/manifest.json`에서 만든다.

<!-- manifest:start -->
| 캐시 파일 (data/cache/ 아래) | 행 수 | 받은 시각 (UTC) | SHA-256 |
|---|---|---|---|
| `allwise_cat/field_mrk573.csv` | 10781 | 2026-11-12T13:47:47Z | `7f6a7b27086dda862067497dfb7edbeb9d26948048cb96af22ed5d7240b94575` |
| `allwise_cat/field_ngc1386.csv` | 11980 | 2026-11-12T13:48:17Z | `e70a6df212b0bd6907e10bdc30ccd08c84c5a85fdfd68972931ea2cea435e501` |
| `allwise_cat/field_ngc2273.csv` | 23606 | 2026-11-12T13:49:03Z | `e6dc24b59740fe9c8fccd32005bc555e486fa1e2e72613d028466e4f5f444f1c` |
| `allwise_cat/field_ngc3081.csv` | 23599 | 2026-11-12T13:49:27Z | `550de6b8b9bda7f6b40143ee9b5dda39ba0efa8278a710c3f00c3ef561898938` |
| `allwise_cat/field_ngc3393.csv` | 20564 | 2026-11-12T13:50:17Z | `a68cd92631a4c593e70493e3594253712af872d24daf370194a30d8d55f0a486` |
| `allwise_cat/field_ngc3786.csv` | 29228 | 2026-11-12T11:12:38Z | `d394d3bf63144d3a958c8cb487332721fb272466c87adcae73460815666344a9` |
| `allwise_cat/field_ngc4941.csv` | 13085 | 2026-11-12T13:51:12Z | `869fbdcbf1008ba81a34997aa296bfaf8cd52cdffcb3822125fea095c4b867fc` |
| `allwise_cat/field_ngc5347.csv` | 10017 | 2026-11-12T13:51:37Z | `7b59e2dc2986e184d3974e93d0ba1338c9e41ea68170fd8e8ff89186464b6074` |
| `allwise_cat/field_ngc5728.csv` | 21231 | 2026-11-12T13:52:10Z | `df5c948b790fada65936ff4a0d0db739901b14fd9624efa96aeed4193a1878b2` |
| `allwise_cat/mrk573_cat.csv` | 1 | 2026-11-12T11:12:46Z | `2b6af9014fc6b3d31a9461889f5ed7efd511e895367ef25690dfbccda082965a` |
| `allwise_cat/ngc1386_cat.csv` | 1 | 2026-11-12T11:13:11Z | `0f4db22afae77f117f11576f5ff921551f9f9b4840a10b481f98569ff72255a6` |
| `allwise_cat/ngc2273_cat.csv` | 1 | 2026-11-12T11:13:31Z | `5c2665cb4876485665df200f7780ebc89818948f9cc1e375eb4f53aac7ff8e14` |
| `allwise_cat/ngc3081_cat.csv` | 1 | 2026-11-12T11:13:51Z | `96abeaa3baff1bf86c5f528ee4e281a8a823ebafa4dd0c055ecd11a384634bea` |
| `allwise_cat/ngc3393_cat.csv` | 1 | 2026-11-12T11:14:08Z | `10f0c9d84ee25349d59272e13caf1468a76056eb2d86dd494ba3770bc425047b` |
| `allwise_cat/ngc3786_cat.csv` | 1 | 2026-11-12T11:12:07Z | `93e6e6f6f2e0f4784daf7fd5bd52d9f8dd8029c7e67f2d5b27d6e6db90544d99` |
| `allwise_cat/ngc4941_cat.csv` | 1 | 2026-11-12T11:14:24Z | `af025199c433b904c951b4474f36554c111b685e629683f5e96fed6b90f2a323` |
| `allwise_cat/ngc5347_cat.csv` | 1 | 2026-11-12T11:14:39Z | `9ea05b68a9e1a1dd162326bd5bf770026fec635924b1a71efc6d8299c1ec89c8` |
| `allwise_cat/ngc5728_cat.csv` | 1 | 2026-11-12T11:14:48Z | `262d58a3c37c7c472f22fdb659076b066b4186371b58f3e791f3aab14dac4723` |
| `allwise_mep/mrk573.csv` | 29 | 2026-11-12T12:18:57Z | `84962570666c568fcadf3dc86e2d9811d0770a0be958a8c4d0a1db12a0649b9f` |
| `allwise_mep/ngc1386.csv` | 64 | 2026-11-12T12:20:53Z | `a7410ac821bb098f666b4967e3195aa649badd43bc66dc3a7cdde373be264e1d` |
| `allwise_mep/ngc2273.csv` | 52 | 2026-11-12T12:23:03Z | `b0d86ce7545a676e8f82012f1c1cf45686c702bb1127bd502e69ddc98ab3c5f3` |
| `allwise_mep/ngc3081.csv` | 32 | 2026-11-12T12:25:15Z | `dad3f0b66b6cc86e4fbede471a74960f76ff6dc4491ba58704c72decf4aff40f` |
| `allwise_mep/ngc3393.csv` | 25 | 2026-11-12T12:27:08Z | `af958fffdde1f1006722f376c6aee292d352bfbb515999008c7446597095d973` |
| `allwise_mep/ngc3786.csv` | 26 | 2026-11-12T11:08:37Z | `e842853f0da95121ef13e22d8c6326ef36c2be9d2a35d28489c77fb22e033b5d` |
| `allwise_mep/ngc4941.csv` | 27 | 2026-11-12T12:28:13Z | `1bb51bb637beaec556d0aee96f591720213b72b36ecd323cc9e7302335c10578` |
| `allwise_mep/ngc5347.csv` | 34 | 2026-11-12T12:29:34Z | `d7ab158568feee2bbb023e71eb0bfdac4e16dd1ade5211e6639ede950c0a2c51` |
| `allwise_mep/ngc5728.csv` | 22 | 2026-11-12T12:31:37Z | `b0a70b4dc3187dacd9654994d60e9b788b6563ba97c51ba8dbda6dc13f9118f5` |
| `allwise_mep/ref01.csv` | 23 | 2026-11-12T11:18:21Z | `6c18b8eb5930641b4c0c6037d7941e56288cf0835b6a91a2d0d62813a5b850e6` |
| `allwise_mep/ref02.csv` | 28 | 2026-11-12T11:20:31Z | `795dd2ca9272de030bb497971a5c43725f92eebb8a2737fbb8345aa7243b71f6` |
| `allwise_mep/ref03.csv` | 25 | 2026-11-12T11:22:22Z | `ceae5f89835c665113140f52db977d3880af4aeeee464ffbcc929576a2cbfaf6` |
| `allwise_mep/ref04.csv` | 29 | 2026-11-12T11:23:22Z | `cd5292fcc76cb1477fbefc4f746f2499b87e5ed94ab933b95f062efe863acdd4` |
| `allwise_mep/ref05.csv` | 27 | 2026-11-12T11:25:07Z | `8b7fc608101710c8b130e438fe4f85bd7c4f166401b786deea1032acb2b09643` |
| `allwise_mep/ref06.csv` | 24 | 2026-11-12T11:27:22Z | `06df8b200e83e61b026df35e0aaa09500f8b237b93b20df02a4fca4f6774eb75` |
| `allwise_mep/ref07.csv` | 58 | 2026-11-12T11:28:57Z | `9cd45b96377e3a4d8249fd728c77faaee49d49d2d3dfa554486fbe3041bfa64f` |
| `allwise_mep/ref08.csv` | 29 | 2026-11-12T11:31:01Z | `bdf1978123f4b57637b44feb123a9f4e82af3ae0415dd50724050607b3ce8273` |
| `allwise_mep/ref09.csv` | 27 | 2026-11-12T11:33:17Z | `0ebc2b89496c2c53f61419e8a1f63b121f034d947203f5cd171e743d35cdae61` |
| `allwise_mep/ref10.csv` | 25 | 2026-11-12T11:34:55Z | `f2dfd54293750daf8e02c86c4279b78a3dd7af297cb8213eef2823695a32de94` |
| `allwise_mep/ref11.csv` | 25 | 2026-11-12T11:36:59Z | `18e55de0d3d80f8114680660839e5511bb13697dd8c4b84f90c11a7c885c42b1` |
| `allwise_mep/ref12.csv` | 30 | 2026-11-12T11:38:43Z | `1773fea30b853ad0d82e8f0b48cab80bccf438cb0841955bbe2bb92db055d530` |
| `allwise_mep/ref13.csv` | 58 | 2026-11-12T11:40:39Z | `5186485d93e27430ae2234dbf7ff0c68e403da26419ae9ee7f08146220800c39` |
| `allwise_mep/ref14.csv` | 29 | 2026-11-12T11:42:49Z | `e2203c7bfe3448cc9df4ed6d875940e1fbe0abb56077740bea1b1d85a5c2a809` |
| `allwise_mep/ref15.csv` | 27 | 2026-11-12T11:44:33Z | `179e2d922f3b628d9adae212e712a2edba7c5e29acf50ac980b700a65dad155b` |
| `allwise_mep/ref16.csv` | 25 | 2026-11-12T11:47:40Z | `4848a0ed2390bd3d22387f7e0745d2ac8a2a7a067985581e18e994d05f3a96e0` |
| `allwise_mep/ref17.csv` | 26 | 2026-11-12T11:49:32Z | `2df0f5bea09332ec59f3cd3192454f2e16241d707d1aed2a7913ca801a45d364` |
| `allwise_mep/ref18.csv` | 26 | 2026-11-12T11:50:44Z | `b55a101ca91e85a3b27189778a56775d0602a0abf3eebd03ad6b6cbbf5b62ef4` |
| `allwise_mep/ref19.csv` | 27 | 2026-11-12T11:52:42Z | `9feeec976ab35d4bbed6b222e4fd2bbe77fd048d6eaabea5891c177e202a0c04` |
| `allwise_mep/ref20.csv` | 27 | 2026-11-12T11:54:39Z | `6e2288bd0402886152b60ccb595f6eb9506b21a18e1a98c8134bfca8f7a364a3` |
| `allwise_mep/ref21.csv` | 25 | 2026-11-12T11:56:29Z | `65c942c3e1cf100ca9df3e083aca760e481c9c03c3be334f988c90748975cc04` |
| `allwise_mep/ref22.csv` | 27 | 2026-11-12T11:58:26Z | `08b8039e4f7c368fcc75972d8d4c066fff59b9f588e8f71174e6a8c87f70fe7d` |
| `allwise_mep/ref23.csv` | 26 | 2026-11-12T12:01:54Z | `1084ccdcfd97290c3cf3a834a38bba6eb780090a8111a3ff07468b25df725942` |
| `allwise_mep/ref24.csv` | 25 | 2026-11-12T12:03:33Z | `bdcf9c1c906cf04529b2bb8badc65d0cb5d157c13c4b31e1a3dcc3fc82e2d1c5` |
| `allwise_mep/ref25.csv` | 48 | 2026-11-12T12:05:19Z | `d49eb1daae7605cef1a51cbff85a2017e9a4a839da0b9549022b602f446e4361` |
| `allwise_mep/ref26.csv` | 24 | 2026-11-12T12:07:01Z | `4ba1b5e3b13de750ea2ccd93bf44b42ebff68809f85e184777c1829872dcb982` |
| `allwise_mep/ref27.csv` | 28 | 2026-11-12T12:09:01Z | `d962518d0c758e61bb3ac61fc8c9730ffa9abcc61cb670e207002c13d2a7bc74` |
| `allwise_mep/ref28.csv` | 50 | 2026-11-12T12:12:44Z | `6d69ab9753db4b6e0bbedb1ae48f9b945f407e95de08aec2b7342c73ec40bad8` |
| `allwise_mep/ref29.csv` | 27 | 2026-11-12T12:14:41Z | `78dc22ad42809839ec63616a3a06e74d74e05e0d0bcf1a0df6529d9c13f3a6fb` |
| `allwise_mep/ref30.csv` | 28 | 2026-11-12T12:16:51Z | `5fb78ba8df2d4f57462ce413135f5d14d4b4c0c952d3983c96faf473cdc75c10` |
| `neowise/mrk573.csv` | 293 | 2026-11-12T12:18:43Z | `d93a35e3135a26cf54531560e371bb5f2a4b27bfe0cef90bfe2777dcecc48599` |
| `neowise/mrk573_s01.csv` | 300 | 2026-11-12T13:53:45Z | `8f89b45aa246bd6d01d75ddaa3f685dd27e9ee85106d2f8ce27232c45e3748f2` |
| `neowise/mrk573_s02.csv` | 301 | 2026-11-12T13:55:50Z | `dd5f273efe952de96f270349d6e0403c30e22e60f272e6580dc2d74be235537d` |
| `neowise/mrk573_s03.csv` | 255 | 2026-11-12T13:57:20Z | `04081edb06dc74ad72455bdf2f248c5a1d5e286ccfcb011668975c31fe018ae6` |
| `neowise/mrk573_s04.csv` | 255 | 2026-11-12T13:59:14Z | `609473a12c45cc65224ddd5bdb557f67e944bae3fe7231e69cad6945fa5c048f` |
| `neowise/mrk573_s05.csv` | 255 | 2026-11-12T14:01:00Z | `2474b4bdf6f81678dabd2acc178267ddeaac97fdca5898c09eaf800e40053d4d` |
| `neowise/mrk573_s06.csv` | 313 | 2026-11-12T14:03:05Z | `c1a44b422cb147a2f9c769cd95d49b7d5e54bc476dc682619fa74f06ed37e1cf` |
| `neowise/mrk573_s07.csv` | 309 | 2026-11-12T14:04:43Z | `f1b666855bd54042f8887aeed7f6b86b1c03d44fc6964ec3940cf931664ff77a` |
| `neowise/mrk573_s08.csv` | 301 | 2026-11-12T14:06:28Z | `0a134596348b66479ca7aa7768fb82355d3dd9feb3a67dc6d92080940eb7fd1a` |
| `neowise/mrk573_s09.csv` | 307 | 2026-11-12T14:08:12Z | `08b67b161479c2a339a8b72dd1a37e981ca79572076d63bc4916e14b684875e3` |
| `neowise/mrk573_s10.csv` | 275 | 2026-11-12T14:09:51Z | `584c2d4d9b401e8a44fbb579ba47024947baea7171338b06ce6e00b28905196b` |
| `neowise/ngc1386.csv` | 394 | 2026-11-12T12:20:38Z | `df676a15e315f6ee26720ebea08939b154e87f2e414743ba3e7817814f078e7f` |
| `neowise/ngc1386_s01.csv` | 441 | 2026-11-12T14:11:58Z | `4b499c75ba1e48fb335604f32ed7740e14a21815ecbbcfa4d30b87d65f55949a` |
| `neowise/ngc1386_s02.csv` | 418 | 2026-11-12T14:15:37Z | `032cf58cfdcec390edf7f033d09c2d01d9281d65f5f9a8c31423cb566ae8c1e9` |
| `neowise/ngc1386_s03.csv` | 429 | 2026-11-12T14:18:02Z | `6110d2fb17282501f0288dc45d6f4b1659af5d096402cdd6bea454e8d76282bd` |
| `neowise/ngc1386_s04.csv` | 372 | 2026-11-12T14:18:51Z | `ffd80946a3cdb8d0b20f2fd16887a85fbffdfdca9ecc7913593e2f8e3a8ed04a` |
| `neowise/ngc1386_s05.csv` | 357 | 2026-11-12T14:19:37Z | `d17f3877301f0a3cb6b92b1ccd530d82de057e36516955388a454c1b4b19d1bc` |
| `neowise/ngc1386_s06.csv` | 435 | 2026-11-12T14:21:27Z | `212e02ab8f28a2da469ff9b384d8e4d29fbd1e477b09401b3d33d7eb76c13b1d` |
| `neowise/ngc1386_s07.csv` | 395 | 2026-11-12T14:22:35Z | `55ba90b108823420da606d33231f0b92a4e0b3118d96cbe734d163eea7a5ad7d` |
| `neowise/ngc1386_s08.csv` | 401 | 2026-11-12T14:24:12Z | `07539d5118de4c9d1b48f4b91ddfadcdc15a77ec7728cbdae4c81a8b1e418f80` |
| `neowise/ngc1386_s09.csv` | 409 | 2026-11-12T14:26:38Z | `2f2b7c2e2cd0b991e347646f0efde2c57427bb54aa90dfccb1827274c79b263b` |
| `neowise/ngc1386_s10.csv` | 382 | 2026-11-12T14:28:30Z | `ad322db5aad47f818da15e46b81f064f6bbfa223ba5b060e7d0159bd04a313a0` |
| `neowise/ngc2273.csv` | 312 | 2026-11-12T12:22:50Z | `992fc9c54c44147a2549321d175ab3b20971e54db4f06f22e91f35f14f2a8b14` |
| `neowise/ngc2273_s01.csv` | 331 | 2026-11-12T14:29:48Z | `849a81ad281a48428335d07ba936ed44293f43511115258063ca6e6e58572af6` |
| `neowise/ngc2273_s02.csv` | 339 | 2026-11-12T14:31:19Z | `0582c6a9cb03f9ea3c3285a922811f488f3de89366bf986ddabf3b67b8c01a5d` |
| `neowise/ngc2273_s03.csv` | 321 | 2026-11-12T14:32:45Z | `a5492572c984d4da96691fc7b4e712c37fd592262a37d14fcea62a9eeb0e4137` |
| `neowise/ngc2273_s04.csv` | 315 | 2026-11-12T14:34:32Z | `0d23a8d618adbeb026f73d73ccea5acf536bb33d3617fd6b95e0e284533c2804` |
| `neowise/ngc2273_s05.csv` | 318 | 2026-11-12T14:36:00Z | `0268d3f55fcf77f93934d61ebce33c3d8e655024bc6bf117abc31ea1be74475a` |
| `neowise/ngc2273_s06.csv` | 326 | 2026-11-12T14:37:38Z | `d5ba57e53cf69669fe4ad0b0611cada4a9136fb0798a79bd302627aae711d8f0` |
| `neowise/ngc2273_s07.csv` | 314 | 2026-11-12T14:39:18Z | `70f53e9014636734e139cb2032c6807fc8d29112f362df018e8115667da00663` |
| `neowise/ngc2273_s08.csv` | 326 | 2026-11-12T14:40:44Z | `b0e7a5d715e3be95ee62885ea3409806c85e83e17178769df8bef435e0bf7905` |
| `neowise/ngc2273_s09.csv` | 316 | 2026-11-12T14:42:10Z | `d7cbd041a6b4e61cdf0b45f0984cbcbc89d6091d1abb32957390d7d9dcc18406` |
| `neowise/ngc2273_s10.csv` | 327 | 2026-11-12T14:44:34Z | `79d90c97a33968efa16e2fd819dfab685af87bc3bd1abf05c10ed475330ecfc4` |
| `neowise/ngc3081.csv` | 333 | 2026-11-12T12:25:03Z | `3f91f49feb7ee14e7b596abe61bba6c10aa313b154ec3ed293a118e2ff946b72` |
| `neowise/ngc3081_s01.csv` | 308 | 2026-11-12T14:46:48Z | `97a2c6bf897a4c6c291d5af65f14ecbf7f96e2460d4599d5679de9535680b441` |
| `neowise/ngc3081_s02.csv` | 352 | 2026-11-12T14:48:15Z | `e4e0954632d622a9597cedaa8a7e6270d1f98669f393e70c2e39a45b5de8e12f` |
| `neowise/ngc3081_s03.csv` | 321 | 2026-11-12T14:50:12Z | `5499e18ddec42b6fd0b5570a2a3a8130f597f7ac9efe68b14dd9f9371a946add` |
| `neowise/ngc3081_s04.csv` | 337 | 2026-11-12T14:51:44Z | `1c7c8c2526a700ad898b7f7bd75e120752f7e75b27f6c9e9c2291198ddcdbf56` |
| `neowise/ngc3081_s05.csv` | 338 | 2026-11-12T14:54:28Z | `20648f0562649d9b7907afe00562b2e6ff7bc6464ebb7c2e556ad32802a1fedb` |
| `neowise/ngc3081_s06.csv` | 313 | 2026-11-12T14:56:08Z | `22a530ded869f6ea075fc4c1900f2d52752caba9d89adcd79ee5e82ba17f0727` |
| `neowise/ngc3081_s07.csv` | 323 | 2026-11-12T14:57:57Z | `5aca59b7ae599da1895fa5953e76085b3f657cc9f8781e548dc6952d7cca89cf` |
| `neowise/ngc3081_s08.csv` | 312 | 2026-11-12T14:59:47Z | `59122e43e381ad8475fbedf3c07e88f2b4cd0fd2cc705c44efb26235c3086f79` |
| `neowise/ngc3081_s09.csv` | 338 | 2026-11-12T15:02:15Z | `30d2035a81ed72b23b19fa5abd8c47a841d6a201a757a4da21a8a6e702ca7722` |
| `neowise/ngc3081_s10.csv` | 336 | 2026-11-12T15:04:42Z | `16a337b3b9aa52e22cfba0aa37fadee7ab7c7da437a6961f34e4b3aa8ff901eb` |
| `neowise/ngc3393.csv` | 340 | 2026-11-12T12:26:52Z | `827e0e4dcc8a5ebf8d4650bcedfd6c06a84de897103949ff26e2484cc2b00fd4` |
| `neowise/ngc3393_s01.csv` | 338 | 2026-11-12T15:06:34Z | `89e7c6c6bd6e9392f0fb53ea48af3e9c16e79ce70f1a1338f68467ea847ed651` |
| `neowise/ngc3393_s02.csv` | 320 | 2026-11-12T15:08:06Z | `820ae2ec63d51fff40ae0ee0a81bb4fbb2bee0a33c385a4bf93af5975af40ad3` |
| `neowise/ngc3393_s03.csv` | 322 | 2026-11-12T15:09:52Z | `be763590ba0f40b5dd8753c3e77a1d39de8c5dade6adf7aa802a0fc07393814b` |
| `neowise/ngc3393_s04.csv` | 322 | 2026-11-12T15:12:09Z | `bbdae645ddfe36676e996906281a06ba5b0f900b5c2e9ded78d98919038934d1` |
| `neowise/ngc3393_s05.csv` | 312 | 2026-11-12T15:14:21Z | `727a33c593836f7c9fa90b6a2d7808217eb212d9bac56661bad8874e58b3ca73` |
| `neowise/ngc3393_s06.csv` | 333 | 2026-11-12T15:16:19Z | `2f01ea6a93dc87ad0e6fb2e40e431fefc257bccf866a02f352d7098cb1efd6d9` |
| `neowise/ngc3393_s07.csv` | 330 | 2026-11-12T15:17:43Z | `62a62c0d79df50d686677295564d4e3535d1c2d1b66d18757ed91a7c33e1e66d` |
| `neowise/ngc3393_s08.csv` | 326 | 2026-11-12T15:19:18Z | `25ccc5921b57998dc2eb92d7cf6a6730ca7cdcedefd50d66345d67d757292cdd` |
| `neowise/ngc3393_s09.csv` | 328 | 2026-11-12T15:20:58Z | `908691a5b58a7d7786317ef604b6dc403f4f3ce192549ae25254f5ae73e71d9e` |
| `neowise/ngc3393_s10.csv` | 314 | 2026-11-12T15:23:00Z | `7bdaede28311ee91248ef52771b1b7b4251c5597ac57738ea9ede68f5e680176` |
| `neowise/ngc3786.csv` | 318 | 2026-11-12T11:08:22Z | `c64a843e7605eea4efcd78074d8f4f850957b54e2d5ffbd39632a0515fef3c1f` |
| `neowise/ngc4941.csv` | 291 | 2026-11-12T12:27:59Z | `9a47edcf3d6415eda35397c34ee95a61905d64bdebc2dfaa4c44b2defde699d8` |
| `neowise/ngc4941_s01.csv` | 290 | 2026-11-12T15:24:54Z | `b19f29975f343e854fb9fa21b79d9371ae4ce1a4f66b19252460d5fc60ed6f46` |
| `neowise/ngc4941_s02.csv` | 292 | 2026-11-12T15:25:51Z | `ade9587c3ea9441f4cef2d539dfd31b96082afa29d13f35076b725ef2b5e002e` |
| `neowise/ngc4941_s03.csv` | 283 | 2026-11-12T15:27:36Z | `d6c932d09ad2fec5f783d28e0bc1bb3e4490726a0f61b5374a4e28ec1209fc71` |
| `neowise/ngc4941_s04.csv` | 287 | 2026-11-12T15:29:01Z | `5c1b9aa8f6af9fd5ee81ad8ceae95d9ab556ab94203b9fc86e10737105abad2e` |
| `neowise/ngc4941_s05.csv` | 287 | 2026-11-12T15:31:09Z | `995f75466e77aa0e7699e38538d3ae648ba774da1435e32808cc9d224450a5ac` |
| `neowise/ngc4941_s06.csv` | 302 | 2026-11-12T15:33:06Z | `6e1833b6125cd678bb309712e2da6fb0a3a3888fdc2ed66a54215d137a577880` |
| `neowise/ngc4941_s07.csv` | 297 | 2026-11-12T15:34:32Z | `d96cdf1d89a4d81329d95d27e3442a771fd0a0feffa8c06d853bc5d5eeb57710` |
| `neowise/ngc4941_s08.csv` | 300 | 2026-11-12T15:36:48Z | `815a455725e71e280a6ab1ac72df03a17728112174234977df754ded03944ed6` |
| `neowise/ngc4941_s09.csv` | 289 | 2026-11-12T15:39:17Z | `ee79248b82ca0ee3ddddd50779b0ce1e4f1c235c8292b322b95fdf0ba2f26857` |
| `neowise/ngc4941_s10.csv` | 292 | 2026-11-12T15:41:32Z | `72dff22e8df4b48e288d629a807d5b9be880a688b2bcff164cc19e7537678b7b` |
| `neowise/ngc5347.csv` | 383 | 2026-11-12T12:29:12Z | `4b00eeec28a9471c95e72600bd79f8bc59c78d3a7029d42b9c010b70cd5ddb01` |
| `neowise/ngc5347_s01.csv` | 378 | 2026-11-12T15:43:37Z | `ce7ec0fde78a977df471eb7e380bda8fe5f4da9c984025f34b9e3067d22cebac` |
| `neowise/ngc5347_s02.csv` | 380 | 2026-11-12T15:45:39Z | `20f8906ee53c8f05d0e08eba007c37b579c9d555688ed52b7cd51bd10973af57` |
| `neowise/ngc5347_s03.csv` | 381 | 2026-11-12T15:47:56Z | `45a9c191937bd5a759edd3fb26ee44553e0dcbbb02fedf2e7c08ca498d29143c` |
| `neowise/ngc5347_s04.csv` | 384 | 2026-11-12T15:52:17Z | `beff13e58e8468cd64c703035092dd48e7183930b6268c4877990e41b4a016f6` |
| `neowise/ngc5347_s05.csv` | 388 | 2026-11-12T15:55:34Z | `d52f98fe88c874f7ecb2c8934b332cf07d23ab6730029397817cae697fa08ad2` |
| `neowise/ngc5347_s06.csv` | 361 | 2026-11-12T15:59:49Z | `05bcdffbff3a5e43e888c97b4dbce585d4e9685ab96625cc295a8feb1174b733` |
| `neowise/ngc5347_s07.csv` | 393 | 2026-11-12T16:02:43Z | `af38dbbf9ca55305cbfff15b1294a338a4a6af7590e54f3385fc042d82395608` |
| `neowise/ngc5347_s08.csv` | 384 | 2026-11-12T16:05:03Z | `973a391dcc9faf06a1a1ac0abe3764cde560ee2cdee7861d8be1315145fd2e0b` |
| `neowise/ngc5347_s09.csv` | 378 | 2026-11-12T16:07:14Z | `33416933cba87395caedbf6f1e78b83c41657fd3698098306dc508d277ba7e9a` |
| `neowise/ngc5347_s10.csv` | 393 | 2026-11-12T16:09:08Z | `681bf23ef19710945bcef09f6723bf895577b5ceab3dd6dfe5fb09ec866ce8ec` |
| `neowise/ngc5728.csv` | 291 | 2026-11-12T12:31:22Z | `92fd70c3fcd0116ad80ce22e1ae9975cde5db11c27020065c8b9d9fb72aac931` |
| `neowise/ngc5728_s01.csv` | 277 | 2026-11-12T16:11:18Z | `f629af0bc218a6c0ad5b25de05f026544a5b4dfd302a85ed965e557a3aecbad9` |
| `neowise/ngc5728_s02.csv` | 277 | 2026-11-12T16:14:10Z | `3a06b60899e503829dd82b5dc21556f59b2d2ca7b5d4395d015ad6044bf3868c` |
| `neowise/ngc5728_s03.csv` | 277 | 2026-11-12T16:16:10Z | `fb413563a75a9c7f7a37a64c3727818f018526ab167b5bcdecf35d278a0a3a1c` |
| `neowise/ngc5728_s04.csv` | 271 | 2026-11-12T16:18:20Z | `450a55072860c2682a31b669fd6dc991dc0a4f7e6d4309240c2757e1988a7550` |
| `neowise/ngc5728_s05.csv` | 276 | 2026-11-12T16:21:11Z | `2331bdf518112093c0231e36e50f34d029597862f94c0b64c3994f1ffb1a134c` |
| `neowise/ngc5728_s06.csv` | 280 | 2026-11-12T16:23:21Z | `87bb986c460dc48651a0470d2683779b18adf9d3adf552f01d7314fbb69a187e` |
| `neowise/ngc5728_s07.csv` | 288 | 2026-11-12T16:26:03Z | `2b25d5fcd49158351c5cba839d0678e26218db9ea69b0613cbdbbb6711aa1ab4` |
| `neowise/ngc5728_s08.csv` | 276 | 2026-11-12T16:28:14Z | `2e2a1e76293e3ea8095e3cbee23d30eac86fdb2799b361a16d8c5bf553210709` |
| `neowise/ngc5728_s09.csv` | 286 | 2026-11-12T16:29:54Z | `9b440f88b11f443265b3cc6eda2d4bcf4506653ef5d7da3f9839b68f90aad96c` |
| `neowise/ngc5728_s10.csv` | 279 | 2026-11-12T16:31:07Z | `11399ff71b28e384a98bf0f74f91ec1d94b92d994876fefc68b970c70b655892` |
| `neowise/ref01.csv` | 308 | 2026-11-12T11:18:03Z | `79387b5b4615fc30ea0099096e71df505e71aa8f9159ac52a55d55ef28d8a380` |
| `neowise/ref02.csv` | 298 | 2026-11-12T11:20:04Z | `55bc91ceb9068e51b80f8c75b2ff0e53a90fbb8c5809f993c13d73d0add84661` |
| `neowise/ref03.csv` | 301 | 2026-11-12T11:22:08Z | `07dff702039020d752bf0aa232ebfed5622316610740768b14576766e5735e7e` |
| `neowise/ref04.csv` | 310 | 2026-11-12T11:23:05Z | `012b4801b060acd98a5ad057e965fa9c4b8c414b62bd0adb97d4ae27f345e3fe` |
| `neowise/ref05.csv` | 309 | 2026-11-12T11:24:54Z | `2f7689c12466d4f83c6380b1fc3d8080f547c87031fa94affec4c5e8e077c9d7` |
| `neowise/ref06.csv` | 319 | 2026-11-12T11:27:09Z | `a5400e5a1777a2ebea487ede9a960081307243c2a3a1276d0421f3724d488d19` |
| `neowise/ref07.csv` | 328 | 2026-11-12T11:28:43Z | `22747cda3d8928579871f41e00ae1d47be88184cd73ce86f3b9513dcd62ecd53` |
| `neowise/ref08.csv` | 317 | 2026-11-12T11:30:48Z | `6a9f7ede541d07a6b06836733034a13a8dddebb4cc526d0172ebbd23cc596501` |
| `neowise/ref09.csv` | 322 | 2026-11-12T11:32:58Z | `0cba9a7914a5aadff665b3ee4efa4b6f9ed8580f2fecc22b87909449e9f6440e` |
| `neowise/ref10.csv` | 324 | 2026-11-12T11:34:43Z | `6ef183f64e2a96390c49ceb33b2361f78a2e06374c97733cad6cede0ec507fb8` |
| `neowise/ref11.csv` | 305 | 2026-11-12T11:36:44Z | `c099d72b17dae71429135ce80feae73acb820cb9630fac25ad73b2e25ced75ab` |
| `neowise/ref12.csv` | 290 | 2026-11-12T11:38:28Z | `7166bea8e903a46b3393dd2dc29dd2c024ec74c2bc0689199447afb9bf875840` |
| `neowise/ref13.csv` | 302 | 2026-11-12T11:40:19Z | `91e6186e6db30a19b4b5504367c59ec175f2e2b7361d8187ad5a0ed7c566ed26` |
| `neowise/ref14.csv` | 319 | 2026-11-12T11:42:24Z | `659ea5e6c5bb713be2cbc2257e4de98d955f4f4b4be9a3d008d1fa67f93b2fa4` |
| `neowise/ref15.csv` | 317 | 2026-11-12T11:44:12Z | `bf113a973c966a79dd80a77723fc6ad7e109e3b6e74097ed07e2ff3ef16aed4a` |
| `neowise/ref16.csv` | 308 | 2026-11-12T11:47:23Z | `255bd1c798bb5c7c5ffea4992a2f85b730373a74cf4c22bf7d2f9d911fd81d6e` |
| `neowise/ref17.csv` | 308 | 2026-11-12T11:49:15Z | `a14aa0b84fcf69cdcac5a1dfdbde41d1d6e8fe615f5548e87cf2da4eec20ea84` |
| `neowise/ref18.csv` | 314 | 2026-11-12T11:50:23Z | `c42f3d9100f46d1ceeed0f4f80349dee177c135894c29841e63c474b80f8debd` |
| `neowise/ref19.csv` | 329 | 2026-11-12T11:52:28Z | `8607db68c88795c1f03a4b74ce14220165ddbc4b42c401be879df406774cd340` |
| `neowise/ref20.csv` | 322 | 2026-11-12T11:54:24Z | `d4e2eb55c47353754686cfc5f2247437543691eeebc9aa85d9eb59610b8c38e5` |
| `neowise/ref21.csv` | 298 | 2026-11-12T11:56:15Z | `7e0202eb3efbb446e47a998f6914c954f2e5d6a2e72699c152cadedcf1b8938a` |
| `neowise/ref22.csv` | 301 | 2026-11-12T11:58:13Z | `c1996933ad143fec6f1189c58c20e039a82e48869f1f4851bc57321d905615b5` |
| `neowise/ref23.csv` | 313 | 2026-11-12T12:01:32Z | `9e4da2e16867efe3c5ec7813b58a4c77714c955add17db561395e7d218b383e7` |
| `neowise/ref24.csv` | 310 | 2026-11-12T12:03:18Z | `1dc734851b261a6d601229999dc4ce8be4f870a37419ef49aa5c713372b40862` |
| `neowise/ref25.csv` | 314 | 2026-11-12T12:05:03Z | `50e98e81e778ba6649259819343a320d238f07ab5122ee592c398074d5c65464` |
| `neowise/ref26.csv` | 324 | 2026-11-12T12:06:43Z | `0f69f0be5c06ecef9113fe71679e91d15b192232140dfac42653f66f04347b31` |
| `neowise/ref27.csv` | 325 | 2026-11-12T12:08:44Z | `1a1f9bc3d3fe2706e48c8470aedf5f4e9a33aa1ad2a589e75b42945834c481b8` |
| `neowise/ref28.csv` | 323 | 2026-11-12T12:12:32Z | `1ff20a42f12851261d3937f2e73e210ae07bb110336a757d3e6d4f2bf7c1e03c` |
| `neowise/ref29.csv` | 314 | 2026-11-12T12:14:23Z | `1929331a112a55125f3efbcbf4b998f865533154fe84739711fec488d358e993` |
| `neowise/ref30.csv` | 318 | 2026-11-12T12:16:28Z | `7c5ab7d8ba2ce5becd15695419ccc5406bcc4600c7c5b56477c4b10a1e333467` |
| `untimely/mrk573.csv` | 34 | 2026-11-12T13:50:09Z | `62a7e732e5e3f67deaa78039ec986b77f28d0753dcb78bfd14d0dfb490fe5884` |
| `untimely/ngc1386.csv` | 34 | 2026-11-12T13:50:33Z | `d9ea81d5ec29c7c05ccd6490b1c4fa785723160137362ffe4c139c3fe5166438` |
| `untimely/ngc2273.csv` | 32 | 2026-11-12T13:50:56Z | `7d4ddb48366bd2a88867a10721554343254896ad33149658c5e11214ba052fbf` |
| `untimely/ngc3081.csv` | 32 | 2026-11-12T13:51:19Z | `2b521b0ecb328b166ce8b71ec753104efcca97f6fcf0bc27367a28949131630d` |
| `untimely/ngc3393.csv` | 32 | 2026-11-12T13:51:39Z | `183dfc22e2196cf2111770a9a7cfa59f40d767a8ae7d0a0a9fe9109fc6b11d95` |
| `untimely/ngc3786.csv` | 32 | 2026-11-12T13:49:55Z | `2918e64fe120f2871636e89d9ea61cdf8e4c929fbffa5696301bf37d0b8ec7a9` |
| `untimely/ngc4941.csv` | 34 | 2026-11-12T13:51:53Z | `3e946e8c7f4d7f32fdf9f262ed7ffc41f58ebcf5f3e26a2dbb86ad16df0e348f` |
| `untimely/ngc5347.csv` | 34 | 2026-11-12T13:52:14Z | `27798b713a69642e7349c492640686725cf2182026a67f29523039bca97d9279` |
| `untimely/ngc5728.csv` | 32 | 2026-11-12T13:52:36Z | `5c941c0efb1c1924ed5373277afd7ab6c758f63aae61ddcd1432ac1390d37eb4` |
| `ztf/ngc2273_g.csv` | 630 | 2026-11-12T14:56:11Z | `ab43fe101872e1e8364c1d0ef759467d4c00cca41d0951af14022488529aa5e1` |
| `ztf/ngc2273_r.csv` | 633 | 2026-11-12T14:56:27Z | `395f5ebef141ba0f86615882021efe4f6503ef1d3133cb4ba1a02b058acf57fd` |
| `ztf/ngc3786_g.csv` | 747 | 2026-11-12T14:52:19Z | `aa3ce98e7166531505007d7f9c2efd8bf5a0123346b275524e1df7893c06df66` |
| `ztf/ngc3786_r.csv` | 933 | 2026-11-12T14:52:29Z | `68aaf20d02f15fa4582a21651d8315677c730b53c44dd654ef1a8abd0258abbe` |

파일 2147개, 모두 9,962,672행, 1935.06 MB.

5a·5b단계의 파일은 수가 많아 묶음별 합계만 적는다. 파일별 조회 조건, 시각, 행 수, SHA-256은 `data/manifest.json`에 있다.

| 묶음 (data/cache/ 아래) | 파일 수 | 행 수 | 크기 (MB) | 받은 시각 (UTC) |
|---|---|---|---|---|
| `allwise_cat/s5_fields_*` | 47 | 3,740,631 | 552.5 | 2026-11-14T10:58:55Z – 2026-11-14T11:32:29Z |
| `allwise_cat/s5_targets_*` | 9 | 136,891 | 19.2 | 2026-11-14T10:53:31Z – 2026-11-14T10:55:45Z |
| `allwise_cat/s6b_fields_*` | 42 | 3,476,515 | 514.0 | 2026-11-16T23:00:32Z – 2026-11-16T23:39:38Z |
| `allwise_cat/s6b_targets_*` | 7 | 107,889 | 15.1 | 2026-11-16T22:55:27Z – 2026-11-16T22:59:29Z |
| `catalogs/*` | 4 | 4,369 | 1.0 | 2026-11-14T10:53:07Z – 2026-11-16T15:57:29Z |
| `s5/allwise_mep/*` | 231 | 92,631 | 14.6 | 2026-11-14T11:33:35Z – 2026-11-14T13:46:09Z |
| `s5/neowise/*` | 231 | 857,888 | 173.1 | 2026-11-14T11:33:16Z – 2026-11-14T13:45:55Z |
| `s5/untimely/*` | 231 | 50,186 | 15.7 | 2026-11-14T11:35:58Z – 2026-11-14T12:51:32Z |
| `s5/ztf/*` | 395 | 300,572 | 76.2 | 2026-11-16T15:26:49Z – 2026-11-16T18:57:56Z |
| `s6/asassn/*` | 4 | 5,143 | 0.7 | 2026-11-16T22:44:57Z – 2026-11-16T22:45:08Z |
| `s6/bat157/*` | 4 | 627 | 0.2 | 2026-11-16T22:45:37Z – 2026-11-16T22:45:42Z |
| `s6b/allwise_mep/*` | 210 | 94,910 | 14.9 | 2026-11-16T23:41:16Z – 2026-11-17T02:16:49Z |
| `s6b/neowise/*` | 210 | 843,129 | 170.1 | 2026-11-16T23:40:54Z – 2026-11-17T02:16:36Z |
| `s6b/untimely/*` | 210 | 44,201 | 13.8 | 2026-11-16T23:40:33Z – 2026-11-17T00:50:31Z |
| `s7/bat157/*` | 3 | 0 | 0.5 | 2026-11-17T04:48:33Z – 2026-11-17T04:48:36Z |
| `s7/nustar/*` | 60 | 0 | 239.3 | 2026-11-17T04:35:07Z – 2026-11-17T04:45:27Z |
| `s7/swift/*` | 60 | 0 | 82.0 | 2026-11-17T04:45:28Z – 2026-11-17T04:48:32Z |
<!-- manifest:end -->
