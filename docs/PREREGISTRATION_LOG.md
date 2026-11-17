# Pre-registration log

This repository starts from a single commit, so its own history does not show the order in which things were done. The table lists the commits of the private development repository that fix that order: predictions, criteria, sample definitions, freezes, and the first time the frozen target light curves were opened. Times are KST (UTC+9), from the commit timestamps; the "opened" rows are from the run logs recorded in the stage reports. The original history (commits, push times on the remote, stage reports) can be shown on request.

| Time (KST) | Commit | What it fixed |
|---|---|---|
| 2026-11-12 20:11:49 | `4f397cba8ea44888b32186c02bdebef1f0da16a4` | Stage 1–2: target selection and flare criteria, before the known-answer run |
| 2026-11-12 22:57:30 | `11b21e4b5fce59f8e04b488cd55700b704ca4766` | Stage 3–4: detector thresholds and injection grid |
| 2026-11-12 23:50:08 | `0aefbe42d0fb345aad27f53da5eb652c5615a9a2` | `PREDICTIONS.md`: hypotheses H1–H3 for NGC 3786, alternative, decision rules |
| 2026-11-14 19:51:20 | `d1ef0dd78a392aad575636649ba7ed3d3691fcc4` | Stage 5a: sample criteria S1–S11, before the sample was built |
| 2026-11-14 19:58:19 | `093fb9af4e390e4f40bea094a2f1324ed01e940c` | `sample_v1` targets (231 objects, z ≤ 0.040) |
| 2026-11-14 20:33:01 | `9fbe19d9d8ebb9f48f4390b840a46f7222c9f13f` | `sample_v1` comparison stars (2,310) |
| 2026-11-14 23:08:30 | `49187d8cad1af6a6318b04838034588f0161b708` | `sample_v1` light curves frozen; SHA-256 in `data/freeze_v1.json` (frozen 2026-11-14 23:03 KST with code at `d4dfb2b`) |
| 2026-11-17 00:25:19 | `8e2768211018a90bec81ea780db685d40659bd12` | `PREDICTIONS_5b.md`: my predictions (Q1–Q3, dust breathing, hysteresis), decision rules, Eddington-ratio threshold |
| 2026-11-17 01:08:25 | `de8d58f001105d653ca841b51b2f3362da85a58e` | Stage 5b estimators validated on toy models; measurability criteria fixed |
| 2026-11-17 01:09:43 | `60eefb5cbe82524f6b4bf52560a43a1e9a0074eb` | Theory note on the slope definition, written before the data were opened |
| 2026-11-17 01:09:49 | (run log) | **Frozen `sample_v1` files first opened**: first run of `scripts/run_stage5b.py`; 231 of 231 hashes matched |
| 2026-11-17 01:49:58 | `db1ad4f98923a628415c4bd15e7a464605ccccc1` | Stage 5b results |
| 2026-11-17 07:25:37 | `8793b03d64fcad26b1eee0085082105055c1144a` | `FLARES_PLAN.md`: what each interpretation of the flare candidates predicts, consistent/inconsistent criteria, before the physical quantities were computed |
| 2026-11-17 07:49:17 | `daa8acce5e570d5d55fdc6ae7efcb95a3e7f872f` | Stage 6A results |
| 2026-11-17 07:54:41 | `4810e08346dccc0c3b0de805eeae7a7ba75db595` | `PREDICTIONS_6.md`: two models, my prediction, decision rule and the definition of the independent sample, before that sample was built |
| 2026-11-17 07:59:41 | `9721046a3d5b979a07375c54f3642114b87b76db` | `sample_v2` targets (210 objects, 0.040 < z ≤ 0.30), before any light curve was downloaded |
| 2026-11-17 08:40:17 | `ec289b1ea0ec48d00f4bd2a7cf9c39a4247fc173` | `sample_v2` comparison stars (2,100) |
| 2026-11-17 11:18:22 | `47663501132161f0d35f08558e53a04f58b72638` | `sample_v2` light curves frozen; SHA-256 in `data/freeze_v2.json` |
| 2026-11-17 11:32:21 | `212cb44a6973f18be48378bd365d056873cf6297` | Stage 6B thresholds and power fixed with the real cadence and errors of `sample_v2`; pushed to the remote before D was computed |
| 2026-11-17 11:32:37 | (run log) | **Frozen `sample_v2` files first opened** to compute D |
| 2026-11-17 11:47:35 | `6d4d16f448ce49627cec9f819cfabd6d93f9b5d3` | Stage 6B results |
| 2026-11-17 13:33:22 | `e34ce17d3b8394104945ed432ff45a5b00fd3c98` | `XRAY_PLAN.md`: predictions and decision rules for CGCG 312-012, before any X-ray data were downloaded |

Notes

- Deviations from the plans were recorded in the stage reports, not hidden. Two that matter for reading the table: in stage 6A the criteria for light-curve shape and colour were written after the shape and colour had already been seen in stage 5b (marked in `FLARES_PLAN.md`); in stage 6B the model-1 threshold was fixed after the download but before D was computed.
- Before each freeze, the only quantities computed for the targets were frame counts, visit counts, rejection counts per cut, saturation fractions, visit gaps, error sizes and the median and scatter of the difference between the two photometries.
- Commit timestamps are set by the committing machine. The push times recorded by the remote are the independent record of the order; they are part of the original history that can be shown on request.
