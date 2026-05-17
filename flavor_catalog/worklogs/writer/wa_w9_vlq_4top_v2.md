# WA Worklog: wa_w9_vlq_4top_v2
**Date**: 2026-05-17
**Cycle**: 2
**Writer agent**: WA-v2-w9-vlq-4top
**Scope**: CR008 (CHK-1), CR010 (CHK-1 + CHK-2)
**Reference checker worklog**: `flavor_catalog/worklogs/checker/ca_w9_vlq_4top.md`

## Summary
- Applied the cycle-2 rework requested by CA cycle-1 for CR008 and CR010 only.
- Followed the L001 / T003 precedent for historical mass-exclusion numerals: removed superseded exclusion numbers from TeX prose rather than promoting them into `pdg_or_equivalent.values`.
- Followed the B023 precedent for CR010 TeX key references: replaced snapshot filename stems with canonical manifest keys.
- Left CR014 and all other catalog entries untouched.

## CR008
- Edited the CMS full-Run-2 broad-envelope paragraph to remove the quoted branching-dependent `1.54 TeV` reach while keeping the citation flow to `cms_2023_arxiv2209_07327.txt` and the canonical `1.48 TeV` envelope already represented in `CR008.yaml`.
- Edited the post-2010 history paragraph to remove historical ATLAS/CMS exclusion numerals `715--950 GeV`, `860 GeV`, and `1.14--1.30 TeV`.
- Preserved qualitative context for Run-1 and early Run-2 searches and kept the current ATLAS 2024 `1.36 TeV` singlet benchmark plus CMS `1.48 TeV` envelope.
- Appended a cycle-2 `WRITER-DONE` transition to `CR008.yaml` and updated `last_updated_at`; left `checker_passed_at` unchanged.

## CR010
- Edited the post-2010 history paragraph to remove the historical `640 GeV` Garberson--Golling reinterpretation numeral and the CMS 2017 `830 GeV` doublet-scenario numeral.
- Kept qualitative wording about the branching-fraction-plane reinterpretation and early CMS boosted-object history.
- Renamed TeX `Key references` tokens from snapshot stems to manifest keys:
  - `pdg2025_tprime_bprime_vlq_limits` -> `PDG2025_TprimeListing`, `PDG2025_BprimeListing`
  - `aguilar_saavedra_2009_arxiv0907_3155` -> `AguilarSaavedra2009_TopPartners`
  - `garberson_golling_2013_arxiv1301_4454` -> `GarbersonGolling2013_ExoticQuarks`
  - `cms_2017_arxiv1706_03408` -> `CMS2017_TBBoosted`
  - `atlas_2018_arxiv1808_02343` -> `ATLAS2018_TBCombination`
  - `cms_2023_arxiv2209_07327` -> `CMS2023_VLQLeptonic`
  - `atlas_2024_arxiv2401_17165` -> `ATLAS2024_VLQLeptonJets`
  - `cms_review_2025_arxiv2405_17605` -> `CMSReview2025_VLQ`
- Appended a cycle-2 `WRITER-DONE` transition to `CR010.yaml` and updated `last_updated_at`; left `checker_passed_at` unchanged.

## Open Issues
None expected.
