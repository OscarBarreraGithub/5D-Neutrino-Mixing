# CA Worklog: ca_w9_vlq_4top
**Date**: 2026-05-17
**Family**: collider_rs
**Batch ID**: Wave-9 VLQ + 4-top
**Cycle**: 1
**Checker agent**: CA-CA-w9-vlq-4top
**Process IDs**: CR008 CR010 CR014

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| CR008 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| CR010 | FAIL | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| CR014 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| CR008 | Singlet VLQ T pair-production mass limit | m_T > 1.36 TeV at 95% CL | ATLAS2024_TsingletPair / PDG2025_TprimeVLQ | 6857ab7bf208 |
| CR008 | CMS T branching-triangle envelope | m_T > 1.48 TeV for all third-generation decays; up to 1.54 TeV depending on branching fractions | CMS2023_VLQPairLeptonic | 54758c30e017 |
| CR008 | ATLAS 8 TeV T branching-plane history | 715--950 GeV mass range | ATLAS2015_VLQPair8TeV | 9c84ca8782aa |
| CR008 | CMS 13 TeV single-lepton history | singlet T below 860 GeV excluded | CMS2017_VLQPairSingleLepton | 560a80326454 |
| CR008 | CMS 35.9 fb^-1 multichannel history | T below 1140--1300 GeV excluded | CMS2018_VLQPairLeptonic | 8faff1fbcef0 |
| CR008 | Low-energy flavor comparison | M_KK^min,p50 = 47.26 TeV at g_* = 3 | local docs/quark_scan_methodology_note.tex:587 | 334d639d4167 |
| CR010 | Weak-isospin (T,B) doublet mass limits | m_T,m_B > 1.37 TeV at 95% CL | PDG2025_TprimeListing / PDG2025_BprimeListing; ATLAS2018_TBCombination | 0b85eb01d887 |
| CR010 | Simplified T endpoint | m_T > 1.70 TeV for B(T -> Wb)=1 | PDG2025_TprimeListing | 0b85eb01d887 |
| CR010 | Simplified B endpoints | m_B > 1.57, 1.56, 1.54 TeV for Hb, Wt, Zb benchmarks | PDG2025_BprimeListing | 0b85eb01d887 |
| CR010 | CMS T branching-triangle envelope | m_T > 1.48 TeV for all third-generation decays | CMS2023_VLQLeptonic | b6b6c2050657 |
| CR010 | Early (T,B)-doublet reinterpretation | 640 GeV exclusion quoted for (T,B) doublet | GarbersonGolling2013_ExoticQuarks | c1fc211e10e7 |
| CR010 | CMS 2017 doublet-scenario T history | T below 830 GeV excluded | CMS2017_TBBoosted | 8c4d2cb370b2 |
| CR010 | Low-energy flavor comparison | M_KK^min(p50,g_s=3) = 47.26 TeV | local docs/quark_scan_methodology_note.tex:587 | 334d639d4167 |
| CR014 | CMS top-philic vector resonance exclusion | m(Z') excluded up to 850 GeV observed at 95% CL, Gamma/m = 50% | CMSB2G25005_PublicResult | a747064ac476 |
| CR014 | ATLAS inclusive four-top cross section | 22.5^{+6.6}_{-5.5} fb | PDG2025_TopFourTop / ATLAS2023_FourTopObservation | 6320d3d26384 |
| CR014 | CMS inclusive four-top cross section | 17.7^{+4.4}_{-4.0} fb in PDG summary | PDG2025_TopFourTop / CMS2023_FourTopObservation | 6320d3d26384 |
| CR014 | ATLAS top-philic observed cross-section upper-limit range | 21--119 fb at 95% CL | ATLASTopPhilic2024 | cb4204207db4 |
| CR014 | CMS resonance scan context | 500 GeV--4 TeV mass scan; 4--50% widths; 1000 GeV expected exclusion for 50% width vector | CMSB2G25005_PublicResult | a747064ac476 |
| CR014 | Low-energy flavor comparison | M_KK^min(p50,g_s^star=3) = 47.26 TeV | QuarkScanMethodology2026 | b2dfbcccbdee |

## Per-process issues
- CR008: CHK-1 fail. `CR008.tex` quotes measured collider mass-exclusion numerals that are not represented in `pdg_or_equivalent.values`: CMS 2023 "up to 1.54 TeV" branching-dependent T reach, ATLAS 8 TeV 715--950 GeV history, CMS 2017 860 GeV singlet history, and CMS 2018 1.14--1.30 TeV history. These are not theory predictions, dataset metadata, EFT translations, or projections, so the L001/B001/B021 carve-out does not apply.
- CR010: CHK-1 fail. `CR010.tex` quotes measured historical mass-exclusion numerals, 640 GeV for the Garberson--Golling (T,B)-doublet reinterpretation and 830 GeV for the CMS 2017 doublet-scenario T search, outside `pdg_or_equivalent.values`. Dataset luminosity and sqrt(s) numerals were correctly left in supporting/context blocks.
- CR010: CHK-2 fail. The `Key references` section at `CR010.tex:131`--`139` lists snapshot filename stems (`pdg2025_tprime_bprime_vlq_limits`, `atlas_2018_arxiv1808_02343`, etc.) rather than the manifest keys (`PDG2025_TprimeListing`, `PDG2025_BprimeListing`, `ATLAS2018_TBCombination`, `CMS2023_VLQLeptonic`, etc.).
- CR014: None.

## Evidence notes
- CHK-1: Applied the L001 / B001_B003 / B021_B023 carve-out without re-litigation. Theory/reference inputs, luminosities, sqrt(s), scan ranges, expected exclusions, future/projection context, and the local 47.26 TeV methodology comparison were not required in `pdg_or_equivalent`. CR014's measured cross sections, observed mass exclusion, and observed cross-section upper-limit range are in `pdg_or_equivalent.values` with strict metadata. CR008 and CR010 fail only on measured collider mass-exclusion values that remain outside `pdg_or_equivalent.values`.
- CHK-2: All `source_key` and `experiment_source_key` fields used by the sidecars resolve to entries in each process `source_manifest.yaml`. Every manifest entry has a non-empty `snapshot_path`, and `git ls-files --error-unmatch` confirmed those text snapshots are tracked. CR008 and CR014 TeX key-reference entries resolve to manifest keys. CR010 fails because its TeX `Key references` section uses snapshot stems instead of manifest keys.
- CHK-3: `find flavor_catalog/references/CR008 flavor_catalog/references/CR010 flavor_catalog/references/CR014 -type f -iname '*.pdf' -o -iname '*.PDF'` returned no files. The reference directories contain text snapshots, manifests, and checksum lists only.
- CHK-4: Before checker updates, all three YAML sidecars showed `DRAFT -> WRITER-INITIATED -> WRITER-DONE` with ISO 8601 `timestamp`/`at` fields and cycle-1 WA completion at `2026-05-17T17:08:10-04:00`.
- CHK-5: All three sidecars claim `code_coverage.status: NO`, which is consistent. Focused VLQ/top-partner/collider-recast greps returned no implementation hits; ATLAS/CMS/LHC greps returned only `tests/test_alpha_s.py:89`, an unrelated CMS/RunDec alpha_s example. CR014 four-top/top-philic greps returned no implementation hits. Cited adjacent evidence lines exist: `scanParams/scan.py:399`, `scanParams/scan.py:523`--`535`, `quarkConstraints/scan.py:359`, `quarkConstraints/scan.py:377`--`441`, `quarkConstraints/validation.py:108`--`145`, `quarkConstraints/couplings.py:96`, and `quarkConstraints/deltaf2.py:449`.
- CHK-6: `HIGH` is consistent for all three collider entries. A faithful live constraint would need an LHC recast or likelihood layer with spectrum/branching inputs, cross sections, detector acceptance, and validation against ATLAS/CMS limits; a hard-coded mass threshold would be misleading.
- CHK-7: The entries do not contradict the rc1.1 / methodology-note flavor result. CR008, CR010, and CR014 explicitly treat the 47.26 TeV low-energy RS-anarchy comparison as stronger than the current direct LHC reach and as context rather than a collider observable.
- CHK-8: No unresolved `\cite`, `\ref`, `\textbf{CHECK}`, `TODO`, `FIXME`, `PLACEHOLDER`, undefined, unresolved, or `??` markers were found. Delimiter sanity checks were balanced: CR008 `\(`/`\)` = 33/33 and `\[`/`\]` = 2/2; CR010 = 54/54 and 2/2; CR014 = 29/29 and 0/0. Dollar math count was zero and brace delta was zero for all three.
