# CA Worklog: ca_w9_vlq_4top_v2
**Date**: 2026-05-17
**Batch**: vlq_4top
**Family**: collider_rs
**Cycle**: 2
**Checker agent**: CA-v2-w9-vlq-4top
**Process IDs**: CR008 CR010
**Reference writer worklog**: `flavor_catalog/worklogs/writer/wa_w9_vlq_4top_v2.md`
**Reference cycle-1 checker worklog**: `flavor_catalog/worklogs/checker/ca_w9_vlq_4top.md`

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| CR008 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| CR010 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| CR008 | Singlet VLQ T pair-production mass limit | m_T > 1.36 TeV at 95% CL | ATLAS2024_TsingletPair / PDG2025_TprimeVLQ | 6857ab7bf208 |
| CR008 | CMS T branching-triangle envelope | m_T > 1.48 TeV for all third-generation bW/tZ/tH decays | CMS2023_VLQPairLeptonic | 54758c30e017 |
| CR008 | Low-energy flavor comparison | M_KK^min,p50 = 47.26 TeV at g_* = 3 | local docs/quark_scan_methodology_note.tex:587 | 334d639d4167 |
| CR010 | Weak-isospin (T,B) doublet mass limits | m_T,m_B > 1.37 TeV at 95% CL | PDG2025_TprimeListing / PDG2025_BprimeListing; ATLAS2018_TBCombination | 0b85eb01d887 |
| CR010 | Simplified T endpoint | m_T > 1.70 TeV for B(T -> Wb)=1 | PDG2025_TprimeListing | 0b85eb01d887 |
| CR010 | Simplified B endpoints | m_B > 1.57, 1.56, 1.54 TeV for Hb, Wt, Zb benchmarks | PDG2025_BprimeListing | 0b85eb01d887 |
| CR010 | CMS T branching-triangle envelope | m_T > 1.48 TeV for all third-generation decays | CMS2023_VLQLeptonic | b6b6c2050657 |
| CR010 | Low-energy flavor comparison | M_KK^min(p50,g_s=3) = 47.26 TeV | local docs/quark_scan_methodology_note.tex:587 | 334d639d4167 |

## Per-process issues
- CR008: None.
- CR010: None.

## Evidence notes
- CHK-1: PASS. CR008 TeX now retains only the current measured collider mass limits represented in `pdg_or_equivalent.values`: ATLAS/PDG `1.36 TeV` and CMS `1.48 TeV`. The cycle-1 historical mass-exclusion numerals were removed from prose per the WA-v2 fix and T003/L001-style carve-out handling, so their removal is not a new issue. CR010 TeX measured mass limits all resolve to `pdg_or_equivalent.values`: the `1.37 TeV` doublet limits, `1.70 TeV` T endpoint, `1.57/1.56/1.54 TeV` B endpoints, and CMS `1.48 TeV` envelope. Dataset luminosities, sqrt(s), charges, branching fractions, local methodology comparisons, and code line numbers are context or model/dataset metadata rather than standalone measured observables under the L001/T003/B021 guidance.
- CHK-2: PASS. CR008 `Key references` entries resolve to the CR008 manifest keys. CR010 `Key references` now use canonical manifest keys per B023 precedent: `PDG2025_TprimeListing`, `PDG2025_BprimeListing`, `AguilarSaavedra2009_TopPartners`, `GarbersonGolling2013_ExoticQuarks`, `CMS2017_TBBoosted`, `ATLAS2018_TBCombination`, `CMS2023_VLQLeptonic`, `ATLAS2024_VLQLeptonJets`, and `CMSReview2025_VLQ`. Sidecar `source_key` and `key` entries match manifest keys, manifest `snapshot_path` fields are non-empty, and `git ls-files flavor_catalog/references/CR008 flavor_catalog/references/CR010` confirms the cited text snapshots are tracked.
- CHK-3: PASS. `find flavor_catalog/references/CR008 flavor_catalog/references/CR010 -type f \( -iname '*.pdf' -o -iname '*.PDF' \) -print` returned no files; the reference directories contain text snapshots, manifests, and checksum lists.
- CHK-4: PASS. Before this checker update, both sidecars had ISO 8601 histories ending in WA-v2 `WRITER-DONE` at `2026-05-17T17:56:57-04:00`, following `DRAFT -> WRITER-INITIATED -> WRITER-DONE -> WRITER-REWORK -> WRITER-DONE`. This cycle-2 PASS appends `CHECKER-DONE` and sets `checker_passed_at` to the same transition timestamp.
- CHK-5: PASS. `code_coverage.status: NO` remains consistent. Focused VLQ/vector-like/top-partner/bottom-partner/collider_rs greps over `quarkConstraints`, `qcd`, `flavorConstraints`, `neutrinos`, `yukawa`, `warpConfig`, `solvers`, `scanParams`, and `tests` produced no implementation hits. The broad ATLAS/CMS/LHC grep found only `tests/test_alpha_s.py:89`, an unrelated CMS/RunDec alpha_s example. The branching/cross-section grep found unrelated mu->e gamma and `ModernPointArtifactHeader` text, not a collider VLQ likelihood or mass-exclusion gate.
- CHK-6: PASS. `implementation_difficulty: HIGH` is still justified for both entries because a live constraint would require spectrum and branching-fraction matching, cross sections, acceptance/efficiency maps or simplified likelihoods, and validation against ATLAS/CMS searches; a hard-coded threshold would not be faithful.
- CHK-7: PASS. Both entries use the local quark-scan `47.26 TeV` low-energy flavor comparison as context and do not contradict the rc1.1/methodology-note result. The direct LHC partner reach is explicitly framed as a weaker collider cross-check for non-anarchic, softened-flavor, or light-partner scenarios.
- CHK-8: PASS. No unresolved `\cite`, `\ref`, `\textbf{CHECK}`, `TODO`, `FIXME`, `PLACEHOLDER`, undefined/unresolved markers, or `??` markers were found in CR008/CR010 TeX/YAML. Delimiter checks were balanced: CR008 `\(`/`\)` = 29/29 and `\[`/`\]` = 2/2; CR010 `\(`/`\)` = 51/51 and `\[`/`\]` = 2/2; dollar math count was zero for both. `git diff --check` was clean on the relevant TeX/YAML paths before status updates.
