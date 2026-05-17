# CA Worklog: ca_w9_custodial
**Date**: 2026-05-17
**Family**: collider_rs
**Batch ID**: CA-w9-custodial
**Cycle**: 1
**Checker agent**: CA-CA-w9-custodial
**Process IDs**: CR002 CR003 CR004

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| CR002 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| CR003 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| CR004 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| CR002 | m(t'(5/3)) lower limit, exclusive B/X -> Wt pair-production interpretation | > 1460 GeV, 95% CL | PDGEncoderQ009_TPrime5Over3 / ATLAS2023_VLQ_MET | 093c2771ed73 |
| CR002 | mass-degenerate vector-like doublet lower limit | > 1590 GeV, 95% CL | ATLAS2023_VLQ_MET | fb8975d1d181 |
| CR002 | m(X_5/3) lower limit, right-handed coupling benchmark | > 1330 GeV, 95% CL | CMS2019_X53_Run2 | e14cb8f94181 |
| CR002 | m(X_5/3) lower limit, left-handed coupling benchmark | > 1300 GeV, 95% CL | CMS2019_X53_Run2 | e14cb8f94181 |
| CR002 | m(T_5/3) lower limit from ATLAS pair production | > 1190 GeV, 95% CL | ATLAS2018_SameCharge_BJets | 700e58bb4698 |
| CR002 | m(T_5/3) lower limit with single production included, unit coupling | > 1600 GeV, 95% CL | PDGEncoderQ009_TPrime5Over3 / ATLAS2018_SameCharge_BJets | 093c2771ed73 |
| CR002 | m(X_5/3) lower limit, early CMS right-handed benchmark | > 1020 GeV, 95% CL | CMS2017_X53_13TeV | 1757a76e808d |
| CR002 | m(X_5/3) lower limit, early CMS left-handed benchmark | > 990 GeV, 95% CL | CMS2017_X53_13TeV | 1757a76e808d |
| CR002 | m(X_5/3) lower limit, CMS 8 TeV same-sign dilepton search | > 800 GeV, 95% CL | CMS2014_X53_8TeV | ffdb33839c2b |
| CR002 | repo-local rc1.1 quark-scan M_KK p50 context, g_* = 3 | 47.26 TeV | docs/quark_scan_methodology_note.tex | 334d639d4167 |
| CR003 | m_T lower limit, B(T -> Wb) = 1 | > 1.70 TeV, 95% CL | PDGLive2026_Q009TPP_TprimePair / ATLAS2024_VLQPairWb_Run2 | 6cc3eac0c3ce |
| CR003 | m_T lower limit, singlet branching pattern Wb:Ht:Zt = 1/2:1/4:1/4 | > 1.36 TeV, 95% CL | PDGLive2026_Q009TPP_TprimePair / ATLAS2024_VLQPairWb_Run2 | 6cc3eac0c3ce |
| CR003 | m_T lower limit, B(T -> Zt) = 1 | > 1.60 TeV, 95% CL | PDGLive2026_Q009TPP_TprimePair / ATLAS2023_VLQPairZt_Run2 | 6cc3eac0c3ce |
| CR003 | m_T lower limit, B(T -> Ht) = 1 | > 1.50 TeV, 95% CL | PDGLive2026_Q009TPP_TprimePair / CMS2023_VLQPairLeptonic_Run2 | 6cc3eac0c3ce |
| CR003 | m_T lower limit for all third-generation Wb/Zt/Ht decay mixtures | > 1.48 TeV, 95% CL | CMS2023_VLQPairLeptonic_Run2 | 8e2d5e1a579b |
| CR003 | repo-local rc1.1 quark-scan M_KK p50 / p95 context, g_* = 3 | 47.26 TeV / 127.13 TeV | QuarkScanMethodology2026_MKKEnvelope | 36da0a9c3d45 |
| CR004 | m_B lower limit, B(B -> H b) = 1 | > 1570 GeV, 95% CL | PDG2025_BPrimeListing / CMS2024_BVLQDilepHad | bca9c8ebe557 |
| CR004 | m_B lower limit, B(B -> Z b) = 1 | > 1540 GeV, 95% CL | PDG2025_BPrimeListing / CMS2024_BVLQDilepHad | bca9c8ebe557 |
| CR004 | m_B lower limit, B(B -> W t) = 1 | > 1560 GeV, 95% CL | PDG2025_BPrimeListing / CMS2023_VLQLeptonic | bca9c8ebe557 |
| CR004 | repo-local rc1.1 quark-scan M_KK p50 / p95 context, g_* = 3 | 47.26 TeV / 127.13 TeV | docs/quark_scan_methodology_note.tex | 334d639d4167 |

## Issues
- CR002: None.
- CR003: None.
- CR004: None.

## Evidence notes
- CHK-1: All measured mass-exclusion limits quoted in the TeX are present in each sidecar's `pdg_or_equivalent.values` with value/display, units, limit type, CL, year, experiment, source URL, access date, snapshot path, and sha256. I applied the L001 / B001_B003 / B021_B023 carve-out: integrated luminosities, sqrt(s), HEPData contour context without a quoted table value, coupling assumptions, theory baselines, and repo-local M_KK comparison values are supporting or auxiliary context and do not require `pdg_or_equivalent` placement.
- CHK-2: Key-reference `texttt` entries in CR002, CR003, and CR004 resolve to the corresponding process `source_manifest.yaml` keys. All manifest entries have non-empty `snapshot_path` fields; a `git ls-files --error-unmatch` loop over every manifest snapshot path returned no missing tracked files. `sha256sum -c sha256sums.txt` passed for CR002, CR003, and CR004.
- CHK-3: `find flavor_catalog/references/CR002 flavor_catalog/references/CR003 flavor_catalog/references/CR004 -type f \( -iname '*.pdf' -o -iname '*.PDF' \) -print` returned no files. The reference directories contain text snapshots, manifests, and sha256 lists only.
- CHK-4: Each sidecar shows `DRAFT -> WRITER-INITIATED -> WRITER-DONE` before checking, with ISO 8601 `timestamp` and `at` fields. I appended `CHECKER-DONE` at `2026-05-17T17:35:23-04:00` with `agent_id: "CA"`, `cycle: 1`, and this worklog reference, and set `checker_passed_at` to the same timestamp.
- CHK-5: The claimed `code_coverage.status: NO` is supported. The CR002 collider-specific grep for `T5/3`, `X5/3`, VLQ, top-partner, collider, and LHC direct-search strings returned no matches. The CR003 VLQ/top-partner/direct-search grep returned no matches. The CR004 collider/VLQ query returned only `tests/test_alpha_s.py:88-89`, a CMS/RunDec alpha_s example unrelated to collider searches. Adjacent evidence lines exist for KK-scale/flavor plumbing: CR002 cites `quarkConstraints/modern/scan.py:1230`, `:1254`, `quarkConstraints/couplings.py:103`, `:116`, `quarkConstraints/modern/matching.py:240`, and `quarkConstraints/deltaf2.py:1`; CR003 cites `quarkConstraints/modern/scan.py:1230`, `:1492`, `:1503`, `quarkConstraints/modern/evaluation.py:643`, and `quarkConstraints/README.md:34-38`; CR004 cites `tests/test_alpha_s.py:88-89`, `scanParams/scan.py:399`, `:523`, and `quarkConstraints/scan.py:359`, `:377`.
- CHK-6: `HIGH` is consistent for all three entries. A real implementation would need an RS custodial-fermion spectrum and branching model, event generation or public likelihood/HEPData contours, detector acceptance, and recast tooling such as CheckMATE, MadAnalysis5, or SModelS; a hard mass veto would not be a faithful collider-RS constraint.
- CHK-7: The entries explicitly state that LHC partner reaches around 1.3-1.7 TeV are weaker than the repo-local anarchic-flavor M_KK reach. The rc1.1 values quoted in the TeX/sidecars match `docs/quark_scan_methodology_note.tex`: `47.26` TeV at 50% acceptance and `127.13` TeV at 95% acceptance for `g_* = 3`. No load-bearing contradiction with the 47 TeV flavor result was found.
- CHK-8: No unresolved `\cite`, `\ref`, `\textbf{CHECK}`, `TODO`, `FIXME`, `PLACEHOLDER`, undefined, or unresolved markers were found. Literal dollar math delimiters are absent. Delimiter sanity checks returned CR002 `\(`/`\)` = 49/49, `\[`/`\]` = 3/3, brace delta 0; CR003 44/44, 2/2, brace delta 0; CR004 54/54, 6/6, brace delta 0. `git diff --check` was clean for the relevant TeX/YAML paths.
