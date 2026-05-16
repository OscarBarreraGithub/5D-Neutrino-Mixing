# CA Worklog: ca_w5b_kaon_v2
**Date**: 2026-05-16
**Family**: kaon
**Process IDs**: K009 K010

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| K009 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| K010 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| K009 | BR(K_L -> pi0 mu+ mu-) canonical limit | <3.8e-10 at 90% CL | PDG 2026 S013.16 / S013R30 | c4ebbf6334ba |
| K009 | KTeV observed/background events and limit | 2 observed; 0.87 +/- 0.15 expected background; <3.8e-10 at 90% CL | KTeV 2000 arXiv hep-ex/0001006 | 0a16b654a73f |
| K009 | BR(K_S -> pi0 mu+ mu-) supporting input | (2.9 +1.5 -1.2 stat +/- 0.2 syst)e-9; 6 events; 0.22 +0.18 -0.11 background | PDG 2026 S012.15 and NA48 2004 | 21f62a439461 / 0cee82c75bac |
| K009 | K_L -> pi0 mu+ mu- theory decomposition | C_mix=(3.7 +/- 0.1)|a_S|^2; C_int=(1.6 +/- 0.1)|a_S|; C_dir=1.0 +/- 0.1; C_CPC=5.2 +/- 1.6; |a_S|=1.2 +/- 0.2; SM constructive=(1.5 +/- 0.3)e-11 | Isidori, Smith, Unterdorfer 2004 | 9634d1af1d93 |
| K010 | BR(K_S -> pi0 e+e-), m_ee > 0.165 GeV | (3.0 +1.5 -1.2 stat +/- 0.2 syst)e-9 | PDG 2026 S012.10 / NA48 | 87a298fc5bb6 |
| K010 | BR(K_S -> pi0 e+e-) full-region extrapolation | (5.8 +2.9 -2.4)e-9 | PDG 2026 S012.10 / NA48 | 87a298fc5bb6 |
| K010 | NA48/1 observed/background events | 7 observed; 0.15 expected background | NA48/1 2003 arXiv hep-ex/0309075 | 406786ac199a |

## Evidence notes
- Required reading completed: plan v1 Section D CA spec, including the WA/CA separation rule, and `flavor_catalog_orchestrator_decisions.md`.
- CHK-1 passed for K009. Grepping `K009.tex` for load-bearing numerical claims found the PDG/KTeV 3.8e-10 90% CL limit, KTeV 2 events and 0.87 +/- 0.15 background, the supporting K_S -> pi0 mu+ mu- 2.9e-9 value with errors plus 6 events and 0.22 +0.18 -0.11 background, and the Isidori-Smith-Unterdorfer 1e-12-normalized coefficient/|a_S|/SM-expectation block. Each is now represented in `pdg_or_equivalent.values` with year, value, uncertainty, source URL, access date, snapshot path, and sha256.
- CHK-1 passed for K010. Grepping `K010.tex` for load-bearing numerical claims found the m_ee > 0.165 GeV partial value, the 5.8e-9 full-region extrapolation, and the NA48/1 7-event and 0.15-background claims. Each is represented in `pdg_or_equivalent.values` with complete metadata.
- CHK-2 passed for both. Every source filename and process-local source key cited in the TeX resolves to `flavor_catalog/references/<process_id>/source_manifest.yaml`; every manifest entry has a non-empty tracked `.txt` snapshot under the matching process directory.
- CHK-3 passed for both. `ls flavor_catalog/references/K009/` and `ls flavor_catalog/references/K010/` showed only text snapshots, manifests, and K010's `sha256sums.txt`; `find ... -name '*.pdf'` returned no files.
- CHK-4 passed for both. `status_history` contains `WRITER-INITIATED` followed by `WRITER-DONE`, then the earlier CA rework cycle, then the WA-v2 `WRITER-DONE`; all timestamps are ISO 8601 and the most recent pre-check entry was `WRITER-DONE`.
- CHK-5 passed for both. The literal requested `rg -l -E "<process keyword>" ...` form fails here because ripgrep parses `-E` as an encoding flag, so the equivalent `rg -l -e ...` searches were used. Targeted searches for K009 (`pi0 mu`, `S013.16`, `S013R30`, `a_S`, `s -> d mu`) and K010 (`pi0 e`, `S012.10`, `a_S`, `s -> d e`) over `quarkConstraints/ qcd/ flavorConstraints/ neutrinos/ yukawa/ warpConfig/ solvers/ scanParams/ tests/` returned no implementation files. The nearby cited Delta-F=2 lines exist at `quarkConstraints/deltaf2.py:209`, `quarkConstraints/deltaf2.py:755`, and `quarkConstraints/modern/phenomenology.py:23`.
- CHK-6 passed for both. HIGH is consistent with the rubric because these are new Delta-S=1 semileptonic rare-kaon observables outside the existing Delta-F=2 SLL/SLR/VLL/VRR/LR1/LR2 basis; K009 additionally needs CPV/interference and CPC two-photon handling.
- CHK-7 passed for both. The TeX does not quote rc1.1 load-bearing rare-decay numbers; CFW2008/arXiv:0804.1954 is used only as RS-flavor baseline context, so no rc1.1 catalog companion contradiction was found.
- CHK-8 passed for both. Greps found no `\textbf{CHECK}`, TODO/FIXME, `\cite`, `\ref`, `\label`, unresolved `??`, or unmatched displayed-equation delimiters in the two TeX entries.
- Snapshot/value spot checks matched local files and sha256s: PDG/KTeV/NA48/Isidori snapshots for K009 and PDG/NA48 snapshots for K010 contain the values listed above.

## Issues (if any)
- None.
