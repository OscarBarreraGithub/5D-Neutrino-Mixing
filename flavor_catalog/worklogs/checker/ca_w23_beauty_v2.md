# CA Worklog: ca_w23_beauty_v2
**Date**: 2026-05-16
**Family**: beauty
**Process IDs**: B017 B018 B032

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| B017 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| B018 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| B032 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| B017 | BR(B0 -> K*(892)0 ell+ ell-) | (9.9 +/- 1.2)e-7 | HFLAV Dec. 2025 | 69c60b8b18cf |
| B017 | BR(B -> X_s ell+ ell-) | (5.84 +/- 0.69)e-6 | HFLAV Dec. 2025 | 06feaa8868e8 |
| B017 | BR(B+ -> K+ ell+ ell-) | (5.76 +/- 0.40)e-7 | HFLAV Dec. 2025 | 08bac0f7f097 |
| B017 | R_K* low q2 | 0.927 +0.093/-0.087 stat +0.036/-0.035 syst | LHCb 2023 | 7f1fc972e825 |
| B017 | R_K* central q2 | 1.027 +0.072/-0.068 stat +0.027/-0.026 syst | LHCb 2023 | 7f1fc972e825 |
| B018 | R_K low q2 | 0.994 +/- 0.090 | HFLAV Dec. 2025 | bcca89f01ff6 |
| B018 | R_K central q2 | 0.947 +/- 0.047 | HFLAV Dec. 2025 | 43fd3133f5ac |
| B018 | R_K full q2, Belle only | 1.08 +/- 0.16 | HFLAV Dec. 2025 | f33098cc54f3 |
| B018 | superseded LHCb 2021 R_K | 0.846 +0.042/-0.039 stat +0.013/-0.012 syst; 3.1 sigma | LHCb 2021 | 6f0e92c95b3 |
| B032 | BR(B+ -> K0 pi+) | (2.392 +/- 0.062)e-5 | HFLAV Dec. 2025 | 522cee5fd679 |
| B032 | BR(B+ -> K+ pi0) | (1.322 +/- 0.043)e-5 | HFLAV Dec. 2025 | 522cee5fd679 |
| B032 | BR(B0 -> K+ pi-) | (2.007 +/- 0.040)e-5 | HFLAV Dec. 2025 | 522cee5fd679 |
| B032 | BR(B0 -> K0 pi0) | (1.012 +/- 0.043)e-5 | HFLAV Dec. 2025 | 522cee5fd679 |
| B032 | A_CP Kpi modes | -2.67 +/- 0.87; +2.7 +/- 1.2; -8.31 +/- 0.31; -1 +/- 13 percent | HFLAV Dec. 2025 | 522cee5fd679 |
| B032 | C(K0 pi0), S(K0 pi0) | 0.00 +/- 0.08; 0.64 +/- 0.13 | PDG 2025 API | dcd2aea4b34a |
| B032 | LHCb A_CP(B+ -> K+ pi0) | 0.025 +/- 0.015 stat +/- 0.006 syst +/- 0.003 external | LHCb 2021 | 3cf10c6a5447 |
| B032 | Belle II Kpi sum-rule test | -0.03 +/- 0.13 stat +/- 0.04 syst | Belle II 2024 | 921fc4417a5f |

## Issues (if any)
- B032: CHK-1 fails under the strict re-check prompt. The new `pdg_or_equivalent.post_2008_measurements` blocks for LHCb `A_CP(B+ -> K+ pi0)` and the Belle II Kpi sum-rule value are complete, but the main HFLAV/PDG value blocks are still incomplete at the per-value level: branching-fraction blocks lack `year` and `access_date`; HFLAV direct-CP blocks lack `year`, `access_date`, and `source_url`; PDG time-dependent-CP blocks lack `year` and `access_date`.

## Evidence notes
- CHK-1: numeric TeX claims were grepped with `rg -n "[0-9]"`. B017 and B018 physical values map to `pdg_or_equivalent` blocks with value, uncertainty, year, source URL, access date, snapshot path, and sha256. B032 source values match snapshots, but metadata completeness fails as noted above.
- CHK-2: all cited process-local reference keys resolve to `flavor_catalog/references/<process_id>/source_manifest.yaml`. All manifest `snapshot_path` fields are non-empty and point to tracked files according to `git ls-files --error-unmatch`.
- CHK-3: `ls flavor_catalog/references/B017/`, `B018/`, and `B032/` showed only `.txt` snapshots and `source_manifest.yaml`; `find ... -iname '*.pdf'` returned no PDFs.
- CHK-4: before this CA update, all three sidecars had `WRITER-INITIATED` followed by WA/WA-v2 `WRITER-DONE` transitions with ISO 8601 timestamps, and the most recent entry was `WRITER-DONE`.
- CHK-5: the literal requested `rg -l -E "<process keyword>" ...` form fails in this installed ripgrep because `-E` is parsed as an encoding flag. Equivalent `rg -l -e` searches found only false positives or adjacent Delta-F=2/notebook matches, not live B017/B018/B032 implementations. The cited lines exist: `quarkConstraints/modern/phenomenology.py:23`, `quarkConstraints/modern/phenomenology.py:165`, `quarkConstraints/deltaf2.py:209`, and `flavorConstraints/muToEGamma.py:75`.
- CHK-6: HIGH is consistent for all three processes because B017/B018 require new Delta-B=1 semileptonic Hamiltonian, Wilson-coefficient, bin, and likelihood machinery, while B032 requires a new nonleptonic b -> s qbar q amplitude/likelihood layer; the existing Delta-F=2 SLL/SLR/VLL/VRR/LR1/LR2 basis does not cover them.
- CHK-7: `git grep` against tag `quarkscan-paper-rc1.1` found plan/audit references and unrelated scan-output numbers, but no rc1.1 paper load-bearing values for these catalog observables; no contradictions found.
- CHK-8: no `\textbf{CHECK}`, TODO/FIXME, `\ref{`, or `\cite{` tokens were found in the three TeX entries, and brace counts are balanced. `chktex` and `lacheck` are not installed in this environment.
