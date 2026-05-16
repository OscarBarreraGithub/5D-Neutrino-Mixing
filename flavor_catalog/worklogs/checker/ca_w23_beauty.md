# CA Worklog: ca_w23_beauty
**Date**: 2026-05-16
**Family**: beauty
**Process IDs**: B002, B005, B017, B018, B025, B032, B033

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| B002 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| B005 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| B017 | PASS | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| B018 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| B025 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| B032 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| B033 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| B002 | sin(2 beta), all charmonium | 0.710 +/- 0.011 | HFLAV Summer 2025 | ea3e6af5b240 |
| B002 | S(B_d -> J/psi K_S) | 0.712 +/- 0.011 | HFLAV Summer 2025 | ea3e6af5b240 |
| B002 | beta physical solution | 22.63 +0.45/-0.44 deg | HFLAV Summer 2025 | ea3e6af5b240 |
| B002 | LHCb S_psi K_S input | 0.717 +/- 0.013 stat +/- 0.008 syst | LHCb 2024 | 6833aa00aa16 |
| B002 | Belle II S input | 0.724 +/- 0.035 stat +/- 0.009 syst | Belle II 2024 | 6ca46cf84f70 |
| B002 | penguin phase bound | abs(Delta phi_d) <= 0.68 deg | Frings-Nierste-Wiebusch 2015 | dcd9b65df65a |
| B005 | BR(B_s0 -> mu+ mu-) | (3.34 +/- 0.27)e-9 | PDG live/API 2026 | c42edfb94b6a |
| B005 | SM BR(B_s -> mu+ mu-) | (3.64 +/- 0.12)e-9 | Czaja-Misiak 2024 | ab871a13e9b0 |
| B005 | HFLAV auxiliary average | (3.45 +/- 0.29)e-9 | HFLAV Apr. 2023 | cb9109a94b5b |
| B005 | CMS input BR(B_s0 -> mu+ mu-) | 3.83e-9 with stat/syst/fsfu errors | CMS 2023 | b821377272dc |
| B005 | LHCb input BR(B_s0 -> mu+ mu-) | 3.09e-9 with stat/syst errors | LHCb 2022 | fc1427e831e0 |
| B005 | ATLAS input BR(B_s0 -> mu+ mu-) | 2.8e-9 +0.8/-0.7 | ATLAS 2019 | 9724967bfdd6 |
| B017 | BR(B0 -> K*(892)0 ell+ ell-) | (9.9 +/- 1.2)e-7 | HFLAV Dec. 2025 | 69c60b8b18cf |
| B017 | BR(B -> X_s ell+ ell-) | (5.84 +/- 0.69)e-6 | HFLAV Dec. 2025 | 06feaa8868e8 |
| B017 | BR(B+ -> K+ ell+ ell-) | (5.76 +/- 0.40)e-7 | HFLAV Dec. 2025 | 08bac0f7f097 |
| B017 | R_K* low q2 | 0.927 +0.093/-0.087 stat +0.036/-0.035 syst | LHCb 2023 | 7f1fc972e825 |
| B017 | R_K* central q2 | 1.027 +0.072/-0.068 stat +0.027/-0.026 syst | LHCb 2023 | 7f1fc972e825 |
| B018 | R_K low q2 | 0.994 +/- 0.090 | HFLAV Dec. 2025 | bcca89f01ff6 |
| B018 | R_K central q2 | 0.947 +/- 0.047 | HFLAV Dec. 2025 | 43fd3133f5ac |
| B018 | R_K full q2, Belle only | 1.08 +/- 0.16 | HFLAV Dec. 2025 | f33098cc54f3 |
| B018 | superseded LHCb 2021 R_K | 0.846 +0.042/-0.039 stat +0.013/-0.012 syst; 3.1 sigma | LHCb 2021 | 6f0e92c95b3b |
| B025 | R(D) | 0.358 +/- 0.024 | HFLAV CKM 2025 | 4da37c4b6e64 |
| B025 | R(D*) joint-fit context | 0.281 +/- 0.011 | HFLAV CKM 2025 | 4da37c4b6e64 |
| B025 | correlation rho(RD,RD*) | -0.374 | HFLAV CKM 2025 | 4da37c4b6e64 |
| B025 | chi2/dof and CL | 16.683/14, CL 0.27 | HFLAV CKM 2025 | 4da37c4b6e64 |
| B025 | SM reference R(D), R(D*) | 0.296 +/- 0.004; 0.254 +/- 0.005 | HFLAV CKM 2025 | 4da37c4b6e64 |
| B025 | combined discrepancy | p=1.48e-4, about 3.8 sigma | HFLAV CKM 2025 | 4da37c4b6e64 |
| B032 | BR(B+ -> K0 pi+) | (2.392 +/- 0.062)e-5 | HFLAV Dec. 2025 | 522cee5fd679 |
| B032 | BR(B+ -> K+ pi0) | (1.322 +/- 0.043)e-5 | HFLAV Dec. 2025 | 522cee5fd679 |
| B032 | BR(B0 -> K+ pi-) | (2.007 +/- 0.040)e-5 | HFLAV Dec. 2025 | 522cee5fd679 |
| B032 | BR(B0 -> K0 pi0) | (1.012 +/- 0.043)e-5 | HFLAV Dec. 2025 | 522cee5fd679 |
| B032 | A_CP Kpi modes | -2.67 +/- 0.87; +2.7 +/- 1.2; -8.31 +/- 0.31; -1 +/- 13 percent | HFLAV Dec. 2025 | 522cee5fd679 |
| B032 | C(K0 pi0), S(K0 pi0) | 0.00 +/- 0.08; 0.64 +/- 0.13 | PDG 2025 API | dcd2aea4b34a |
| B033 | sin(2 beta_eff)(phi K0) | 0.74 +/- 0.12 | HFLAV Summer 2025 | f02d219f3508 |
| B033 | C_CP(phi K0) | -0.09 +/- 0.12 | HFLAV Summer 2025 | f02d219f3508 |
| B033 | comparison sin(2 beta) | 0.710 +/- 0.011; J/psi K_S 0.712 +/- 0.011 | HFLAV Summer 2025 | f02d219f3508 |
| B033 | Belle II S and C input | S=0.54 +/- 0.26 +0.06/-0.08; C=-0.31 +/- 0.20 +/- 0.05 | Belle II 2023 | 46d751b86a3d |

## Issues (if any)
- B017: CHK-2 fails. `flavor_catalog/references/B017/source_manifest.yaml` and every B017 snapshot are present only as untracked files. `git ls-files --error-unmatch flavor_catalog/references/B017/source_manifest.yaml flavor_catalog/references/B017/*.txt` reports that none are known to Git, so cited reference keys do not resolve to tracked text snapshots.
- B018: CHK-1 fails. The `pdg_or_equivalent` numeric entries for the `R_K` values include values, uncertainties, snapshot paths, and hashes, but the per-value blocks omit the required `source_url`, `year`, and `access_date`; the top `pdg_or_equivalent` block also lacks `source_url`.
- B032: CHK-1 fails. `B032.tex` includes post-2008 numerical claims for LHCb `A_CP(B+ -> K+ pi0)=0.025 +/- 0.015 +/- 0.006 +/- 0.003` and the Belle II `Kpi` sum-rule value `-0.03 +/- 0.13 +/- 0.04`; these values appear in local snapshots but not in corresponding `pdg_or_equivalent` YAML value blocks.

## Evidence notes
- CHK-2/CHK-3: all non-B017 manifest entries have non-empty `snapshot_path` fields pointing to tracked `.txt` files. `find ... -name '*.pdf'` returned no PDFs for the seven process reference directories.
- CHK-4: all seven sidecars had `WRITER-INITIATED` followed by `WRITER-DONE` with ISO 8601 timestamps before CA status updates.
- CHK-5: line citations in `quarkConstraints/deltaf2.py`, `quarkConstraints/qcd_running.py`, `quarkConstraints/modern/evaluation.py`, `quarkConstraints/modern/phenomenology.py`, and `flavorConstraints/muToEGamma.py` exist and match the described adjacent coverage. The literal requested `rg -l -E` form is invalid for this installed ripgrep because `-E` is parsed as encoding; the equivalent `rg -l -e` searches were used. NO-coverage hits were false positives or adjacent Delta-F=2 / notebook matches, not live implementations of the listed processes.
- CHK-6: B002 LOW is consistent with reusing the existing Delta-F=2 operator basis; B005 MEDIUM is consistent with adding a standard Delta-B=1 leptonic rare-decay operator layer; B017, B018, B025, B032, and B033 HIGH are consistent with new semileptonic/nonleptonic mode calculations or new likelihood machinery.
- CHK-7: `git grep` against tag `quarkscan-paper-rc1.1` found only plan/audit references for this batch, not rc1.1 paper load-bearing process values; no contradictions found.
- CHK-8: no `\textbf{CHECK}`, TODO/FIXME, `\ref`, or `\cite` tokens were found in the seven TeX entries, and brace counts are balanced. `chktex` and `lacheck` are not installed in this environment.
