# CA Worklog: ca_w23_charged_lepton
**Date**: 2026-05-16
**Family**: charged_lepton
**Process IDs**: L001 L002 L007

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| L001 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| L002 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| L007 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| L001 | BR(mu+ -> e+ gamma), primary limit | < 1.5 x 10^-13 at 90% CL; sensitivity 2.2 x 10^-13 | MEG II 2025 arXiv:2504.15711 | 3560d6556552 |
| L001 | BR(mu+ -> e+ gamma), PDG listing | < 3.1 x 10^-13 at 90% CL | PDG 2025 muon listing S004 | 4baeea03e1ea |
| L001 | BR(mu+ -> e+ gamma), MEG 2016 | < 4.2 x 10^-13 at 90% CL | MEG arXiv:1605.05081 | 883288d06a72 |
| L001 | MEG II 2021-only and MEG+MEG II combination | < 7.5 x 10^-13; < 3.1 x 10^-13 at 90% CL | MEG II arXiv:2310.12614 | bf3dddb7b5fd |
| L001 | Perez-Randall paper-era BR/prefactor/C | BR <= 1.2 x 10^-11; prefactor 4 x 10^-8; C = 0.02 | Perez-Randall arXiv:0805.4652 | 850dccbe5532 |
| L002 | BR(mu+ -> e+e-e+), PDG current limit | < 1.0 x 10^-12 at 90% CL | PDG 2026 pdgLive S004R4 | 72b7c9181196 |
| L002 | Gamma(mu -> 3e) / Gamma(mu -> e 2nu) | < 1.0e-12 at 90% CL | SINDRUM/Bellgardt INSPIRE metadata | 87086d504690 |
| L002 | Mu3e phase-I sensitivity and rate | 2e-15 single-event sensitivity; 10^8 muon decays/s | Mu3e TDR arXiv:2009.11690 | 4123a6876f35 |
| L002 | Mu3e status/prospects and SM rate | O(10^-15), O(10^-16), O(10^-54) | Mu3e status arXiv:2501.14667 | 7d4cc8303d28 |
| L007 | BR(tau- -> mu- gamma), PDG current limit | < 4.2 x 10^-8 at 90% CL | PDG 2025 pdgLive S035.31 | b94475dc141a |
| L007 | BR(tau+- -> mu+- gamma), Belle | <= 4.2 x 10^-8 at 90% CL; 988 fb^-1 | Belle arXiv:2103.12994 | 99810fa89b36 |
| L007 | BR(tau -> mu gamma), BaBar | < 4.4 x 10^-8 at 90% CL; (963 +- 7) x 10^6 tau decays | BaBar arXiv:0908.2381 | 2528b26a9cba |
| L007 | Belle II projection | roughly < 5 x 10^-9 with 50 ab^-1 | Belle II Physics Book arXiv:1808.10567 | d3789bb9ccd9 |

## Check Evidence
- CHK-1: Grepped each TeX for numerical claims and matched load-bearing values to sidecar `pdg_or_equivalent` entries and local snapshots. L001 and L002 fail the strict metadata rule because some nested numerical entries lack per-entry `source_url`/`access_date` fields; L007 satisfies the rule for its current, prior, and projected numerical values.
- CHK-2: All process-local reference keys listed in the TeX resolve to entries in `flavor_catalog/references/<process_id>/source_manifest.yaml`; each manifest entry has a non-empty `snapshot_path`, and `git ls-files` shows the snapshots are tracked.
- CHK-3: `ls flavor_catalog/references/{L001,L002,L007}/` shows text snapshots and manifests only; `find ... -name '*.pdf'` returned no publisher PDFs.
- CHK-4: Each sidecar contains `WRITER-INITIATED` before `WRITER-DONE` with ISO 8601 timestamps, and the most recent pre-CA entry was `WRITER-DONE`.
- CHK-5: The literal prompt form `rg -l -E "<process keyword>" ...` fails because this ripgrep treats `-E` as an encoding flag, so I reran the semantic checks with `-e`. L001 has live `mu->e gamma` coverage in `flavorConstraints/muToEGamma.py`, `scanParams/scan.py`, and tests. L002 and L007 targeted searches returned no implementation hits in required code paths.
- CHK-6: L001 LOW is consistent with plan v1 because the `mu->e gamma` pattern already exists. L002 MEDIUM is consistent because it needs a new four-lepton/contact observable without lattice inputs. L007 MEDIUM is consistent because it can reuse the dipole pattern but needs tau-mu indices and tau-rate normalization.
- CHK-7: No load-bearing rc1.1 paper number is silently contradicted. L001 preserves the paper-era Perez-Randall `C=0.02` context while separately documenting the current repo default; L002 and L007 are companion catalog entries.
- CHK-8: No `\textbf{CHECK}`, TODO/FIXME, `\cite`, `\ref`, or unresolved-reference markers were found in the three TeX files. Spot inspection found no blocking math-mode or brace issues.

## Issues (if any)
- L001: CHK-1 fails. `L001.tex` cites prior MEG/MEG II limits and the live repo-default `C\simeq1.936\times10^{-3}`; corresponding `pdg_or_equivalent.prior_experimental_limits` and `repo_default` entries do not carry complete strict metadata fields (`year`/`source_url`/`access_date`/`sha256` as applicable to each numerical claim).
- L002: CHK-1 fails. `L002.tex` cites Mu3e sensitivity/prospect values, the `10^8` muon-decay-rate context, and the `O(10^-54)` SM rate; corresponding `pdg_or_equivalent.prospects` entries lack per-entry `source_url` and `access_date` metadata.
