# CA Worklog: ca_w23_charged_lepton_v2
**Date**: 2026-05-16
**Family**: charged_lepton
**Process IDs**: L001 L002

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| L001 | PASS | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| L002 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| L001 | BR(mu+ -> e+ gamma), primary limit | < 1.5 x 10^-13 at 90% C.L.; sensitivity 2.2 x 10^-13 | MEG II 2025 arXiv:2504.15711 | 3560d6556552 |
| L001 | BR(mu+ -> e+ gamma), PDG listing | < 3.1 x 10^-13 at 90% C.L. | PDG 2025 muon listing S004 | 4baeea03e1ea |
| L001 | BR(mu+ -> e+ gamma), MEG II first dataset / combination | < 7.5 x 10^-13 and < 3.1 x 10^-13 at 90% C.L. | MEG II arXiv:2310.12614 | bf3dddb7b5fd |
| L001 | BR(mu+ -> e+ gamma), MEG full dataset | < 4.2 x 10^-13 at 90% C.L. | MEG arXiv:1605.05081 | 883288d06a72 |
| L001 | Perez-Randall paper-era BR/prefactor/C/reference scale | BR <= 1.2 x 10^-11; prefactor 4 x 10^-8; C = 0.02; 3 TeV | Perez-Randall arXiv:0805.4652 | 850dccbe5532 |
| L001 | Repo scan default BR limit | 1.5 x 10^-13 | scanParams/scan.py local source | 4aab2bec426a |
| L001 | Repo LFV prefactor and coefficients | prefactor 4.0 x 10^-8; C ~= 1.936 x 10^-3; C_PAPER = 0.02 | flavorConstraints/muToEGamma.py local source | bf4e2ee1ea1c |
| L002 | BR(mu+ -> e+e-e+), PDG current limit | < 1.0 x 10^-12 at 90% C.L. | PDG 2026 pdgLive S004R4 | 72b7c9181196 |
| L002 | Gamma(mu -> 3e) / Gamma(mu -> e 2nu) | < 1.0 x 10^-12 at 90% C.L. | SINDRUM/Bellgardt INSPIRE metadata | 87086d504690 |
| L002 | Mu3e phase-I sensitivity and rate | 2 x 10^-15 single-event sensitivity; 10^8 muon decays/s | Mu3e TDR arXiv:2009.11690 | 4123a6876f35 |
| L002 | Mu3e status/prospects and SM rate | O(10^-15), O(10^-16), O(10^-54) | Mu3e status arXiv:2501.14667 | 7d4cc8303d28 |

## Check evidence
- CHK-1: Grepped both TeX files for numerical claims. The WA-v2 metadata additions cover the prior MEG/MEG II limits, repo-default L001 numbers, Mu3e prospect numbers, 10^8 muon-rate context, and O(10^-54) SM rate with year/source URL/access date/sha256 metadata in `pdg_or_equivalent`.
- CHK-2: L002 passes: all seven TeX source keys resolve to `flavor_catalog/references/L002/source_manifest.yaml`, and each manifest entry has a tracked text snapshot. L001 fails: `L001.tex:51` cites arXiv:0804.1954, but `flavor_catalog/references/L001/source_manifest.yaml` and `flavor_catalog/references/L001/` contain no matching CFW/0804.1954 entry.
- CHK-3: `ls flavor_catalog/references/{L001,L002}/` shows text snapshots, manifests, and L001 `sha256sums.txt`; `find ... -name '*.pdf'` returned no PDFs.
- CHK-4: Both sidecars contain `WRITER-INITIATED` before `WRITER-DONE` with ISO 8601 timestamps, and the most recent pre-CA entry was the WA-v2 `WRITER-DONE`.
- CHK-5: The prompt's literal `rg -l -E` form fails in this ripgrep because `-E` is parsed as an encoding flag. Rerunning the semantic checks with `-e` found L001 coverage in `flavorConstraints/muToEGamma.py`, `scanParams/scan.py`, and tests, while the L002 mu->3e/Mu3e/contact search returned no hits in the required code paths.
- CHK-6: L001 LOW is consistent for cataloging because the mu->e gamma LFV dipole check is already implemented and tested. L002 MEDIUM is consistent because it needs a new four-lepton/contact observable and matching hook but no lattice or hadronic long-distance input for a first implementation.
- CHK-7: No load-bearing rc1.1 quark-paper number is silently contradicted. L001 preserves the paper-era Perez-Randall `C=0.02` context while separately documenting the current repo default; L002 is a companion charged-lepton entry.
- CHK-8: Fixed-string greps found no `\textbf{CHECK}`, TODO/FIXME, `\cite`, `\ref`, `\label`, or `??` markers. Spot inspection found no blocking math-mode or brace issue.

## Issues (if any)
- L001: CHK-2 fails. `L001.tex:51` cites `arXiv:0804.1954` as the RS-flavor baseline, but `flavor_catalog/references/L001/source_manifest.yaml` has no corresponding key and no tracked process-local snapshot. Rework should either add a L001 process-local source-manifest entry plus tracked text snapshot for CFW arXiv:0804.1954, or remove/replace that explicit TeX citation.
