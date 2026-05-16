# CA Worklog: ca_w23_charged_lepton_v3
**Date**: 2026-05-16
**Family**: charged_lepton
**Process IDs**: L001

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| L001 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| L001 | BR(mu+ -> e+ gamma), primary current limit | < 1.5 x 10^-13 at 90% C.L.; sensitivity 2.2 x 10^-13 | MEG II 2025 arXiv:2504.15711 | 3560d6556552 |
| L001 | BR(mu+ -> e+ gamma), PDG listing | < 3.1 x 10^-13 at 90% C.L. | PDG 2025 muon listing S004 | 4baeea03e1ea |
| L001 | BR(mu+ -> e+ gamma), MEG II first dataset / combination | < 7.5 x 10^-13 and < 3.1 x 10^-13 at 90% C.L. | MEG II arXiv:2310.12614 | bf3dddb7b5fd |
| L001 | BR(mu+ -> e+ gamma), MEG full dataset | < 4.2 x 10^-13 at 90% C.L. | MEG arXiv:1605.05081 | 883288d06a72 |
| L001 | Perez-Randall paper-era formula context | BR <= 1.2 x 10^-11; prefactor 4 x 10^-8; C = 0.02; reference scale 3 TeV | Perez-Randall arXiv:0805.4652 | 850dccbe5532 |
| L001 | Repo scan default BR limit | 1.5 x 10^-13 | scanParams/scan.py local source | 4aab2bec426a |
| L001 | Repo LFV prefactor and coefficients | prefactor 4.0 x 10^-8; C ~= 1.936 x 10^-3; C_PAPER = 0.02; reference_scale = 3000 GeV in code | flavorConstraints/muToEGamma.py local source | bf4e2ee1ea1c |

## Check evidence
- CHK-1: Grepped `L001.tex` for numerical claims. The BR limits, sensitivity, prefactor, `C` values, and repo defaults trace to `pdg_or_equivalent` source blocks with hashes. Strict failure remains because `L001.tex:42-44` includes the formula factor `(3 TeV / M_KK)^4`, while `L001.yaml:pdg_or_equivalent` has no explicit value-bearing `reference_scale` entry for `3 TeV` / `3000 GeV` with year, value, uncertainty, source URL, access date, and sha256.
- CHK-2: PASS. `L001.tex` no longer cites `arXiv:0804.1954`. The TeX source keys `MEGII2025_MuEGamma`, `MEGII2024_MuEGammaFirst`, `MEG2016_MuEGammaFull`, `PDG2025_MuonListing`, and `PerezRandall2008_WarpedNeutrinos` all appear in `flavor_catalog/references/L001/source_manifest.yaml`, and each manifest `snapshot_path` is a tracked text snapshot.
- CHK-3: PASS. `ls flavor_catalog/references/L001/` shows only text snapshots, `sha256sums.txt`, and `source_manifest.yaml`; `git ls-files flavor_catalog/references/L001/*.pdf` returned no tracked PDFs.
- CHK-4: PASS on the pre-CA state. `status_history` contains `WRITER-INITIATED` before the WA `WRITER-DONE` transitions, all with ISO 8601 timestamps, and the most recent pre-check state was `WRITER-DONE` at `2026-05-16T13:11:19-04:00`.
- CHK-5: PASS. The prompt's literal `rg -l -E "mu.?e.?gamma" ... 2>&1` form was run and fails in this ripgrep because `-E` is parsed as an encoding flag. The semantic rerun with `rg -l -e` found live coverage in `flavorConstraints/muToEGamma.py`, `scanParams/scan.py`, `tests/test_mu_to_e_gamma.py`, and `tests/test_scan.py`. The cited lines exist: `flavorConstraints/muToEGamma.py:1`, `:20`, `:31`, `:61`, `:75`, `:136`; `scanParams/scan.py:32`, `:33`, `:264`; `tests/test_mu_to_e_gamma.py:40`; `tests/test_scan.py:50`.
- CHK-6: PASS. `LOW` is consistent for the current catalog entry because the `mu -> e gamma` LFV dipole/NDA check and scan default are already implemented. The sidecar's `MEDIUM` note for a fresh missing implementation is also consistent because a new lepton-dipole observable would be needed but not new lattice inputs or new RG.
- CHK-7: PASS. `git grep` against `quarkscan-paper-rc1.1` found the same L001 LFV code/default context and no silently contradicted load-bearing rc1.1 quark-paper number.
- CHK-8: PASS. No `\textbf{CHECK}`, `\cite{}`, `\ref{}`, TODO/FIXME, or unresolved-reference markers were found in `L001.tex`; manual inspection found no cosmetic math-mode or brace issue. `chktex` is not installed in this environment.

## Issues (if any)
- L001: CHK-1 fails under the strict sidecar traceability rule. Add an explicit `pdg_or_equivalent` value block for the `3 TeV` / `3000 GeV` reference scale used in the Perez-Randall NDA formula, with year, value, units, `uncertainty: null`, source URL, access date, snapshot path, and sha256. If the scalar `paper_era_reference` fields remain as-is, also consider normalizing them into value-bearing subentries before Opus arbitration.
