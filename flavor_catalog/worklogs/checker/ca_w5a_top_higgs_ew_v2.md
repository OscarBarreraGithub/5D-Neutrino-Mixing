# CA Worklog: ca_w5a_top_higgs_ew_v2
**Date**: 2026-05-16
**Family**: top_higgs_ew
**Process IDs**: T005, T019

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| T005 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| T019 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| T005 | ATLAS dataset for \(cg\to t\) | 139 fb^-1 at sqrt(s) = 13 TeV | PDG 2025 / ATLAS 2022 | b58775f55944 |
| T005 | sigma(cg -> t) B(t -> bW) B(W -> l nu) | < 4.7 pb, 95% CL | PDG 2025 / ATLAS 2022 | b58775f55944 |
| T005 | \|C_uG^{ct}\| / Lambda^2 | < 0.14 TeV^-2, 95% CL | PDG 2025 / ATLAS 2022 | b58775f55944 |
| T005 | B(t -> c g), ATLAS EFT interpretation | < 3.7e-4, 95% CL | PDG 2025 / ATLAS 2022 | b58775f55944 |
| T005 | CMS dataset | 5.0 fb^-1 at 7 TeV and 19.7 fb^-1 at 8 TeV | CMS 2017 | 75352c02b40e |
| T005 | \|kappa_tcg\| / Lambda | < 1.8e-2 TeV^-1, 95% CL | CMS 2017 | 75352c02b40e |
| T005 | B(t -> c g), CMS interpretation | < 4.1e-4, 95% CL | CMS 2017 | 75352c02b40e |
| T005 | SM B(t -> c g) benchmark | 4.6e-12 | Aguilar-Saavedra 2004 | 098fd502e902 |
| T005 | CFW ordinary RS KK-gluon context | about 21 TeV | Csaki-Falkowski-Weiler 2008 | 1ecc5891da03 |
| T005 | CFW composite-pNGB Higgs KK-gluon context | about 33 TeV | Csaki-Falkowski-Weiler 2008 | 1ecc5891da03 |
| T019 | B(h -> e tau), ATLAS Run-2 | < 0.20%; expected < 0.12%; 95% CL | PDG 2025 / ATLAS 2023 | 67ea71f88aaa |
| T019 | ATLAS Run-2 dataset | 138 fb^-1 at sqrt(s) = 13 TeV | ATLAS 2023 | a430814893db |
| T019 | B(h -> e tau), CMS Run-2 | < 0.22%; expected < 0.16%; 95% CL | PDG 2025 / CMS 2021 | 67ea71f88aaa |
| T019 | CMS Run-2 dataset | 137 fb^-1 at sqrt(s) = 13 TeV | CMS 2021 | 6b4f4c02719f |
| T019 | strongest observed single-experiment h -> e tau limit | < 2.0e-3 | ATLAS 2023 | a430814893db |
| T019 | generic LFV-Higgs EFT allowance | order 10% | Harnik-Kopp-Zupan 2012 | 50662238edb6 |
| T019 | B(H -> e tau), ATLAS partial Run-2 | < 0.47%; expected < 0.34% (+0.13%/-0.10%); 95% CL | ATLAS 2019 | 726a7510795c |
| T019 | ATLAS partial Run-2 dataset | 36.1 fb^-1 at sqrt(s) = 13 TeV | ATLAS 2019 | 726a7510795c |

## Evidence notes
- CHK-1: Grepped both TeX files for numbers. After WA-v2, all load-bearing physics numbers in the TeX are present in `pdg_or_equivalent.values` with year, value, uncertainty or null uncertainty, source URL, access date, and sha256. Non-physics dates, catalog ids, source ids, and code line numbers were checked under CHK-2, CHK-4, and CHK-5 rather than treated as physics observables.
- CHK-2: Every cited source key in the TeX resolves to `flavor_catalog/references/<process_id>/source_manifest.yaml`. Value-id texttt tags such as `PDG2025:T019:atlas_run2` resolve to YAML `pdg_or_equivalent.values`. Each manifest entry has a non-empty process-local `snapshot_path`.
- CHK-3: `ls flavor_catalog/references/T005` and `ls flavor_catalog/references/T019` show only `.txt` snapshots, `source_manifest.yaml`, and T005's `sha256sums.txt`; `find ... -name '*.pdf'` returned no PDF files.
- CHK-4: Before this CA update, both sidecars contained `WRITER-INITIATED` followed by `WRITER-DONE` with ISO 8601 timestamps, and the latest pre-checker entry was the cycle-2 `WRITER-DONE` from WA-v2. This re-check appended `CHECKER-DONE`, set `checker_passed_at`, and updated `last_updated_at`.
- CHK-5: The prompt's literal `rg -l -E` form errors in this ripgrep because `-E` is parsed as an encoding option. Re-ran the process searches with `rg -l -e` and `-g '!*.ipynb'`. Targeted searches for `tcg`, `tqg`, `C_uG`, `cg -> t`, `t -> c g`, `h -> e tau`, `etau`, `Y_e_tau`, and `Y_tau_e` found no observable implementation in the required code directories. The cited nearby/generic code lines exist at `quarkConstraints/PAPER_0710_1869.md:35`, `tests/test_paper_couplings.py:258`, `quarkConstraints/modern/phenomenology.py:23`, `quarkConstraints/modern/phenomenology.py:165`, `flavorConstraints/muToEGamma.py:1`, `flavorConstraints/muToEGamma.py:75`, `scanParams/scan.py:523`, `scanParams/scan.py:524`, `yukawa/charged_lepton.py:1`, and `yukawa/charged_lepton.py:36`.
- CHK-6: `HIGH` is consistent for both processes. T005 requires a new top chromomagnetic FCNC convention, RS-to-SMEFT matching, decay/production normalization, collider/PDF recast choices, and QCD-running choices. T019 requires a new charged-lepton Higgs-Yukawa misalignment calculation, \(Y_{e\tau}\)/\(Y_{\tau e}\) conventions, and Higgs production/width interpretation.
- CHK-7: Process-specific searches outside `flavor_catalog/` found no rc1.1 paper numbers for T005 or T019 to contradict; these entries remain companion-mode catalog additions.
- CHK-8: No `CHECK`, `TODO`, `\ref{}`, `\cite{}`, `\label{}`, or `??` markers were found in the TeX files; manual inspection found no obvious math-mode or brace issues.

## Issues (if any)
- None.
