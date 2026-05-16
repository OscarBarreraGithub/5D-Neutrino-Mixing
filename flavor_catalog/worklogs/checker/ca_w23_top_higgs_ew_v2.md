# CA Worklog: ca_w23_top_higgs_ew_v2
**Date**: 2026-05-16
**Family**: top_higgs_ew
**Process IDs**: T002, T018

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| T002 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| T018 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| T002 | B(t -> Z u), left-handed t-u-Z benchmark | < 6.2e-5 | PDG 2025 top listing / ATLAS 2023 | 94cf070f2608 |
| T002 | B(t -> Z u), right-handed t-u-Z benchmark | < 6.6e-5 | PDG 2025 top listing / ATLAS 2023 | 94cf070f2608 |
| T002 | ATLAS dataset | 139 fb^-1 at sqrt(s) = 13 TeV | ATLAS 2023 | 4f5e39fa0d55 |
| T002 | B(t -> Z u), CMS comparison limit | < 0.022% | CMS 2017 | 7c60c3092af2 |
| T002 | CMS dataset | 19.7 fb^-1 at sqrt(s) = 8 TeV | CMS 2017 | 7c60c3092af2 |
| T002 | Standard Model B(t -> u Z) estimate | 8e-17 | Aguilar-Saavedra 2004 | b50f50126f63 |
| T002 | Minimal-RS t -> u Z suppression relative to t -> c Z | approximately two orders of magnitude below t -> c Z | Casagrande et al. 2008 | 22d2db0df79d |
| T018 | B(h -> mu tau), CMS Run-2 direct search | < 0.15%, expected < 0.15% | PDG 2025 Higgs LFV review / CMS 2021 | 789633a52f05 |
| T018 | CMS dataset | 137 fb^-1 at sqrt(s) = 13 TeV | CMS 2021 | 1f45cde9d68a |
| T018 | B(h -> mu tau), ATLAS Run-2 direct search | < 0.18%, expected < 0.09% | PDG 2025 Higgs LFV review / ATLAS 2023 | 789633a52f05 |
| T018 | ATLAS dataset | 138 fb^-1 at sqrt(s) = 13 TeV | ATLAS 2023 | 9e03404c2aa3 |
| T018 | CMS 8 TeV local excess and best fit | 2.4 sigma, p = 0.010, B(H -> mu tau) = (0.84 +0.39 -0.37)% | CMS 2015 | 3e0638a894a0 |
| T018 | CMS 8 TeV limit and dataset | < 1.51%; 19.7 fb^-1 at sqrt(s) = 8 TeV | CMS 2015 | 3e0638a894a0 |
| T018 | Generic LFV-Higgs EFT allowance | order 10% | Harnik, Kopp, Zupan 2012 | 3cdc2b58e96 |

## Evidence notes
- CHK-1: Grepped both TeX files for numeric claims. T002's `pdg_or_equivalent.values` now contains the PDG/ATLAS 6.2e-5 and 6.6e-5 limits, CMS 0.022% comparison, 19.7 fb^-1 and 139 fb^-1 datasets, 8e-17 SM estimate, and RS two-orders suppression statement, each with year, value, uncertainty, source URL, access date, and sha256. T018's block contains the CMS/ATLAS 0.15%, 0.18%, and 0.09% limits, 137/138/19.7 fb^-1 datasets, CMS 2015 2.4 sigma, p = 0.010, best-fit and 1.51% claims, plus the Harnik-Kopp-Zupan order-10% context with the same metadata fields.
- CHK-2: Every source key in each Key references section resolves to `flavor_catalog/references/<process_id>/source_manifest.yaml`; each manifest entry has a non-empty `snapshot_path` under the same process reference directory. `git ls-files flavor_catalog/references/T002 flavor_catalog/references/T018` confirms all snapshots and manifests are tracked.
- CHK-3: `ls -la flavor_catalog/references/T002` and `ls -la flavor_catalog/references/T018` showed only `.txt` snapshots plus `source_manifest.yaml`; `find ... -name '*.pdf'` returned no files.
- CHK-4: Both YAML sidecars contain `WRITER-INITIATED` followed by `WRITER-DONE`, `WRITER-REWORK`, and the cycle-2 `WRITER-DONE`, all with ISO 8601 timestamps. Before this CA update, the most recent entry was the WA-v2 `WRITER-DONE` at 2026-05-16T12:43:15-04:00.
- CHK-5: The prompt's literal `rg -l -E "<process keyword>" ... 2>&1` form was run for both processes and failed because this ripgrep interprets `-E` as an encoding flag. Re-running the same searches with `rg -n -e` found only generic FCNC comments for T002 at `quarkConstraints/PAPER_0710_1869.md:35` and `tests/test_paper_couplings.py:258`, and no process-specific T018 implementation hits. T018's cited existing LFV support lines exist at `flavorConstraints/muToEGamma.py:1`, `flavorConstraints/muToEGamma.py:75`, and `scanParams/scan.py:523`; `yukawa/charged_lepton.py:1` and `:36` are diagonal helper code, not an h -> mu tau observable.
- CHK-6: HIGH is consistent for both entries. Neither uses the existing Delta F=2 SLL/SLR/VLL/VRR/LR1/LR2 operator basis; T002 needs new top electroweak FCNC matching/decay-width interpretation, and T018 needs a new charged-lepton Higgs-Yukawa misalignment and Higgs-rate/width interpretation.
- CHK-7: Repo-wide grep outside `flavor_catalog/` found only plan/phase-log mentions of T002/T018 and no paper-number reuse for these processes, so no rc1.1 load-bearing number is contradicted.
- CHK-8: No `CHECK`, `TODO`, unresolved `\ref{}`, `\cite{}`, or `\label{}` markers were found. Brace counts are balanced for both TeX files. `chktex`, `lacheck`, and `latexmk` are not installed in this environment.

## Issues (if any)
- None.
