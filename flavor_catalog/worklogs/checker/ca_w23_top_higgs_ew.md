# CA Worklog: ca_w23_top_higgs_ew
**Date**: 2026-05-16
**Family**: top_higgs_ew
**Process IDs**: T002, T007, T018

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| T002 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| T007 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| T018 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| T002 | B(t -> Z u), left-handed tuZ | < 6.2e-5 | PDG 2025 top listing / ATLAS 2023 | 94cf070f2608 |
| T002 | B(t -> Z u), right-handed tuZ | < 6.6e-5 | PDG 2025 top listing / ATLAS 2023 | 94cf070f2608 |
| T002 | B(t -> Z u), CMS comparison | < 0.022% | CMS 2017 | 7c60c3092af2 |
| T002 | CMS dataset | 19.7 fb^-1 at sqrt(s) = 8 TeV | CMS 2017 | 7c60c3092af2 |
| T002 | ATLAS dataset | 139 fb^-1 at sqrt(s) = 13 TeV | ATLAS 2023 | 4f5e39fa0d55 |
| T002 | SM B(t -> u Z) estimate | 8e-17 | Aguilar-Saavedra 2004 | b50f50126f63 |
| T002 | minimal-RS t -> u Z scaling | about two orders below t -> c Z | Casagrande et al. 2008 | 22d2db0df79d |
| T007 | B(t -> H c), PDG headline | < 3.4e-4 | PDG 2026 pdgLive top datablock | 6b6d2c99c13e |
| T007 | ATLAS dataset | 140 fb^-1 at sqrt(s) = 13 TeV | ATLAS 2024 | b0373e82980a |
| T007 | B(t -> H c), ATLAS multilepton-only | < 3.3e-4, expected < 3.8e-4 | ATLAS 2024 | b0373e82980a |
| T007 | B(t -> H c), ATLAS combination | < 3.4e-4, expected < 2.3e-4 | ATLAS 2024 | b0373e82980a |
| T007 | CMS dataset | 138 fb^-1 at sqrt(s) = 13 TeV, 2016-2018 | CMS 2025 | 4d2d78674a2e |
| T007 | B(t -> H c), CMS combination | < 0.037%, expected < 0.035% | CMS 2025 | 4d2d78674a2e |
| T018 | B(h -> mu tau), CMS Run-2 | < 0.15%, expected < 0.15% | PDG 2025 Higgs LFV review / CMS 2021 | 789633a52f05 |
| T018 | CMS dataset | 137 fb^-1 at sqrt(s) = 13 TeV | CMS 2021 | 1f45cde9d68a |
| T018 | B(h -> mu tau), ATLAS Run-2 | < 0.18%, expected < 0.09% | PDG 2025 Higgs LFV review / ATLAS 2023 | 789633a52f05 |
| T018 | ATLAS dataset | 138 fb^-1 at sqrt(s) = 13 TeV | ATLAS 2023 | 9e03404c2aa3 |
| T018 | CMS 8 TeV hint | 2.4 sigma, p = 0.010, best fit (0.84 +0.39 -0.37)% | CMS 2015 | 3e0638a894a0 |
| T018 | CMS 8 TeV limit and dataset | < 1.51%; 19.7 fb^-1 at sqrt(s) = 8 TeV | CMS 2015 | 3e0638a894a0 |
| T018 | generic LFV Higgs EFT allowance | order 10% | Harnik, Kopp, Zupan 2012 | 3cdc2b58e96d |

## Evidence notes
- CHK-1: Grepped the TeX for numerical claims. T007 maps its numerical claims to `pdg_or_equivalent.values`, with year, value, uncertainty, source URL, access date, and sha256. T002 maps several numerical claims only to `auxiliary_values` (`CMS2017:T002:*`, `ATLAS2023:T002:*`, `AguilarSaavedra2004:T002:SM`, `Casagrande2008:T002:RS_suppression`), and those blocks lack explicit `year` fields. T018's `pdg_or_equivalent.values` omit an explicit `uncertainty` field, and several TeX numerical claims map only to `auxiliary_values` without explicit year/uncertainty fields.
- CHK-2: Every source key in the Key references sections resolves to `flavor_catalog/references/<process_id>/source_manifest.yaml`; each manifest entry has a non-empty `snapshot_path` under the process reference directory. `git ls-files` confirms the snapshot text files are tracked.
- CHK-3: `ls flavor_catalog/references/T002`, `T007`, and `T018` shows only `.txt` snapshots plus `source_manifest.yaml`; `find ... -name '*.pdf'` returned no publisher PDFs.
- CHK-4: All three YAML sidecars contain `WRITER-INITIATED` followed by `WRITER-DONE`, with ISO 8601 timestamps; before CA edits, `WRITER-DONE` was the latest status.
- CHK-5: The prompt's literal `rg -l -E` form errors with this installed ripgrep because `-E` is the encoding flag. Re-ran the same searches with `rg -l -e`. With `-g '!*.ipynb'`, process-specific searches for T002, T007, and T018 returned no implementation hits; generic FCNC terms only hit `quarkConstraints/PAPER_0710_1869.md:35` and `tests/test_paper_couplings.py:258`. T018's cited LFV support lines exist at `flavorConstraints/muToEGamma.py:1`, `flavorConstraints/muToEGamma.py:75`, `scanParams/scan.py:523`, `yukawa/charged_lepton.py:1`, and `yukawa/charged_lepton.py:36`.
- CHK-6: HIGH is consistent for all three processes because each needs a new top/Higgs/LFV matching or decay-mode calculation, not the existing Delta F=2 operator basis.
- CHK-7: Repo-wide grep outside `flavor_catalog/` found only catalog-plan mentions for these processes, so no load-bearing rc1.1 paper numbers are silently revised.
- CHK-8: No `CHECK`, `TODO`, `??`, `\ref{}`, or `\cite{}` markers were found; brace counts are balanced for all three TeX files. `chktex` and `lacheck` are not installed.

## Issues (if any)
- T002: CHK-1 fails. Move every TeX numerical claim into `pdg_or_equivalent.values` or extend the accepted numeric-value schema so each auxiliary block carries explicit `year`, `value`, `uncertainty`, `source_url`, `access_date`, and `sha256`.
- T018: CHK-1 fails. Add explicit `uncertainty` fields to the headline `pdg_or_equivalent.values`, and move or normalize auxiliary numerical claims so each TeX number is represented in a block with `year`, `value`, `uncertainty`, `source_url`, `access_date`, and `sha256`.
