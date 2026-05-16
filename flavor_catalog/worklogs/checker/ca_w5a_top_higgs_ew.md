# CA Worklog: ca_w5a_top_higgs_ew
**Date**: 2026-05-16
**Family**: top_higgs_ew
**Process IDs**: T005, T015, T019

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| T005 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| T015 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| T019 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| T005 | ATLAS dataset for cg -> t | 139 fb^-1 at sqrt(s) = 13 TeV | PDG 2025 / ATLAS 2022 | b58775f55944 |
| T005 | sigma(cg -> t) B(t -> bW) B(W -> l nu) | < 4.7 pb | PDG 2025 / ATLAS 2022 | b58775f55944 |
| T005 | \|C_uG^{ct}\| / Lambda^2 | < 0.14 TeV^-2 | PDG 2025 / ATLAS 2022 | b58775f55944 |
| T005 | B(t -> c g), ATLAS EFT interpretation | < 3.7e-4 | PDG 2025 / ATLAS 2022 | b58775f55944 |
| T005 | CMS dataset | 5.0 fb^-1 at 7 TeV and 19.7 fb^-1 at 8 TeV | CMS 2017 | 75352c02b40e |
| T005 | \|kappa_tcg\| / Lambda | < 1.8e-2 TeV^-1 | CMS 2017 | 75352c02b40e |
| T005 | B(t -> c g), CMS interpretation | < 4.1e-4 | CMS 2017 | 75352c02b40e |
| T005 | SM B(t -> c g) benchmark | 4.6e-12 | Aguilar-Saavedra 2004 | 098fd502e902 |
| T005 | CFW ordinary/composite-pNGB RS KK-gluon context | about 21 TeV / about 33 TeV | Csaki-Falkowski-Weiler 2008 | 1ecc5891da03 |
| T015 | B(Z -> e mu), CMS observed expected | < 1.9e-7; expected < 2.0e-7 | CMS 2025 | 48418c00b815 |
| T015 | CMS dataset | 138 fb^-1 at sqrt(s) = 13 TeV | CMS 2025 | 48418c00b815 |
| T015 | B(Z -> e mu), PDG/ATLAS Run-2 | < 2.62e-7 | PDG 2025 / ATLAS 2023 | ac32ce7dd18e |
| T015 | ATLAS Run-2 dataset | 139 fb^-1 at sqrt(s) = 13 TeV | ATLAS 2023 | df940028eec5 |
| T015 | B(Z -> e mu), ATLAS Run-1 | < 7.5e-7 | ATLAS 2014 | 62c1fd32c929 |
| T015 | ATLAS Run-1 dataset | 20.3 fb^-1 at sqrt(s) = 8 TeV | ATLAS 2014 | 62c1fd32c929 |
| T015 | Tera-Z scale for LFV Z studies | O(1e12) Z decays | Calibbi-Marcano-Roy 2021 | e422d713e693 |
| T019 | B(h -> e tau), ATLAS Run-2 | < 0.20%; expected < 0.12% | PDG 2025 / ATLAS 2023 | 67ea71f88aaa |
| T019 | ATLAS Run-2 dataset | 138 fb^-1 at sqrt(s) = 13 TeV | ATLAS 2023 | a430814893db |
| T019 | B(h -> e tau), CMS Run-2 | < 0.22%; expected < 0.16% | PDG 2025 / CMS 2021 | 67ea71f88aaa |
| T019 | CMS Run-2 dataset | 137 fb^-1 at sqrt(s) = 13 TeV | CMS 2021 | 6b4f4c02719f |
| T019 | ATLAS strongest observed single-experiment limit | < 2.0e-3 | ATLAS 2023 | a430814893db |
| T019 | B(H -> e tau), ATLAS partial Run-2 | < 0.47%; expected < 0.34% (+0.13%/-0.10%) | ATLAS 2019 | 726a7510795c |
| T019 | ATLAS partial Run-2 dataset | 36.1 fb^-1 at sqrt(s) = 13 TeV | ATLAS 2019 | 726a7510795c |
| T019 | generic LFV-Higgs EFT allowance | order 10% | Harnik-Kopp-Zupan 2012 | 50662238edb6 |

## Evidence notes
- CHK-1: Grepped each TeX for numerical claims. T015 maps its physics numerical claims to `pdg_or_equivalent.values` with year, value, uncertainty/null uncertainty, source URL, access date, and sha256. T005's measured ATLAS/CMS values satisfy this, but the TeX also cites the SM `4.6e-12` benchmark and CFW `21 TeV`/`33 TeV` context outside `pdg_or_equivalent`. T019's ATLAS/CMS values satisfy this, but the TeX also cites the Harnik-Kopp-Zupan order-`10%` theory context outside `pdg_or_equivalent`.
- CHK-2: Source-key references in each TeX resolve to `flavor_catalog/references/<process_id>/source_manifest.yaml`. Value-id tags such as `CMS2025:T015:zemu_limit` resolve to YAML `pdg_or_equivalent.values` and were checked under CHK-1. Each manifest entry has a non-empty process-local tracked `.txt` snapshot path.
- CHK-3: `ls flavor_catalog/references/T005`, `T015`, and `T019` shows only `.txt` snapshots, optional `sha256sums.txt`, and `source_manifest.yaml`; `find ... -name '*.pdf'` returned no PDFs.
- CHK-4: Before CA edits, all three YAML sidecars contained `WRITER-INITIATED` followed by `WRITER-DONE`, both with ISO 8601 timestamps, and the latest pre-checker state was `WRITER-DONE`.
- CHK-5: The prompt's literal `rg -l -E` form errors with this installed ripgrep because `-E` is the encoding flag. Re-ran the same searches with `rg -l -e`. With `-g '!*.ipynb'`, process-specific searches for `tcg`, `Zemu`, and `etau` and broader process regexes returned no implementation hits. The cited generic/nearby code evidence exists at `quarkConstraints/PAPER_0710_1869.md:35`, `tests/test_paper_couplings.py:258`, `flavorConstraints/muToEGamma.py:1`, `flavorConstraints/muToEGamma.py:75`, `scanParams/scan.py:523`, `scanParams/scan.py:524`, `yukawa/charged_lepton.py:1`, and `yukawa/charged_lepton.py:36`.
- CHK-6: HIGH is consistent for all three processes. T005 needs new top chromomagnetic FCNC matching and collider reinterpretation; T015 needs LFV electroweak Z-coupling matching and a combination policy with low-energy muon constraints; T019 needs off-diagonal charged-lepton Higgs-Yukawa matching and Higgs-rate/width interpretation.
- CHK-7: No contradictions with rc1.1 load-bearing numbers were found; the batch is companion-mode and the relevant numerical values are process-local catalog additions or explicitly marked theory context.
- CHK-8: No `CHECK`, `TODO`, `\ref{}`, `\cite{}`, or `\label{}` markers were found in the three TeX files; no obvious math-mode or brace issues were found. `chktex`/`lacheck` were not used.

## Issues (if any)
- T005: CHK-1 fails under the strict CA placement rule. Move or duplicate the SM `B(t -> c g)=4.6e-12` benchmark and CFW `21 TeV`/`33 TeV` context into `pdg_or_equivalent.values`, with explicit `year`, `value`, `uncertainty`, `source_url`, `access_date`, and `sha256`, or remove those numerical claims from the TeX.
- T019: CHK-1 fails under the strict CA placement rule. Move or duplicate the Harnik-Kopp-Zupan order-`10%` LFV-Higgs EFT context into `pdg_or_equivalent.values`, with explicit `year`, `value`, `uncertainty`, `source_url`, `access_date`, and `sha256`, or remove that numerical claim from the TeX.
