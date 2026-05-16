# CA Worklog: ca_w6_charged_lepton_top_v2
**Date**: 2026-05-16
**Family**: charged_lepton_top
**Process IDs**: L023 T020

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| L023 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| T020 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| L023 | CHARM-II sigma_exp/sigma_SM | 1.58 +/- 0.64 | DUNE trident compilation, arXiv:1902.06765 | 19acfb9bc4c3 |
| L023 | CCFR sigma_exp/sigma_SM | 0.82 +/- 0.28 | CCFR PRL/PubMed snapshot plus DUNE compilation | 903c8a531035 |
| L023 | CCFR corrected trident events vs SM | 37.0 +/- 12.4 vs 45.3 +/- 2.3 | CCFR PRL/PubMed snapshot | 903c8a531035 |
| L023 | NuTeV sigma_exp/sigma_SM | 0.72 +1.73 -0.72 | DUNE compilation and NuTeV arXiv snapshot | 19acfb9bc4c3 |
| L023 | Altmannshofer 2014 CCFR contours | 95% CL | arXiv:1406.2332 snapshot | b360bc6bd857 |
| L023 | Belle-II auxiliary reach inputs | sqrt(s) = 10.58 GeV; 50 ab^-1 | Kaneta-Shimomura arXiv:1701.00156 | 9f1614fe9aa0 |
| L023 | DUNE projection wording | about 25% after about 3 years per beam mode; 40% baseline and 25% improved contour | DUNE trident compilation, arXiv:1902.06765 | 19acfb9bc4c3 |
| T020 | Higgs mass hypothesis in direct searches | 125 GeV | PDG 2025 Higgs LFV review | 8178cb0c1055 |
| T020 | PDG/CMS B(h -> e mu) | <4.4 x 10^-5, expected <4.7 x 10^-5, 95% CL | PDG 2025 Higgs LFV review | 8178cb0c1055 |
| T020 | PDG/ATLAS B(h -> e mu) | <6.2 x 10^-5, expected <5.9 x 10^-5, 95% CL | PDG 2025 Higgs LFV review | 8178cb0c1055 |
| T020 | CMS direct search dataset and scan | 138 fb^-1 at 13 TeV; 110--160 GeV scan; largest excess near 146 GeV with 3.8/2.8 sigma | CMS HIG-22-002 public page | 072412def18c |
| T020 | ATLAS primary B(H -> e mu) | <6.1 x 10^-5, expected <5.8 x 10^-5; 139 fb^-1 at 13 TeV | ATLAS arXiv:1909.10235 snapshot | e1e4c131be4c |

## Issues (if any)
- None.

## Checklist evidence
- CHK-1: Numeric greps of both TeX files were compared to `pdg_or_equivalent`. WA-v2 added the prior missing L023 auxiliary values and T020 125 GeV mass-hypothesis value under `pdg_or_equivalent.values[]` with year, value, source URL, access date, snapshot path, and sha256. The original measured-value rows also retain source URL, access date, snapshot path, and sha256 metadata.
- CHK-2: All bare process-local source keys cited in `L023.tex` and `T020.tex` resolve to the corresponding `source_manifest.yaml`. Every manifest entry has a non-empty `snapshot_path`, and `git ls-files` confirms the paths are tracked.
- CHK-3: `ls flavor_catalog/references/L023/` and `ls flavor_catalog/references/T020/` showed text snapshots plus manifest/checksum files only. `find flavor_catalog/references/L023 flavor_catalog/references/T020 -type f -name '*.pdf'` returned no PDFs.
- CHK-4: Before this CA update, both sidecars had ordered ISO-8601 `WRITER-INITIATED` and terminal `WRITER-DONE` entries. This work appended `CHECKER-DONE`, set `checker_passed_at`, updated `last_updated_at`, and recorded `checker_agent_id: "CA"`.
- CHK-5: The literal requested `rg -l -E "<process keyword>" ...` form errors in this installed ripgrep because `-E` is parsed as an encoding flag. The equivalent `rg -l -e` searches over the required directories found no L023 trident implementation and no T020 h -> e mu or `Y_e_mu`/`Y_mu_e` implementation. Broad L023 `CHARM` hits are only QCD charm-mass constants/tests. T020's nearby cited code lines exist at `flavorConstraints/muToEGamma.py:1`, `flavorConstraints/muToEGamma.py:75`, `scanParams/scan.py:523`, `scanParams/scan.py:524`, `yukawa/charged_lepton.py:1`, and `yukawa/charged_lepton.py:36`.
- CHK-6: HIGH is consistent for both entries. L023 needs a new neutrino-trident mode calculation with lepton-current operators, target/form-factor handling, flux integration, and experiment-specific rates. T020 needs new charged-lepton Higgs-Yukawa misalignment matching plus Higgs production/width interpretation.
- CHK-7: Searches outside `flavor_catalog/` found only plan rows for L023 and T020, not rc1.1 load-bearing values to contradict.
- CHK-8: `rg` found no `CHECK`, `TODO`, `FIXME`, `??`, `\cite`, or `\ref` markers in either TeX file. `chktex` and `lacheck` are not installed; visual inspection found no obvious math-mode or bracing defects.
