# CA Worklog: ca_w23_beauty_v3
**Date**: 2026-05-16
**Family**: beauty
**Process IDs**: B032

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| B032 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| B032 | BR(B+ -> K0 pi+) | (2.392 +/- 0.062)e-5 | HFLAV Dec. 2025 | 522cee5fd679 |
| B032 | BR(B+ -> K+ pi0) | (1.322 +/- 0.043)e-5 | HFLAV Dec. 2025 | 522cee5fd679 |
| B032 | BR(B0 -> K+ pi-) | (2.007 +/- 0.040)e-5 | HFLAV Dec. 2025 | 522cee5fd679 |
| B032 | BR(B0 -> K0 pi0) | (1.012 +/- 0.043)e-5 | HFLAV Dec. 2025 | 522cee5fd679 |
| B032 | A_CP(B+ -> K0_S pi+) | -2.67 +/- 0.87 percent | HFLAV Dec. 2025 | 522cee5fd679 |
| B032 | A_CP(B+ -> K+ pi0) | +2.7 +/- 1.2 percent | HFLAV Dec. 2025 | 522cee5fd679 |
| B032 | A_CP(B0 -> K+ pi-) | -8.31 +/- 0.31 percent | HFLAV Dec. 2025 | 522cee5fd679 |
| B032 | A_CP(B0 -> K0 pi0) | -1 +/- 13 percent | HFLAV Dec. 2025 | 522cee5fd679 |
| B032 | C(K0 pi0) in B0 -> K0 pi0 | 0.00 +/- 0.08 | PDG 2025 API | dcd2aea4b34 |
| B032 | S(K0 pi0) in B0 -> K0 pi0 | 0.64 +/- 0.13 | PDG 2025 API | dcd2aea4b34 |
| B032 | LHCb A_CP(B+ -> K+ pi0) | 0.025 +/- 0.015 stat +/- 0.006 syst +/- 0.003 external | LHCb 2021 | 3cf10c6a5447 |
| B032 | Belle II Kpi sum-rule test | -0.03 +/- 0.13 stat +/- 0.04 syst | Belle II 2024 | 921fc4417a5f |

## Issues (if any)
- None.

## Evidence notes
- CHK-1: `rg -n "[0-9]" flavor_catalog/processes/beauty/B032.tex` found the load-bearing numerical observable claims at lines 28-31, 36-41, 46, 62, and 70. Each maps to a `pdg_or_equivalent` value item in `B032.yaml` with `year`, `value`, `uncertainty`, `source_url`, `access_date`, `snapshot_path`, and `sha256`. The HFLAV, PDG, LHCb, and Belle II source snapshots contain the same values, and `sha256sum flavor_catalog/references/B032/*.txt` matches `source_shas` and the per-value sha256 metadata.
- CHK-2: TeX reference keys `hflav_dec2025_btopik`, `pdg2025_k0pi0_time_dependent_cp_api`, `belleii2024_btokpi`, `lhcb2021_bplus_kplus_pi0_acp`, `belleii2023_b0_ks_pi0_cp`, and `cfw2008_rs_flavor` all resolve in `flavor_catalog/references/B032/source_manifest.yaml`. Each manifest entry has a non-empty tracked `.txt` `snapshot_path` under `flavor_catalog/references/B032/`.
- CHK-3: `ls -la flavor_catalog/references/B032` shows only `.txt` snapshots plus `source_manifest.yaml`; `find flavor_catalog/references/B032 -maxdepth 1 -type f -iname '*.pdf' -print` returned no PDFs.
- CHK-4: Before this CA update, `status_history` contained `WRITER-INITIATED` at `2026-05-16T11:58:56-04:00` followed by `WRITER-DONE` transitions, with the most recent pre-check entry `WRITER-DONE` at `2026-05-16T13:11:30-04:00`. Timestamps are ISO 8601. This worklog appends `CHECKER-DONE` for cycle 3.
- CHK-5: The literal requested `rg -l -E "<process keyword>" ...` form fails in this installed ripgrep because `-E` is parsed as an encoding flag. The equivalent `rg -l -e "B.?->.?K.?pi|B.?to.?K.?pi|B.?->.?pi.?K|B.?to.?pi.?K|BtoKpi|Kpi_puzzle|Kpi|K_pi|pi_K|charmless|nonleptonic|electroweak.?penguin" ... -g '!*.ipynb'` returned no files. Without the notebook exclusion, the only hit is a base64 PNG blob in `quarkConstraints/quark_benchmark_validation.ipynb`, not an implementation. The cited lines `quarkConstraints/modern/phenomenology.py:23`, `quarkConstraints/modern/phenomenology.py:165`, `quarkConstraints/deltaf2.py:209`, and `flavorConstraints/muToEGamma.py:75` exist and support the `NO` coverage claim.
- CHK-6: `HIGH` is consistent with the rubric: this process needs a new nonleptonic `b -> s qbar q` amplitude/likelihood layer with electroweak-penguin matching and hadronic strong phases, not only the existing Delta-F=2 SLL/SLR/VLL/VRR/LR1/LR2 operator basis.
- CHK-7: Targeted greps for B -> K pi / pi K / charmless nonleptonic language in the rc1.1 documentation surface outside the catalog found no load-bearing rc1.1 B032 values to contradict. This entry is companion-mode and new to the catalog.
- CHK-8: `B032.tex` has no `\textbf{CHECK}`, `\ref{}`, `\cite{}`, or `undefined` markers. Cosmetic LaTeX scan found no missing math-mode delimiters or brace issues in the touched process entry.
