# CA Worklog: ca_w7_new_v2
**Date**: 2026-05-16
**Family**: mixed
**Process IDs**: T003 T004 B012

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| T003 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| T004 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| B012 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| T003 | CMS dataset for tqgamma search | 138 fb^-1 at sqrt(s) = 13 TeV | CMS-TOP-21-013 / PRD 109 (2024) 072004 | 5f094c0620a2 |
| T003 | B(t -> c gamma), CMS | <1.51e-5 observed, <1.54e-5 expected at 95% CL | CMS-TOP-21-013 / PRD 109 (2024) 072004 | 5f094c0620a2 |
| T003 | ATLAS dataset for tqgamma search | 139 fb^-1 at sqrt(s) = 13 TeV | ATLAS arXiv:2205.02537 | 82494cd375fd |
| T003 | B(t -> c gamma), ATLAS LH/RH | <4.2e-5 LH, <4.5e-5 RH at 95% CL | ATLAS arXiv:2205.02537 | 82494cd375fd |
| T003 | Gamma(t -> gamma q)/Gamma(t -> Wb), q = u,c | <9.5e-6 at 95% CL | PDG 2026 pdgLive top summary | e05906a43ebc |
| T003 | SM B(t -> c gamma) context | 4.6e-14 central value, TeX now says 10^-14 level | Aguilar-Saavedra 2004, stored outside pdg_or_equivalent | 01040933ce18 |
| T004 | B(t -> gamma q), q = u,c | <9.5e-6 at 95% CL | PDG 2026 pdgLive t-quark listing | 5e72da33cd1b |
| T004 | B(t -> u gamma), ATLAS LH/RH | <0.85e-5 LH, <1.2e-5 RH at 95% CL | ATLAS arXiv:2205.02537 | 1b92967f8b52 |
| T004 | B(t -> u gamma), CMS | <0.95e-5 observed, <1.20e-5 expected at 95% CL | CMS-TOP-21-013 | 1c729a489ac5 |
| B012 | B(B0 -> K*0 gamma) | (41.63 +/- 0.92)e-6 | HFLAV rare decays, end-of-year 2024 | dc826a3f6d24 |
| B012 | B(B+ -> K*+ gamma) | (39.5 +/- 1.0)e-6 | HFLAV rare decays, end-of-year 2024 | dc826a3f6d24 |
| B012 | Delta_0+(B -> K* gamma) | 0.059 +/- 0.014 | HFLAV rare decays, end-of-year 2024 | dc826a3f6d24 |
| B012 | A_CP(B -> K* gamma) | -0.004 +/- 0.011 | HFLAV rare decays, end-of-year 2024 | dc826a3f6d24 |
| B012 | S_{K*gamma} | -0.08 +/- 0.17 | PDG pdgLive 2025 | 3ba65881b19e |
| B012 | C_{K*gamma} | 0.03 +/- 0.10 | PDG pdgLive 2025 | 3ba65881b19e |

## Issues (if any)
- T003: CHK-1 fail. `T003.tex` still states that the Standard Model expectation remains at the `10^-14` level. The backing source/value is present only under `paper_era_reference` (`AguilarSaavedra2004:T003:SM`, value `4.6e-14`), not under `pdg_or_equivalent`, so the TeX numerical physics claim does not meet the checklist rule.
- T004: no issues.
- B012: no issues.

## Evidence notes
- CHK-1: I grepped each TeX for numeric observable claims and checked the sidecar `pdg_or_equivalent` blocks for year, value, uncertainty or CL, source URL, access date, and sha256. T004 and B012 now satisfy the v2 re-check. T003 still has the SM `10^-14`-level claim outside `pdg_or_equivalent`.
- CHK-2: TeX key references resolve to `flavor_catalog/references/<process_id>/source_manifest.yaml`. Every manifest entry has a non-empty `snapshot_path` under the process reference directory, and `git ls-files flavor_catalog/references/{T003,T004,B012}` shows the snapshots are tracked.
- CHK-3: `ls flavor_catalog/references/{T003,T004,B012}` shows only `.txt` snapshots, manifests, and B012 `sha256sums.txt`; no `.pdf` files are present.
- CHK-4: Before checker updates, all three YAML sidecars had `WRITER-INITIATED` followed by `WRITER-DONE` in order with ISO 8601 timestamps, and the most recent entry was the WA-v2 `WRITER-DONE`.
- CHK-5: In this environment, `rg -l -E "<pattern>"` fails because ripgrep 15.1.0 uses `-E` for encoding, so I reran the checks with `rg -l -e "<pattern>"` across the required directories. T003 and T004 returned no source-code hits. B012 returned notebook false positives only for a broad `C7` pattern; excluding notebooks and using the process-specific B -> K* gamma/photon-polarization pattern returned no implementation hits. The cited code lines exist at `quarkConstraints/deltaf2.py:209`, `quarkConstraints/modern/phenomenology.py:23`, and `flavorConstraints/muToEGamma.py:75`.
- CHK-6: The `HIGH` rubric is consistent for all three entries. T003 and T004 need new top-FCNC photon/dipole conventions and width or production recast machinery; B012 needs b -> s gamma dipole matching/RG plus exclusive form-factor and helicity/time-dependent CP treatment. None uses the existing Delta F = 2 SLL/SLR/VLL/VRR/LR1/LR2 basis.
- CHK-7: No TeX entry revises rc1.1 scan numbers; references to existing repo code are only code-coverage context.
- CHK-8: Manual and grep-based checks found no `CHECK`, `TODO`, `??`, `\ref{...}`, or `\cite{...}` markers in the three TeX files.
