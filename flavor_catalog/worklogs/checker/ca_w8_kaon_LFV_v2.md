# CA Worklog: ca_w8_kaon_LFV_v2
**Date**: 2026-05-17
**Family**: kaon, Wave-8 SECONDARY
**Cycle**: 2
**Process IDs**: K020
**Agent**: CA-v2-K020

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | CHK-W8 | overall |
|---|---|---|---|---|---|---|---|---|---|---|
| K020 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| K020 | BR(K+ -> pi+ mu+ e-), Sher/BNL E865 final limit newly promoted in cycle 2 | < 2.1e-11 at 90% CL | Sher et al. / BNL E865, arXiv:hep-ex/0502020 | 0ae93baf10b7 |
| K020 | BR(K+ -> pi+ mu+ e-), PDG combined E865/E777 limit | < 1.3e-11 at 90% CL | PDG Live / PDG API listing S010.29, 2025 edition | 8ddac950390a |
| K020 | BR(K+ -> pi+ mu- e+), PDG current limit | < 6.6e-11 at 90% CL | PDG Live / PDG API listing S010.25 / NA62 2021 | 8ddac950390a |
| K020 | BR(K+ -> pi+ mu- e+), NA62 direct abstract spot-check | < 6.6e-11 at 90% CL | NA62 2021, arXiv:2105.06759 | b1bae8275544 |

## Issues
- K020: None.

## Evidence notes
- CHK-1: PASS. The cycle-1 failure is fixed: `pdg_or_equivalent.values` now contains `Sher2005:K020:Kplus_piplus_mup_em_e865_only_limit` with value `2.1e-11`, `uncertainty: null`, units, upper-limit type, `90% CL`, year 2005, source URL, access date, snapshot path, and sha256. The TeX PDG `1.3e-11` and `6.6e-11` current limits also map to structured value blocks. The L001 / B001_B003 carve-out was applied in the intended direction: this Sher/E865 number is a measured observable, so promotion is appropriate and no remaining theory/dataset-context number is being misclassified.
- CHK-2: PASS. The TeX reference keys `PDG2025_Kplus_LFV_semileptonic`, `Sher2005_Kplus_piplus_mup_em`, `NA622021_Kplus_piplus_mum_ep`, `CFW2008_RS_flavor`, `BenekeMochRohrwild2015_RS_LFV`, and `AngelescuFaroughySumensari2020_semileptonic_LFV` all resolve in `flavor_catalog/references/K020/source_manifest.yaml`. `git ls-files --error-unmatch` confirmed all manifest snapshot paths are tracked text files.
- CHK-3: PASS. `find flavor_catalog/references/K020` found no publisher PDF or office-document files; the directory contains text snapshots, `source_manifest.yaml`, and `sha256sums.txt`.
- CHK-4: PASS. `K020.yaml` has the legal `DRAFT -> WRITER-INITIATED -> WRITER-DONE -> WRITER-REWORK -> WRITER-DONE` cycle-2 path with ISO 8601 `timestamp`/`at` fields. The pre-check state was `WRITER-DONE`.
- CHK-5: PASS. Exact K020/S010/charged-kaon LFV searches over `quarkConstraints`, `qcd`, `flavorConstraints`, `neutrinos`, `yukawa`, `warpConfig`, `solvers`, `scanParams`, and `tests` found only generic LFV metadata and the existing `mu -> e gamma` dipole path. The cited nearby implementation lines exist: `flavorConstraints/muToEGamma.py:75`, `scanParams/scan.py:524`, `quarkConstraints/deltaf2.py:755`, and `quarkConstraints/modern/phenomenology.py:23`.
- CHK-6: PASS. `HIGH` remains consistent with the rubric because K020 needs a new Delta S = 1 semileptonic LFV operator basis, K -> pi e mu rate/form-factor treatment, and cross-sector quark-lepton RS matching, not just a new threshold on an existing observable.
- CHK-7: PASS. Searches in `flavor_catalog` and `docs` found only the plan/deferred-scope K020 row, Wave-8 worklogs, and general lepton-sector deferral notes; no upstream or methodology load-bearing number contradicts the K020 limits, coverage classification, or difficulty rating.
- CHK-8: PASS. `rg` found no unresolved `\cite`, `\ref`, TODO/FIXME/MISSING/CHECK, or `??` markers in `K020.tex`. Delimiter counts were balanced: `\(`/`\)` = 40/40, `\[`/`\]` = 2/2, dollar count = 0, and braces = 50/50. `git diff --check` was clean before the checker status edit.
- CHK-W8: PASS. `K020.tex` and `K020.yaml` live under `flavor_catalog/processes/secondary/kaon/`; the sidecar has `priority_tier: SECONDARY`, a non-empty `priority_rationale`, and `promoted_in_wave: 8`; the TeX prose explicitly identifies the entry as SECONDARY Wave-8 and points to `flavor_catalog/PRIORITY_TIERS.md`.
