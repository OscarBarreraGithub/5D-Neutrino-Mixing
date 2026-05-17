# CA Worklog: ca_w8_B_radiative
**Date**: 2026-05-17
**Family**: beauty
**Batch ID**: Wave-8 B-radiative
**Cycle**: 1
**Checker agent**: CA-CA-w8-B-radiative
**Process IDs**: B013 B014

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | CHK-W8 | overall |
|---|---|---|---|---|---|---|---|---|---|---|
| B013 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| B014 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| B013 | B(B_s0 -> phi gamma), PDG average | (3.4 +/- 0.4) x 10^-5 | PDG 2025 B_s0 listing | 59578d332213 |
| B013 | B(B_s0 -> phi gamma), HFLAV average | (34.0 +/- 3.2) x 10^-6 | HFLAV rare decays Dec. 2024 | 0572d6a8d33f |
| B013 | S_phi gamma | 0.43 +/- 0.30(stat) +/- 0.11(syst) | LHCb 2019, arXiv:1905.06284 | a8f9cea7df59 |
| B013 | C_phi gamma | 0.11 +/- 0.29(stat) +/- 0.11(syst) | LHCb 2019, arXiv:1905.06284 | a8f9cea7df59 |
| B013 | A_Delta_phi gamma | -0.67 +0.37 -0.41(stat) +/- 0.17(syst) | LHCb 2019, arXiv:1905.06284 | a8f9cea7df59 |
| B014 | B(B+ -> rho+ gamma) | (9.8 +/- 2.5) x 10^-7 | PDG 2025 B+/- listing | 8041b539fb0b |
| B014 | B(B0 -> rho0 gamma) | (8.6 +/- 1.5) x 10^-7 | PDG 2025 B0 listing | 9a004583a300 |
| B014 | B(B0 -> omega gamma) | (4.4 +1.8 -1.6) x 10^-7 | PDG 2025 B0 listing | 9a004583a300 |
| B014 | B(B -> rho gamma) | (1.40 +/- 0.22) x 10^-6 | HFLAV rare decays Dec. 2024 | e3d392743924 |
| B014 | B(B -> rho/omega gamma) | (1.30 +/- 0.18) x 10^-6 | HFLAV rare decays Dec. 2024 | e3d392743924 |
| B014 | Gamma(B+ -> rho+ gamma)/(2 Gamma(B0 -> rho0 gamma)) - 1 | -0.46 +/- 0.17 | HFLAV rare decays Dec. 2024 | e3d392743924 |
| B014 | B(B -> Xd gamma), related normalization | (9.2 +/- 3.0) x 10^-6 | HFLAV rare decays Dec. 2024 | e3d392743924 |

## Issues
- B013: None.
- B014: None.

## Evidence notes
- CHK-1: B013's load-bearing measured observables in the TeX are the PDG/HFLAV branching fractions and LHCb 2019 S, C, and A_Delta observables; each appears in `pdg_or_equivalent` with year, value, uncertainty fields, units, source URL, access date, snapshot path, and sha256. B014's TeX measured observables are the PDG component branching fractions, HFLAV combined branching fractions, HFLAV naive isospin average, and related HFLAV B -> Xd gamma average; each likewise has complete `pdg_or_equivalent` metadata. L001/B001/B021 carve-out applied: dataset descriptors, theory context, arXiv IDs, code line numbers, and source years were not treated as measured-observable CHK-1 requirements.
- CHK-2: Manifest/key validation found 8 B013 source keys and 12 B014 source keys. Every process-local key listed in each TeX resolves to the corresponding `source_manifest.yaml`, and every manifest entry has a non-empty `snapshot_path` tracked by `git ls-files`.
- CHK-3: `find flavor_catalog/references/B013 flavor_catalog/references/B014 -type f -iname '*.pdf'` returned no files. Reference directories contain text snapshots, manifests, and sha256 lists only.
- CHK-4: Both sidecars show `DRAFT -> WRITER-INITIATED -> WRITER-DONE` with ISO 8601 `timestamp`/`at` fields; the latest pre-check state was `WRITER-DONE` cycle 1.
- CHK-5: B013 and B014 both claim `code_coverage.status: NO`. Focused greps for B013-specific `B_s -> phi gamma`, S_phi, A_Delta, photon-polarization, and b -> s gamma dipole patterns returned no implementation hits. Focused greps for B014-specific `B -> rho gamma`, `B -> omega gamma`, b -> d gamma, Xd gamma, and B014 returned no matches. The cited evidence lines exist: `quarkConstraints/modern/phenomenology.py:23` lists only neutral-meson systems, `quarkConstraints/deltaf2.py:903`, `:922`, and `:941` are Delta F=2 mixing evaluators, and `flavorConstraints/muToEGamma.py:75` / `flavorConstraints/README.md:3` are lepton-sector dipole code/docs.
- CHK-6: `HIGH` is consistent with the rubric for both entries. B013 requires new b -> s gamma dipole matching/RG plus exclusive B_s -> phi form-factor and time-dependent helicity likelihood machinery. B014 requires analogous b -> d gamma dipole matching/RG plus exclusive B -> rho/omega form factors, annihilation/isospin treatment, CP inputs, and CKM-ratio handling.
- CHK-7: Searches of the methodology note, rc1.1 phase logs, and neighboring radiative beauty entries B011/B012 found no contradictory load-bearing B013/B014 numerical claim. B011/B012 mention these modes only as qualitative radiative companions; the old round-001 signoff one-line B011 wording is not a load-bearing numerical conflict.
- CHK-8: No `\textbf{CHECK}`, `TODO`, `FIXME`, `\cite`, `\ref`, unresolved-reference marker, or undefined-reference marker was found in either TeX. A lightweight balance check found B013 `\(`/`\)` = 50/50, `\[`/`\]` = 4/4, brace delta 0; B014 `\(`/`\)` = 58/58, `\[`/`\]` = 4/4, brace delta 0. `git diff --check` on the relevant TeX, sidecar, and worklog paths was clean.
- CHK-W8: Both sidecars contain `priority_tier: SECONDARY`, non-empty `priority_rationale`, and `promoted_in_wave: 8`. Both TeX files contain a SECONDARY-tier note pointing to `flavor_catalog/PRIORITY_TIERS.md`, and both TeX/YAML pairs live under `flavor_catalog/processes/secondary/beauty/`.
