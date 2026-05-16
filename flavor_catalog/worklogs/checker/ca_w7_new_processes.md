# CA Worklog: ca_w7_new_processes
**Date**: 2026-05-16
**Family**: mixed
**Process IDs**: T003 T004 T008 T012 B012

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| T003 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| T004 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| T008 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| T012 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| B012 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| T003 | B(t -> c gamma), CMS observed expected | <1.51e-5 observed, <1.54e-5 expected at 95% CL | CMS-TOP-21-013 / PRD 109 (2024) 072004 | 5f094c0620a2 |
| T003 | B(t -> c gamma), ATLAS LH/RH | <4.2e-5 LH, <4.5e-5 RH at 95% CL | ATLAS arXiv:2205.02537 | 82494cd375fd |
| T003 | Gamma(t -> gamma q)/Gamma(t -> Wb), q = u,c | <9.5e-6 at 95% CL | PDG 2026 pdgLive top summary | e05906a43ebc |
| T003 | SM B(t -> c gamma) | 4.6e-14 central value | Aguilar-Saavedra 2004 | 01040933ce18 |
| T004 | B(t -> gamma q), q = u,c | <9.5e-6 at 95% CL | PDG 2026 pdgLive t-quark listing | 5e72da33cd1b |
| T004 | B(t -> u gamma), ATLAS LH/RH | <0.85e-5 LH, <1.2e-5 RH at 95% CL | ATLAS arXiv:2205.02537 | 1b92967f8b52 |
| T004 | B(t -> u gamma), CMS observed expected | <0.95e-5 observed, 1.20e-5 expected at 95% CL | CMS-TOP-21-013 | 1c729a489ac5 |
| T004 | HL-LHC/FCC-hh projected reach | 1e-5 / 1e-6 level | Aguilar-Saavedra arXiv:1709.03975 | 21341a8358a9 |
| T008 | B(t -> H u), PDG headline | <1.9e-4 at 95% CL | PDG 2026 pdgLive t -> H u datablock | 33554bc8abef |
| T008 | B(t -> H u), CMS diphoton | <0.019% observed, <0.031% expected at 95% CL | CMS arXiv:2111.02219 | 6fa74c1e8b51 |
| T008 | B(t -> H u), ATLAS multilepton/combination | <2.8(3.0)e-4 and <2.6(1.8)e-4 observed(expected) | ATLAS arXiv:2404.02123 | 358dbdb2af5b |
| T008 | B(t -> H u), CMS same-sign/combination | <0.072%(0.059%) and <0.019%(0.027%) observed(expected) | CMS arXiv:2407.15172 | 6de117876a96 |
| T012 | R_c^0 | 0.1721 +/- 0.0030 | PDG 2025 Z-boson listing | 4394ce4d5f6f |
| T012 | A_FB^{0,c} | 0.0707 +/- 0.0035, converted from 7.07 +/- 0.35% | PDG 2025 Z-boson listing | 4394ce4d5f6f |
| T012 | A_c | 0.670 +/- 0.027 | PDG 2025 Z-boson listing | 4394ce4d5f6f |
| T012 | LEP/SLC charm cross-check | Same R_c^0, A_FB^{0,c}, A_c central values and uncertainties | LEP/SLC final Z-resonance combination | 6f373a37f4cc |
| B012 | B(B0 -> K*0 gamma) | (41.63 +/- 0.92)e-6 | HFLAV rare decays, end-of-year 2024 | dc826a3f6d24 |
| B012 | B(B+ -> K*+ gamma) | (39.5 +/- 1.0)e-6 | HFLAV rare decays, end-of-year 2024 | dc826a3f6d24 |
| B012 | Delta_0+(B -> K* gamma) | 0.059 +/- 0.014 | HFLAV rare decays, end-of-year 2024 | dc826a3f6d24 |
| B012 | A_CP(B -> K* gamma) | -0.004 +/- 0.011 | HFLAV rare decays, end-of-year 2024 | dc826a3f6d24 |
| B012 | S_{K*gamma}, C_{K*gamma} | -0.08 +/- 0.17, 0.03 +/- 0.10 | PDG pdgLive 2025 | 3ba65881b19e |
| B012 | Belle/LHCb post-2008 significance claims | 3.1 sigma isospin evidence, 5.2 sigma photon-polarization observation | Belle arXiv:1707.00394; LHCb arXiv:1402.6852 | 826984017c41; 664bac221a55 |

## Evidence notes
- CHK-1: Source snapshots were reopened with `rg` against the quoted values. T008 and T012 numerical observable claims trace to `pdg_or_equivalent` entries with source metadata. T003, T004, and B012 have verified source values, but some TeX numerical claims do not reside in `pdg_or_equivalent` as required.
- CHK-2: All process-local key references in the TeX `Key references` sections resolve to entries in `flavor_catalog/references/<process_id>/source_manifest.yaml`. TeX value-id markers such as `CMS2024:T003:tcgamma` were checked against the YAML sidecar. Each manifest entry has a non-empty `snapshot_path` under the process reference directory and every snapshot path is tracked by `git ls-files`.
- CHK-3: `ls flavor_catalog/references/{T003,T004,T008,T012,B012}` showed only `.txt` snapshots, manifests, and B012 `sha256sums.txt`; no `.pdf` files are tracked.
- CHK-4: Before checker updates, each YAML sidecar had `WRITER-INITIATED` followed by `WRITER-DONE` with ISO 8601 timestamps, and `WRITER-DONE` was the most recent entry.
- CHK-5: The requested `rg -l -E "<pattern>" ...` form fails in this ripgrep because `-E` is the encoding flag, so the checks were rerun with `rg -l -e "<pattern>" ...`. T003, T004, T008, and B012 returned no source-code implementation hits after excluding notebooks where needed. T012 returned only unrelated `R_chi` Delta-F=2 hadronic helpers and bridge-matching class names; the cited generic Z-mass support lines exist at `qcd/running.py:3`, `qcd/constants.py:11`, and `quarkConstraints/qcd_running.py:100`.
- CHK-6: `HIGH` is consistent for all five entries: each would need a new top-FCNC, electroweak-pole, or b -> s gamma dipole/RG/mode calculation path rather than the existing Delta-F=2 operator basis.
- CHK-7: No entry revises rc1.1 scan numbers; the TeX only cites existing Delta-F=2 or lepton-dipole code as coverage context.
- CHK-8: No `CHECK`, `TODO`, `??`, `\ref{...}`, or `\cite{...}` markers were found in the five TeX files. `chktex` and `lacheck` are not installed in this environment, so the cosmetic check was manual plus grep-based.

## Issues (if any)
- T003: CHK-1 fail. `T003.tex` quotes `B(t -> c gamma)=4.6e-14` for the SM expectation; the YAML records this under `paper_era_reference`, not under `pdg_or_equivalent`.
- T004: CHK-1 fail. `T004.tex` quotes 139 fb^-1, 138 fb^-1, and 10^-5/10^-6 reach claims that are in `auxiliary_experimental_context` or `paper_era_reference`, not in `pdg_or_equivalent`.
- B012: CHK-1 fail. The canonical `pdg_or_equivalent` observable blocks lack explicit `year` keys and use `sha256_of_text_snapshot` rather than `sha256`; the 3.1 sigma, 5.2 sigma, and 2019--2022 TeX claims are under `post_2008_sources`, not `pdg_or_equivalent`.
