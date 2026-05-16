# CA Worklog: ca_w5a_kaon_charm_v2
**Date**: 2026-05-16
**Family**: kaon_charm
**Process IDs**: K008

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| K008 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| K008 | BR(K_L -> pi0 e+e-) canonical limit | <2.8e-10 at 90% CL | PDG 2026 S013.20 / KTeV combined | 6578aa4b77ae |
| K008 | KTeV standalone BR(K_L -> pi0 e+e-) | <3.5e-10 at 90% CL | KTeV 2004 arXiv hep-ex/0309072 | 4133e47fd482 |
| K008 | KTeV observed/background events | 1 observed; 0.99 +/- 0.35 expected background | KTeV 2004 arXiv hep-ex/0309072 | 4133e47fd482 |
| K008 | CPC part from K_L -> pi0 gamma gamma | (0.0047 +0.0022 -0.0018)e-10 | PDG 2026 S013.20 | 6578aa4b77ae |
| K008 | BR(K_S -> pi0 e+e-), m_ee > 0.165 GeV | (3.0 +1.5 -1.2 stat +/-0.2 syst)e-9 | PDG 2026 S012.10 / NA48 | 27b6bbf006fc |
| K008 | BR(K_S -> pi0 e+e-) full-region extrapolation | (5.8 +2.9 -2.4)e-9 | PDG 2026 S012.10 / NA48 | 27b6bbf006fc |
| K008 | KTeV theory context for CPV/direct-CP scale | SM order 3e-11; CPV 8-45e-12; direct CP 3-6e-12; interference 40% | KTeV 2004 arXiv hep-ex/0309072 | 4133e47fd482 |
| K008 | Isidori-Smith-Unterdorfer electron coefficients | C_mix=(15.7 +/-0.3) abs(a_S)^2; C_int=(6.2 +/-0.3) abs(a_S); C_dir=2.4 +/-0.2; C_CPC approximately 0; abs(a_S)=1.2 +/-0.2; Im lambda_t/1e-4=1.36 +/-0.12 | Isidori, Smith, Unterdorfer 2004 | 4cf0225eda40 |
| K008 | Post-2008 long-distance context | Exploratory lattice K -> pi l+l- calculation; no superseding K008 limit | Christ et al. 2016 | 1af95b56cdff |
| K008 | Post-2008 rare-kaon context | K_{L,S} -> pi0 l+l- branching ratios and q^2 spectra remain in rare-kaon program | Aebischer, Buras, Kumar 2022 | 20cacd8c2c54 |

## Evidence notes
- Required reading completed: plan v1 Section D CA spec including writer/checker separation and worklog requirement, plus orchestrator decisions.
- CHK-1 PASS: numeric physics claims in `K008.tex` now have corresponding `pdg_or_equivalent.values` entries with year, value, uncertainty or null, source URL, access date, snapshot path, and sha256. This covers the WA-v1 failure items: KTeV event/background and standalone limit, CPC supporting value, K_S supporting values including the 0.165 GeV selection, KTeV theory-context numbers, and the Isidori-Smith-Unterdorfer coefficient block.
- CHK-2 PASS: every process-local key or snapshot filename cited in the TeX resolves to `flavor_catalog/references/K008/source_manifest.yaml`. Each manifest entry has a non-empty `snapshot_path`, and `git ls-files flavor_catalog/references/K008` shows the referenced text snapshots are tracked.
- CHK-3 PASS: `ls flavor_catalog/references/K008/` shows only `.txt` snapshots plus `source_manifest.yaml`; `find flavor_catalog/references/K008 -maxdepth 1 -type f -name '*.pdf'` returned no files.
- CHK-4 PASS: `status_history` contains `WRITER-INITIATED` followed by cycle-1 `WRITER-DONE`, `WRITER-REWORK`, and the cycle-2 `WRITER-DONE` from WA-v2, all with ISO 8601 timestamps. The latest pre-CA state was `WRITER-DONE`.
- CHK-5 PASS: the requested literal `rg -l -E "<process keyword>" ...` form fails in this environment because ripgrep treats `-E` as an encoding option. Equivalent `rg -l -e` and `rg -n -e` searches over `quarkConstraints/`, `qcd/`, `flavorConstraints/`, `neutrinos/`, `yukawa/`, `warpConfig/`, `solvers/`, `scanParams/`, and `tests/` found no `K_L/K0L -> pi0 e+e-`, `K_S/K0S -> pi0 e+e-`, `S013.20`, `a_S`, `C_mix/C_int/C_dir`, or `s -> d e+e-` implementation. Existing kaon hits are Delta F=2 only, e.g. `quarkConstraints/deltaf2.py:209` and `quarkConstraints/modern/phenomenology.py:23`.
- CHK-6 PASS: HIGH is consistent with the rubric because K008 needs a new Delta S=1 semileptonic rare-kaon observable, short-distance `s -> d e+e-` matching, CPV/CPC decomposition with interference, and long-distance inputs; the existing Delta F=2 SLL/SLR/VLL/VRR/LR1/LR2 basis is insufficient.
- CHK-7 PASS: targeted rc1.1/paper-facing searches found no load-bearing K008 rare-decay numbers outside the catalog process, reference, or worklog material; no contradiction with sealed rc1.1 paper numbers was found.
- CHK-8 PASS: no `\textbf{CHECK}`, TODO/FIXME, unresolved `??`, `\cite`, `\ref`, or unbalanced braces were found. Inline and display math delimiter counts are balanced. `chktex` is not installed in this environment.
- External/source spot checks were performed against PDG API S013.20 and S012.10 plus arXiv pages for hep-ex/0309072, hep-ex/0309075, hep-ph/0404127, 0804.1954, 1608.07585, and 2203.09524. The checked public source text matches the local snapshots and YAML metadata.

## Issues (if any)
- None.
