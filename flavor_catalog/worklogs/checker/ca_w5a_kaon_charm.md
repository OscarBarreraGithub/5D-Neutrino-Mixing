# CA Worklog: ca_w5a_kaon_charm
**Date**: 2026-05-16
**Family**: kaon_charm
**Process IDs**: K008 C005

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| K008 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| C005 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| K008 | BR(K_L -> pi0 e+e-) | <2.8e-10 at 90% CL | PDG 2026 S013.20 / KTeV combined | 6578aa4b77ae |
| K008 | KTeV standalone BR(K_L -> pi0 e+e-) | <3.5e-10 at 90% CL | KTeV 2004 arXiv hep-ex/0309072 | 4133e47fd482 |
| K008 | KTeV observed/background events | 1 observed; 0.99 +/- 0.35 expected background | KTeV 2004 arXiv hep-ex/0309072 | 4133e47fd482 |
| K008 | CPC part from K_L -> pi0 gamma gamma | (0.0047 +0.0022 -0.0018)e-10 | PDG 2026 S013.20 | 6578aa4b77ae |
| K008 | BR(K_S -> pi0 e+e-), m_ee > 0.165 GeV | (3.0 +1.5 -1.2 stat +/-0.2 syst)e-9 | PDG 2026 S012.10 / NA48 | 27b6bbf006fc |
| K008 | BR(K_S -> pi0 e+e-) full-region extrapolation | (5.8 +2.9 -2.4)e-9 | PDG 2026 S012.10 / NA48 | 27b6bbf006fc |
| K008 | KTeV theory context for CPV/direct-CP scale | SM order 3e-11; CPV 8-45e-12; direct CP 3-6e-12; interference 40% | KTeV 2004 arXiv hep-ex/0309072 | 4133e47fd482 |
| K008 | Isidori-Smith-Unterdorfer electron coefficients | C_mix=(15.7 +/-0.3) abs(a_S)^2; C_int=(6.2 +/-0.3) abs(a_S); C_dir=2.4 +/-0.2; abs(a_S)=1.2 +/-0.2; Im lambda_t/1e-4=1.36 +/-0.12 | Isidori, Smith, Unterdorfer 2004 | 4cf0225eda40 |
| C005 | BR(D0 -> e+e-) canonical limit | <7.9e-8 at 90% CL | PDG 2026 S032.39 | 11656d94c38f |
| C005 | BR(D0 -> e+e-) Belle limit | <7.9e-8 at 90% CL; 660 fb^-1 | Belle 2010 arXiv:1003.2345 | cb1da7b14457 |
| C005 | BR(D0 -> e+e-) BaBar limit | <1.7e-7 at 90% CL; 468 fb^-1 | BaBar 2012 arXiv:1206.5419 | baa274320714 |
| C005 | Rare-charm model-dependence context | c -> u l+l- short-distance and long-distance context | Burdman et al. 2002 | 25f176716b7e |
| C005 | RS/anarchic-flavor baseline context | CFW 2008 source used only for RS framing | Csaki-Falkowski-Weiler 2008 | 91056272adfb |

## Evidence notes
- Required reading completed: plan v1 CA spec and writer/checker separation rule, plus orchestrator decisions.
- K008 CHK-1 failed under the strict batch rule. The main PDG limit is under `pdg_or_equivalent`, and all checked numbers match local snapshots and `source_shas`, but several numerical claims in `K008.tex` are only in `supporting_experimental_values`, `ktev_theory_context`, or `theory_decomposition`, not in `pdg_or_equivalent` entries with year, value, uncertainty, source URL, access date, and sha256. Examples are the KTeV one-candidate/background count, the standalone 3.5e-10 limit, the KTeV theory scale/range/interference numbers, and the Isidori-Smith-Unterdorfer coefficient block.
- C005 CHK-1 passed: the PDG, Belle, and BaBar numerical claims in `C005.tex` are represented in `pdg_or_equivalent` subblocks with year, value, confidence level, source URL, access date, snapshot path, and sha256.
- CHK-2 passed for both: every TeX source key or snapshot filename resolves to a `source_manifest.yaml` entry, and each manifest entry has a non-empty tracked `.txt` snapshot path under the matching process reference directory.
- CHK-3 passed for both: `find flavor_catalog/references/K008 flavor_catalog/references/C005 -maxdepth 1 -type f -name '*.pdf'` returned no PDF files; the directories contain only tracked `.txt` snapshots plus `source_manifest.yaml`.
- CHK-4 passed for both: sidecar `status_history` contained `WRITER-INITIATED` followed by `WRITER-DONE`, with ISO 8601 timestamps, and `WRITER-DONE` was the most recent pre-CA entry.
- CHK-5 passed for both. The installed `rg` treats `-E` as an encoding option, so the literal `rg -l -E` form errors before searching; equivalent ripgrep regex searches with `rg -l -e ... -g '!*.ipynb'` returned no K008 or C005 implementation matches in `quarkConstraints/`, `qcd/`, `flavorConstraints/`, `neutrinos/`, `yukawa/`, `warpConfig/`, `solvers/`, `scanParams/`, or `tests/`. Nearby C005 hits are D0 mixing or `muToEGamma`, not `D0 -> e+e-`.
- CHK-6 passed for both: HIGH is consistent with the rubric because K008 needs a new Delta S=1 semileptonic rare-kaon mode with CP decomposition/long-distance inputs, and C005 needs a new Delta C=1 rare leptonic mode calculation rather than the existing Delta F=2 operator lane.
- CHK-7 passed for both: no load-bearing K008 or C005 rare-decay numbers appear in the rc1.1 paper-facing files; repo hits are catalog planning rows or existing Delta F=2 D0/K mixing machinery.
- CHK-8 passed for both: no `\textbf{CHECK}`, `\cite`, `\ref`, TODO/FIXME markers, unresolved `??`, or unbalanced braces were found. `chktex` is not installed in this environment, so this was a manual/grep check.
- External spot checks were performed against PDG Live/API S013.20, S012.10, and S032.39 and arXiv abstract pages for hep-ex/0309072, 1003.2345, 1206.5419, and hep-ph/0112235; the checked values match the local snapshots.
- K010 exists as its own process draft carrying the K_S -> pi0 e+e- measurement used as K008 supporting input.

## Issues (if any)
- K008: CHK-1 rework required. Move or duplicate every load-bearing numerical claim in `K008.tex` into YAML `pdg_or_equivalent` entries with complete metadata, or reduce the TeX numeric claims to only those already represented there. At minimum cover the KTeV standalone/event/background numbers, the PDG CPC supporting value, the K_S supporting values if retained in K008 text, the KTeV theory-context numbers, and the Isidori-Smith-Unterdorfer coefficient block.
