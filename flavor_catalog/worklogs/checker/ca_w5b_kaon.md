# CA Worklog: ca_w5b_kaon
**Date**: 2026-05-16
**Family**: kaon
**Process IDs**: K009 K010

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| K009 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| K010 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| K009 | BR(K_L -> pi0 mu+ mu-) canonical limit | <3.8e-10 at 90% CL | PDG 2026 S013.16 / S013R30 | c4ebbf6334ba |
| K009 | KTeV observed/background events and limit | 2 observed; 0.87 +/- 0.15 expected background; <3.8e-10 at 90% CL | KTeV 2000 arXiv hep-ex/0001006 | 0a16b654a73f |
| K009 | BR(K_S -> pi0 mu+ mu-) supporting input | (2.9 +1.5 -1.2 stat +/-0.2 syst)e-9; 6 events; 0.22 +0.18 -0.11 background | PDG 2026 S012.15 / NA48 | 21f62a439461 |
| K009 | K_L -> pi0 mu+ mu- theory decomposition | C_mix=(3.7 +/-0.1)|a_S|^2; C_int=(1.6 +/-0.1)|a_S|; C_dir=1.0 +/-0.1; C_CPC=5.2 +/-1.6; |a_S|=1.2 +/-0.2; SM constructive=(1.5 +/-0.3)e-11 | Isidori, Smith, Unterdorfer 2004 | 9634d1af1d93 |
| K010 | BR(K_S -> pi0 e+e-), m_ee > 0.165 GeV | (3.0 +1.5 -1.2 stat +/-0.2 syst)e-9 | PDG 2026 S012.10 / NA48 | 87a298fc5bb6 |
| K010 | BR(K_S -> pi0 e+e-) full-region extrapolation | (5.8 +2.9 -2.4)e-9 | PDG 2026 S012.10 / NA48 | 87a298fc5bb6 |
| K010 | NA48/1 observed/background events | 7 observed; 0.15 expected background | NA48/1 2003 arXiv hep-ex/0309075 | 406786ac199a |

## Evidence notes
- Required reading completed: plan v1 CA spec and writer/checker separation rule, plus orchestrator decisions.
- CHK-1 failed for K009 under the strict batch rule. The canonical PDG limit is in `pdg_or_equivalent`, and the supporting/theory values match local snapshots and live spot checks, but several TeX numerical claims are only in `supporting_experimental_values` or `theory_decomposition`, not `pdg_or_equivalent` entries with full metadata.
- CHK-1 failed for K010 under the strict batch rule. The PDG partial-rate and extrapolated-total values are in `pdg_or_equivalent`, but the TeX claim of seven observed events and 0.15 expected background is only in `supporting_experimental_values`.
- CHK-2 passed for both: every TeX source key or source filename resolves in `flavor_catalog/references/<process_id>/source_manifest.yaml`, and every manifest entry has a non-empty tracked `.txt` snapshot path under the matching process reference directory.
- CHK-3 passed for both: `ls flavor_catalog/references/K009/` and `ls flavor_catalog/references/K010/` showed only text snapshots, manifests, and K010's `sha256sums.txt`; `find ... -iname '*.pdf'` returned no PDFs.
- CHK-4 passed for both: pre-CA `status_history` contained `WRITER-INITIATED` followed by `WRITER-DONE` with ISO 8601 timestamps, and `WRITER-DONE` was the latest pre-check entry.
- CHK-5 passed for both. The literal requested `rg -l -E "<process keyword>" ...` form fails in this environment because `rg -E` is an encoding flag; equivalent `rg -l -e ...` searches found only existing Delta-F=2 kaon-mixing identifiers, not K009/K010 rare semileptonic implementations. The nearby cited lines `quarkConstraints/deltaf2.py:209`, `quarkConstraints/deltaf2.py:755`, and `quarkConstraints/modern/phenomenology.py:23` exist.
- CHK-6 passed for both: HIGH is consistent with the rubric because both processes need new Delta-S=1 semileptonic rare-kaon observable handling beyond the existing Delta-F=2 SLL/SLR/VLL/VRR/LR machinery; K009 additionally needs CPV/interference and CPC two-photon bookkeeping.
- CHK-7 passed for both: focused repo searches outside `flavor_catalog/` found only catalog planning rows for K009/K010 and existing Delta-F=2 kaon-mixing code; no rc1.1 load-bearing rare-decay number is contradicted.
- CHK-8 passed for both: no `\textbf{CHECK}`, TODO/FIXME, `\cite`, `\ref`, unresolved `??`, or simple brace imbalance was found. `chktex` is not installed in this environment, so this was a manual/grep check.
- External spot checks against PDG Live/API S013.16, S012.15, S012.10 and arXiv abstract pages hep-ex/0001006 and hep-ex/0309075 matched the local snapshots.

## Issues (if any)
- K009: CHK-1 rework required. Move or duplicate every load-bearing numerical claim in `K009.tex` into YAML `pdg_or_equivalent` entries with complete metadata, or reduce the TeX numeric claims to only those already represented there. At minimum cover the KTeV event/background numbers, the supporting K_S -> pi0 mu+ mu- value/background numbers if retained in K009 text, and the Isidori-Smith-Unterdorfer coefficient/SM-expectation block.
- K010: CHK-1 rework required. Move or duplicate the NA48/1 seven-observed-events and 0.15-background claims into YAML `pdg_or_equivalent` entries with complete metadata, or remove those event/background numbers from the TeX.
