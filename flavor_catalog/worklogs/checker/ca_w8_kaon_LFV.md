# CA Worklog: ca_w8_kaon_LFV
**Date**: 2026-05-17
**Family**: kaon LFV, Wave-8 SECONDARY
**Process IDs**: K019 K020 K021
**Cycle**: 1

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | CHK-W8 | overall |
|---|---|---|---|---|---|---|---|---|---|---|
| K019 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| K020 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| K021 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| K019 | BR(K_L -> e^+/- mu^-/+), PDG canonical limit | < 4.7e-12 at 90% CL | PDG 2025 K_L0 listing | b5797d539dbe |
| K019 | BR(K_L -> mu^+/- e^-/+), direct experiment | < 4.7e-12 at 90% CL | Ambrose et al. / BNL E871, arXiv:hep-ex/9811038 | 98edf5c4d18c |
| K020 | BR(K+ -> pi+ mu+ e-), PDG combined limit | < 1.3e-11 at 90% CL | PDG Live / PDG API S010.29 | 8ddac950390a |
| K020 | BR(K+ -> pi+ mu- e+), PDG current limit | < 6.6e-11 at 90% CL | PDG Live / PDG API S010.25 / NA62 2021 | 8ddac950390a |
| K020 | BR(K+ -> pi+ mu+ e-), E865-only limit quoted in TeX | < 2.1e-11 at 90% CL | Sher et al. / BNL E865, arXiv:hep-ex/0502020 | 0ae93baf10b7 |
| K021 | BR(K_L -> pi0 e^+- mu^-+), PDG canonical limit | < 7.6e-11 at 90% CL | PDG 2025 K_L0 listing | a6e37f9494c9 |
| K021 | BR(K_L -> pi0 mu e), exact KTeV limit | < 7.56e-11 at 90% CL | Abouzaid et al. / KTeV, arXiv:0711.3472 | e17b4f9fced0 |
| K021 | BR(K+ -> pi+ mu- e+), supporting NA62 companion | < 6.6e-11 at 90% CL | NA62 2021, arXiv:2105.06759 | 14998087e135 |
| K021 | BR(K+ -> pi- mu+ e+), supporting NA62 companion | < 4.2e-11 at 90% CL | NA62 2021, arXiv:2105.06759 | 14998087e135 |
| K021 | BR(pi0 -> mu- e+), supporting NA62 companion | < 3.2e-10 at 90% CL | NA62 2021, arXiv:2105.06759 | 14998087e135 |

## Issues
- K019: None.
- K020: CHK-1 failure. `K020.tex` quotes the Sher/E865-only experimental limit `BR(K+ -> pi+ mu+ e-) < 2.1 x 10^-11` in the Post-2008 developments paragraph. The source snapshots verify this number, and the PDG snapshot also contains the same E865-only subentry, but `K020.yaml` has no dedicated `pdg_or_equivalent.values` entry carrying that observable with value, uncertainty, year, units, source URL, access date, snapshot path, and sha256. This is a measured branching-fraction upper limit, not a theory normalization scale, dataset metadata, or lattice/theory average covered by the L001/B001/B021 carve-out.
- K021: None.

## Evidence notes
- CHK-1: K019's TeX numerical experimental observable is the PDG/E871 `4.7e-12` upper limit, represented in `pdg_or_equivalent.values` with complete metadata. K021's PDG `7.6e-11`, KTeV `7.56e-11`, and NA62 companion limits quoted/ranged in prose are represented in `pdg_or_equivalent.values`; the KTeV zero-event statement is dataset/event-yield context and is also kept under `supporting_measurements`, so it is not a CHK-1 failure under the L001/B001/B021 policy. K020's current PDG `1.3e-11` and `6.6e-11` limits are structured, but the TeX-retained Sher/E865-only `2.1e-11` measured limit is not.
- CHK-2: All process-local keys listed in the TeX `Key references` sections resolve to entries in `flavor_catalog/references/K019/source_manifest.yaml`, `K020/source_manifest.yaml`, or `K021/source_manifest.yaml`. `rg` found no empty `snapshot_path`, and `git ls-files --error-unmatch` confirmed every manifest snapshot path is tracked.
- CHK-3: `find flavor_catalog/references/K019 K020 K021 -type f` for publisher-document extensions found no PDFs or office documents. The reference directories contain manifests, checksum files, and text snapshots only.
- CHK-4: Each sidecar has `DRAFT -> WRITER-INITIATED -> WRITER-DONE` with ISO 8601 `timestamp`/`at` fields. The most recent pre-check state was `WRITER-DONE` for all three.
- CHK-5: The `NO` code-coverage claims are supported. K019 exact LFV-kaon searches find no implementation; cited nearby lines exist at `flavorConstraints/muToEGamma.py:1`, `:15-21`, `scanParams/scan.py:523-535`, and `quarkConstraints/deltaf2.py:610-624`. K020 exact K020/S010/K+ pi+ e mu searches return only generic LFV metadata and the mu->e gamma lane; cited lines exist at `flavorConstraints/muToEGamma.py:75`, `scanParams/scan.py:524`, `quarkConstraints/deltaf2.py:755`, and `quarkConstraints/modern/phenomenology.py:23`. K021 exact pi0 e mu searches return no hits, and semileptonic LFV kaon searches return no implementation; cited nearby lines exist at `quarkConstraints/deltaf2.py:1`, `:615`, and `flavorConstraints/muToEGamma.py:1`.
- CHK-6: K019 `MEDIUM` matches the rubric because it needs a new Delta S = 1 charged-LFV semileptonic operator/matching interface with standard clean two-body inputs, but no new lattice/RG/long-distance calculation. K020 and K021 `HIGH` match the rubric because production support needs new semileptonic LFV operator plumbing plus K -> pi e mu mode-rate/form-factor calculations and cross-sector RS matching.
- CHK-7: Searches outside the batch found only plan/deferred-scope rows, nearby already-approved kaon channels, writer/orchestration logs, and generic lepton-sector scope notes. No rc1.1 or methodology-note load-bearing number contradicts the K019/K020/K021 limits or difficulty/coverage classifications.
- CHK-8: `rg` found no unresolved `\cite`, `\ref`, `\textbf{CHECK}`, TODO/FIXME/MISSING, or `??` markers in the three TeX files. `git diff --check` was clean for the relevant TeX/YAML paths. TeX delimiter spot checks showed balanced `{}` counts, no dollar math markers, and even `\(`/`\)` plus `\[`/`\]` counts.
- CHK-W8: All three `.tex` and `.yaml` files live under `flavor_catalog/processes/secondary/kaon/`. Each sidecar has `priority_tier: SECONDARY`, a non-empty `priority_rationale`, and `promoted_in_wave: 8`. Each TeX file carries a SECONDARY Wave-8 note pointing at `flavor_catalog/PRIORITY_TIERS.md`.
