# CA Worklog: ca_w4_charm
**Date**: 2026-05-16
**Family**: charm
**Process IDs**: C002, C003, C004, C007

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| C002 | PASS | PASS | PASS | FAIL | PASS | PASS | PASS | PASS | WRITER-REWORK |
| C003 | PASS | PASS | PASS | FAIL | PASS | PASS | PASS | PASS | WRITER-REWORK |
| C004 | PASS | PASS | PASS | FAIL | PASS | FAIL | PASS | PASS | WRITER-REWORK |
| C007 | FAIL | PASS | PASS | FAIL | PASS | PASS | PASS | PASS | WRITER-REWORK |

## Evidence notes
- Batch handoff: no WA batch worklog for C002, C003, C004, or C007 was found under `flavor_catalog/worklogs/writer/`. The C002 and C003 sidecars stop at `PKA-DONE`; C004 and C007 stop at `WRITER-INITIATED`. None has `WRITER-DONE`, so CHK-4 fails for every process.
- CHK-1: C002, C003, and C004 TeX numerical physics claims are represented in `pdg_or_equivalent` entries with year, source URL, access date, snapshot path, and local sha256. C007 fails because the TeX claims "25 rare or forbidden" modes and short-distance SM branching fractions of order `10^-12`, both visible in `lhcb2021_arxiv2011_00217.txt`, but neither is recorded in a `pdg_or_equivalent` value block.
- CHK-2/3: all TeX source keys or filenames resolve to entries in each process-local `source_manifest.yaml`; every manifest snapshot path is non-empty and tracked by `git ls-files`. `find flavor_catalog/references/C002 ... C007 -name '*.pdf'` returned no PDFs.
- CHK-5: the exact prompt form `rg -l -E "<process keyword>" ...` fails in this ripgrep build because `-E` is parsed as an encoding flag; I reran the equivalent searches with `-e`. C002 PARTIAL is supported by the cited D0 mixing and paper-mode complex `M12_D0_NP_GeV` lines, with no `|q/p|` or `phi_D` evaluator. C003, C004, and C007 `NO` claims are supported: focused searches found only unrelated neutrino `delta_CP`, D0 mixing/hadronic variables, or `muToEGamma` hits, not the requested charm observables.
- CHK-6: C002 LOW is consistent with reuse of the existing Delta F = 2 basis plus a missing CPV wrapper; C003 HIGH and C007 HIGH are consistent with new Delta C = 1 operator/mode work. C004 MEDIUM is inconsistent with the stated rubric because a live `D0 -> mu+ mu-` constraint needs a new Delta C = 1 rare-leptonic mode calculation and Wilson normalization; by the prompt rubric this is HIGH.
- CHK-7: targeted non-catalog grep for these process keywords found only flavor-catalog plan/signoff rows, not rc1.1 paper values, so there is no detected contradiction with load-bearing rc1.1 numbers.
- CHK-8: no `\textbf{CHECK}`, `TODO`, unresolved `\ref`/`\cite`, or `??` markers were found; visual inspection did not find malformed math mode or missing braces.

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| C002 | `|q/p|` | `0.983 +0.015/-0.014`, 95% C.L. `[0.955, 1.012]` | HFLAV CKM25 all-CPV global fit | `5908d28e0388` |
| C002 | `phi_D` | `-1.51 +1.03/-1.06 deg`, 95% C.L. `[-3.63, 0.51]` | HFLAV CKM25 all-CPV global fit | `5908d28e0388` |
| C002 | no-indirect-CPV test | `Delta chi2 = 2.16`, `0.95 sigma` | HFLAV CKM25 all-CPV global fit | `5908d28e0388` |
| C002 | PDG cross-check | `|q/p| = 0.983 +/- 0.015`; `phi_D = -1.51 +/- 1.04 deg` | PDG 2025 D0-D0bar mixing review excerpt | `4ffedbd5dcda` |
| C003 | `Delta a_CP^dir` | `(-0.159 +/- 0.029)%` | HFLAV direct/indirect CPV combination | `50781ba7963b` |
| C003 | `a_CP^ind` | `(-0.010 +/- 0.012)%` | HFLAV direct/indirect CPV combination | `50781ba7963b` |
| C003 | no-CPV discrepancy | `Delta chi2 = 32.3`; `5.3 sigma` | HFLAV direct/indirect CPV combination | `50781ba7963b` |
| C003 | `Delta A_CP` | `(-15.4 +/- 2.9) x 10^-4 = (-0.154 +/- 0.029)%` | LHCb 2019 | `c0d7f9a85048` |
| C003 | `A_CP^K - A_CP^pi` | `(-0.154 +/- 0.029)%` | HFLAV CKM25 DCPV input | `ca6fa2bc4a1` |
| C003 | `A_pi`, `A_K` | `A_pi = (0.225 +/- 0.057)%`; `A_K = (0.068 +/- 0.051)%` | HFLAV CKM25 table extraction | `ca6fa2bc4a1` |
| C004 | `BR(D0 -> mu+ mu-)` | `<2.1 x 10^-9` at 90% C.L. | PDG Live/API S032.28 | `34cff61f22c2` |
| C004 | CMS companion limit | `<2.4 x 10^-9` at 95% C.L.; `64.5 fb^-1`, `13.6 TeV`, `2022-2023` | CMS 2025 public result page | `46aaaf8dc7fb` |
| C004 | LHCb predecessor limit | `<3.1 x 10^-9` at 90% C.L.; companion `<3.5 x 10^-9` at 95% C.L.; `9 fb^-1` at `7, 8, 13 TeV` | LHCb 2023 | `a94f1dbcfa8d` |
| C004 | long-distance two-photon relation | `2.7 x 10^-5 * BR(D0 -> gamma gamma)`; floor at least `3 x 10^-13` | Burdman et al. 2002 | `ec7271413986` |
| C007 | `BR(D+ -> pi+ mu+ mu-)` | `<6.7 x 10^-8` at 90% C.L. | PDG Live/API S031.42 | `1bfa609f88a7` |
| C007 | LHCb 2021 limit | `<6.7 x 10^-8` at 90% C.L.; companion `<7.4 x 10^-8` at 95% C.L.; `1.6 fb^-1` in `2016` | LHCb 2021 | `7e293a7a298e` |
| C007 | LHCb 2013 predecessor | `<7.3(8.3) x 10^-8` at 90% (95%) C.L.; `1.0 fb^-1` at `7 TeV` | LHCb 2013 | `64daf4ae376` |
| C007 | LHCb search scope | `25` rare/forbidden modes | LHCb 2021 snapshot; missing YAML `pdg_or_equivalent` value | `7e293a7a298e` |
| C007 | short-distance SM branching scale | order `10^-12` | LHCb 2021 snapshot; missing YAML `pdg_or_equivalent` value | `7e293a7a298e` |

## Issues (if any)
- C002: CHK-4 fail. Append or rerun WA so the sidecar has an ordered `WRITER-INITIATED` -> `WRITER-DONE` history with most recent pre-CA status `WRITER-DONE`, and add a WA batch worklog that includes C002.
- C003: CHK-4 fail. Append or rerun WA so the sidecar has an ordered `WRITER-INITIATED` -> `WRITER-DONE` history with most recent pre-CA status `WRITER-DONE`, and add a WA batch worklog that includes C003.
- C004: CHK-4 fail. Append or rerun WA so the sidecar has an ordered `WRITER-INITIATED` -> `WRITER-DONE` history with most recent pre-CA status `WRITER-DONE`, and add a WA batch worklog that includes C004.
- C004: CHK-6 fail. Change the implementation difficulty to HIGH or justify why this does not count as a new mode calculation under the prompt rubric.
- C007: CHK-1 fail. Add `pdg_or_equivalent` entries for the TeX claims that the LHCb 2021 search covers 25 modes and that short-distance SM branching fractions are of order `10^-12`, or remove those numerical claims in the next WA pass.
- C007: CHK-4 fail. Append or rerun WA so the sidecar has an ordered `WRITER-INITIATED` -> `WRITER-DONE` history with most recent pre-CA status `WRITER-DONE`, and add a WA batch worklog that includes C007.
