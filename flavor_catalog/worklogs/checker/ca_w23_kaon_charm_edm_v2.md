# CA Worklog: ca_w23_kaon_charm_edm_v2
**Date**: 2026-05-16
**Family**: kaon_charm_edm
**Process IDs**: K002 C001 E001

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| K002 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| C001 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| E001 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| K002 | Delta m_K, CPT-assuming PDG fit | (0.5293 +/- 0.0009) x 10^10 hbar s^-1, scale factor 1.3 | PDG pdgLive S013D 2026 | b19ce40faf7d |
| K002 | Delta m_K, no-CPT companion fit | (0.5289 +/- 0.0010) x 10^10 hbar s^-1 | PDG pdgLive S013D 2026 | b19ce40faf7d |
| K002 | Exploratory lattice Delta M_K | 3.19(41)(96) x 10^-12 MeV | Bai et al. 2014 | 94522a76f8bf |
| K002 | Physical-mass lattice statistical uncertainty | near 9% | Wang 2023 | 26d8bab47c57 |
| C001 | x_D | (0.405 +/- 0.043)%, 95% C.L. [0.320, 0.489]% | HFLAV CKM25 all-CPV fit | 85254267f434 |
| C001 | y_D | (0.636 +/- 0.024)%, 95% C.L. [0.590, 0.682]% | HFLAV CKM25 all-CPV fit | 85254267f434 |
| C001 | Delta m_D | (0.997 +/- 0.116) x 10^10 hbar s^-1 = (6.56 +/- 0.76) x 10^-15 GeV | PDG Live S032D / HFLAV 2025 | 398d3420f232 |
| C001 | Delta y | (0.031 +/- 0.035 +/- 0.013)% | HFLAV CKM25 input table | 024569bebcc2 |
| C001 | DeltaGamma_D/Gamma_D | 2 y_D = (1.272 +/- 0.048)% | HFLAV CKM25 all-CPV fit | 85254267f434 |
| C001 | Generic RS KK-gluon scale | around 21 TeV | Csaki-Falkowski-Weiler 2008 | 224e28f70091 |
| C001 | Composite pseudo-Goldstone KK-gluon scale | around 33 TeV | Csaki-Falkowski-Weiler 2008 | 224e28f70091 |
| E001 | Current electron EDM limit | \|d_e\| < 4.1 x 10^-30 e cm at 90% CL | PDG Live S003EDM 2026 | aab24570becf |
| E001 | Roussy et al. measurement | (-1.3 +/- 2.0_stat +/- 0.6_syst) x 10^-30 e cm; limit < 4.1 x 10^-30 e cm at 90% confidence | Roussy et al. 2023 | 85dfc60b010c |
| E001 | Roussy improvement factor | about 2.4 over previous best upper bound | Roussy et al. 2023 | 85dfc60b010c |
| E001 | ACME/Andreev 2018 benchmark | \|d_e\| < 1.1 x 10^-29 e cm at 90% CL; (4.3 +/- 3.1 +/- 2.6) x 10^-30 e cm | ACME / Andreev et al. 2018 | c81e0a1936ec |

## Issues (if any)
- None.

## Verification notes
- CHK-1: Grepped the three TeX files for numerical claims and checked them against the current `pdg_or_equivalent` blocks. The K002 promoted no-CPT PDG fit, Bai 2014 lattice value, and Wang 2023 near-9% claim are now present with year, URL, access date, snapshot path, and sha256. C001 now has the two CFW 21 TeV and 33 TeV contextual scale blocks. E001 has explicit years plus structured value, uncertainty, limit, confidence-level, source URL, access date, snapshot path, and sha256 fields for the PDG, Roussy, and ACME numerical claims.
- CHK-2: Every process-local source key cited in the TeX resolves to `flavor_catalog/references/<process_id>/source_manifest.yaml`. Each manifest entry has a non-empty `snapshot_path`; `git ls-files --error-unmatch` and `test -s` passed for every snapshot path.
- CHK-3: `ls flavor_catalog/references/{K002,C001,E001}/` showed only text snapshots, source manifests, and C001 `sha256sums.txt`; `find ... -iname '*.pdf'` returned no PDFs.
- CHK-4: Each sidecar contains `WRITER-INITIATED` before `WRITER-DONE`, followed by the cycle-2 `WRITER-DONE` from WA-v2 with ISO-8601 timestamps. This CA pass appends `CHECKER-DONE` as the most recent state.
- CHK-5: The literal requested `rg -l -E "<pattern>" ...` form fails on this installed ripgrep because `-E` is parsed as an encoding flag. Equivalent `rg -l -e` sweeps confirmed K002 hits `DELTA_M_K/evaluate_delta_mk`, C001 hits `D0/d_mix/evaluate_d0_mixing`, and a focused E001 grep for electron-EDM implementation terms returned no hits in the required implementation/test directories. The cited file:line references exist.
- CHK-6: Difficulty labels are consistent with the rubric. K002 and C001 are LOW because they use the existing Delta F = 2 VLL/VRR/SLL/SLR/LR machinery for conservative neutral-meson bounds. E001 is HIGH because it needs a new flavor-diagonal CP-odd lepton dipole observable, matching, and likely running.
- CHK-7: No rc1.1 load-bearing number is silently revised by these entries. K002 is consistent with the existing Delta m_K budget, C001 quotes current catalog source values without changing the repo's code constant, and E001 is catalog-only/new to this scope.
- CHK-8: No unresolved `\textbf{CHECK}`, `\cite`, or `\ref` tokens were found in the three TeX files. Math mode/bracing looked consistent in the edited entries; the CA made no TeX changes.
