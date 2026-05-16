# CA Worklog: ca_wave1_beauty
**Date**: 2026-05-16
**Family**: beauty
**Process IDs**: B009, B011, B015

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| B009 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| B011 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| B015 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| B009 | BR(B+ -> tau+ nu), HFLAV average | (1.12 +/- 0.19) x 10^-4 | HFLAV End Dec. 2025 | 7a4e4079a822 |
| B009 | BR(B+ -> tau+ nu), PDG average | (1.09 +0.25/-0.24) x 10^-4 | PDG 2025 REST S041.184 | 7e0c8745233b |
| B009 | BR(B+ -> tau+ nu), Belle hadronic tag | (0.72 +0.27/-0.25 +/- 0.11) x 10^-4 | Belle 2013 / HFLAV | 34fba2f4ff09 |
| B009 | BR(B+ -> tau+ nu), Belle semileptonic tag | (1.25 +/- 0.28 +/- 0.27) x 10^-4 | Belle 2015 / HFLAV | aa5550faa3fc |
| B009 | BR(B+ -> tau+ nu), Belle II hadronic tag | (1.24 +/- 0.41 +/- 0.19) x 10^-4 | Belle II 2025 / HFLAV | 8aed99257179 |
| B009 | BR(B+ -> tau+ nu), BaBar hadronic tag | (1.83 +0.53/-0.49 +/- 0.24) x 10^-4 | BaBar 2013 / HFLAV | eb7d51abe66f |
| B009 | BR(B+ -> tau+ nu), BaBar semileptonic tag | (1.70 +/- 0.80 +/- 0.20) x 10^-4 | BaBar 2010 / HFLAV | 0e3378f2fc0e |
| B009 | BR(B -> tau nu), SM prediction | (0.865 +/- 0.041) x 10^-4 | UTfit Summer 2024 | c025b53f3641 |
| B011 | BR(Bbar -> X_s gamma), exp | (3.49 +/- 0.19) x 10^-4, E_gamma > 1.6 GeV | HFLAV 2024/2026 snapshot | 81dde0df2037 |
| B011 | HFLAV table-unit average | 349 +/- 19 in 10^-6 units | HFLAV Dec. 2024 Table 67 snapshot | 43581f22065d |
| B011 | BR(Bbar -> X_s gamma), SM | (3.40 +/- 0.17) x 10^-4, E_gamma > 1.6 GeV | Misiak-Rehman-Steinhauser 2020 | 8a3ec66a22c3 |
| B011 | BR(b -> s gamma), PDG-reviewed SM quote | (3.36 +/- 0.23) x 10^-4, E_gamma >= 1.6 GeV | Misiak et al. 2015 / PDG 2024 | c252dbcd2b42 |
| B015 | BR(B -> X_s ell+ ell-), HFLAV average | (5.84 +/- 0.69) x 10^-6 | HFLAV Dec. 2024 | 94a3f65a6240 |
| B015 | BR(B -> X_s ell+ ell-), PDG listed | (5.84 +1.31/-1.23) x 10^-6 | HFLAV Dec. 2024 table | 94a3f65a6240 |
| B015 | BR(B -> X_s ell+ ell-), BaBar input | (6.73 +0.70/-0.64 +0.60/-0.56) x 10^-6 | BaBar 2014 / HFLAV | fcaba867afcd |
| B015 | BR(B -> X_s ell+ ell-), Belle input | (4.11 +/- 0.83 +0.85/-0.81) x 10^-6 | Belle 2005 / HFLAV | 1f51dd799054 |
| B015 | Low-q2 weighted experimental average | (1.58 +/- 0.37) x 10^-6 for 1 < q^2 < 6 GeV^2 | Huber-Hurth-Lunghi 2015 | 735177a8007f |
| B015 | High-q2 weighted experimental average | (0.48 +/- 0.10) x 10^-6 | Huber-Hurth-Lunghi 2015 | 735177a8007f |
| B015 | Belle II projection luminosity | 50 ab^-1 | Huber-Hurth-Lunghi 2015 | 735177a8007f |

## Issues (if any)
- B009: CHK-1 fails. `B009.tex` lines 28-32 list five post-2008 experimental-input BR values. The corresponding `pdg_or_equivalent.experimental_inputs` blocks in `B009.yaml` lines 77-122 have year, value, uncertainty, snapshot, and sha256, but do not carry the required `source_url` and `access_date` fields for each numerical claim.
- B015: CHK-1 fails. `B015.tex` lines 24-35 list the HFLAV/PDG average, BaBar/Belle inputs, and low/high-q2 bin values. The corresponding `B015.yaml` observable blocks in lines 44-87 have value, uncertainty, snapshot, and sha256, but not per-observable year/source_url/access_date metadata; the top-level block also has no `source_url`.
- B015: CHK-1 also fails for the Belle II `50 ab^-1` projection in `B015.tex` line 51. The value appears in the Huber-Hurth-Lunghi snapshot, but no `pdg_or_equivalent` observable records that projection/luminosity.
- Coverage-grep note: the literal requested `rg -l -E "<process keyword>" ...` form fails in this environment because ripgrep uses `-E` as an encoding flag. I reran the process searches with `rg -l -e`; only notebooks matched, and rerunning with `-g '!*.ipynb'` returned no implementation hits in the required code/test directories.
