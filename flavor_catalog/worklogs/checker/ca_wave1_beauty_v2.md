# CA Worklog: ca_wave1_beauty_v2
**Date**: 2026-05-16
**Family**: beauty
**Process IDs**: B009, B015

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| B009 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| B015 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| B009 | HFLAV average BR(B+ -> tau+ nu_tau) | (1.12 +/- 0.19) x 10^-4 | HFLAV End Dec. 2025; access 2026-05-16; `flavor_catalog/references/B009/hflav_dec2025_btaunu.txt` | 7a4e4079a822 |
| B009 | PDG REST average BR(B+ -> tau+ nu_tau) | (1.09 +0.25/-0.24) x 10^-4 | PDG 2025 REST S041.184; access 2026-05-16; `flavor_catalog/references/B009/pdg2025_btaunu_api.txt` | 7e0c8745233b |
| B009 | Belle hadronic-tag input | (0.72 +0.27/-0.25 +/- 0.11) x 10^-4 | Belle 2013 arXiv:1208.4678; access 2026-05-16; `flavor_catalog/references/B009/belle2013_arxiv1208_4678.txt` | 34fba2f4ff09 |
| B009 | Belle semileptonic-tag input | (1.25 +/- 0.28 +/- 0.27) x 10^-4 | Belle 2015 arXiv:1503.05613; access 2026-05-16; `flavor_catalog/references/B009/belle2015_arxiv1503_05613.txt` | aa5550faa3fc |
| B009 | Belle II hadronic-tag input | (1.24 +/- 0.41 +/- 0.19) x 10^-4 | Belle II 2025 arXiv:2502.04885; access 2026-05-16; `flavor_catalog/references/B009/belleii2025_arxiv2502_04885.txt` | 8aed99257179 |
| B009 | BaBar hadronic-tag input | (1.83 +0.53/-0.49 +/- 0.24) x 10^-4 | BaBar 2013 arXiv:1207.0698; access 2026-05-16; `flavor_catalog/references/B009/babar2013_arxiv1207_0698.txt` | eb7d51abe66f |
| B009 | BaBar semileptonic-tag input | (1.70 +/- 0.80 +/- 0.20) x 10^-4 | BaBar 2010 arXiv:0912.2453; access 2026-05-16; `flavor_catalog/references/B009/babar2010_arxiv0912_2453.txt` | 0e3378f2fc0e |
| B009 | UTfit SM prediction | (0.865 +/- 0.041) x 10^-4 | UTfit Summer 2024 SM; access 2026-05-16; `flavor_catalog/references/B009/utfit_summer2024_btaunu.txt` | c025b53f3641 |
| B015 | HFLAV average BR(B -> X_s ell+ ell-) | (5.84 +/- 0.69) x 10^-6 | HFLAV End Dec. 2024; access 2026-05-16; `flavor_catalog/references/B015/hflav_dec2024_B_to_Xsll.txt` | 94a3f65a6240 |
| B015 | PDG-listed BR(B -> X_s ell+ ell-) | (5.84 +1.31/-1.23) x 10^-6 | HFLAV End Dec. 2024 table; access 2026-05-16; `flavor_catalog/references/B015/hflav_dec2024_B_to_Xsll.txt` | 94a3f65a6240 |
| B015 | BaBar input BR(B -> X_s ell+ ell-) | (6.73 +0.70/-0.64 +0.60/-0.56) x 10^-6 | BaBar 2014 arXiv:1312.5364 and HFLAV table; access 2026-05-16; `flavor_catalog/references/B015/babar_1312_5364_arxiv.txt` | fcaba867afcd |
| B015 | Belle input BR(B -> X_s ell+ ell-) | (4.11 +/- 0.83 +0.85/-0.81) x 10^-6 | Belle 2005 arXiv:hep-ex/0503044; access 2026-05-16; `flavor_catalog/references/B015/belle_hepex_0503044_arxiv.txt` | 1f51dd799054 |
| B015 | Experimental low-q2 weighted average | (1.58 +/- 0.37) x 10^-6 for 1 < q^2 < 6 GeV^2 | Huber-Hurth-Lunghi 2015 arXiv:1503.04849; access 2026-05-16; `flavor_catalog/references/B015/huber_hurth_lunghi_1503_04849_arxiv.txt` | 735177a8007f |
| B015 | Experimental high-q2 weighted average | (0.48 +/- 0.10) x 10^-6 | Huber-Hurth-Lunghi 2015 arXiv:1503.04849; access 2026-05-16; `flavor_catalog/references/B015/huber_hurth_lunghi_1503_04849_arxiv.txt` | 735177a8007f |
| B015 | Belle II projection luminosity | 50 ab^-1 | Huber-Hurth-Lunghi 2015 arXiv:1503.04849; access 2026-05-16; `flavor_catalog/references/B015/huber_hurth_lunghi_1503_04849_arxiv.txt` | 735177a8007f |

## Evidence notes
- B009 CHK-1: all branching-fraction values and the UTfit SM prediction quoted in `B009.tex` have matching `pdg_or_equivalent` blocks with year, value, uncertainty, source URL, access date, snapshot path, and sha256. Local `sha256sum` output matched `source_shas`.
- B009 CHK-2: every process-local snapshot named in `B009.tex` resolves to `flavor_catalog/references/B009/source_manifest.yaml`; all manifest `snapshot_path` entries are tracked text snapshots.
- B009 CHK-3: `ls flavor_catalog/references/B009/` showed only `.txt` snapshots and `source_manifest.yaml`; no `.pdf` files are present or tracked.
- B009 CHK-4: pre-CA `status_history` contained `WRITER-INITIATED` followed by cycle-2 `WRITER-DONE`, with ISO 8601 timestamps; cycle-2 `WRITER-DONE` was the latest pre-CA entry.
- B009 CHK-5: `B009.tex` says coverage is NO. The prompt's `rg -E` form is not valid for ripgrep 15.1.0 because `-E` is the encoding flag, so I reran the equivalent `rg -l -e` search over the required directories. Without notebook exclusion, broad terms hit only `.ipynb` embedded output/noise; with `-g '!*.ipynb'`, the B+ -> tau nu / b -> u / charged-current search returned no source hits.
- B009 CHK-6: MEDIUM is consistent with the rubric because promotion needs a new charged-current observable and model matching, but not new lattice/RG machinery at catalog level.
- B009 CHK-7: no load-bearing rc1.1 paper number is revised; the CFW 2008 snapshot is used only as baseline context.
- B009 CHK-8: no `\textbf{CHECK}`, `TODO`, `\cite`, or `\ref` markers; visual TeX inspection found balanced display math and acceptable math mode.
- B015 CHK-1: all branching-fraction values, q2-bin weighted averages, and the Belle II `50 ab^-1` projection quoted in `B015.tex` have matching `pdg_or_equivalent.observables` blocks with year, value, uncertainty, source URL, access date, snapshot path, and sha256. Local `sha256sum` output matched `source_shas`.
- B015 CHK-2: every process-local source key or snapshot named in `B015.tex` resolves to `flavor_catalog/references/B015/source_manifest.yaml`; all manifest `snapshot_path` entries are tracked text snapshots.
- B015 CHK-3: `ls flavor_catalog/references/B015/` showed only `.txt` snapshots, `sha256sums.txt`, and `source_manifest.yaml`; no `.pdf` files are present or tracked.
- B015 CHK-4: pre-CA `status_history` contained `WRITER-INITIATED` followed by cycle-2 `WRITER-DONE`, with ISO 8601 timestamps; cycle-2 `WRITER-DONE` was the latest pre-CA entry.
- B015 CHK-5: `B015.tex` says coverage is NO. The equivalent `rg -l -e` search for b -> s / X_s / ell+ ell- / C7/C9/C10 over the required directories returned only notebook matches before exclusion and no source hits with `-g '!*.ipynb'`; cited lines `quarkConstraints/modern/phenomenology.py:23`, `quarkConstraints/modern/phenomenology.py:166`, and `flavorConstraints/muToEGamma.py:75` exist and match the stated context.
- B015 CHK-6: HIGH is consistent with the rubric because production use needs a new Delta B = 1 semileptonic operator basis, matching/running, bin-aware predictions, charm/QED treatment, and a likelihood/covariance model.
- B015 CHK-7: no load-bearing rc1.1 paper number is revised; the CFW 2008 snapshot is used only as baseline context.
- B015 CHK-8: no `\textbf{CHECK}`, `TODO`, `\cite`, or `\ref` markers; visual TeX inspection found balanced display math and acceptable math mode.

## Issues (if any)
- None.
