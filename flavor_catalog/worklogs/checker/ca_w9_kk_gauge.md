# CA Worklog: ca_w9_kk_gauge
**Date**: 2026-05-17
**Family**: collider_rs, Wave-9 PRIMARY
**Batch ID**: CA-w9-kk-gauge
**Cycle**: 1
**Checker agent**: CA-CA-w9-kk-gauge
**Process IDs**: CR001 CR005 CR006 CR007

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| CR001 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| CR005 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| CR006 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| CR007 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| CR001 | RS KK-gluon excluded mass interval, pp -> g_KK -> t tbar | 0.5 < m(g_KK) < 5.5 TeV excluded at 95% CL | CMS-B2G-25-009, arXiv:2603.23454 | c8991eda1916 |
| CR001 | Historical CMS RS KK-gluon excluded mass interval | 0.50 < m(g_KK) < 4.55 TeV excluded at 95% CL | CMS-B2G-17-017, arXiv:1810.05905 | 93a713d209e4 |
| CR001 | PDG Live S071KKG ATLAS review-table lower bound | m(g_KK) > 3.8 TeV at 95% CL | PDG Live S071KKG | b0d08fb66567 |
| CR001 | ATLAS 7 TeV lepton+jets KK-gluon exclusion | m(g_KK) < 2.07 TeV excluded at 95% credibility level | ATLAS TOPQ-2012-14, arXiv:1305.2756 | 0f2049c8120 |
| CR001 | ATLAS boosted 7 TeV KK-gluon exclusion | m(g_KK) < 1.5 TeV excluded at 95% credibility level | ATLAS boosted ttbar, arXiv:1207.2409 | 79a0e8ade86 |
| CR005 | Z'_SSM -> ee + mumu observed lower mass limit | m(Z'_SSM) > 5.15 TeV at 95% CL | CMS high-mass dilepton, arXiv:2103.02708 | c33742062e29 |
| CR005 | Z'_psi -> ee + mumu observed lower mass limit | m(Z'_psi) > 4.56 TeV at 95% CL | CMS high-mass dilepton, arXiv:2103.02708 | c33742062e29 |
| CR005 | ATLAS Z'_SSM combined dilepton lower mass limit | m(Z'_SSM) > 5.1 TeV at 95% CL | ATLAS high-mass dilepton, arXiv:1903.06248 | eacbb1bc7c23 |
| CR005 | ATLAS fiducial zero-width dilepton limit at 6 TeV | sigma_fid x B < about 0.014 fb at 95% CL | ATLAS high-mass dilepton, arXiv:1903.06248 | eacbb1bc7c23 |
| CR005 | PDG prompt-dilepton Z' cross-section summary limit | cross-section upper limits as low as 0.02 fb at 95% CL | PDG 2024 Z'-Boson Searches review | f12a1d992f09 |
| CR006 | W'_SSM -> e nu lower mass limit | M(W'_SSM) > 6.0 TeV at 95% CL | PDG 2025 W' searches, quoting ATLAS arXiv:1906.05609 | 78c2ad8e5581 |
| CR006 | W'_SSM -> mu nu lower mass limit | M(W'_SSM) > 5.6 TeV at 95% CL | PDG 2025 W' searches, quoting CMS arXiv:2202.06075 | 78c2ad8e5581 |
| CR006 | CMS combined e/mu W'_SSM lower mass limit | M(W'_SSM) > 5.7 TeV at 95% CL | CMS lepton plus missing-pT, arXiv:2202.06075 | 7f341586575 |
| CR006 | CMS W'_R -> tb lower mass limit | M(W'_R) > 4.3 TeV at 95% CL | CMS W' -> tb leptonic, arXiv:2310.19893 | b699d14ea827 |
| CR006 | CMS W'_L -> tb lower mass limit | M(W'_L) > 3.9 TeV at 95% CL | CMS W' -> tb leptonic, arXiv:2310.19893 | b699d14ea827 |
| CR007 | Generic RS graviton diphoton lower mass limit | m(G_KK^(1)) > 4.8 TeV at 95% CL | PDG 2025 Extra Dimensions listing / CMS 2024 diphoton | 9d9619662dbd |
| CR007 | Generic RS graviton dilepton lower mass limit | m(G_KK^(1)) > 4.78 TeV at 95% CL | PDG 2025 Extra Dimensions listing / CMS 2021 dilepton | 9d9619662dbd |
| CR007 | Bulk-RS graviton all-jets WW/ZZ lower mass limit | m(G_bulk) > 1.4 TeV at 95% CL | PDG 2025 Extra Dimensions listing / CMS 2023 all-jets | 9d9619662dbd |
| CR007 | ATLAS combined bulk-RS graviton lower mass limit | m(G_KK) > 2.3 TeV at 95% CL | ATLAS combined bosonic/leptonic search, arXiv:1808.02380 | 1ae94abf395 |

## Issues
- CR001: None.
- CR005: None.
- CR006: None.
- CR007: None.

## Evidence notes
- CHK-1: PASS. Measured collider limits in the TeX are present in `pdg_or_equivalent.values` with strict metadata: value, uncertainty/null limit form, units, CL, year, source URL, access date, snapshot path, and sha256. I applied the L001/B001_B003/B021_B023 carve-out: dataset luminosities, sqrt(s), event yields, search mass ranges, RS theory reach/projections, EFT translations, and local quark-scan comparison scales remain in `supporting_measurements`, `paper_era_reference`, or `auxiliary_theory_inputs`. The CR005 ATLAS `0.014 fb` and PDG `0.02 fb` measured cross-section upper limits are correctly promoted.
- CHK-2: PASS. Process-local TeX key lists resolve to the corresponding `source_manifest.yaml` entries for CR001, CR005, CR006, and CR007. Every manifest entry has a non-empty `snapshot_path`, and `git ls-files --error-unmatch` plus non-empty file checks found no untracked or missing manifest snapshots. `sha256sum -c` passed for all four reference directories.
- CHK-3: PASS. `find flavor_catalog/references/CR001 flavor_catalog/references/CR005 flavor_catalog/references/CR006 flavor_catalog/references/CR007 -type f \( -iname '*.pdf' -o -iname '*.ps' \)` returned no files. The directories contain text/JSON snapshots, manifests, and checksum lists only.
- CHK-4: PASS. Before checker edits, all four sidecars had legal `DRAFT -> WRITER-INITIATED -> WRITER-DONE` transitions with ISO 8601 `timestamp`/`at` fields. The latest pre-check state was `WRITER-DONE`, `agent_id: "WA"`, `cycle: 1`, at `2026-05-17T17:08:57-04:00`.
- CHK-5: PASS. All four sidecars claim `code_coverage.status: "NO"`, which matches the grep evidence. Broad collider/recast searches over `quarkConstraints`, `qcd`, `flavorConstraints`, `neutrinos`, `yukawa`, `warpConfig`, `solvers`, `scanParams`, and `tests` found no direct ATLAS/CMS limit ingestion, likelihood, HEPData constraint layer, or CheckMATE/MadAnalysis5/SModelS path. KK-gluon hits are low-energy Delta F=2/coupling infrastructure (`quarkConstraints/couplings.py:1`, `quarkConstraints/deltaf2.py:1`, `quarkConstraints/modern/matching.py:241-243`), not CR001 collider logic. The collider-looking false positive is `tests/test_alpha_s.py:89`, a CMS/RunDec alpha_s calibration comment. Cited adjacent scan lines exist: `quarkConstraints/scan.py:359-380`, `quarkConstraints/scan.py:418-439`, `scanParams/scan.py:524-535`, `quarkConstraints/modern/scan.py:1229-1255`, and `quarkConstraints/modern/evaluation.py:643-666`.
- CHK-6: PASS. `HIGH` is consistent for each collider-RS entry. CR001 needs ttbar resonance recasting with top tagging, widths, branching fractions, and interference; CR005 needs neutral electroweak-KK production/couplings and dilepton likelihood/HEPData reinterpretation; CR006 needs charged-vector production, tb/lnu branching, finite-width/interference, and detector acceptance; CR007 needs spin-2 mass/width/cross-section/branching predictions and diboson/diphoton/dilepton recasts. These are collider reinterpretations, not simple additions to an existing low-energy operator basis.
- CHK-7: PASS. The entries do not contradict rc1.1 or methodology-note load-bearing numbers. `docs/quark_scan_methodology_note.tex:587-588`, `:703-705`, and `:1170-1174` support the quoted `47.26 TeV` p50 and `127.13 TeV` p95 low-energy flavor comparison scales at `g_s^* = 3`. The TeX correctly states that current direct LHC limits are weaker than the anarchic-flavor result and are collider cross-checks or non-anarchic diagnostics.
- CHK-8: PASS. `rg` found no unresolved `\cite`, `\ref`, `\textbf{CHECK}`, `CHECK`, `TODO`, `FIXME`, `MISSING`, or `??` markers in the four TeX files. Lightweight delimiter counts were balanced: CR001 `\(`/`\)` = 41/41, `\[`/`\]` = 0/0, braces = 91/91; CR005 35/35, 2/2, braces = 77/77; CR006 58/58, 4/4, braces = 74/74; CR007 29/29, 3/3, braces = 80/80. Dollar count was zero in all four.
