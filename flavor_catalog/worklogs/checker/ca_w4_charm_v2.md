# CA Worklog: ca_w4_charm_v2
**Date**: 2026-05-16
**Family**: charm
**Process IDs**: C002 C003 C004 C007

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| C002 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| C003 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| C004 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| C007 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| C002 | \(|q/p|\) | \(0.983^{+0.015}_{-0.014}\), 95% C.L. [0.955, 1.012] | HFLAV CKM25 all-CPV fit | 5908d28e0388 |
| C002 | \(\phi_D\) | \(-1.51^{+1.03}_{-1.06}\) deg, 95% C.L. [-3.63, 0.51] | HFLAV CKM25 all-CPV fit | 5908d28e0388 |
| C002 | no-indirect-CPV point | \(\Delta\chi^2=2.16\), 0.95 sigma | HFLAV CKM25 all-CPV fit | 5908d28e0388 |
| C002 | PDG cross-check | \(|q/p|=0.983\pm0.015\), \(\phi_D=-1.51\pm1.04\) deg | PDG 2025 D0-D0bar mixing review excerpt | 4ffedbd5dcda |
| C003 | \(\Delta a_{CP}^{\rm dir}\) | \((-0.159\pm0.029)\)% | HFLAV direct/indirect CPV combination | 50781ba7963b |
| C003 | \(a_{CP}^{\rm ind}\) | \((-0.010\pm0.012)\)% | HFLAV direct/indirect CPV combination | 50781ba7963b |
| C003 | no-CPV discrepancy | \(\Delta\chi^2=32.3\), 5.3 sigma | HFLAV direct/indirect CPV combination | 50781ba7963b |
| C003 | \(\Delta A_{CP}\) | \((-15.4\pm2.9)\times10^{-4}=(-0.154\pm0.029)\)% | LHCb 2019 arXiv:1903.08726 | c0d7f9a85048 |
| C003 | \(A_{CP}^{K}-A_{CP}^{\pi}\) | \((-0.154\pm0.029)\)% | HFLAV CKM25 DCPV input | ca6fa2bc4a1c |
| C003 | \(A_{\pi}\) | \((0.225\pm0.057)\)% | HFLAV CKM25 table-results extract | ca6fa2bc4a1c |
| C003 | \(A_K\) | \((0.068\pm0.051)\)% | HFLAV CKM25 table-results extract | ca6fa2bc4a1c |
| C004 | \(\mathcal B(D^0\to\mu^+\mu^-)\) | \(<2.1\times10^{-9}\) at 90% C.L. | PDG Live/API S032.28 | 34cff61f22c2 |
| C004 | \(\mathcal B(D^0\to\mu^+\mu^-)\) | \(<2.4\times10^{-9}\) at 95% C.L.; 64.5 fb\(^{-1}\), 13.6 TeV, 2022-2023 | CMS 2025 public result | 46aaaf8dc7fb |
| C004 | LHCb predecessor limit | \(<3.1\times10^{-9}\) at 90% C.L.; companion \(<3.5\times10^{-9}\) at 95% C.L.; 9 fb\(^{-1}\) at 7, 8, 13 TeV | LHCb 2023 arXiv:2212.11203 | a94f1dbcfa8d |
| C004 | long-distance context | \(2.7\times10^{-5}\mathcal B(D^0\to\gamma\gamma)\); floor \(3\times10^{-13}\); VMD \(D^0\to\gamma\gamma=(3.5^{+4.0}_{-2.6})\times10^{-8}\) | Burdman et al. 2002 | ec7271413986 |
| C004 | rare-charm EFT context | \(C_{7,9,10}^{(\prime)}\) Wilson coefficients discussed | Gisbert, Hiller, Suelmann 2024 | 7432c6e611db |
| C007 | \(\mathcal B(D^+\to\pi^+\mu^+\mu^-)\) | \(<6.7\times10^{-8}\) at 90% C.L. | PDG Live/API S031.42 | 1bfa609f88a7 |
| C007 | LHCb current result | \(<6.7\times10^{-8}\) at 90% C.L.; companion \(<7.4\times10^{-8}\) at 95% C.L.; 1.6 fb\(^{-1}\), 2016 | LHCb 2021 arXiv:2011.00217 | 7e293a7a298e |
| C007 | LHCb 2021 search scope | 25 rare or forbidden \(D^+\) and \(D_s^+\) modes | LHCb 2021 arXiv:2011.00217 | 7e293a7a298e |
| C007 | short-distance SM scale | order \(10^{-12}\) branching fractions | LHCb 2021 arXiv:2011.00217 | 7e293a7a298e |
| C007 | LHCb predecessor limit | \(<7.3(8.3)\times10^{-8}\) at 90% (95%) C.L.; 1.0 fb\(^{-1}\) at 7 TeV | LHCb 2013 arXiv:1304.6365 | 64daf4ae3760 |

## Evidence notes
- C002: All numerical TeX claims are covered by `q_over_p`, `phi_D`, `no_indirect_cpv_test`, and `pdg_cross_check` YAML blocks. Manifest keys `hflav_ckm25_dmixing_cpv_fit`, `pdg2025_dmixing_review_cpv`, `lhcb2023_arxiv2208_06512`, `belle_belleii2025_arxiv2410_22961`, and `cfw2008_arxiv0804_1954` all have tracked `.txt` snapshots and matching local hashes. Code-coverage citations resolve at `quarkConstraints/deltaf2.py:252`, `quarkConstraints/deltaf2.py:941`, `quarkConstraints/scan.py:103`, `quarkConstraints/scan.py:439`, `quarkConstraints/modern/phenomenology.py:46`, and `quarkConstraints/paper_0710_1869/eft_deltaf2/observables.py:1403`; focused grep found no exported \(|q/p|\) or \(\phi_D\) evaluator.
- C003: All numerical TeX claims are covered by `hflav_direct_cpv_average`, `hflav_indirect_cpv_companion`, `hflav_no_cpv_test`, `lhcb2019_discovery_average`, `hflav_ckm25_delta_acp_input`, `hflav_ckm25_Api_fit`, and `hflav_ckm25_AK_fit`. Manifest keys all resolve to tracked `.txt` snapshots with matching hashes. Focused direct-CP grep returned only unrelated neutrino `delta_CP` symbols and no C003 observable; the `NO` coverage claim is consistent. HIGH difficulty is consistent with the need for new \(\Delta C=1\) nonleptonic/RG/hadronic treatment.
- C004: All numerical TeX claims are covered by `canonical_current_limit`, `cms_current_95cl_limit`, `lhcb_previous_limit`, `standard_model_long_distance_context`, and `post_2008_eft_context`. Manifest keys all resolve to tracked `.txt` snapshots with matching hashes. Focused D0-dimuon/rare-charm grep returned no live \(D^0\to\mu^+\mu^-\) implementation; nearby hits are D0 mixing or unrelated \(\mu\to e\gamma\). HIGH difficulty is now consistent with the rubric because the mode needs a new \(\Delta C=1\) rare-leptonic calculation and Wilson normalization.
- C007: WA-v2 added the missing `lhcb_2021_search_scope` and `lhcb_2021_short_distance_sm_scale` blocks, so all numerical TeX claims now trace to YAML. Manifest keys all resolve to tracked `.txt` snapshots with matching hashes. Focused rare-semileptonic charm grep returned no \(D^+\to\pi^+\mu^+\mu^-\) implementation; only unrelated notebook/QCD README noise appeared in broad searches. HIGH difficulty is consistent with a new semileptonic \(\Delta C=1\) operator/mode calculation.
- For all four processes: `ls flavor_catalog/references/<process_id>/` showed no `.pdf` files. `status_history` contains `WRITER-INITIATED` before `WRITER-DONE` with ISO 8601 timestamps, and the pre-check most recent entry was `WRITER-DONE`. No `\textbf{CHECK}`, `TODO`, `FIXME`, `\ref`, `\eqref`, `\cite`, or unresolved-reference marker was present in the TeX files. Targeted rc1.1/repo searches found only D0-mixing companion-code material and catalog plan rows, not contradictory load-bearing catalog numbers.

## Verification command note
- The literal form `rg -l -E "<process keyword>" ...` fails with this installed ripgrep because `-E` is the encoding flag. I ran that form once and then used the equivalent ripgrep regex form `rg -l -e "<process keyword>" ... 2>&1`, plus focused `rg -n -g '!*.ipynb'` checks and line-level reads for cited code locations.

## Issues (if any)
- None.
