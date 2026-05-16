# CA Worklog: ca_w6_beauty_v2
**Date**: 2026-05-16
**Family**: beauty
**Process IDs**: B001, B003, B016

Cycle-2 interpretation note: CHK-1 is enforced strictly as written in the
dispatch prompt. Numerical TeX claims for physics/source values must trace to
`pdg_or_equivalent` entries, not only to other sidecar sections.

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| B001 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| B003 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| B016 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| B001 | Delta m_d | 0.5069 +/- 0.0019 ps^-1 | HFLAV/PDG 2025, accessed 2026-05-16, `hflav_pdg2025_bd_mixing.txt` | 5427ddc287e4 |
| B001 | x_d | 0.7697 +/- 0.0035 | HFLAV/PDG 2025, accessed 2026-05-16, `hflav_pdg2025_bd_mixing.txt` | 5427ddc287e4 |
| B001 | chi_d | 0.1860 +/- 0.0011 | HFLAV/PDG 2025, accessed 2026-05-16, `hflav_pdg2025_bd_mixing.txt` | 5427ddc287e4 |
| B001 | CFW RS-flavor context scales | 21 TeV standard RS; 33 TeV pseudo-Goldstone scenario | Csaki-Falkowski-Weiler 2008 / arXiv:0804.1954, accessed 2026-05-16, `cfw_2008_rs_flavor_arxiv.txt` | cf93d80c47eb |
| B001 | Belle II Delta m_d and dataset | 0.516 +/- 0.008(stat) +/- 0.005(syst) ps^-1; 190 fb^-1 | Belle II 2023 / arXiv:2302.12791, accessed 2026-05-16, `belleii_2023_bd_mixing_arxiv.txt` | 10ad4816deaa |
| B003 | Delta m_s, HFLAV recommended | 17.766 +/- 0.006 ps^-1 | HFLAV Fall 2024, accessed 2026-05-16, `hflav_2024_dms.txt` | 21d1748c46c3 |
| B003 | HFLAV companion quantities | x_s = 26.94 +/- 0.10; chi_s = 0.499314 +/- 0.000005; 1/Gamma_s = 1.516 +/- 0.006 ps; DeltaGamma_s = +0.0781 +/- 0.0035 ps^-1 | HFLAV Fall 2024, accessed 2026-05-16, `hflav_2024_dms.txt` | 21d1748c46c3 |
| B003 | Delta m_s, HFLAV/PDG 2025 cross-check | 17.765 +/- 0.006 ps^-1 | HFLAV/PDG 2025, accessed 2026-05-16, `hflav_pdg2025_dms.txt` | 7b6e88fed6d4 |
| B003 | LHCb Delta m_s and dataset | 17.7683 +/- 0.0051(stat) +/- 0.0032(syst) ps^-1; LHCb combination 17.7656 +/- 0.0057 ps^-1; 6 fb^-1 | LHCb 2022 / arXiv:2104.04421, accessed 2026-05-16, `lhcb_2021_deltams_arxiv.txt` | 2827716498d9 |
| B003 | FLAG B_s mixing inputs | f_Bs sqrt(Bhat_Bs) = 256.1 +/- 5.7 MeV; xi = 1.216 +/- 0.016 | FLAG Review 2024 / arXiv:2411.04268, accessed 2026-05-16, `flag_2024_bmixing.txt` | 0d34c342c10b |
| B016 | BR(B+ -> K+ ell+ ell-) | (5.76 +/- 0.40) x 10^-7 | HFLAV rare B Dec. 2025, accessed 2026-05-16, `hflav_dec2025_Bplus_to_Kplus_ll.txt` | d57204eb4634 |
| B016 | BR(B0 -> K0 ell+ ell-) | (3.28 +/- 0.32) x 10^-7 | HFLAV rare B Dec. 2025, accessed 2026-05-16, `hflav_dec2025_B0_to_K0_ll.txt` | 44551ba39055 |

## Issues (if any)
- B001: CHK-1 fails. `B001.tex` quotes the CFW `21 TeV` and `33 TeV` RS-flavor context values plus the Belle II `Delta m_d = 0.516 +/- 0.008(stat) +/- 0.005(syst) ps^-1` and `190 fb^-1` dataset. WA-v2 added complete metadata for these values under `paper_era_reference` and `recent_experimental_inputs`, but the strict CA prompt requires a `pdg_or_equivalent` entry for every numerical TeX claim.
- B003: CHK-1 fails. The HFLAV companion quantities and LHCb Run-2 values now live under `pdg_or_equivalent` with complete metadata, but `B003.tex` still quotes the FLAG `f_Bs sqrt(Bhat_Bs) = 256.1 +/- 5.7 MeV` and `xi = 1.216 +/- 0.016` auxiliary theory inputs outside `pdg_or_equivalent`.

## Verification notes
- CHK-1: `rg -n "[0-9]"` over the three TeX files identified the value-bearing numerical claims in the table. Snapshot searches verified the quoted values, and `sha256sum` matched sidecar `source_shas`; B001 and B003 still fail only because some quoted numbers are not in `pdg_or_equivalent`.
- CHK-2: every process-local source key named in the TeX files resolves to `flavor_catalog/references/<process_id>/source_manifest.yaml`; each manifest entry has a non-empty `snapshot_path`, and `git ls-files` confirms those snapshots are tracked.
- CHK-3: `ls` and `find ... -name '*.pdf'` under `references/B001`, `references/B003`, and `references/B016` found no PDFs. The directories contain text snapshots, manifests, and B003's local `sha256sums.txt`.
- CHK-4: before this CA-v2 update, all three sidecars contained `WRITER-INITIATED` followed by `WRITER-DONE` with ISO 8601 timestamps, and the most recent status was the WA-v2 `WRITER-DONE`.
- CHK-5: the prompt's literal `rg -l -E "<process keyword>" ...` form errors in this ripgrep because `-E` is parsed as an encoding flag, so I reran equivalent `rg -l -e` searches over the required directories. B001 and B003 `YES` claims match live Delta-F=2 B_d/B_s code and cited lines exist. B016 `NO` matches targeted non-notebook searches for `B_to_K`, `BtoK`, `B -> K`, `B\\s*to\\s*K`, `b.?->.?s.?ell`, `C9`, `C10`, word-boundary `R_K`, `Bplus`, and `Kplus`, which returned no implementation.
- CHK-6: B001 and B003 `LOW` are consistent because existing Delta-F=2 SLL/SLR/VLL/VRR/LR1/LR2 machinery covers neutral B mixing. B016 `HIGH` is consistent because it requires a new Delta-B=1 semileptonic Hamiltonian, q2-bin/form-factor/nonlocal-charm treatment, and covariance support.
- CHK-7: no load-bearing rc1.1 contradiction found. The catalog entries are companion-mode; B001/B003 align with existing B_d/B_s implementation context, and B016 remains catalog-only.
- CHK-8: no `\textbf{CHECK}`, `TODO`, `FIXME`, `\ref{}`, `\cite{}`, or `??` markers were found; brace-balance checks over all three TeX files returned zero.
