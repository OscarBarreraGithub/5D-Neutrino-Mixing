# CA Worklog: ca_w5b_charm_edm
**Date**: 2026-05-16
**Family**: charm_edm
**Process IDs**: C006, C008, E007

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| C006 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| C008 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| E007 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| C006 | BR(D0 -> e+- mu-+) | < 1.3e-8 at 90% CL | PDG Live/API S032.40; LHCb 2016 | 87cbe3dec3cf |
| C006 | BR(D0 -> e+- mu-+) | < 1.3e-8 at 90% CL, 3.0 fb^-1 at 7 and 8 TeV | LHCb 2016 arXiv:1512.00322 | caf3d6e3bae8 |
| C006 | BR(D0 -> e+ mu-) + BR(D0 -> mu+ e-) | < 2.6e-7 at 90% CL | Belle 2010 arXiv:1003.2345 | da7d0fe5b963 |
| C006 | BR(D0 -> e mu) | < 3.3e-7 at 90% CL | BaBar 2012 arXiv:1206.5419 | e5042790a01b |
| C008 | BR(D+ -> pi+ e+ mu-) | < 2.1e-7 at 90% CL | PDG Live/API S031.110; LHCb 2021 | fda49cf83fef |
| C008 | BR(D+ -> pi+ e- mu+) | < 2.2e-7 at 90% CL | PDG Live/API S031.111; LHCb 2021 | 118afebb925f |
| C008 | LHCb 2021 search scope and companion limits | 25 modes; 1.6 fb^-1; 95% CL limits 2.3e-7 and 2.2e-7 | LHCb 2021 arXiv:2011.00217 | ccf385471c70 |
| C008 | BaBar predecessor limits | 384 fb^-1; < 2.9e-6 and < 3.6e-6 at 90% CL | BaBar 2011 arXiv:1107.4465 | 7615bb1fbf36 |
| E007 | abs(d(225Ra)) | < 1.4e-23 e cm at 95% CL | Bishof et al. 2016 arXiv:1606.04931 | a9042e38b7b3 |
| E007 | abs(d(225Ra)) first measurement | < 5.0e-22 e cm at 95% CL | Parker et al. 2015 arXiv:1504.07477 | ff5adb5cdce9 |
| E007 | d_A(129Xe) | (1.4 +/- 6.6_stat +/- 2.0_syst)e-28 e cm; < 1.4e-27 e cm at 95% CL | Sachdeva et al. 2019 arXiv:1909.12800 | abce59966ada |
| E007 | d_Xe independent measurement | (-4.7 +/- 6.4)e-28 e cm; < 1.5e-27 e cm at 95% CL | Allmendinger et al. 2019 arXiv:1904.12295 | 94e4abc0ed9f |

## Check evidence
- C006: CHK-1 passed for all observable/source numbers in the TeX: 1.3e-8, 3.0 fb^-1, 7 and 8 TeV, 2.6e-7, and 3.3e-7 all appear in `pdg_or_equivalent` blocks with source URL, access date, snapshot path, and sha256. CHK-2 passed: every process-local source key resolves in `flavor_catalog/references/C006/source_manifest.yaml`; `git ls-files` confirms all snapshot paths are tracked. CHK-3 passed: `ls flavor_catalog/references/C006/` shows only `.txt` snapshots plus `source_manifest.yaml`. CHK-4 passed: `WRITER-INITIATED` precedes `WRITER-DONE`, and `WRITER-DONE` was terminal before this CA update. CHK-5 passed: the prompt's literal `rg -l -E` form errors in this ripgrep because `-E` is an encoding flag, so I reran the intended regex search with `rg -l -e`; non-notebook hits are D0-mixing `mu_had_GeV` false positives, and cited nearby lines `deltaf2.py:252,941`, `modern/phenomenology.py:23,668`, `muToEGamma.py:75`, and `scan.py:524` exist. CHK-6 passed: MEDIUM matches a new Delta C = 1 LFV operator/rate wrapper without new lattice or RG. CHK-7 passed: no rc1.1 load-bearing number is revised. CHK-8 passed: no `\cite`, `\ref`, `TODO`, `CHECK`, or `??` markers found.
- C008: CHK-1 passed for 2.1e-7, 2.2e-7, 1.6 fb^-1, 2016 data, companion 95% CL limits 2.3e-7 and 2.2e-7, 384 fb^-1, 2.9e-6, 3.6e-6, and 25 modes; each is present in `pdg_or_equivalent` with source URL, access date, snapshot path, and sha256. CHK-2 passed: all cited keys resolve in the C008 source manifest with tracked snapshot paths. CHK-3 passed: the reference directory has only `.txt`, `sha256sums.txt`, and `source_manifest.yaml`. CHK-4 passed: `WRITER-INITIATED` precedes terminal `WRITER-DONE` before this CA update. CHK-5 passed: the intended `rg -l -e` coverage search finds no non-notebook C008 implementation, and cited nearby D0-mixing / mu->e gamma lines exist. CHK-6 passed: HIGH matches the new semileptonic mode calculation/form-factor requirement. CHK-7 passed: no rc1.1 load-bearing number is revised. CHK-8 passed: no unresolved LaTeX markers found.
- E007: CHK-1 failed. The measured direct limits in the PDG-equivalent section are traced in `pdg_or_equivalent`, but the TeX also contains numerical projected/context claims at lines 70-84: `1e-28 e cm`, `4e-29 e cm`, `15` day half-life, spin `1/2`, `2--3` orders, `2021--2026`, `1.5e-27 e cm` as a proposal baseline, more than `2` orders, `2023--2027`, and several orders. These trace to YAML `post_2008_context`, not to `pdg_or_equivalent`, so the strict CHK-1 rule is not met. CHK-2 passed: all cited source keys resolve in `flavor_catalog/references/E007/source_manifest.yaml` with tracked snapshot paths. CHK-3 passed: no PDF files are tracked under `references/E007`. CHK-4 passed: `WRITER-INITIATED`, `PKA-DONE`, and `WRITER-DONE` are ordered with ISO timestamps, and `WRITER-DONE` was terminal before this CA update. CHK-5 passed: the intended broad search returns only generic Randall/dipole false positives; a focused Ra/Xe/Schiff/diamagnetic/atomic-EDM search returns no implementation, and the cited `muToEGamma.py:3,21,81` nearby lines exist. CHK-6 passed: HIGH matches the new CP-odd matching/RG and nuclear/atomic response requirement. CHK-7 passed: no rc1.1 load-bearing number is revised. CHK-8 passed: no unresolved LaTeX markers found.

## Issues (if any)
- E007: CHK-1 rework required. Either move the numerical future/context trace blocks that support TeX lines 70-84 into `pdg_or_equivalent` with year/value/source URL/access date/snapshot/sha256 fields, or remove/soften those numeric claims from the TeX in the next WA cycle. I did not edit the TeX per CA hard rules.
