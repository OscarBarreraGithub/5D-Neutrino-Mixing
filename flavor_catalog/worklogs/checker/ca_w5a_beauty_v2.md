# CA Worklog: ca_w5a_beauty_v2
**Date**: 2026-05-16
**Family**: beauty
**Process IDs**: B021, B022, B023

Cycle-2 re-check after `wa_w5a_beauty_v2`. I enforced the same CHK-1
interpretation as `ca_w5a_beauty.md`: numerical TeX claims must trace to
`pdg_or_equivalent`; values kept only in supporting/theory blocks are rework.

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| B021 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| B022 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| B023 | PASS | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| B021 | BR(Lambda_b0 -> Lambda mu+ mu-) | (1.08 +/- 0.28)e-6 | PDG Live 2025 | 73959783ca3c |
| B021 | dBR/dq2, 15 < q2 < 20 GeV^2/c^4 | (1.18 +0.09/-0.08 +/-0.03 +/-0.27)e-7 (GeV^2/c^4)^-1 | LHCb 2015 | 2feb8efe212b |
| B021 | A_FB^l, A_FB^h | -0.05 +/-0.09 +/-0.03; -0.29 +/-0.07 +/-0.03 | LHCb 2015 | 2feb8efe212b |
| B021 | LHCb dataset size quoted in TeX | 3.0 fb^-1 | LHCb 2015 | 2feb8efe212b |
| B021 | CDF observation claims | 24 events; 5.8 sigma; BR=(1.73 +/-0.42 +/-0.55)e-6 | CDF 2011 | 83084eabd1d4 |
| B022 | BR(B+ -> K+ nu nubar), HFLAV average | (1.38 +/- 0.35)e-5 | HFLAV Dec. 2025 | 4c562ccd073f |
| B022 | BR(B+ -> K+ nu nubar), PDG display | (2.30 +0.71/-0.64)e-5 | PDG Live/API 2026 | af37fb365951 |
| B022 | Belle II evidence measurement | (2.3 +/-0.5 stat +0.5/-0.4 syst)e-5; 362 fb^-1; 3.5 sigma; 2.7 sigma above SM | Belle II 2024 | 5e8e54e24cbc |
| B022 | BaBar hadronic-tag limit | <3.7e-5 at 90% CL | BaBar 2013 | 8366c09006d7 |
| B022 | SM BR(B+ -> K+ nu nubar) | (5.58 +/- 0.37)e-6 | HPQCD 2023 | 66d8ba9efc52 |
| B023 | BR(B0 -> K*0 nu nubar) | <1.8e-5 at 90% CL | PDG/HFLAV 2025 | f8d7014f1e7f |
| B023 | BR(B+ -> K*+ nu nubar) | <4.0e-5 at 90% CL | PDG/HFLAV 2025 | d9aa9db63226 |
| B023 | Combined BR(B -> K* nu nubar) | <2.7e-5 at 90% CL | Belle 2017 | 2a51560a4570 |
| B023 | SM BR(B0 -> K*0 nu nubar) | (9.2 +/- 1.0)e-6 | Buras et al. 2015 | 20932455faf2 |

## Issues (if any)
- B021: CHK-1 fail. WA-v2 promoted the CA-flagged CDF values and added years
  to the existing observable blocks, but the TeX still says LHCb used
  `3.0 fb^-1`; that number is verified in the LHCb snapshot and remains only
  under `supporting_measurements`, not `pdg_or_equivalent`.
- B023: CHK-2 fail. The TeX post-2008 paragraph cites
  `BaBar2013:BKstarNuNu`, but `references/B023/source_manifest.yaml` contains
  `BaBar2013:BKstarNunu`; all manifest entries otherwise have non-empty
  `snapshot_path` fields pointing to tracked text snapshots.
- B021/B022/B023 CHK-3: `ls` and `find ... -name '*.pdf'` showed no tracked
  PDFs; the reference directories contain text snapshots, manifests, and
  B021's process-local `sha256sums.txt`.
- B021/B022/B023 CHK-4: `status_history` contains `WRITER-INITIATED` followed
  by `WRITER-DONE` with ISO 8601 timestamps; before this CA update the most
  recent entry for each process was the WA-v2 `WRITER-DONE` transition.
- B021/B022/B023 CHK-5: the prompt's literal `rg -l -E "<process keyword>"`
  form is invalid for this installed ripgrep because `-E` is an encoding flag.
  I recorded that failure and ran the equivalent `rg -l -e` searches. B021
  returned only `Lambda_IR` / unrelated `muToEGamma` false positives after
  excluding notebooks; B022 and B023 had no code hits after excluding notebook
  base64 false positives. Cited file:line evidence exists.
- B021/B022/B023 CHK-6: `HIGH` is consistent with the rubric because these
  entries require new Delta B = 1 semileptonic or invisible-mode Hamiltonians,
  form-factor treatment, and likelihood/CL handling beyond the existing
  Delta F = 2 SLL/SLR/VLL/VRR/LR1/LR2 machinery.
- B021/B022/B023 CHK-7: non-catalog TeX/Markdown grep found only plan rows or
  unrelated gauge-parameter text, not rc1.1 load-bearing values for these
  catalog entries.
- B021/B022/B023 CHK-8: no `\textbf{CHECK}`, `TODO`, unresolved `\ref{...}` or
  `\cite{...}` macros were found; brace counts are balanced in all three TeX
  files.
