# CA Worklog: ca_w5a_beauty
**Date**: 2026-05-16
**Family**: beauty
**Process IDs**: B021, B022, B023, B034

Interpretation note: CHK-1 was enforced exactly as written. Numerical claims in
the TeX must trace to the YAML `pdg_or_equivalent` block with source metadata;
values stored only in `supporting_measurements`, `theory_inputs`, or
`theory_context` are treated as rework even when the local snapshot verifies
the number.

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| B021 | FAIL | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| B022 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| B023 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| B034 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| B021 | BR(Lambda_b0 -> Lambda mu+ mu-) | (1.08 +/- 0.28)e-6 | PDG Live 2025 | 73959783ca3c |
| B021 | dBR/dq2, 15 < q2 < 20 GeV^2/c^4 | (1.18 +0.09/-0.08 +/-0.03 +/-0.27)e-7 (GeV^2/c^4)^-1 | LHCb 2015 | 2feb8efe212b |
| B021 | A_FB^l, A_FB^h | -0.05 +/-0.09 +/-0.03; -0.29 +/-0.07 +/-0.03 | LHCb 2015 | 2feb8efe212b |
| B021 | CDF observation claims | 24 events; 5.8 sigma; 6.8 fb^-1 at 1.96 TeV; BR=(1.73 +/-0.42 +/-0.55)e-6 | CDF 2011 | 83084eabd1d4 |
| B022 | BR(B+ -> K+ nu nubar), HFLAV average | (1.38 +/- 0.35)e-5 | HFLAV Dec. 2025 | 4c562ccd073f |
| B022 | BR(B+ -> K+ nu nubar), PDG display | (2.30 +0.71/-0.64)e-5 | PDG Live/API 2026 | af37fb365951 |
| B022 | Belle II evidence measurement | (2.3 +/-0.5 stat +0.5/-0.4 syst)e-5; 362 fb^-1; 3.5 sigma; 2.7 sigma above SM | Belle II 2024 | 5e8e54e24cbc |
| B022 | BaBar hadronic-tag limit | <3.7e-5 at 90% CL | BaBar 2013 | 8366c09006d7 |
| B022 | SM BR(B+ -> K+ nu nubar) | (5.58 +/- 0.37)e-6 | HPQCD 2023 | 66d8ba9efc52 |
| B023 | BR(B0 -> K*0 nu nubar) | <1.8e-5 at 90% CL | PDG/HFLAV 2025 | f8d7014f1e7f |
| B023 | BR(B+ -> K*+ nu nubar) | <4.0e-5 at 90% CL | PDG/HFLAV 2025 | d9aa9db63226 |
| B023 | Combined BR(B -> K* nu nubar) | <2.7e-5 at 90% CL | Belle 2017 | 2a51560a4570 |
| B023 | SM BR(B0 -> K*0 nu nubar) | (9.2 +/- 1.0)e-6 | Buras et al. 2015 | 20932455faf2 |
| B034 | beta_s(b -> s sbar s) | (3.7 +/- 3.5)e-2 rad | PDG Live 2025 | 4c88ea2de1ba |
| B034 | phi_s^{s sbar s}, LHCb combined | -0.074 +/- 0.069 rad | LHCb 2023 / PDG note | b003cf5fbbda |
| B034 | \|lambda\|, LHCb combined | 1.009 +/- 0.030 | LHCb 2023 / PDG note | b003cf5fbbda |
| B034 | BR(B_s0 -> phi phi) | (1.85 +/- 0.14)e-5 | HFLAV May 2024 | 0109068c253d |
| B034 | LHCb 2014 phi_s milestone | 3.0 fb^-1; about 4000 decays; phi_s=-0.17 +/-0.15 +/-0.03 rad | LHCb 2014 | 21b1b1b7fa25 |
| B034 | LHCb 2019 phi_s, \|lambda\| milestone | 5.0 fb^-1; about 9000 decays; phi_s=-0.073 +/-0.115 +/-0.027 rad; \|lambda\|=0.99 +/-0.05 +/-0.01 | LHCb 2019 | 37b596af8f67 |
| B034 | LHCb 2023 precision fit | 6 fb^-1 at 13 TeV; phi_s=-0.042 +/-0.075 +/-0.009 rad; \|lambda\|=1.004 +/-0.030 +/-0.009 | LHCb 2023 | b003cf5fbbda |

## Issues (if any)
- B021: CHK-1 fail. TeX post-2008 CDF numerical claims are verified in
  `supporting_measurements`, but not under `pdg_or_equivalent`. The
  `pdg_or_equivalent` observables also lack an explicit year field in the
  block, even though the manifest supplies years.
- B021: CHK-2 fail. All TeX source keys resolve, all snapshots exist and are
  tracked, and hashes match, but `source_manifest.yaml` uses
  `local_snapshot_path` rather than the required `snapshot_path` key for every
  entry.
- B022: CHK-1 fail. The HPQCD SM prediction quoted in TeX is verified in
  `theory_inputs`, but not under `pdg_or_equivalent`.
- B023: CHK-1 fail. The Belle 2017 combined vector limit and Buras 2015 SM
  prediction quoted in TeX are verified in `experiment_inputs` /
  `theory_context`, but not under `pdg_or_equivalent`.
- Code coverage checks: the prompt's literal `rg -l -E "<process keyword>"`
  form is invalid for this installed ripgrep because `-E` is an encoding flag.
  I recorded that failure and ran the intended equivalent with `rg -l -e`.
  B021 returned only `Lambda_IR` / unrelated `muToEGamma` false positives
  after excluding notebooks; B022, B023, and B034 had no code hits after
  excluding notebook base64 false positives. Cited file:line evidence exists.
- CHK-3: `find ... -name '*.pdf'` found no PDFs in B021, B022, B023, or B034
  reference directories. `ls` showed only `.txt` snapshots, `source_manifest.yaml`,
  and process-local `sha256sums.txt` where present.
- CHK-7: non-catalog TeX grep found no rc1.1 paper references to these
  processes or load-bearing numbers, so no contradiction was found.
- CHK-8: no `CHECK`, `TODO`, unresolved `\ref{...}` or `\cite{...}` markers
  were found; a brace-count check across the four TeX files was balanced.
