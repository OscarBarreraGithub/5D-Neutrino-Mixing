# CA Worklog: ca_w4_ew
**Date**: 2026-05-16
**Family**: top_higgs_ew
**Process IDs**: EW001 EW002 EW003

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| EW001 | PASS | PASS | PASS | FAIL | PASS | PASS | PASS | PASS | WRITER-REWORK |
| EW002 | PASS | PASS | PASS | FAIL | PASS | PASS | PASS | PASS | WRITER-REWORK |
| EW003 | FAIL | PASS | PASS | FAIL | PASS | PASS | PASS | PASS | WRITER-REWORK |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| EW001 | PDG U=0 fit S | 0.026 +/- 0.075 | PDG 2025 electroweak review | c902522a3c60 |
| EW001 | PDG U=0 fit T | 0.047 +/- 0.066 | PDG 2025 electroweak review | c902522a3c60 |
| EW001 | PDG U=0 rho(S,T) | 0.90 | PDG 2025 electroweak review | c902522a3c60 |
| EW001 | PDG floating-U fit S,T,U | 0.021 +/- 0.096; 0.04 +/- 0.12; 0.008 +/- 0.092 | PDG 2025 electroweak review | c902522a3c60 |
| EW001 | Warped S estimate and bound | S approx 30 v^2/M_KK^2; S <= 0.18 at 95%; M_KK > about 3.2 TeV | PDG 2025 electroweak review | c902522a3c60 |
| EW001 | Gfitter S,T,U | 0.05 +/- 0.11; 0.09 +/- 0.13; 0.01 +/- 0.11 | Gfitter oblique page | c9fde1fad4e1 |
| EW001 | CDF W mass | 80.4335 +/- 0.0094 GeV | de Blas et al. 2022 | f5155be1508e |
| EW001 | HEPfit/de Blas ST and U context | S=0.100 +/- 0.073; T=0.202 +/- 0.056; U=0.134 +/- 0.087 | de Blas et al. 2022 | f5155be1508e |
| EW002 | PDG superallowed V_ud | 0.97367 +/- 0.00032 | PDG Vud/Vus/CKM unitarity review | 9ad375863344 |
| EW002 | PDG kaon-sector V_us average | 0.22431 +/- 0.00085; scale factor 2.5 | PDG Vud/Vus/CKM unitarity review | 9ad375863344 |
| EW002 | PDG first-row sum | 0.9983(6)(4), about 2 sigma; about 3 sigma without nuclear-structure uncertainties | PDG Vud/Vus/CKM unitarity review | 9ad375863344 |
| EW002 | FLAG K_l3 product and f_+(0) | 0.21654(41); 0.9698(17) for N_f=2+1+1 | FLAG Review 2024 | b9436b897ef8 |
| EW002 | FLAG V_us from f_+(0) | 0.22328(58) | FLAG Review 2024 | b9436b897ef8 |
| EW002 | FLAG first-row sums | 0.99802(66), 3.0 sigma; 0.99888(67), 1.7 sigma | FLAG Review 2024 | b9436b897ef8 |
| EW002 | Hardy-Towner V_ud | 0.97373 +/- 0.00031 | Hardy and Towner 2020 | 913ac2c95b82 |
| EW002 | CFW RS context | about 21 TeV and about 33 TeV | Csaki-Falkowski-Weiler 2008 | d22463e7e9e0 |
| EW003 | PDG V_cb inclusive/exclusive/average | 42.2 +/- 0.5; 39.8 +/- 0.6; 41.1 +/- 1.2, all x 10^-3 | PDG 2024 Vcb/Vub review | 095c00ac3711 |
| EW003 | PDG V_ub inclusive/exclusive/average | 4.13 with 0.12 exp, +0.13/-0.14 theory, 0.18 model; 3.70 +/- 0.10 +/- 0.12; 3.82 +/- 0.20, all x 10^-3 | PDG 2024 Vcb/Vub review | 095c00ac3711 |
| EW003 | PDG V_cb marginal consistency | about 3.0 sigma | PDG 2024 Vcb/Vub review | 095c00ac3711 |
| EW003 | HFLAV inclusive V_cb | 41.97 +/- 0.48 x 10^-3 | HFLAV 2023 semileptonic report | e063d36ade2f |
| EW003 | HFLAV exclusive V_ub,V_cb | 3.43 +/- 0.12; 39.77 +/- 0.46, all x 10^-3 | HFLAV 2023 semileptonic report | e063d36ade2f |
| EW003 | FLAG exclusive V_cb | 39.23(65) x 10^-3 for B -> D* l nu, N_f=2+1 | FLAG Review 2024 | 6f400dee4a10 |
| EW003 | FLAG exclusive V_ub | 3.61 +/- 0.16 x 10^-3 for B -> pi l nu | FLAG Review 2024 | 6f400dee4a10 |

## Check evidence
- CHK-1: Grepped the TeX files for numerical claims. EW001 and EW002 numerical physics claims trace to `pdg_or_equivalent` entries with year, value, uncertainty or CL/significance, source URL, access date, snapshot path, and sha256. EW003 fails because its `pdg_or_equivalent` entries lack `year`, and the PDG `3.0 sigma` consistency statement is outside `pdg_or_equivalent`.
- CHK-2: Every key listed in the TeX "Key references" sections appears in `flavor_catalog/references/<process_id>/source_manifest.yaml`. Each manifest entry has a non-empty `snapshot_path`, and `git ls-files flavor_catalog/references/EW00{1,2,3}` shows those snapshots are tracked.
- CHK-3: `ls flavor_catalog/references/EW001`, `EW002`, and `EW003` shows only `.txt` snapshots plus `source_manifest.yaml`; no `.pdf` files are tracked.
- CHK-4: All three sidecars contain `WRITER-INITIATED`, but none contains `WRITER-DONE`. `rg` over `flavor_catalog/worklogs/writer/*top_higgs_ew*.md` found WA logs only for T001, T010, T002, T007, and T018, not EW001-EW003.
- CHK-5: The requested literal `rg -l -E "<process keyword>" ... 2>&1` fails under this ripgrep because `-E` is the encoding flag. I reran the equivalent `rg -l -e` checks. EW001 has no oblique/STU implementation hit, aside from a generic `m_W` comment under the detailed grep. EW002/EW003 have CKM target/tolerance hits at `quarkConstraints/scan.py:42`, `quarkConstraints/fit.py:47`, `quarkConstraints/modern/inputs.py:631`, and `tests/test_quark_fit.py:575`, but no first-row unitarity or semileptonic-B likelihood implementation.
- CHK-6: HIGH is consistent for all three. EW001 needs a new EW precision likelihood and RS electroweak matching; EW002 needs new charged-current/G_F/CKM-unitarity likelihood machinery; EW003 needs new semileptonic likelihoods, lattice form-factor inputs, inclusive scheme handling, and optional charged-current WET operators. None reuses the existing Delta F = 2 SLL/SLR/VLL/VRR/LR1/LR2 basis as a live constraint.
- CHK-7: No rc1.1 load-bearing number is contradicted. EW002 preserves the CFW contextual 21 TeV / 33 TeV values as literature context; EW001 and EW003 are companion-catalog additions and do not revise the sealed rc1.1 scan headline.
- CHK-8: Fixed-string greps for `\textbf{CHECK}`, `TODO`, `\ref{`, `\cite{`, `MISSING`, and `??` returned no hits. Manual scan found no blocking math-mode or brace issues.

## Issues (if any)
- EW001: CHK-4 fails. Sidecar needs a WA `WRITER-DONE` transition after writer polish, and the WA batch worklog for EW001-EW003 needs to exist or be identified.
- EW002: CHK-4 fails. Sidecar needs a WA `WRITER-DONE` transition after writer polish, and the WA batch worklog for EW001-EW003 needs to exist or be identified.
- EW003: CHK-1 fails. Add `year` to every `pdg_or_equivalent` value block and promote the `3.0 sigma` marginal-consistency claim from `additional_numerical_context` into `pdg_or_equivalent` with full metadata.
- EW003: CHK-4 fails. Sidecar needs a WA `WRITER-DONE` transition after writer polish, and the WA batch worklog for EW001-EW003 needs to exist or be identified.
