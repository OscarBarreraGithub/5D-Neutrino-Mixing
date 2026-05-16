# Charm Family Fact-Check Report

Fact-checker: factcheck-codex-charm
Date: 2026-05-16
Family: charm
Processes: C001 C002 C003 C004 C005 C006 C007 C008

Summary: all 8 charm processes are VERIFIED. I found 0 mismatches and 0 unresolvable source URLs. CFW RS scale numerals and similar normalization/context values were treated as policy-level theory context, not measured observables, per the L001 signoff policy.

## C001

Verdict: VERIFIED

| Claim | source_url | verified? | Note |
|---|---|---:|---|
| x_D = (0.405 +/- 0.043)% with 95% C.L. [0.320, 0.489] | https://hflav-eos.web.cern.ch/hflav-eos/charm/CKM25/results_mix_cpv.html | Y | HFLAV CKM25 table_results.pdf lists the all-CPV fit value and interval. |
| y_D = (0.636 +/- 0.024)% with 95% C.L. [0.590, 0.682] | https://hflav-eos.web.cern.ch/hflav-eos/charm/CKM25/results_mix_cpv.html | Y | HFLAV CKM25 table_results.pdf lists the all-CPV fit value and interval. |
| Delta m_D = (0.997 +/- 0.116) x 10^10 hbar s^-1 = (6.562 +/- 0.764) x 10^-15 GeV | https://pdgprod.lbl.gov/pdgprod/pdgLiveJson/DataBlock.action?node=S032D | Y | PDG Live S032D shows 0.997 +/- 0.116 as the HFLAV-produced evaluation; the GeV conversion uses hbar = 6.582119569e-25 GeV s. |
| Delta y = (0.031 +/- 0.035 +/- 0.013)% | https://hflav-eos.web.cern.ch/hflav-eos/charm/CKM25/results_mix_cpv.html | Y | HFLAV CKM25 measured-input table contains the LHCb Delta y row with this value. |
| Delta Gamma_D / Gamma_D = (1.272 +/- 0.048)% | https://hflav-eos.web.cern.ch/hflav-eos/charm/CKM25/results_mix_cpv.html | Y | Derived directly as 2*y_D from the verified HFLAV y_D value. |
| CFW 21 TeV / 33 TeV contextual RS scales | https://arxiv.org/abs/0804.1954 | N/A | Skipped as theory-normalization/context under policy; arXiv metadata matches Csaki, Falkowski, Weiler, "The Flavor of the Composite Pseudo-Goldstone Higgs." |

Key reference metadata checked: arXiv:2106.03744 is the LHCb paper "Observation of the mass difference between neutral charm-meson eigenstates"; arXiv:0804.1954 metadata matches CFW 2008.

## C002

Verdict: VERIFIED

| Claim | source_url | verified? | Note |
|---|---|---:|---|
| |q/p| = 0.983 +0.015/-0.014 with 95% C.L. [0.955, 1.012] | https://hflav-eos.web.cern.ch/hflav-eos/charm/CKM25/results_mix_cpv.html | Y | HFLAV CKM25 table_results.pdf lists this all-CPV fit value. |
| phi_D = -1.51 +1.03/-1.06 deg with 95% C.L. [-3.63, 0.51] | https://hflav-eos.web.cern.ch/hflav-eos/charm/CKM25/results_mix_cpv.html | Y | HFLAV CKM25 table_results.pdf lists this all-CPV fit value. |
| No-indirect-CPV point has Delta chi2 = 2.16 and 0.95 sigma | https://hflav-eos.web.cern.ch/hflav-eos/charm/CKM25/results_mix_cpv.html | Y | HFLAV page states the (|q/p|, phi) = (1,0) point has Delta chi2 = 2.16 and is consistent at 0.95 sigma. |
| PDG cross-check: |q/p| = 0.983 +/- 0.015, phi = -1.51 +/- 1.04 deg | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-d-dbar-mixing.pdf | Y | PDG 2025 Table 70.7 lists the rounded all-CPV values and intervals. |
| CFW 2008 context predates modern CPV inputs | https://arxiv.org/abs/0804.1954 | N/A | Context only; arXiv metadata matches the manifest. |

Key reference metadata checked: arXiv:2208.06512 matches the LHCb charm-mixing measurement; arXiv:2410.22961 matches the Belle/Belle II model-independent D0 mixing measurement; arXiv:0804.1954 matches CFW 2008.

## C003

Verdict: VERIFIED

| Claim | source_url | verified? | Note |
|---|---|---:|---|
| Delta a_CP^dir = (-0.159 +/- 0.029)% | https://hflav-eos.web.cern.ch/hflav-eos/charm/CKM23/DCPV/direct_indirect_cpv.html | Y | HFLAV direct/indirect CPV page updated 2025 lists this fit value. |
| a_CP^ind = (-0.010 +/- 0.012)% | https://hflav-eos.web.cern.ch/hflav-eos/charm/CKM23/DCPV/direct_indirect_cpv.html | Y | HFLAV direct/indirect CPV page lists this fit value. |
| No-CPV test Delta chi2 = 32.3, CL table 6.9e-8, CL text 9.7e-8, significance 5.3 sigma | https://hflav-eos.web.cern.ch/hflav-eos/charm/CKM23/DCPV/direct_indirect_cpv.html | Y | HFLAV page gives the table CL, text CL, Delta chi2, and significance. |
| LHCb 2019 Delta A_CP = (-15.4 +/- 2.9) x 10^-4 = (-0.154 +/- 0.029)% | https://arxiv.org/abs/1903.08726 | Y | arXiv metadata matches the discovery paper; the local snapshot records the exact value. |
| HFLAV CKM25 input A_CP^K - A_CP^pi = (-0.154 +/- 0.029)% | https://hflav-eos.web.cern.ch/hflav-eos/charm/CKM25/results_mix_cpv.html | Y | HFLAV CKM25 measured-input table lists this LHCb input. |
| A_pi = (0.225 +/- 0.057)% with 95% C.L. [0.11, 0.34] | https://hflav-eos.web.cern.ch/hflav-eos/charm/CKM25/figures/table_results.pdf | Y | HFLAV CKM25 table_results.pdf lists this all-CPV fit parameter. |
| A_K = (0.068 +/- 0.051)% with 95% C.L. [-0.03, 0.17] | https://hflav-eos.web.cern.ch/hflav-eos/charm/CKM25/figures/table_results.pdf | Y | HFLAV CKM25 table_results.pdf lists this all-CPV fit parameter. |
| CFW context with 21/33 TeV scale wording | https://arxiv.org/abs/0804.1954 | N/A | Skipped as theory-normalization/context under policy; metadata matches the manifest. |

Key reference metadata checked: arXiv:1903.08726 matches LHCb "Observation of CP violation in charm decays"; arXiv:0804.1954 matches CFW 2008.

## C004

Verdict: VERIFIED

| Claim | source_url | verified? | Note |
|---|---|---:|---|
| PDG BR(D0 -> mu+ mu-) < 2.1e-9 at 90% C.L. | https://pdgprod.lbl.gov/pdgprod/pdgLiveJson/api/listings/S032.28 | Y | PDG 2026 API listing S032.28 gives CHEKHOVSKY 2025F value <2.1E-9. |
| CMS BR(D0 -> mu+ mu-) < 2.4e-9 at 95% C.L.; 64.5 fb^-1 at 13.6 TeV in 2022-2023 | https://cms-results.web.cern.ch/cms-results/public-results/publications/BPH-23-008/ | Y | CMS public page and arXiv:2506.06152 abstract contain the limit, luminosity, energy, and years. |
| LHCb predecessor BR(D0 -> mu+ mu-) < 3.1e-9 at 90% C.L.; companion 95% limit <3.5e-9; 9 fb^-1 at 7, 8, 13 TeV | https://arxiv.org/abs/2212.11203 | Y | arXiv abstract contains the 3.1e-9, 9 fb^-1, and energies; arXiv source defines 3.1(3.5)e-9 at 90(95)% C.L.; PDG footnote also confirms 3.5e-9. |
| Burdman et al. relation Br(D0 -> mu+ mu-)_(gamma gamma) ~= 2.7e-5 Br(D0 -> gamma gamma), minimum >=3e-13, VMD Br(D0 -> gamma gamma) = (3.5 +4.0 -2.6)e-8 | https://arxiv.org/abs/hep-ph/0112235 | Y | arXiv metadata matches the theory paper; local TeX-source snapshot records the quoted relation and values. |
| Gisbert-Hiller-Suelmann rare-charm EFT context and C7,C9,C10 coefficient discussion | https://link.springer.com/article/10.1007/JHEP12(2024)102 | N/A | Context source has no measured value claim; title/authors/abstract match the manifest. |
| CFW 21/33 TeV RS scales | https://arxiv.org/abs/0804.1954 | N/A | Skipped as theory-normalization/context under policy; metadata matches the manifest. |

Key reference metadata checked: CMS arXiv:2506.06152 matches the CMS rare D0 dimuon paper; arXiv:2212.11203 matches LHCb "Search for rare decays of D0 mesons into two muons"; arXiv:hep-ph/0112235 and the Springer JHEP page match the theory/context citations.

## C005

Verdict: VERIFIED

| Claim | source_url | verified? | Note |
|---|---|---:|---|
| PDG BR(D0 -> e+ e-) < 7.9e-8 at 90% C.L. | https://pdgprod.lbl.gov/pdgprod/pdgLiveJson/api/listings/S032.39 | Y | PDG 2026 API listing S032.39 gives the Belle 2010 value <7.9E-8 and marks it used. |
| Belle BR(D0 -> e+ e-) < 7.9e-8 at 90% C.L.; 660 fb^-1 | https://arxiv.org/abs/1003.2345 | Y | arXiv abstract gives the e+e- limit and 660 fb^-1 dataset. |
| BABAR BR(D0 -> e+ e-) < 1.7e-7 at 90% C.L.; 468 fb^-1 | https://arxiv.org/abs/1206.5419 | Y | arXiv abstract gives the e+e- limit and 468 fb^-1 dataset. |
| CFW 21/33 TeV RS scales | https://arxiv.org/abs/0804.1954 | N/A | Skipped as theory-normalization/context under policy; metadata matches the manifest. |
| Burdman rare-charm model-dependence context | https://arxiv.org/abs/hep-ph/0112235 | N/A | Context only for C005; metadata matches the manifest. |

Key reference metadata checked: arXiv:1003.2345 matches Belle "Search for leptonic decays of D0 mesons"; arXiv:1206.5419 matches BABAR "Search for the decay modes D0 -> e+e-, D0 -> mu+mu-, and D0 -> e mu"; arXiv:hep-ph/0112235 and arXiv:0804.1954 match the theory references.

## C006

Verdict: VERIFIED

| Claim | source_url | verified? | Note |
|---|---|---:|---|
| PDG BR(D0 -> mu+- e-+) < 1.3e-8 at 90% C.L. | https://pdgprod.lbl.gov/pdgprod/pdgLiveJson/api/listings/S032.40 | Y | PDG 2026 API listing S032.40 gives the LHCb 2016 value <1.3E-8 and marks it used. |
| LHCb BR(D0 -> e+- mu-+) < 1.3e-8 at 90% C.L.; 3.0 fb^-1 at 7 and 8 TeV; normalized to D0 -> K- pi+ | https://arxiv.org/abs/1512.00322 | Y | arXiv abstract contains the limit, luminosity, collision energies, and normalization channel. |
| Belle combined D0 -> e+mu- plus mu+e- limit <2.6e-7 at 90% C.L.; 660 fb^-1 | https://arxiv.org/abs/1003.2345 | Y | arXiv abstract gives the combined LFV limit and 660 fb^-1 dataset. |
| BABAR D0 -> e mu limit <3.3e-7 at 90% C.L.; 468 fb^-1 | https://arxiv.org/abs/1206.5419 | Y | arXiv abstract gives the LFV limit and 468 fb^-1 dataset. |
| CFW 21/33 TeV RS scales | https://arxiv.org/abs/0804.1954 | N/A | Skipped as theory-normalization/context under policy; metadata matches the manifest. |

Key reference metadata checked: arXiv:1512.00322 matches LHCb "Search for the lepton-flavour violating decay D0 -> e+- mu-+"; arXiv:1003.2345, arXiv:1206.5419, and arXiv:0804.1954 match their manifest metadata.

## C007

Verdict: VERIFIED

| Claim | source_url | verified? | Note |
|---|---|---:|---|
| PDG BR(D+ -> pi+ mu+ mu-) < 6.7e-8 at 90% C.L. | https://pdgprod.lbl.gov/pdgprod/pdgLiveJson/api/listings/S031.42 | Y | PDG 2026 API listing S031.42 gives the LHCb 2021 value <6.7E-8 and marks it used. |
| LHCb 2021 BR(D+ -> pi+ mu+ mu-) <6.7e-8 at 90% C.L.; companion 95% limit <7.4e-8; 1.6 fb^-1 in 2016 | https://arxiv.org/abs/2011.00217 | Y | arXiv abstract gives the 25-mode search scope and 1.6 fb^-1; arXiv source and local snapshot give the 67/74 x 10^-9 table values. |
| LHCb 2021 search covered 25 rare/forbidden D+ and Ds+ modes | https://arxiv.org/abs/2011.00217 | Y | arXiv abstract states no evidence for the 25 investigated modes. |
| LHCb 2021 quotes short-distance SM branching fractions of order 10^-12 for FCNC D modes | https://arxiv.org/abs/2011.00217 | Y | arXiv source introduction contains the order(10^-12) statement; local snapshot records it. |
| LHCb 2013 nonresonant BR(D+ -> pi+ mu+ mu-) <7.3(8.3)e-8 at 90(95)% C.L.; 1.0 fb^-1 at 7 TeV | https://arxiv.org/abs/1304.6365 | Y | arXiv abstract contains the 7.3(8.3)e-8 limit and 1.0 fb^-1 at sqrt(s)=7 TeV. |
| de Boer-Hiller rare-charm EFT context | https://arxiv.org/abs/1510.00311 | N/A | Context source has no measured C007 value claim; metadata and abstract match. |
| CFW 21/33 TeV RS scales | https://arxiv.org/abs/0804.1954 | N/A | Skipped as theory-normalization/context under policy; metadata matches the manifest. |

Key reference metadata checked: arXiv:2011.00217, arXiv:1304.6365, arXiv:1510.00311, and arXiv:0804.1954 match the manifest authors/titles.

## C008

Verdict: VERIFIED

| Claim | source_url | verified? | Note |
|---|---|---:|---|
| PDG BR(D+ -> pi+ e+ mu-) < 2.1e-7 at 90% C.L.; companion 95% limit <2.3e-7 | https://pdgprod.lbl.gov/pdgprod/pdgLiveJson/api/listings/S031.110 | Y | PDG 2026 API listing S031.110 gives <2.1E-7; LHCb 2021 arXiv source/local snapshot gives the 210/230 x 10^-9 table values. |
| PDG BR(D+ -> pi+ e- mu+) < 2.2e-7 at 90% C.L.; companion 95% limit <2.2e-7 | https://pdgprod.lbl.gov/pdgprod/pdgLiveJson/api/listings/S031.111 | Y | PDG 2026 API listing S031.111 gives <2.2E-7; LHCb 2021 arXiv source/local snapshot gives the 220/220 x 10^-9 table values. |
| LHCb 2021 search scope: 25 modes, 1.6 fb^-1 in 2016, no significant deviation | https://arxiv.org/abs/2011.00217 | Y | arXiv abstract contains the 25-mode search, 1.6 fb^-1 dataset, and no-evidence statement. |
| BABAR predecessor: 384 fb^-1, D+ -> pi+ e+ mu- <2.9e-6 and D+ -> pi+ e- mu+ <3.6e-6 at 90% C.L. | https://arxiv.org/abs/1107.4465 | Y | arXiv metadata matches the BABAR search; local snapshot records the exact two C008 charge-mode limits. |
| CFW historical RS baseline | https://arxiv.org/abs/0804.1954 | N/A | Context only; no measured C008 observable is taken from this source. |

Key reference metadata checked: arXiv:2011.00217 matches LHCb "Searches for 25 rare and forbidden decays of D+ and Ds+ mesons"; arXiv:1107.4465 matches BABAR "Searches for Rare or Forbidden Semileptonic Charm Decays"; arXiv:0804.1954 matches CFW 2008.

## Unresolvable Sources

None.

## Mismatches

None.
