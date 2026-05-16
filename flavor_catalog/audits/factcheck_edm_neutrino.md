# Flavor Catalog Fact Check: edm_neutrino

Agent: `factcheck-codex-edm_neutrino`  
Date: 2026-05-16  
Scope: E001, E002, E004, E006, E007, E008, E009

Method: read each process TeX sidecar and source manifest, opened the cited PDG/arXiv/collaboration URLs, and checked exact numerical claims against either the fetched abstract/datablock or the local source snapshot under `flavor_catalog/references/<process_id>/`. Theory normalization constants and broad policy-scale estimates were not treated as measured observables, but claimed-to-be-quoted theory benchmark formulae were checked.

Summary: 6 VERIFIED, 1 PARTIAL, 0 MISMATCH. The only PARTIAL item is E009's INSPIRE key-reference URL for Weinberg 1989, which WebFetch renders as a JavaScript-only page; the key-reference metadata was cross-checked against APS and the local snapshot, and no value-bearing catalog claim depends on unresolved INSPIRE text.

## E001

Verdict: VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| PDG electron EDM limit `|d_e| < 4.1e-30 e cm` at 90% CL; table value `<0.041` in `10^-28 e cm` | https://pdgprod.lbl.gov/pdgprod/pdgLive/DataBlock.action?node=S003EDM | Y | Current PDG page/search and a direct curl of the PDG datablock show S003EDM, `<0.041`, 90% CL, Roussy 2023, and the local snapshot records the equivalent `4.1e-30 e cm`. |
| Roussy 2023 measurement `d_e = (-1.3 +/- 2.0_stat +/- 0.6_syst)e-30 e cm`; limit `4.1e-30`; improvement factor about 2.4 | https://arxiv.org/abs/2212.11841 | Y | arXiv metadata matches Tanya S. Roussy et al., "A new bound on the electron's electric dipole moment"; abstract gives the factor `~2.4`, and the local snapshot records the signed measurement and limit. |
| ACME/Andreev 2018 superseded limit `|d_e| < 1.1e-29 e cm`; measurement note `(4.3 +/- 3.1 +/- 2.6)e-30 e cm` | https://authors.library.caltech.edu/records/6swh6-6ya42 | Y | CaltechAUTHORS metadata matches "Improved Limit on the Electric Dipole Moment of the Electron", V. Andreev / ACME Collaboration; local snapshot records the quoted ACME/PDG numerical values. |
| Panico-Pomarol-Riembau 2019 key reference metadata and ACME-bound context | https://arxiv.org/abs/1810.09413 | Y | arXiv metadata matches title/authors; abstract quotes the ACME `1.1e-29 e cm` bound used only as context. |
| CFW 2008 key reference metadata | https://arxiv.org/abs/0804.1954 | Y | arXiv metadata matches Csaki, Falkowski, Weiler, "The Flavor of the Composite Pseudo-Goldstone Higgs". |

Mismatches: none.

## E002

Verdict: VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| PDG muon EDM direct limit `|d_mu| < 1.8e-19 e cm` at 95% CL | https://pdgprod.lbl.gov/pdgprod/pdgLive/DataBlock.action?cl=T&desig=&init=0&node=S004EDM&parCode=S004EDM | Y | PDG Live S004EDM lists `<1.8` in units of `10^-19 e cm`, CL 95%, Bennett 2009. |
| PDG combined measurement `(0.0 +/- 0.9)e-19 e cm` | https://pdgprod.lbl.gov/pdgprod/pdgLive/DataBlock.action?cl=T&desig=&init=0&node=S004EDM&parCode=S004EDM | Y | PDG note states the combined Bennett result is also presented as `(0.0 +/- 0.9) x 10^-19 e cm`. |
| Bennett 2009 primary measurement `d = -0.1(0.9)e-19 e cm`; limit `<1.9e-19 e cm` at 95% CL | https://arxiv.org/abs/0811.1207 | Y | arXiv metadata matches G.W. Bennett et al., "An Improved Limit on the Muon Electric Dipole Moment"; abstract quotes both values. |
| PSI proposal: `p = 125 MeV/c`, `|B| = 3 T`, rest-frame `E = 1 GV/m`, sensitivity `<= 6e-23 e cm`, and current direct limit `1.8e-19 e cm` | https://arxiv.org/abs/2102.08838 | Y | arXiv metadata matches Adelmann et al.; abstract quotes all projection values. |
| Renga 2024 status: first phase by 2026; ultimate factor 100 improvement by early 2030s | https://arxiv.org/abs/2409.20050 | Y | arXiv metadata matches Francesco Renga for the muEDM collaboration; abstract quotes the 2026 and early-2030s projection claims. |
| CFW 2008 key reference metadata | https://arxiv.org/abs/0804.1954 | Y | arXiv title/authors match manifest. |

Mismatches: none.

## E004

Verdict: VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| PDG neutron EDM current limit `|d_n| < 1.8e-26 e cm` at 90% CL; table `<0.18` in `10^-25 e cm` | https://pdgprod.lbl.gov/pdgprod/pdgLive/DataBlock.action?node=S017EDM | Y | PDG Live S017EDM lists `<0.18`, CL 90%, Abel 2020. |
| Abel 2020 measurement `d_n = (0.0 +/- 1.1_stat +/- 0.2_sys)e-26 e cm`; PDG limit `<1.8e-26` | https://arxiv.org/abs/2001.11966 | Y | arXiv metadata matches C. Abel et al., "Measurement of the permanent electric dipole moment of the neutron"; abstract quotes the signed measurement, and PDG/local snapshot gives the corresponding limit. |
| Superseded Pendlebury 2015 PDG row `<0.30` in `10^-25 e cm`, equivalent `<3.0e-26 e cm` at 90% CL | https://pdgprod.lbl.gov/pdgprod/pdgLive/DataBlock.action?node=S017EDM | Y | PDG Live S017EDM lists `<0.30`, 90% CL, Pendlebury 2015. |
| Koenig-Neubert-Straub 2014 key reference metadata and neutron-EDM/composite-dipole context | https://arxiv.org/abs/1403.2756 | Y | arXiv metadata matches title/authors; abstract confirms composite/warped dipole context. |
| Bhattacharya et al. 2022 key reference metadata and lattice target | https://arxiv.org/abs/2203.03746 | Y | arXiv metadata matches title/authors and abstract names dimension-4/dimension-6 gluonic CP violation plus isovector qCEDM. |
| CFW 2008 key reference metadata | https://arxiv.org/abs/0804.1954 | Y | arXiv title/authors match manifest. |

Mismatches: none.

## E006

Verdict: VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| Graner 2016 direct Hg-199 EDM result `(2.20 +/- 2.75_stat +/- 1.48_syst)e-30 e cm`; limit `<7.4e-30 e cm` at 95% CL | https://arxiv.org/abs/1601.04339 | Y | arXiv metadata matches B. Graner et al., "Reduced Limit on the Permanent Electric Dipole Moment of 199Hg"; abstract quotes the measurement and limit. |
| Graner 2016 improvement factor `4` | https://arxiv.org/abs/1601.04339 | Y | arXiv abstract says the new limit improved the previous limit by a factor of 4. |
| PDG neutron EDM cross-reference rows from Hg+theory: Graner `<0.16`, Sahoo `<0.22`, units `10^-25 e cm`, 95% CL | https://pdgprod.lbl.gov/pdgprod/pdgLive/DataBlock.action?node=S017EDM | Y | PDG Live S017EDM lists the two not-used rows with the stated comments and confidence levels. |
| Sahoo 2017 interpretation limits `d_n < 2.2e-26`, `d_p < 2.1e-25`, `|theta_bar| < 1.1e-10`, `|tilde d_u - tilde d_d| < 5.5e-27` | https://arxiv.org/abs/1612.09371 | Y | arXiv metadata matches B.K. Sahoo; abstract quotes all four limits. |
| CFW 2008 key reference metadata | https://arxiv.org/abs/0804.1954 | Y | arXiv title/authors match manifest. |

Mismatches: none.

## E007

Verdict: VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| Bishof 2016 Ra-225 direct limit `<1.4e-23 e cm` at 95% CL; improvement factor 36 | https://arxiv.org/abs/1606.04931 | Y | arXiv metadata matches "Improved limit on the 225Ra electric dipole moment"; abstract quotes the limit and factor 36. |
| Parker 2015 first Ra-225 measurement `<5.0e-22 e cm` at 95% confidence | https://arxiv.org/abs/1504.07477 | Y | arXiv metadata matches R.H. Parker et al.; abstract quotes the limit. |
| Sachdeva 2019 Xe-129 result `(1.4 +/- 6.6_stat +/- 2.0_syst)e-28 e cm`; limit `<1.4e-27 e cm` at 95% CL; factor five improvement | https://arxiv.org/abs/1909.12800 | Y | arXiv metadata matches N. Sachdeva et al.; abstract quotes all values. |
| Allmendinger 2019 Xe-129 result `(-4.7 +/- 6.4)e-28 e cm`; limit `<1.5e-27 e cm` at 95% CL | https://arxiv.org/abs/1904.12295 | Y | arXiv metadata matches F. Allmendinger et al.; abstract quotes the measurement and limit. |
| Argonne Ra-225 program context: half-life 15 d, spin 1/2, 2-3 orders sensitivity vs Hg | https://www.phy.anl.gov/mep/atta/research/radiumedm.html | Y | Argonne page text contains all three contextual values. |
| Bishof future sensitivity `1e-28 e cm` statistical and `4e-29 e cm` systematic after upgrades | https://arxiv.org/abs/1606.04931 | Y | arXiv abstract quotes both projection values. |
| DFG Xe-129 upgrade: 2021-2026 term, starting limit `1.5E-27 ecm`, improvement by more than 2 orders, 3He/129Xe clock comparison, SQUID readout | https://gepris.dfg.de/gepris/projekt/455449148?language=en | Y | DFG page contains the project term and all stated upgrade/projection details. |
| KU Leuven RaF/AcF context: 2023-2027 period, Ra/Ac isotopes in dipolar molecules, several-orders Schiff-moment sensitivity enhancement | https://research.kuleuven.be/portal/en/project/3E230783 | Y | KU Leuven page contains the project dates and RaF/AcF Schiff-moment context. |
| CFW 2008 key reference metadata | https://arxiv.org/abs/0804.1954 | Y | arXiv title/authors match manifest. |

Mismatches: none.

## E008

Verdict: VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| Neutron EDM anchor `|d_n| < 1.8e-26 e cm` at 90% CL; table `<0.18` in `10^-25 e cm` | https://pdgprod.lbl.gov/pdgprod/pdgLive/DataBlock.action?node=S017EDM | Y | PDG Live S017EDM lists the Abel 2020 row as `<0.18`, 90% CL. |
| Hg-199 anchor `(2.20 +/- 2.75_stat +/- 1.48_syst)e-30 e cm`; limit `<7.4e-30`; improvement factor 4 | https://arxiv.org/abs/1601.04339 | Y | Graner arXiv abstract quotes the signed result, limit, and factor 4. |
| Pospelov-Ritz qCEDM neutron formula with `(1 +/- 0.5)`, coefficient `1.1`, and combination `tilde d_d + 0.5 tilde d_u` | https://arxiv.org/abs/hep-ph/0010037 | Y | arXiv metadata matches Pospelov and Ritz; abstract quotes the formula. |
| Derived neutron qCEDM bound `<1.6e-26 cm` | https://pdgprod.lbl.gov/pdgprod/pdgLive/DataBlock.action?node=S017EDM; https://arxiv.org/abs/hep-ph/0010037 | Y | Arithmetic checks: `1.8e-26 / 1.1 = 1.636e-26 cm`, rounded to `1.6e-26 cm`. |
| Olive-Pospelov-Ritz-Santoso Hg qCEDM formula with `7e-3 e(tilde d_u - tilde d_d) + 1e-2 d_e + ...` | https://arxiv.org/abs/hep-ph/0506106 | Y | arXiv metadata matches title/authors; local snapshot contains the formula used for the translation. |
| Derived Hg isovector qCEDM bound `<1.1e-27 cm` | https://arxiv.org/abs/1601.04339; https://arxiv.org/abs/hep-ph/0506106 | Y | Arithmetic checks: `7.4e-30 / 7.0e-3 = 1.057e-27 cm`, rounded to `1.1e-27 cm`. |
| CFW context: about 21 TeV generic KK-gluon mass and about 33 TeV pseudo-Goldstone scenario | https://arxiv.org/abs/0804.1954 | Y | arXiv abstract quotes both values; context only, not a measured observable. |
| Koenig-Neubert-Straub context: older neutron limit `2.9e-26 e cm`, scale bounds `Q_uug 1.45 TeV`, `Q_ddg 3.79 TeV` | https://arxiv.org/abs/1403.2756 | Y | arXiv metadata matches title/authors; local snapshot records the table/formula values. |
| Bhattacharya et al. 2024 isovector qCEDM lattice context: four 2+1+1-flavor HISQ ensembles | https://arxiv.org/abs/2404.13516 | Y | arXiv metadata matches title/authors; abstract says the calculation used four ensembles. |

Mismatches: none.

## E009

Verdict: PARTIAL

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| Neutron EDM anchor `|d_n| < 1.8e-26 e cm`, table `<0.18`, measurement `(0.0 +/- 1.1_stat +/- 0.2_sys)e-26 e cm` | https://pdgprod.lbl.gov/pdgprod/pdgLive/DataBlock.action?node=S017EDM | Y | PDG Live S017EDM lists `<0.18`, 90% CL, and the Abel measurement note. |
| Abel 2020 primary experiment measurement and corresponding limit | https://arxiv.org/abs/2001.11966 | Y | arXiv metadata matches C. Abel et al.; abstract quotes the signed measurement, and local snapshot records the `1.8e-26` limit. |
| Pospelov-Ritz Weinberg response `|d_n(w)| approximately e * 22 MeV * w(1 GeV)` | https://arxiv.org/abs/hep-ph/0504231 | Y | arXiv metadata matches Pospelov and Ritz; local snapshot records the quoted formula. |
| Pospelov-Ritz reference scale `w(1 GeV)` | https://arxiv.org/abs/hep-ph/0504231 | Y | Local snapshot records the scale as part of the quoted benchmark formula. |
| Derived `|w(1 GeV)| < 4.1e-11 GeV^-2` | https://pdgprod.lbl.gov/pdgprod/pdgLive/DataBlock.action?node=S017EDM; https://arxiv.org/abs/hep-ph/0504231 | Y | Arithmetic checks: `1.8e-26 cm / 1.97327e-14 = 9.12e-13 GeV^-1`; `9.12e-13 / 0.022 = 4.15e-11 GeV^-2`, rounded to `4.1e-11`. |
| Haisch-Hala O6 response `(d_n/e)_O6 = 74(1 +/- 0.5) MeV` | https://arxiv.org/abs/1909.08955 | Y | arXiv metadata matches Haisch and Hala; local snapshot records the quoted equation. |
| Derived `|C_6| < 1.2e-11 GeV^-2` | https://pdgprod.lbl.gov/pdgprod/pdgLive/DataBlock.action?node=S017EDM; https://arxiv.org/abs/1909.08955 | Y | Arithmetic checks: `9.12e-13 / 0.074 = 1.23e-11 GeV^-2`, rounded to `1.2e-11`. |
| Weinberg 1989 key reference metadata | https://inspirehep.net/literature/26600 | N/A | Direct WebFetch of INSPIRE renders a JavaScript-only app. The same key-reference title/year was cross-checked through the APS page `https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.63.2333` and the local snapshot. |
| Chupp-Ramsey-Musolf global-analysis key reference metadata | https://arxiv.org/abs/1407.1064 | Y | arXiv metadata matches title/authors; abstract confirms global EDM analysis context. |
| Koenig-Neubert-Straub composite/warped context including older `2.9e-26 e cm` neutron limit and about `4 TeV` average fermion-resonance scale | https://arxiv.org/abs/1403.2756 | Y | arXiv metadata matches title/authors; local snapshot records the numerical context. |
| Bhattacharya et al. 2022 lattice context: dimension-4 and dimension-6 gluonic CPV, Weinberg operator mixing | https://arxiv.org/abs/2203.03746 | Y | arXiv metadata and abstract match the manifest and context. |
| CFW 2008 key reference metadata | https://arxiv.org/abs/0804.1954 | Y | arXiv title/authors match manifest. |

Unresolvable URL:

- `https://inspirehep.net/literature/26600`: WebFetch returned only "You need to enable JavaScript to run this app." Alternate APS record and the local snapshot confirm Weinberg 1989 metadata. This is a key-reference metadata caveat, not a numerical mismatch.

Mismatches: none.
