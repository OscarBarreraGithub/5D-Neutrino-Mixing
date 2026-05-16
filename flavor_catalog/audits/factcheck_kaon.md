# Flavor Catalog Fact Check: Kaon Family

Date: 2026-05-16

Agent: factcheck-codex-kaon

Family: kaon

Scope: K001 K002 K003 K004 K005 K006 K008 K009 K010 K012 K013 K017 K018.

Method: read each TeX process file, YAML sidecar, and source manifest; live-fetched every cited process/manifest URL; checked arXiv abs pages for title/author consistency; checked the cited numerical values against the live source when directly visible and against the local source snapshots under `flavor_catalog/references/<id>/`. All cited URLs fetched successfully with HTTP 200. Additional uncited probes of `https://pdg.lbl.gov/2026/listings/rpp2026-list-K-zero-L.pdf` and `https://pdg.lbl.gov/2026/listings/rpp2026-list-K-zero-S.pdf` returned 404; those are not cited catalog URLs, and the cited 2025 listing PDFs / pdgLive API pages remain the resolvable anchors.

## Summary

| process_id | verdict | notes |
|---|---|---|
| K001 | VERIFIED | PDG 2026 epsilon_K, BGS 2020 SM prediction, and FLAG bag values match. |
| K002 | VERIFIED | PDG Delta m_K fits and lattice/supporting values match. |
| K003 | VERIFIED | PDG, KTeV, NA48, RBC/UKQCD, and Aebischer/Buras values match. |
| K004 | VERIFIED | NA62 2026/2025 and SM prediction values match. |
| K005 | VERIFIED | PDG/KOTO limit and rare-kaon program values match. |
| K006 | VERIFIED | PDG K_L -> mu mu value and supporting theory/measurement values match. |
| K008 | VERIFIED | PDG/KTeV/NA48 values and ISU coefficient block match; RS context skipped. |
| K009 | VERIFIED | PDG/KTeV/NA48 values and ISU coefficient block match; RS context skipped. |
| K010 | VERIFIED | PDG/NA48 K_S -> pi0 e e values and event counts match. |
| K012 | VERIFIED | PDG/LHCb limits and displayed theory-context numbers match. |
| K013 | VERIFIED | PDG average, KTeV/NA48 inputs, and NA48/2 related values match. |
| K017 | VERIFIED | PDG API, NA62, KLOE, and SM prediction values match. |
| K018 | VERIFIED | PDG K_l3, FlaviaNet, and FLAG values match. |

## K001

Verdict: VERIFIED

| claim | source_url | verified? | note |
|---|---|---|---|
| `|epsilon| = (2.228 +/- 0.011) x 10^-3`; related `phi_epsilon = (43.5 +/- 0.5) deg`, `A_L = (3.32 +/- 0.06) x 10^-3` | https://pdg.lbl.gov/2026/reviews/rpp2026-rev-cp-viol-kl-decays.pdf | Y | Live PDF and local `pdg2026_epsilon_k.txt` contain these values. |
| `|epsilon_K|_SM = 2.161e-3`, grouped as `2.16(18) x 10^-3` | https://arxiv.org/abs/1911.06822 | Y | arXiv page title/authors match Brod, Gorbahn, Stamou; local extraction contains `2.161` and `2.16(18)`. |
| FLAG `Bhat_K(Nf=2+1)=0.7533(91)`, `B_K(2 GeV)=0.5503(66)`, `B_K(3 GeV)=0.5324(64)`, `Bhat_K(Nf=2+1+1)=0.717(18)(16)` | https://arxiv.org/abs/2411.04268 | Y | arXiv page title/authors match FLAG Review 2024; local FLAG snapshot contains the quoted bag values. |
| CFW 21 TeV / 33 TeV RS context | https://arxiv.org/abs/0804.1954 | N/A | Source exists and metadata matches; treated as contextual theory normalization, not a measured-observable fact-check target. |

Mismatch: none.

Unresolvable: none.

## K002

Verdict: VERIFIED

| claim | source_url | verified? | note |
|---|---|---|---|
| PDG fit assuming CPT: `0.5293 +/- 0.0009` in `10^10 hbar s^-1` | https://pdgprod.lbl.gov/pdgprod/pdgLive/DataBlock.action?home=sumtabM&node=S013D | Y | Live pdgLive DataBlock and local snapshot contain `0.5293 +/- 0.0009`. |
| PDG fit not assuming CPT: `0.5289 +/- 0.0010` | https://pdgprod.lbl.gov/pdgprod/pdgLive/DataBlock.action?home=sumtabM&node=S013D | Y | Same live DataBlock/local snapshot contain `0.5289 +/- 0.0010`. |
| Bai 2014 lattice `Delta M_K = 3.19(41)(96) x 10^-12 MeV` | https://arxiv.org/abs/1406.0916 | Y | arXiv metadata matches Bai et al.; local snapshot contains `3.19(41)(96)`. |
| Wang 2023 statistical error approaches `9%` | https://arxiv.org/abs/2301.01387 | Y | arXiv metadata matches Bigeng Wang; local snapshot contains the `9%` statement. |
| KTeV 2011 `Delta m = (5270 +/- 12) x 10^6 hbar/s` | https://arxiv.org/abs/1011.0127 | Y | arXiv metadata matches KTeV; local snapshot contains the value. |
| Christ 2013 range `6.58(30)` to `11.89(81) x 10^-12 MeV` | https://arxiv.org/abs/1212.5931 | Y | arXiv metadata matches Christ et al.; local snapshot contains both endpoints. |
| Boyle 2024 BSM bag parameters `B1=0.5240(17)(54)`, `B2=0.4794(25)(35)`, `B3=0.746(13)(17)`, `B4=0.897(02)(10)`, `B5=0.6882(78)(94)` | https://arxiv.org/abs/2404.02297 | Y | arXiv metadata matches Boyle et al.; local snapshot contains the displayed bag values. |
| CFW RS context | https://arxiv.org/abs/0804.1954 | N/A | Context only; source metadata matches. |

Mismatch: none.

Unresolvable: none.

## K003

Verdict: VERIFIED

| claim | source_url | verified? | note |
|---|---|---|---|
| PDG listing `Re(epsilon'/epsilon) = 0.00166 +/- 0.00023` | https://pdglive.lbl.gov/Particle.action?home=MXXX020&init=0&node=S013 | Y | Live PDG page fetched; local snapshot contains the exact listing value. |
| PDG DataBlock `OUR FIT: 1.67 +/- 0.23`, `OUR AVERAGE: 1.68 +/- 0.20` in `10^-3` | https://pdgprod.lbl.gov/pdgprod/pdgLiveJson/DataBlock.action?cl=CPV&desig=&init=0&node=S013EPS&parCode=S013EPS | Y | Live DataBlock fetched; local snapshot contains both numbers. |
| KTeV 2011 `(19.2 +/- 2.1) x 10^-4` | https://arxiv.org/abs/1011.0127 | Y | arXiv metadata matches KTeV; local snapshot contains the value. |
| NA48 2002 `(14.7 +/- 2.2) x 10^-4` | https://arxiv.org/abs/hep-ex/0208009 | Y | arXiv metadata matches NA48; local snapshot contains the value. |
| RBC/UKQCD 2020 `21.7(2.6)(6.2)(5.0) x 10^-4` | https://arxiv.org/abs/2004.09440 | Y | arXiv metadata matches RBC/UKQCD; local snapshot contains the value. |
| Aebischer/Buras 2020 octet/nonet SM values `(17.4 +/- 6.1)` and `(13.9 +/- 5.2) x 10^-4` | https://arxiv.org/abs/2005.05978 | Y | arXiv metadata matches; local snapshot contains both values. |
| BSM master/anatomy references | https://arxiv.org/abs/1807.02520 and https://arxiv.org/abs/1808.00466 | Y | arXiv pages fetched and title/authors match the manifest; no catalog value claim from these was mismatched. |
| CFW RS context | https://arxiv.org/abs/0804.1954 | N/A | Context only; source metadata matches. |

Mismatch: none.

Unresolvable: none.

## K004

Verdict: VERIFIED

| claim | source_url | verified? | note |
|---|---|---|---|
| NA62 2016-2024 preliminary `B(K+ -> pi+ nu nubar) = (9.6 +1.9/-1.8) x 10^-11` | https://arxiv.org/abs/2604.12649 | Y | arXiv metadata matches NA62 Moriond contribution; abstract/local snapshot contain the value. |
| NA62 2016-2022 observation `(13.0 +3.3/-3.0) x 10^-11`, 51 candidates, background `18 +3/-2` | https://arxiv.org/abs/2412.12015 | Y | arXiv metadata matches NA62; local snapshot contains branching ratio and event counts. |
| Buras/Venturini SM `(8.60 +/- 0.42) x 10^-11` | https://arxiv.org/abs/2203.10099 | Y | arXiv metadata matches; local snapshot contains the value. |
| Brod/Gorbahn/Stamou SM `7.73(61) x 10^-11` | https://arxiv.org/abs/2105.02868 | Y | arXiv metadata matches; local snapshot contains the value. |
| CFW 21 TeV / 33 TeV RS context | https://arxiv.org/abs/0804.1954 | N/A | Context only; source metadata matches. |

Mismatch: none.

Unresolvable: none.

## K005

Verdict: VERIFIED

| claim | source_url | verified? | note |
|---|---|---|---|
| PDG/KOTO `B(K_L -> pi0 nu nubar) < 2.2e-9` at 90% CL | https://pdgprod.lbl.gov/pdgprod/pdgLive/DataBlock.action?node=S013R40 | Y | Live pdgLive DataBlock and local snapshot contain `<0.22` in `10^-8`, i.e. `<2.2e-9`. |
| KOTO 2025 direct limit `<2.2e-9`, expected background `0.252`, single-event sensitivity `9.33e-10`, zero observed signal-region events | https://arxiv.org/abs/2411.11237 | Y | arXiv metadata matches KOTO PRL; local snapshot contains the limit and supporting numbers. |
| NA62 related charged mode `9.6e-11` with `+1.9/-1.8e-11` | https://arxiv.org/abs/2604.12649 | Y | Same NA62 live arXiv/local snapshot as K004 contains the value. |
| Brod/Gorbahn/Stamou SM `B(K_L -> pi0 nu nubar) = 2.59(29)e-11` | https://arxiv.org/abs/2105.02868 | Y | arXiv metadata matches; local snapshot contains the value. |
| Buras/Venturini SM `(2.94 +/- 0.15)e-11` | https://arxiv.org/abs/2203.10099 | Y | arXiv metadata matches; local snapshot contains the value. |
| Buchalla/Buras 1996 reference | https://arxiv.org/abs/hep-ph/9607447 | Y | arXiv page fetched and metadata matches; no catalog value mismatch found. |

Mismatch: none.

Unresolvable: none.

## K006

Verdict: VERIFIED

| claim | source_url | verified? | note |
|---|---|---|---|
| PDG `BR(K_L -> mu+ mu-) = (6.84 +/- 0.11) x 10^-9` | https://pdg.lbl.gov/2024/listings/rpp2024-list-K-zero-L.pdf | Y | Live 2024 PDF and 2025 K_L listing PDF both contain `6.84 +/- 0.11`; local snapshot matches. |
| BNL E871 `7.18 +/- 0.17 x 10^-9` and normalization ratio `3.474 +/- 0.057 x 10^-6` | https://pubmed.ncbi.nlm.nih.gov/11017525/ | Y | PubMed page fetched; local snapshot contains both abstract values. |
| Isidori/Unterdorfer short-distance bound `<2.5e-9` | https://arxiv.org/abs/hep-ph/0311084 | Y | arXiv metadata matches; local snapshot contains the bound. |
| Dery/Ghosh/Grossman/Schacht hadronic uncertainty below `1%` | https://arxiv.org/abs/2104.06427 | Y | arXiv metadata matches; local snapshot contains the statement. |
| Chao/Christ target lattice accuracy `10%` | https://arxiv.org/abs/2406.07447 | Y | arXiv metadata matches; local snapshot contains the target accuracy. |
| D'Ambrosio/Kitahara 2017 reference | https://arxiv.org/abs/1707.06999 | Y | arXiv page fetched and title/authors match manifest; no catalog value mismatch found. |

Mismatch: none.

Unresolvable: none.

## K008

Verdict: VERIFIED

| claim | source_url | verified? | note |
|---|---|---|---|
| PDG/API `BR(K_L -> pi0 e+ e-) < 2.8e-10` | https://pdgprod.lbl.gov/pdgprod/pdgLiveJson/api/listings/S013.20 | Y | Live API fetched; local snapshot contains `<2.8 E-10`. |
| KTeV standalone `<3.5e-10`, one observed event, background `0.99 +/- 0.35` | https://arxiv.org/abs/hep-ex/0309072 | Y | arXiv abstract contains these values; local snapshot matches. |
| PDG-listed CP-conserving part `(0.0047 +0.0022 -0.0018) E-10` | https://pdgprod.lbl.gov/pdgprod/pdgLiveJson/api/listings/S013.20 | Y | Live API/local snapshot contain the PDG note. |
| Supporting `K_S -> pi0 e+ e-`: `(3.0 +1.5 -1.2 stat +/- 0.2 syst) x 10^-9`, `m_ee > 0.165 GeV`, extrapolated `5.8 +2.9 -2.4 x 10^-9` | https://pdgprod.lbl.gov/pdgprod/pdgLiveJson/api/listings/S012.10 | Y | Live API fetched; local snapshot contains all three values. |
| KTeV theory-context values `3e-11`, `8e-12 to 45e-12`, `3e-12 to 6e-12`, `40%` | https://arxiv.org/abs/hep-ex/0309072 | Y | Local text extraction contains these introduction values. |
| ISU electron coefficients `15.7 +/- 0.3`, `6.2 +/- 0.3`, `2.4 +/- 0.2`, `C_CPC approx 0`, `|a_S| = 1.2 +/- 0.2`, `Im lambda_t/1e-4 = 1.36 +/- 0.12` | https://arxiv.org/abs/hep-ph/0404127 | Y | arXiv metadata matches Isidori/Smith/Unterdorfer; local extraction contains the coefficient block. |
| CFW RS context | https://arxiv.org/abs/0804.1954 | N/A | Context only; source metadata matches. |
| Christ et al. 2016 and Aebischer/Buras/Kumar 2022 references | https://arxiv.org/abs/1608.07585 and https://arxiv.org/abs/2203.09524 | Y | arXiv pages fetched and metadata matches manifest. |

Mismatch: none.

Unresolvable: none.

## K009

Verdict: VERIFIED

| claim | source_url | verified? | note |
|---|---|---|---|
| PDG/API `BR(K_L -> pi0 mu+ mu-) < 3.8e-10` | https://pdgprod.lbl.gov/pdgprod/pdgLiveJson/api/listings/S013.16 | Y | Live API fetched; local snapshot contains `<3.8 x 10^-10`. |
| KTeV 1997-data limit `<3.8e-10`, two candidates, background `0.87 +/- 0.15` | https://arxiv.org/abs/hep-ex/0001006 | Y | arXiv abstract/local snapshot contain these values. |
| Supporting `K_S -> pi0 mu+ mu-` value `(2.9 +1.5 -1.2 stat +/- 0.2 syst) x 10^-9` | https://pdgprod.lbl.gov/pdgprod/pdgLiveJson/api/listings/S012.15 | Y | Live API fetched; local snapshot contains the value. |
| NA48 `K_S -> pi0 mu+ mu-`: six candidates, background `0.22 +0.18 -0.11` | https://arxiv.org/abs/hep-ex/0409011 | Y | arXiv metadata matches NA48/1; local snapshot contains both event counts. |
| ISU muon coefficients `3.7 |a_S|^2`, `1.6 |a_S|`, `1.0`, `5.2`, `|a_S|=1.2`, constructive SM `(1.5 +/- 0.3) x 10^-11` | https://arxiv.org/abs/hep-ph/0404127 | Y | Local extraction contains the formula and coefficient/prediction block. |
| CFW RS context | https://arxiv.org/abs/0804.1954 | N/A | Context only; source metadata matches. |
| Christ et al. 2016 and Aebischer/Buras/Kumar 2022 references | https://arxiv.org/abs/1608.07585 and https://arxiv.org/abs/2203.09524 | Y | arXiv pages fetched and metadata matches manifest. |

Mismatch: none.

Unresolvable: none.

## K010

Verdict: VERIFIED

| claim | source_url | verified? | note |
|---|---|---|---|
| PDG/NA48 partial `BR(K_S -> pi0 e+ e-, m_ee > 0.165 GeV) = (3.0 +1.5 -1.2 stat +/- 0.2 syst) x 10^-9` | https://pdgprod.lbl.gov/pdgprod/pdgLiveJson/api/listings/S012.10 | Y | Live API/local snapshot contain the value. |
| PDG/NA48 full-region extrapolation `5.8 +2.9 -2.4 x 10^-9` | https://pdgprod.lbl.gov/pdgprod/pdgLiveJson/api/listings/S012.10 | Y | Live API/local snapshot contain the footnote value. |
| NA48/1 seven signal candidates and background `0.15` events | https://arxiv.org/abs/hep-ex/0309075 | Y | arXiv abstract/local snapshot contain both numbers and the branching fractions. |
| ISU rare-kaon theory reference | https://arxiv.org/abs/hep-ph/0404127 | Y | arXiv page fetched and metadata matches manifest; no K010 measured value mismatch found. |
| Aebischer/Buras/Kumar 2022 and CFW 2008 references | https://arxiv.org/abs/2203.09524 and https://arxiv.org/abs/0804.1954 | Y / N/A | Rare-kaon white paper metadata matches; CFW is context only. |

Mismatch: none.

Unresolvable: none.

## K012

Verdict: VERIFIED

| claim | source_url | verified? | note |
|---|---|---|---|
| PDG headline `BR(K_S^0 -> mu+ mu-) < 2.1e-10` at 90% CL | https://pdg.lbl.gov/2025/listings/rpp2025-list-K-zero-S.pdf | Y | Live 2025 PDF and local snapshot contain `<2.1 x 10^-10`. |
| LHCb 2020 standalone `<2.2e-10` and combined `<2.1e-10` | https://arxiv.org/abs/2001.10354 | Y | arXiv metadata matches LHCb; local snapshot contains both limits. |
| LHCb 2017 run-1 `<0.8e-9` at 90% and `<1.0e-9` at 95% | https://arxiv.org/abs/1706.00758 | Y | arXiv metadata matches LHCb; local snapshot contains both limits. |
| Chobanova 2018 SM estimate approximately `5e-12` | https://arxiv.org/abs/1711.11030 | Y | arXiv metadata matches; local snapshot contains the SM estimate. |
| Dery/Ghosh/Grossman/Schacht hadronic uncertainty below `1%` | https://arxiv.org/abs/2104.06427 | Y | arXiv metadata matches; local snapshot contains the statement. |

Mismatch: none.

Unresolvable: none.

## K013

Verdict: VERIFIED

| claim | source_url | verified? | note |
|---|---|---|---|
| PDG `BR(K_L -> pi0 gamma gamma) = (1.273 +/- 0.033) x 10^-6` | https://pdg.lbl.gov/2025/listings/rpp2025-list-K-zero-L.pdf | Y | Live 2025 PDF and local snapshot contain the average. |
| PDG-rescaled KTeV input `1.28 +/- 0.06 +/- 0.01 x 10^-6`; original KTeV `(1.29 +/- 0.03 +/- 0.05) x 10^-6`, `a_V=-0.31 +/- 0.05 +/- 0.07` | https://arxiv.org/abs/0805.0031 | Y | arXiv metadata matches KTeV; local snapshots contain PDG-rescaled and original values. |
| PDG-rescaled NA48 input `1.27 +/- 0.04 +/- 0.01 x 10^-6`; original NA48 `(1.36 +/- 0.03 +/- 0.03 +/- 0.03) x 10^-6`, `a_V=-0.46 +/- 0.03 +/- 0.04` | https://iris.unito.it/handle/2318/110341 | Y | IRIS page fetched; local snapshot contains title/authors and original values; PDG local snapshot contains rescaled value. |
| NA48/2 related charged-mode values `0.877 +/- 0.089 x 10^-6` and `0.910 +/- 0.075 x 10^-6` | https://arxiv.org/abs/1310.5499 | Y | arXiv metadata matches NA48/2; local snapshot contains both values. |
| Cappiello/Cata/D'Ambrosio 2018 and Gabbiani/Valencia 2001 references | https://link.springer.com/article/10.1140/epjc/s10052-018-5748-6 and https://arxiv.org/abs/hep-ph/0105006 | Y | Live pages fetched; metadata matches manifest; no catalog value mismatch found. |

Mismatch: none.

Unresolvable: none.

## K017

Verdict: VERIFIED

| claim | source_url | verified? | note |
|---|---|---|---|
| PDG API `R_K = (2.488 +/- 0.009) x 10^-5` | https://pdgapi.lbl.gov/summaries/S010R28 | Y | Live API and local snapshot contain `value_text: (2.488+-0.009)E-5`. |
| NA62 2013 `R_K = (2.488 +/- 0.007_stat +/- 0.007_syst) x 10^-5` | https://arxiv.org/abs/1212.4012 | Y | arXiv metadata matches NA62; local snapshot contains the value. |
| KLOE 2009 `R_K = (2.493 +/- 0.025_stat +/- 0.019_syst) x 10^-5` | https://arxiv.org/abs/0907.3594 | Y | arXiv metadata matches KLOE; local snapshot contains the value. |
| Cirigliano/Rosell SM `R_K^SM = (2.477 +/- 0.001) x 10^-5` | https://arxiv.org/abs/0707.4464 | Y | arXiv metadata matches; local snapshot contains the value. |
| FlaviaNet 2010 and CFW 2008 references | https://arxiv.org/abs/1005.2323 and https://arxiv.org/abs/0804.1954 | Y / N/A | FlaviaNet metadata matches; CFW is context only. |

Mismatch: none.

Unresolvable: none.

## K018

Verdict: VERIFIED

| claim | source_url | verified? | note |
|---|---|---|---|
| PDG `|V_us| f_+(0) = 0.21656(35)` average | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-vud-vus.pdf | Y | Live 2025 and 2026 PDFs and local snapshot contain the average. |
| PDG Table 67.1 mode values: `K+- e3 0.2169(6)`, `K+- mu3 0.2168(10)`, `K_L e3 0.2162(5)`, `K_L mu3 0.2165(6)`, `K_S e3 0.2169(8)`, `K_S mu3 0.2125(47)` | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-vud-vus.pdf | Y | Live PDF text and local snapshot contain all mode values. |
| FlaviaNet 2010 `|V_us| f_+(0) = 0.2163(5)`, fit quality `chi2/ndf = 0.77/4`, `r_mu e = 1.002(5)` | https://arxiv.org/abs/1005.2323 | Y | arXiv metadata matches FlaviaNet; local snapshot contains the values. |
| PDG `|V_us| = 0.22330(53)` from K_l3 and FLAG `f_+(0)` | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-vud-vus.pdf | Y | Live PDF/local snapshot contain the equation and uncertainty breakdown. |
| FLAG 2024 `f_+(0)=0.9698(17)` and derived `|V_us|=0.22328(58)` | https://arxiv.org/abs/2411.04268 | Y | arXiv metadata matches FLAG Review 2024; local snapshot contains both values. |
| CFW 21 TeV / 33 TeV RS context | https://arxiv.org/abs/0804.1954 | N/A | Context only; source metadata matches. |

Mismatch: none.

Unresolvable: none.
