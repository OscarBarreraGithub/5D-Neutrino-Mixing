# Beauty Family Fact-Check Audit

Agent: `factcheck-codex-beauty`
Date: 2026-05-16
Branch: `flavor-catalog/2026q2`
Family: beauty
Processes: B001 B002 B003 B004 B005 B006 B009 B011 B012 B015 B016 B017 B018 B019 B021 B022 B023 B025 B026 B032 B033 B034

## Summary

All 22 beauty-family process entries were independently checked against the cited web source or the local source snapshot under `flavor_catalog/references/<process_id>/`. I re-fetched all cited arXiv abstract pages used as key references and checked that the author/title metadata matched the catalog citation at the level needed for identification. Measured observables in `pdg_or_equivalent` and measured-observable supporting blocks were checked for the quoted numerical values. Quoted SM predictions or theory paper numerical predictions were checked where the catalog claims to quote them. Theory normalization constants and RS-flavor policy-scale estimates were skipped per the L001 policy carve-out.

Verdict: 22 VERIFIED / 0 MISMATCH / 0 UNRESOLVABLE.

## B001 - VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| HFLAV/PDG 2025 `Delta m_d = 0.5069 +/- 0.0019 ps^-1` | https://hflav-eos.web.cern.ch/hflav-eos/osc/PDG_2025/ | Y | Fetched page contains the quoted `Delta m_d` value. |
| Companion `x_d = 0.7697 +/- 0.0035`, `chi_d = 0.1860 +/- 0.0011` | https://hflav-eos.web.cern.ch/hflav-eos/osc/PDG_2025/ | Y | Same HFLAV table contains both companion values. |
| Belle II `Delta m_d = 0.516 +/- 0.008 +/- 0.005 ps^-1`, `190 fb^-1` | https://arxiv.org/abs/2302.12791 | Y | arXiv abstract and local snapshot contain the value and luminosity. |
| Key references: Belle II 2023, FLAG 2024, CFW 2008 | https://arxiv.org/abs/2302.12791 ; https://arxiv.org/abs/2411.04268 ; https://arxiv.org/abs/0804.1954 | Y | Authors/titles match the cited works; CFW scale values treated as theory-normalization policy context. |

## B002 - VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| HFLAV summer 2025 all-charmonium `sin(2 beta) = 0.710 +/- 0.011` | https://hflav-eos.web.cern.ch/hflav-eos/triangle/summer2025/ | Y | Fetched HFLAV page contains the quoted average. |
| HFLAV `J/psi K_S` mode `0.712 +/- 0.011`; beta ambiguity `22.63 deg`/`67.37 deg` | https://hflav-eos.web.cern.ch/hflav-eos/triangle/summer2025/ | Y | Values appear in the same HFLAV fit table. |
| LHCb Run 2 `S = 0.717 +/- 0.013 +/- 0.008`, `C = 0.008 ...` | https://arxiv.org/abs/2309.09728 | Y | arXiv metadata matched; local snapshot contains the quoted values. |
| Belle II `S = 0.724 +/- 0.035 +/- 0.009`, `C = -0.035 ...` | https://arxiv.org/abs/2402.17260 | Y | arXiv metadata matched; local snapshot contains the quoted values. |
| Frings et al. penguin-shift bound `0.68 deg` | https://arxiv.org/abs/1503.00859 | Y | Value found in the local snapshot for the cited paper. |
| Key references | https://arxiv.org/abs/0809.0842 ; https://arxiv.org/abs/0804.1954 | Y | Authors/titles match; CFW scale context skipped as theory normalization. |

## B003 - VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| HFLAV 2024 `Delta m_s = 17.766 +/- 0.006 ps^-1` | https://hflav-eos.web.cern.ch/hflav-eos/osc/HFLAV_2024/ | Y | Fetched HFLAV page contains the quoted value. |
| HFLAV companions `x_s = 26.94 +/- 0.10`, `chi_s = 0.499314 +/- 0.000005` | https://hflav-eos.web.cern.ch/hflav-eos/osc/HFLAV_2024/ | Y | Both companion values appear in the fetched HFLAV table. |
| HFLAV companions `1/Gamma_s = 1.516 +/- 0.006 ps`, `DeltaGamma_s = 0.0781 +/- 0.0035 ps^-1` | https://hflav-eos.web.cern.ch/hflav-eos/osc/HFLAV_2024/ | Y | Values appear in the fetched HFLAV table. |
| PDG 2025 `Delta m_s = 17.765 +/- 0.006 ps^-1` | https://hflav-eos.web.cern.ch/hflav-eos/osc/PDG_2025/ | Y | Current HFLAV/PDG page contains the close updated value. |
| LHCb 2021 `17.7683 +/- 0.0051 +/- 0.0032 ps^-1`; LHCb combination `17.7656 +/- 0.0057` | https://arxiv.org/abs/2104.04421 | Y | arXiv abstract/local snapshot contain the measurement and combination. |
| Key references: LHCb, FLAG, CFW | https://arxiv.org/abs/2104.04421 ; https://arxiv.org/abs/2411.04268 ; https://arxiv.org/abs/0804.1954 | Y | Authors/titles match; CFW theory scales skipped. |

## B004 - VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| HFLAV `phi_s = -0.041 +/- 0.016 rad` all-combined | https://hflav-eos.web.cern.ch/hflav-eos/osc/PDG_2025/HFLAV_phis_inputs/HFLAV_phis_inputs.pdf | Y | Fetched PDF text contains all-combined `-0.041 +/- 0.016`. |
| HFLAV `DeltaGamma_s = +0.079 +/- 0.006 ps^-1`; `Bs -> J/psi K K` component `-0.050 +/- 0.017` | https://hflav-eos.web.cern.ch/hflav-eos/osc/PDG_2025/HFLAV_phis_inputs/HFLAV_phis_inputs.pdf | Y | Fetched PDF text contains these companion values. |
| LHCb Run 2 value and SM reference `phi_s^SM = -0.0368...` | https://arxiv.org/abs/2308.01468 | Y | arXiv metadata matched; local snapshot contains quoted measurement/SM reference. |
| Key references | https://arxiv.org/abs/2308.01468 ; https://arxiv.org/abs/0804.1954 | Y | Authors/titles match; CFW theory scale context skipped. |

## B005 - VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| PDG `B(B_s -> mu+ mu-) = (3.34 +/- 0.27) x 10^-9` | https://pdgprod.lbl.gov/pdgprod/pdgLiveJson/api/listings/S086.15 | Y | API/listing source and local snapshot contain the quoted average. |
| HFLAV 2023 average `3.45 +/- 0.29 x 10^-9` | https://hflav-eos.web.cern.ch/hflav-eos/rare/Apr2023/html/Bs/BR_Bs0_mu%2B_mu-.html | Y | Fetched/source snapshot contains the quoted HFLAV average. |
| SM prediction `3.64 +/- 0.12 x 10^-9` | https://arxiv.org/abs/2407.03810 | Y | Czaja/Misiak/Rehman/Steinhauser source snapshot contains the quoted SM value. |
| CMS, LHCb, ATLAS individual inputs | https://arxiv.org/abs/2212.10311 ; https://arxiv.org/abs/2108.09283 ; https://arxiv.org/abs/1812.03017 | Y | Each arXiv metadata page matched; local snapshots contain the claimed input values. |
| Key references | https://arxiv.org/abs/0804.1954 | Y | CFW citation metadata matched; theory scale context skipped. |

## B006 - VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| PDG `B(B0 -> mu+ mu-) < 1.5 x 10^-10` and 95% note `<1.9 x 10^-10` | https://pdgprod.lbl.gov/pdgprod/pdgLiveJson/api/listings/S042.7 | Y | API/listing source contains the quoted limits. |
| HFLAV ratio `B(B0 -> mu mu)/B(B_s -> mu mu) < 0.081` at 90% CL; 95% `<0.095` | https://hflav-eos.web.cern.ch/hflav-eos/rare/Apr2023/html/radll/RBR/BR_B0_mu%2B_mu-..BR_Bs0_mu%2B_mu-.html | Y | Fetched HFLAV ratio page/snapshot contains both limits. |
| SM prediction `1.06 +/- 0.09 x 10^-10` | https://arxiv.org/abs/1311.0903 | Y | Bobeth et al. snapshot contains the quoted SM value. |
| CMS/LHCb/ATLAS supporting inputs | https://arxiv.org/abs/2212.10311 ; https://arxiv.org/abs/2108.09283 ; https://arxiv.org/abs/1812.03017 | Y | arXiv metadata matched and numerical claims were found in snapshots. |

## B009 - VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| HFLAV direct average `B(B+ -> tau+ nu) = (1.12 +/- 0.19) x 10^-4` | https://hflav-eos.web.cern.ch/hflav-eos/rare/Dec2025/html/radll/Leptonic/BR_B+_tau+_nutau.html | Y | Fetched HFLAV page contains the average and component measurements. |
| PDG summary `1.09 +0.25 -0.24 x 10^-4` | https://pdgapi.lbl.gov/summaries/S041.184 | Y | PDG summary endpoint contains the quoted value. |
| UTfit SM indirect `0.865 +/- 0.041 x 10^-4` | https://www.utfit.org/foswiki/bin/view/UTfit/ResultsSummer2024SM | Y | Source fetched with TLS workaround; page contains the quoted UTfit value. |
| Belle, Belle II, BaBar supporting measurements | https://arxiv.org/abs/1208.4678 ; https://arxiv.org/abs/1503.05613 ; https://arxiv.org/abs/2502.04885 ; https://arxiv.org/abs/1207.0698 ; https://arxiv.org/abs/0912.2453 | Y | arXiv metadata matched and supporting values were present in snapshots/HFLAV input table. |

## B011 - VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| HFLAV inclusive `B -> X_s gamma = (3.49 +/- 0.19) x 10^-4` | https://hflav.web.cern.ch/results | Y | Fetched HFLAV results page contains the quoted average. |
| Supporting HFLAV/PDG radiative-decay sources | https://hflav-eos.web.cern.ch/hflav-eos/rare/Dec2024/hflav_rare_decays_Dec24.pdf ; https://pdg.lbl.gov/2024/reviews/rpp2024-rev-b-meson-prod-decay.pdf | Y | Local snapshots match the catalog's quoted average context. |
| SM prediction `3.40 +/- 0.17 x 10^-4` | https://arxiv.org/abs/2002.01548 | Y | Misiak/Rehman/Steinhauser snapshot contains the quoted SM value. |
| Older SM comparison `3.36 +/- 0.23 x 10^-4`; Belle II projection/source | https://arxiv.org/abs/1503.01789 ; https://arxiv.org/abs/2210.10220 ; https://arxiv.org/abs/1608.02344 | Y | Metadata matched; numerical values/projection context verified where claimed. |

## B012 - VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| HFLAV Dec. 2024 `B(B0 -> K*(892)0 gamma) = (41.63 +/- 0.92) x 10^-6` | https://hflav-eos.web.cern.ch/hflav-eos/rare/Dec2024/hflav_rare_decays_Dec24.pdf | Y | Fetched HFLAV PDF Table 62 contains the quoted average; local snapshot agrees. |
| HFLAV Dec. 2024 `B(B+ -> K*(892)+ gamma) = (39.5 +/- 1.0) x 10^-6` | https://hflav-eos.web.cern.ch/hflav-eos/rare/Dec2024/hflav_rare_decays_Dec24.pdf | Y | Fetched HFLAV PDF Table 57 contains the quoted average; local snapshot agrees. |
| HFLAV Dec. 2024 `Delta_0+(B -> K* gamma) = 0.059 +/- 0.014` | https://hflav-eos.web.cern.ch/hflav-eos/rare/Dec2024/hflav_rare_decays_Dec24.pdf | Y | Fetched HFLAV PDF text and local snapshot contain the quoted value and the HFLAV caution that Table 76 isospin averages are naive. |
| HFLAV Dec. 2024 `A_CP(B -> K* gamma) = -0.004 +/- 0.011` | https://hflav-eos.web.cern.ch/hflav-eos/rare/Dec2024/hflav_rare_decays_Dec24.pdf | Y | Fetched HFLAV PDF Table 96 contains Belle/BaBar inputs and the `-0.004 +/- 0.011` average. |
| PDG pdgLive 2025 `S_{K*(892)0 gamma} = -0.08 +/- 0.17` | https://pdglive.lbl.gov/DataBlock.action?node=S042SX4 | Y | Fetched pdgLive HTML contains S042SX4 and `OUR AVERAGE -0.08 +/- 0.17`; local snapshot agrees. |
| PDG pdgLive 2025 `C_{K*(892)0 gamma} = 0.03 +/- 0.10` | https://pdglive.lbl.gov/Particle.action?node=S042 | Y | Fetched current B0 listing contains `C_{K*(892)0 gamma}` with `0.03 +/- 0.10`; local snapshot agrees. |
| Belle II 2024 supporting measurement: 2019-2022 data, `365 fb^-1`, `(387 +/- 6) x 10^6 Upsilon(4S)` events, and quoted B -> K* gamma values | https://arxiv.org/abs/2411.10127 | Y | arXiv metadata matches Belle II/I. Adachi; the abstract contains the dataset and the quoted branching-fraction/CP/isospin values. |
| Belle 2017 support: evidence for isospin violation, `3.1 sigma`, `772 x 10^6 BBbar` pairs | https://arxiv.org/abs/1707.00394 | Y | arXiv metadata matches T. Horiguchi et al. (Belle); abstract contains the significance, dataset, and observable statement. |
| BaBar 2009 support: `383 million BBbar`, `BF(B0)=4.47... x 10^-5`, `BF(B+)=4.22... x 10^-5`, CP/isospin intervals | https://arxiv.org/abs/0906.2177 | Y | arXiv metadata matches B. Aubert et al. (BaBar); abstract contains the quoted values. |
| LHCb 2014 support: photon-polarization observation in `b -> s gamma`, `5.2 sigma` | https://arxiv.org/abs/1402.6852 | Y | arXiv metadata matches R. Aaij et al. (LHCb); abstract contains the `5.2 sigma` observation statement. |
| CFW 2008 key reference metadata | https://arxiv.org/abs/0804.1954 | Y | arXiv title/authors match Csaki, Falkowski, and Weiler; RS flavor-scale numerals treated as theory context under the L001 carve-out. |

## B015 - VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| HFLAV inclusive `B -> X_s l+l- = (5.84 +/- 0.69) x 10^-6` | https://hflav-eos.web.cern.ch/hflav-eos/rare/Dec2024/html/radll/Badmix/BR_Badmix_Xs_l+_l-.html | Y | Fetched HFLAV page contains the average and PDG/input values. |
| BaBar `6.73... x 10^-6`; Belle `4.11... x 10^-6` supporting inputs | https://hflav-eos.web.cern.ch/hflav-eos/rare/Dec2024/html/radll/Badmix/BR_Badmix_Xs_l+_l-.html | Y | Values appear in the fetched HFLAV input table. |
| Theory low/high-q2 predictions `1.58 +/- 0.37`, `0.48 +/- 0.10`, `[1,6] ee = 1.67 +/- 0.10`, `[1,6] mumu = 1.62 +/- 0.09` | https://arxiv.org/abs/1312.5364 ; https://arxiv.org/abs/hep-ex/0503044 ; https://arxiv.org/abs/1503.04849 | Y | arXiv metadata matched; local snapshots contain the quoted numerical predictions. |
| Belle II prospect `50 ab^-1` | https://arxiv.org/abs/1503.04849 | Y | Projection source exists and snapshot contains the luminosity/prospect context. |

## B016 - VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| HFLAV `B+ -> K+ l+l- = (5.76 +/- 0.40) x 10^-7` | https://hflav-eos.web.cern.ch/hflav-eos/rare/Dec2025/html/radll/Bu/BR_B+_K+_l+_l-.html | Y | Fetched HFLAV page contains the quoted average. |
| HFLAV `B0 -> K0 l+l- = (3.28 +/- 0.32) x 10^-7` | https://hflav-eos.web.cern.ch/hflav-eos/rare/Dec2025/html/radll/Bd/BR_B0_K0_l+_l-.html | Y | Fetched HFLAV page contains the quoted average. |
| Key references | https://arxiv.org/abs/1908.01848 ; https://arxiv.org/abs/1403.8044 ; https://arxiv.org/abs/0807.4119 ; https://arxiv.org/abs/0804.1954 | Y | arXiv metadata matched; CFW theory scale context skipped. |

## B017 - VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| HFLAV inclusive `B -> X_s l+l- = (5.84 +/- 0.69) x 10^-6` | https://hflav-eos.web.cern.ch/hflav-eos/rare/Dec2025/html/radll/Badmix/BR_Badmix_Xs_l+_l-.html | Y | Fetched HFLAV page contains the quoted average. |
| HFLAV `B+ -> K+ l+l- = (5.76 +/- 0.40) x 10^-7` | https://hflav-eos.web.cern.ch/hflav-eos/rare/Dec2025/html/radll/Bu/BR_B+_K+_l+_l-.html | Y | Fetched HFLAV page contains the quoted average. |
| HFLAV `B0 -> K*0 l+l- = (9.9 +/- 1.2) x 10^-7` | https://hflav-eos.web.cern.ch/hflav-eos/rare/Dec2025/html/radll/Bd/BR_B0_K*0_l+_l-.html | Y | Fetched HFLAV page contains the quoted average. |
| LHCb LFU ratios `R_K`, `R_K*`, and SM-compatibility statement | https://arxiv.org/abs/2212.09153 | Y | arXiv metadata matched; local snapshot contains the quoted ratios and `0.2 sigma` statement. |
| Key references | https://arxiv.org/abs/2003.04831 ; https://arxiv.org/abs/2305.19038 ; https://arxiv.org/abs/0804.1954 | Y | Metadata matched; CFW theory scale context skipped. |

## B018 - VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| HFLAV low-bin `R_K = 0.994 +/- 0.090` | https://hflav-eos.web.cern.ch/hflav-eos/rare/Dec2025/html/radll/RBR/BR_B%2B_K%2B_mu%2B_mu-..B%2B_K%2B_e%2B_e-__bin1.html | Y | Fetched HFLAV page contains the quoted low-bin value. |
| HFLAV central-bin `R_K = 0.947 +/- 0.047` | https://hflav-eos.web.cern.ch/hflav-eos/rare/Dec2025/html/radll/RBR/BR_B%2B_K%2B_mu%2B_mu-..B%2B_K%2B_e%2B_e-__bin2b.html | Y | Fetched HFLAV page contains the quoted central-bin value. |
| HFLAV full/Belle value `1.08 +/- 0.16` and input values | https://hflav-eos.web.cern.ch/hflav-eos/rare/Dec2025/html/radll/RBR/BR_B%2B_K%2B_mu%2B_mu-..BR_B%2B_K%2B_e%2B_e-.html | Y | Fetched HFLAV page contains the quoted value and input table. |
| LHCb 2021 `R_K = 0.846 ...` and `3.1 sigma` | https://arxiv.org/abs/2105.10303 | Y | Value appears in source/snapshot. Metadata note: arXiv author line is Davide Lancierini, but title/value identify the LHCb `Updated measurement of R_K` source used by the catalog. |
| Key references | https://arxiv.org/abs/2212.09153 ; https://arxiv.org/abs/1605.07633 ; https://arxiv.org/abs/0804.1954 | Y | Metadata matched; CFW theory scale context skipped. |

## B019 - VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| HFLAV low-bin `R_K* = 0.927 +/- 0.097` | https://hflav-eos.web.cern.ch/hflav-eos/rare/Dec2025/html/radll/RBR/BR_B0_K*0_mu+_mu-..BR_B0_K*0_e+_e-__bin1b.html | Y | Fetched HFLAV page contains the quoted low-bin value. |
| HFLAV central-bin `R_K* = 1.028 +/- 0.074` | https://hflav-eos.web.cern.ch/hflav-eos/rare/Dec2025/html/radll/RBR/BR_B0_K*0_mu+_mu-..BR_B0_K*0_e+_e-__bin2.html | Y | Fetched HFLAV page contains the quoted central-bin value. |
| LHCb 2017 low/central values `0.66...`, `0.69...` and significance ranges | https://arxiv.org/abs/1705.05802 | Y | arXiv metadata matched; local snapshot contains the quoted numbers. |
| 2022 LHCb update and SM theory values from BIP 2016 | https://arxiv.org/abs/2212.09153 ; https://arxiv.org/abs/1605.07633 | Y | arXiv metadata matched; source snapshots contain the quoted ratios/predictions. |
| Key references | https://arxiv.org/abs/0804.1954 | Y | CFW metadata matched; theory scale context skipped. |

## B021 - VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| PDG summary `B(Lambda_b -> Lambda mu+ mu-) = (1.08 +/- 0.28) x 10^-6` | https://pdgapi.lbl.gov/summaries/S040.26 | Y | PDG summary endpoint contains the average; listing URL also resolves. |
| LHCb `dB/dq2 = (1.18... +/- 0.27) x 10^-7`, `A_FB^l = -0.05 ...`, `A_FB^h = -0.29 ...`, `3.0 fb^-1` | https://arxiv.org/abs/1503.07138 | Y | arXiv abstract contains the quoted measurement values. |
| CDF 2011 `24` signal events, `5.8 sigma`, `6.8 fb^-1`, `1.73 +/- 0.42 +/- 0.55 x 10^-6` | https://arxiv.org/abs/1107.3753 | Y | arXiv metadata matched; local snapshot contains the quoted observation values. |
| Key references | https://arxiv.org/abs/1212.4827 ; https://arxiv.org/abs/0804.1954 | Y | Metadata matched; CFW theory scale context skipped. |

## B022 - VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| HFLAV `B(B+ -> K+ nu nubar) = (1.38 +/- 0.35) x 10^-5` | https://hflav-eos.web.cern.ch/hflav-eos/rare/Dec2025/html/radll/Bu/BR_B+_K+_nu_nubar.html | Y | Fetched HFLAV page contains the quoted average. |
| PDG/Belle II `2.3 +/- 0.5 +0.5 -0.4 x 10^-5` | https://pdgprod.lbl.gov/pdgprod/pdgLiveJson/api/listings/S041.273 ; https://docs.belle2.org/pub_data/documents/4404/ | Y | PDG API and Belle II document page contain the quoted value. |
| Belle II `362 fb^-1`, inclusive `2.7 +/- 0.5 +/- 0.5`, evidence `3.5 sigma`, SM tension `2.7 sigma` | https://docs.belle2.org/pub_data/documents/4404/ | Y | Belle II page contains the listed dataset, values, and significance statements. |
| BaBar limit `<3.7 x 10^-5`; HPQCD SM `5.58(37) x 10^-6` | https://export.arxiv.org/api/query?id_list=1303.7465 ; https://arxiv.org/abs/2207.13371 | Y | Metadata and local snapshots contain the quoted limit/prediction. |
| Belle II prospect/reference | https://docs.belle2.org/pub_data/publications/4277/ | Y | Collaboration page resolves and matches the cited projection context. |

## B023 - VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| PDG/HFLAV `B0 -> K*0 nu nubar < 1.8 x 10^-5` at 90% CL | https://pdgapi.lbl.gov/summaries/S042.152 ; https://hflav-eos.web.cern.ch/hflav-eos/rare/Dec2025/html/radll/Bd/BR_B0_K%2A0_nu_nubar.html | Y | PDG summary and HFLAV page contain the quoted limit. |
| PDG/HFLAV `B+ -> K*+ nu nubar < 4.0 x 10^-5` at 90% CL | https://pdgapi.lbl.gov/summaries/S041.490 ; https://hflav-eos.web.cern.ch/hflav-eos/rare/Dec2025/html/radll/Bu/BR_B%2B_K%2A%2B_nu_nubar.html | Y | PDG summary and HFLAV page contain the quoted limit. |
| Belle 2017 combined vector limit `<2.7 x 10^-5` | https://arxiv.org/abs/1702.03224 | Y | arXiv metadata matched; local snapshot contains the quoted limit. |
| Buras et al. SM `B0 -> K*0 nu nubar = 9.2 +/- 1.0 x 10^-6` | https://arxiv.org/abs/1409.4557 | Y | Source snapshot contains the quoted SM prediction. |
| Key references | https://arxiv.org/abs/1303.3719 ; https://arxiv.org/abs/1303.7465 ; https://arxiv.org/abs/0804.1954 | Y | Metadata matched; CFW theory scale context skipped. |

## B025 - VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| HFLAV CKM2025 `R(D) = 0.358 +/- 0.024`; companion `R(D*) = 0.281 +/- 0.011` | https://hflav-eos.web.cern.ch/hflav-eos/semi/ckm25/html/RDsDsstar/RDRDs.html | Y | Fetched HFLAV page/snapshot contains both averages. |
| HFLAV SM references `R(D) = 0.296 +/- 0.004`, `R(D*) = 0.254 +/- 0.005`; tension `p = 1.48e-4`, about `3.8 sigma` | https://hflav-eos.web.cern.ch/hflav-eos/semi/ckm25/html/RDsDsstar/RDRDs.html | Y | Same HFLAV page/snapshot contains the quoted SM reference and tension. |
| HFLAV fit detail `chi2/dof = 16.683/14` | https://hflav-eos.web.cern.ch/hflav-eos/semi/ckm25/r_dtaunu/logfile_ckm25_new.txt | Y | HFLAV logfile source contains the fit detail. |
| Selected BaBar, Belle, LHCb, Belle II inputs | https://export.arxiv.org/api/query?id_list=1205.5442 ; https://export.arxiv.org/api/query?id_list=1303.0571 ; https://export.arxiv.org/api/query?id_list=1507.03233 ; https://export.arxiv.org/api/query?id_list=1910.05864 ; https://export.arxiv.org/api/query?id_list=2302.02886 ; https://export.arxiv.org/api/query?id_list=2406.03387 ; https://export.arxiv.org/api/query?id_list=2504.11220 | Y | Metadata matched and source snapshots/HFLAV table contain the claimed input values. |
| Key references: FLAG/CFW | https://export.arxiv.org/api/query?id_list=2411.04268 ; https://export.arxiv.org/api/query?id_list=0804.1954 | Y | Metadata matched; CFW theory scale context skipped. |

## B026 - VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| HFLAV CKM2025 `R(D*) = 0.281 +/- 0.011`; companion `R(D) = 0.358 +/- 0.024` | https://hflav-eos.web.cern.ch/hflav-eos/semi/ckm25/html/RDsDsstar/RDRDs.html | Y | Fetched HFLAV page/snapshot contains both averages. |
| HFLAV SM references `R(D*) = 0.254 +/- 0.005`, `R(D) = 0.296 +/- 0.004`; tension `p = 1.48e-4` | https://hflav-eos.web.cern.ch/hflav-eos/semi/ckm25/html/RDsDsstar/RDRDs.html | Y | Same HFLAV page/snapshot contains the quoted SM reference and tension. |
| HFLAV fit detail `chi2/dof = 16.683/14` | https://hflav-eos.web.cern.ch/hflav-eos/semi/ckm25/r_dtaunu/logfile_ckm25_new.txt | Y | HFLAV logfile source contains the fit detail. |
| Selected BaBar, Belle, LHCb, Belle II inputs | https://arxiv.org/abs/1205.5442 ; https://arxiv.org/abs/1303.0571 ; https://arxiv.org/abs/1507.03233 ; https://arxiv.org/abs/1612.00529 ; https://arxiv.org/abs/1709.00129 ; https://arxiv.org/abs/1910.05864 ; https://export.arxiv.org/api/query?id_list=2302.02886 ; https://arxiv.org/abs/2305.01463 ; https://arxiv.org/abs/2406.03387 ; https://arxiv.org/abs/2401.02840 ; https://arxiv.org/abs/2504.11220 | Y | Metadata matched and source snapshots/HFLAV table contain the claimed input values. |
| Key references: FLAG/CFW | https://arxiv.org/abs/2411.04268 ; https://arxiv.org/abs/0804.1954 | Y | Metadata matched; CFW theory scale context skipped. |

## B032 - VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| HFLAV branching fractions: `2.392 +/- 0.062`, `1.322 +/- 0.043`, `2.007 +/- 0.040`, `1.012 +/- 0.043` in units of `10^-5` | https://hflav-eos.web.cern.ch/hflav-eos/rare/Dec2025/html/charmless/Bu/BR_B%2B_K0_pi%2B.html ; https://hflav-eos.web.cern.ch/hflav-eos/rare/Dec2025/html/charmless/Bu/BR_B%2B_K%2B_pi0.html ; https://hflav-eos.web.cern.ch/hflav-eos/rare/Dec2025/html/charmless/Bd/BR_B0_K%2B_pi-.html ; https://hflav-eos.web.cern.ch/hflav-eos/rare/Dec2025/html/charmless/Bd/BR_B0_K0_pi0.html | Y | Each fetched HFLAV page contains the quoted branching-fraction average. |
| HFLAV direct CP values `-2.67 +/- 0.87%`, `+2.7 +/- 1.2%`, `-8.31 +/- 0.31%`, `-1 +/- 13%` | https://hflav.web.cern.ch/content/rare-b-decays | Y | Verified from the local HFLAV rare-decays snapshot for the cited general rare-decays source. |
| PDG `C_{K pi} = 0.00 +/- 0.08`, `S_{K pi} = 0.64 +/- 0.13` | https://pdgapi.lbl.gov/summaries/S042CKP ; https://pdgapi.lbl.gov/summaries/S042SKP | Y | PDG summary endpoints contain the quoted CP parameters. |
| LHCb/Belle II sum-rule supporting values | https://arxiv.org/abs/2012.12789 ; https://arxiv.org/abs/2310.06381 | Y | arXiv metadata matched; snapshots contain `0.025 +/- 0.015 +/- 0.006 +/- 0.003` and sum-rule `-0.03 +/- 0.13 +/- 0.04`. |
| Key references | https://arxiv.org/abs/2305.07555 ; https://arxiv.org/abs/0804.1954 | Y | Metadata matched; CFW theory scale context skipped. |

## B033 - VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| HFLAV `S(phi K0) = 0.74 +/- 0.12` and `C(phi K0) = -0.09 +/- 0.12` | https://hflav-eos.web.cern.ch/hflav-eos/triangle/summer2025/ | Y | Fetched HFLAV page contains the `phi K0` averages. |
| HFLAV all-charmonium comparator `sin(2 beta) = 0.710 +/- 0.011` | https://hflav-eos.web.cern.ch/hflav-eos/triangle/summer2025/ | Y | Same HFLAV source contains the quoted benchmark average. |
| Belle II `S = 0.54...`, `A = 0.31...`, `162 +/- 17` events, `387 +/- 6` million B pairs | https://arxiv.org/abs/2307.02802 | Y | arXiv metadata matched; source snapshot and HFLAV input table contain the quoted values. |
| Key references | https://arxiv.org/abs/1201.5897 ; https://arxiv.org/abs/1007.3848 ; https://arxiv.org/abs/0804.1954 | Y | Metadata matched; CFW theory scale context skipped. |

## B034 - VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| PDG/LHCb CP values `phi_s = -0.074 +/- 0.069`, `|lambda| = 1.009 +/- 0.030` | https://pdg.lbl.gov/encoder_listings/s086.pdf ; https://cds.cern.ch/record/2857424 | Y | CDS record contains the values directly; local PDG snapshot records the PDG listing values from the cited PDF. |
| LHCb 2023 component `phi_s = -0.042 +/- 0.075 +/- 0.009`, `|lambda| = 1.004 +/- 0.030 +/- 0.009`, `6 fb^-1` | https://cds.cern.ch/record/2857424 | Y | Fetched CDS record contains the quoted measurement, combination, and dataset. |
| HFLAV `B(B_s -> phi phi) = (1.85 +/- 0.14) x 10^-5` | https://hflav-eos.web.cern.ch/hflav-eos/rare/May2024/html/Bs/BR_Bs0_phi_phi.html | Y | Fetched HFLAV page contains the average and input values. |
| Earlier LHCb/CDF supporting values and `B0 -> phi phi < 2.7 x 10^-8` | https://arxiv.org/abs/1407.2222 ; https://cds.cern.ch/record/2684085 | Y | Metadata/source snapshots contain the quoted earlier measurements and limit. |
| Key references | https://arxiv.org/abs/0804.1954 | Y | CFW metadata matched; theory scale context skipped. |

## Discrepancies

None.

## Unresolvable Sources

None. UTfit required a TLS/certificate workaround (`curl -k`) but the page content was fetched and the quoted value was visible. The PDG `s086.pdf` value for B034 was verified against the local PDG snapshot and the CDS record carrying the same LHCb combined values.
