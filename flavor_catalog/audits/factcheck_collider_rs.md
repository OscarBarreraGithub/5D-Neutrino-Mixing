# Fact-check report: collider_rs

Date: 2026-05-17T18:44:20-04:00
Agent: factcheck-codex-collider_rs-w9
Family: collider_rs (PRIMARY tier; Wave-9 additions)

## Summary

I checked the 14 Wave-9 `collider_rs` entries against their TeX files, YAML sidecars, source manifests, live cited URLs, and local snapshots under `flavor_catalog/references/<process_id>/`. Live URL checks returned HTTP 200 for all cited arXiv, PDG/PdgLive, HEPData, CMS public-result, and CMS table URLs. For arXiv pages I verified title, author/collaboration, and year metadata against the manifests; for PDG PDFs and experiment public pages I used the fetched URL plus the tracked text snapshot for value-level checking. L001 carve-outs were applied to dataset metadata, local quark-scan comparison scales, theory normalizations, future/projection reach, and EFT translations; those are marked `N/A` rather than mismatches.

| Process | Verdict | Mismatches | Note |
|---|---:|---:|---|
| CR001 | VERIFIED | 0 | CMS 2026 `g_KK -> ttbar` 5.5 TeV exclusion, CMS 2019, PDG Live, and ATLAS historical limits matched. |
| CR002 | VERIFIED | 0 | PDG/ATLAS/CMS `T_{5/3}` pair-production mass limits matched; projection and scan-context items skipped. |
| CR003 | VERIFIED | 0 | PDG Live and ATLAS/CMS charge-2/3 top-partner mass-limit rows matched. |
| CR004 | VERIFIED | 0 | PDG 2025 bottom-partner `bH`, `bZ`, and `tW` mass limits matched. |
| CR005 | VERIFIED | 0 | CMS/ATLAS/PDG dilepton mass and fiducial cross-section limits matched. |
| CR006 | VERIFIED | 0 | PDG/CMS charged-resonance `lnu` and `tb` mass limits matched. |
| CR007 | VERIFIED | 0 | PDG RS/bulk-graviton and ATLAS bulk-graviton limits matched. |
| CR008 | VERIFIED | 0 | ATLAS singlet-T and CMS all-branching-envelope mass limits matched. |
| CR009 | VERIFIED | 0 | PDG ATLAS/CMS compositeness-scale endpoints matched; EFT translation marked N/A. |
| CR010 | VERIFIED | 0 | PDG VLQ doublet and simplified `T/B` endpoint limits matched. |
| CR011 | VERIFIED | 0 | ATLAS/CMS longitudinal same-sign `WW` fiducial cross-section limits matched. |
| CR012 | VERIFIED | 0 | PDG/CMS HVT diboson mass limits and CMS Table 2 intervals matched. |
| CR013 | VERIFIED | 0 | PDG/CMS/ATLAS diphoton RS-graviton mass limits matched. |
| CR014 | VERIFIED | 0 | PDG four-top cross sections, CMS top-philic vector limit, and ATLAS cross-section-limit range matched. |

## Fetch exceptions

| URL | Status | Impact |
|---|---:|---|
| None | All cited URLs returned HTTP 200 | No fetch exception affected the audit. PDF value checks used the local text snapshots after live URL/metadata verification. |

### CR001 - VERIFIED

Key references checked: arXiv:2603.23454, HEPData ins3134005, PDG Live S071KKG, arXiv:1810.05905, arXiv:1305.2756, arXiv:1207.2409, arXiv:2512.17856, arXiv:2005.05138, arXiv:1804.10823, arXiv:hep-ph/0701166.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| `CMS2026:CR001:gkk_ttbar_mass_exclusion`: `0.5 < m(g_KK) < 5.5 TeV` excluded at 95% CL; TeX also uses `138 fb^-1`, 13 TeV, 0/1/2-lepton channels. | https://arxiv.org/abs/2603.23454 | Y | Live arXiv metadata matches CMS 2026. Local snapshot contains the dataset/topology text and the 0.5-5.5 TeV KK-gluon exclusion. |
| HEPData record for CMS-B2G-25-009, DOI `10.17182/hepdata.168403.v1`. | https://www.hepdata.net/record/ins3134005 | Y | Live HEPData record metadata matches the CMS result; local JSON snapshot is present and matches the sidecar. |
| `CMS2019:CR001:gkk_ttbar_mass_exclusion`: KK gluon excluded up to `4.55 TeV` in the 2016 CMS ttbar resonance combination. | https://arxiv.org/abs/1810.05905 | Y | Live arXiv metadata matches CMS 2019; local snapshot states KK gluons are excluded up to 4.55 TeV. |
| `PDGLive2026:CR001:s071kkg_aaboud2018bi`: `m(g_KK) > 3.8 TeV`, `Gamma/m = 15.3%`. | https://pdglive.lbl.gov/DataBlock.action?node=S071KKG | Y | Live PDG Live page and local snapshot contain the 3.8 TeV row and width benchmark. |
| `ATLAS2013:CR001:gkk_ttbar_mass_exclusion`: `m(g_KK) < 2.07 TeV` excluded. | https://arxiv.org/abs/1305.2756 | Y | Live arXiv metadata matches ATLAS; local snapshot contains the 2.07 TeV exclusion. |
| `ATLAS2012:CR001:gkk_ttbar_mass_exclusion`: `m(g_KK) < 1.5 TeV` excluded with `2.05 fb^-1` at 7 TeV. | https://arxiv.org/abs/1207.2409 | Y | Live arXiv metadata matches ATLAS; local snapshot contains the 1.5 TeV exclusion and dataset. |
| ATLAS 2025/2020/2018 supporting dataset and search-range metadata: `140 fb^-1`, `139 fb^-1`, `36.1 fb^-1`, and `0.4-5.0 TeV` search range. | https://arxiv.org/abs/2512.17856 | N/A | Dataset/search metadata are policy carve-outs; live metadata and local snapshots were still checked and support the TeX statements. |
| Local RUNA comparison scales `47.26 TeV` and `127.13 TeV`, and Lillie-Randall-Wang `5 TeV` discovery-reach context. | docs/quark_scan_methodology_note.tex and https://arxiv.org/abs/hep-ph/0701166 | N/A | Local scan context and future-reach/theory baseline, not measured collider facts. |

### CR002 - VERIFIED

Key references checked: PDG encoder `q009.pdf`, arXiv:2212.05263, arXiv:1810.03188, arXiv:1807.11883, arXiv:1705.10967, arXiv:1312.2391, arXiv:0801.1679, arXiv:0909.3977.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| `PDG2026:CR002:ATLAS2023:X53_pair_Wt`: `m_{t'(5/3)} > 1460 GeV` at 95% CL. | https://pdg.lbl.gov/encoder_listings/q009.pdf | Y | Live PDG PDF returned 200; local PDG encoder snapshot contains the `>1460` ATLAS 23AG row. |
| `ATLAS2023:CR002:mass_degenerate_doublet`: stronger `1.59 TeV` mass-degenerate doublet limit, with `139 fb^-1` at 13 TeV. | https://arxiv.org/abs/2212.05263 | Y | Live arXiv metadata matches ATLAS; local snapshot contains `1.59 TeV` and the Run-2 dataset. |
| `CMS2019:CR002:X53_RH` and `CMS2019:CR002:X53_LH`: `1.33 TeV` and `1.30 TeV`. | https://arxiv.org/abs/1810.03188 | Y | Live arXiv metadata matches CMS; local snapshot contains both right- and left-handed limits. |
| `ATLAS2018:CR002:T53_pair_only`: `1.19 TeV`; pair-plus-single interpretation reaches `1.6 TeV` for unit coupling. | https://arxiv.org/abs/1807.11883 | Y | Live arXiv metadata matches ATLAS; local snapshot and PDG encoder row contain the 1.19 TeV and 1.6 TeV values. |
| `CMS2017:CR002:X53_RH/LH`: `1.02 TeV` and `0.99 TeV`. | https://arxiv.org/abs/1705.10967 | Y | Live arXiv metadata matches CMS; local snapshot supports both values. |
| `CMS2014:CR002:X53_same_sign`: `800 GeV`. | https://arxiv.org/abs/1312.2391 | Y | Live arXiv metadata matches CMS; local snapshot contains the 800 GeV exclusion. |
| Dataset metadata for ATLAS/CMS searches: `139`, `35.9`, `36.1`, `2.3`, and `19.5 fb^-1`. | multiple arXiv URLs above | N/A | Dataset metadata are policy carve-outs; snapshots support the values. |
| Quark-scan `47.26 TeV`, Contino-Servant theory context, and Mrazek-Wulzer `1.5 TeV` reach context. | https://arxiv.org/abs/0801.1679 and https://arxiv.org/abs/0909.3977 | N/A | Theory/projection and local-scan context only. |

### CR003 - VERIFIED

Key references checked: PDG encoder `q009.pdf`, arXiv:2401.17165, arXiv:2210.15413, arXiv:2209.07327, arXiv:1707.03347, arXiv:1803.09678, arXiv:1808.01771, arXiv:0801.1679.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| PDG Live/encoder charge-2/3 `T` entries: `m_T > 1.70 TeV`, `1.36 TeV`, `1.60 TeV`, and `1.50 TeV` for the stated Wb/singlet/Zt/Ht benchmarks. | https://pdg.lbl.gov/encoder_listings/q009.pdf | Y | Live PDF returned 200; local PDG snapshot contains all four mass-limit rows. |
| `CMS2023:CR003:T_all_third_generation_decay_mixtures`: `m_T > 1.48 TeV`. | https://arxiv.org/abs/2209.07327 | Y | Live arXiv metadata matches CMS; local snapshot contains the all-decay-mixture 1.48 TeV statement. |
| ATLAS 2024, ATLAS 2023, and CMS 2023 supporting datasets: `140`, `139`, and `138 fb^-1` at 13 TeV. | https://arxiv.org/abs/2401.17165 | N/A | Dataset metadata are policy carve-outs; values are visible in local snapshots. |
| Post-2010 search-history references and their arXiv IDs. | https://arxiv.org/abs/1707.03347 | Y | Live arXiv metadata for ATLAS/CMS historical references match manifest entries; no additional value mismatch found. |
| Local low-energy scan values `47.26 TeV` and `127.13 TeV`; Contino-Servant top-partner baseline. | docs/quark_scan_methodology_note.tex and https://arxiv.org/abs/0801.1679 | N/A | Local-scan and theory-baseline context only. |

### CR004 - VERIFIED

Key references checked: PDG 2025 `b-prime` listing, arXiv:2402.13808, HEPData ins2760468, arXiv:1808.02343, arXiv:2008.09835, arXiv:2209.07327, arXiv:2210.15413, arXiv:2212.05263, arXiv:2405.17605, arXiv:0801.1679.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| `PDG2025:CR004:B_bH_pair_mass_limit`: `m_B > 1570 GeV`. | https://pdg.lbl.gov/2025/listings/rpp2025-list-b-prime-quark.pdf | Y | Live PDG PDF returned 200; local PDG snapshot contains the `bH` 1570 GeV row. |
| `PDG2025:CR004:B_bZ_pair_mass_limit`: `m_B > 1540 GeV`. | https://pdg.lbl.gov/2025/listings/rpp2025-list-b-prime-quark.pdf | Y | Local PDG snapshot contains the `bZ` 1540 GeV row. |
| `PDG2025:CR004:B_tW_pair_mass_limit`: `m_B > 1560 GeV`. | https://pdg.lbl.gov/2025/listings/rpp2025-list-b-prime-quark.pdf | Y | Local PDG snapshot contains the `tW` 1560 GeV row. |
| CMS-B2G-20-014 dataset context and HEPData contour record. | https://arxiv.org/abs/2402.13808 and https://www.hepdata.net/record/ins2760468 | N/A | Dataset and contour metadata are policy carve-outs; live URL metadata and snapshots matched. |
| Post-2010 ATLAS/CMS VLQ search-history references. | multiple arXiv URLs above | Y | Live arXiv metadata matches manifest records; no unsupported numerical claim found. |
| Local RS-flavor comparison `47.26 TeV` and `127.13 TeV`; CMS review synthesis. | repo-local:docs/quark_scan_methodology_note.tex and https://arxiv.org/abs/2405.17605 | N/A | Local-scan and review-context material, not a measured CR004 observable. |

### CR005 - VERIFIED

Key references checked: arXiv:2103.02708, HEPData ins1849964, arXiv:1903.06248, HEPData ins1725190, PDG 2024 Z-prime review, arXiv:1209.2535, arXiv:1707.02424, arXiv:1803.06292, arXiv:hep-ph/9911262, arXiv:hep-ph/0308036, arXiv:hep-ph/0411254, arXiv:1004.2432.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| `CMS2021:CR005:ZSSM_ll_mass_lower_limit`: `m(Z'_SSM) > 5.15 TeV`. | https://arxiv.org/abs/2103.02708 | Y | Live arXiv metadata matches CMS; local snapshot contains the 5.15 TeV CMS Table 5 value. |
| `CMS2021:CR005:Zpsi_ll_mass_lower_limit`: `m(Z'_psi) > 4.56 TeV`. | https://arxiv.org/abs/2103.02708 | Y | Same source/snapshot contains the 4.56 TeV value. |
| CMS HEPData structured dilepton record. | https://www.hepdata.net/record/ins1849964?version=2 | Y | Live HEPData record metadata matches CMS 2021; used as structured cross-check. |
| `ATLAS2019:CR005:ZSSM_ll_mass_lower_limit`: `m(Z'_SSM) > 5.1 TeV`. | https://arxiv.org/abs/1903.06248 | Y | Live arXiv metadata matches ATLAS; local snapshot contains the rounded 5.1 TeV value. |
| `ATLAS2019:CR005:fiducial_xsec_limit_6tev_zero_width`: `sigma_fid x B < about 0.014 fb` at `m_X = 6 TeV`. | https://arxiv.org/abs/1903.06248 | Y | Local ATLAS snapshot and HEPData cross-check contain the 6 TeV zero-width fiducial limit context. |
| `PDG2024:CR005:dilepton_zprime_xsec_summary_limit`: cross-section limits as low as `0.02 fb`. | https://pdg.lbl.gov/2024/reviews/rpp2024-rev-zprime-searches.pdf | Y | Live PDG PDF returned 200; local PDG text snapshot contains the 0.02 fb summary. |
| Dataset metadata for CMS/ATLAS full and early Run-2/Run-1 inputs. | multiple arXiv URLs above | N/A | Dataset metadata are policy carve-outs; local snapshots support the values. |
| Bulk-gauge, custodial-RS, warped-GUT, and KK electroweak theory references plus local `47.26` and `127.13 TeV` scan context. | multiple arXiv URLs above | N/A | Theory and local-scan context only. |

### CR006 - VERIFIED

Key references checked: PDG 2025 W-prime review, arXiv:1906.05609, arXiv:2202.06075, arXiv:2310.19893, arXiv:2104.04831, arXiv:1807.10473, arXiv:1801.07893, arXiv:0810.1497, arXiv:hep-ph/0403143.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| `PDG2025:CR006:Wprime_SSM_enu_mass_lower_bound`: `M(W'_SSM) > 6.0 TeV`. | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-wprime-searches.pdf | Y | Live PDG PDF returned 200; local snapshot contains the 6.0 TeV electron-channel limit. |
| `PDG2025:CR006:Wprime_SSM_munu_mass_lower_bound`: `M(W'_SSM) > 5.6 TeV`. | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-wprime-searches.pdf | Y | Local PDG snapshot contains the 5.6 TeV muon-channel value. |
| `CMS2022:CR006:Wprime_SSM_enu_munu_combined_mass_lower_bound`: combined `e + mu` limit below `5.7 TeV`. | https://arxiv.org/abs/2202.06075 | Y | Live arXiv metadata matches CMS; local snapshot contains the 5.7 TeV combined statement. |
| `CMS2024:CR006:Wprime_R_tb_mass_lower_bound`: `M(W'_R) > 4.3 TeV`. | https://arxiv.org/abs/2310.19893 | Y | Live arXiv metadata matches CMS; local snapshot contains the 4.3 TeV right-handed `tb` bound. |
| `CMS2024:CR006:Wprime_L_tb_mass_lower_bound`: `M(W'_L) > 3.9 TeV`. | https://arxiv.org/abs/2310.19893 | Y | Same source/snapshot contains the 3.9 TeV left-handed `tb` bound. |
| ATLAS/CMS `lnu` and `tb` dataset metadata: `139` and `138 fb^-1` at 13 TeV. | multiple arXiv URLs above | N/A | Dataset metadata are policy carve-outs; snapshots support the values. |
| Agashe warped charged-gauge `2-3 TeV` reach context and local `47.26/127.13 TeV` scan context. | https://arxiv.org/abs/0810.1497 | N/A | Future-reach/theory and local-scan context only. |

### CR007 - VERIFIED

Key references checked: PDG 2025 extra-dimensions listing, arXiv:2405.09320, arXiv:2210.00043, arXiv:2103.02708, arXiv:2004.14636, arXiv:1906.08589, arXiv:1808.02380, arXiv:hep-ph/9905221, arXiv:hep-ph/0701186, arXiv:0709.0007.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| `PDG2025:CR007:RS_graviton_diphoton_mass_kMpl_0p1`: `m(G_KK^(1)) > 4.8 TeV` for `k/Mbar_Pl = 0.1`. | https://pdg.lbl.gov/2025/listings/rpp2025-list-extra-dimensions.pdf | Y | Live PDG PDF returned 200; local PDG snapshot contains the 4.8 TeV diphoton row. |
| `PDG2025:CR007:RS_graviton_dilepton_mass_kMpl_0p1`: `m(G_KK^(1)) > 4.78 TeV`. | https://pdg.lbl.gov/2025/listings/rpp2025-list-extra-dimensions.pdf | Y | Local PDG snapshot contains the 4.78 TeV dilepton row. |
| `PDG2025:CR007:bulk_graviton_diboson_mass_kMpl_0p5`: `m(G_bulk) > 1.4 TeV` for `k/Mbar_Pl = 0.5`. | https://pdg.lbl.gov/2025/listings/rpp2025-list-extra-dimensions.pdf | Y | Local PDG snapshot contains the bulk-graviton diboson 1.4 TeV row. |
| `ATLAS2018:CR007:bulk_graviton_combined_mass_kMpl_1`: `m(G_KK) > 2.3 TeV`. | https://arxiv.org/abs/1808.02380 | Y | Live arXiv metadata matches ATLAS; local snapshot contains the 2.3 TeV combined limit. |
| CMS/ATLAS supporting datasets: `138 fb^-1`, `36.1 fb^-1`, 13 TeV. | multiple arXiv URLs above | N/A | Dataset metadata are policy carve-outs. |
| Original RS and warped-graviton/gauge theory references; local `47.26 TeV` scan comparison. | multiple arXiv URLs above | N/A | Theory and local-scan context only. |

### CR008 - VERIFIED

Key references checked: PDG 2025 T-prime listing, arXiv:2401.17165, arXiv:2209.07327, arXiv:1505.04306, arXiv:1706.03408, arXiv:1805.04758, arXiv:1710.01539, arXiv:0907.3155, arXiv:0801.1800, arXiv:1306.0572.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| `ATLAS2024:CR008:T_singlet_pair_mass_limit`: singlet `T` limit `m_T > 1.36 TeV`. | https://arxiv.org/abs/2401.17165 | Y | Live arXiv metadata matches ATLAS; local snapshot contains the 1.36 TeV singlet benchmark. |
| PDG cross-check for the ATLAS singlet benchmark. | https://pdg.lbl.gov/2025/listings/rpp2025-list-t-prime-quark.pdf | Y | Live PDG PDF returned 200; local PDG snapshot contains the same singlet limit. |
| `CMS2023:CR008:T_pair_all_third_generation_decays_envelope`: `m_T > 1.48 TeV` for all third-generation branching patterns. | https://arxiv.org/abs/2209.07327 | Y | Live arXiv metadata matches CMS; local snapshot contains the 1.48 TeV all-branching envelope. |
| ATLAS and CMS Run-2 dataset metadata: `140` and `138 fb^-1`. | https://arxiv.org/abs/2401.17165 | N/A | Dataset metadata are policy carve-outs; snapshots support the values. |
| VLQ theory/review references and local `47.26 TeV` scan comparison. | multiple arXiv URLs above | N/A | Theory and local-scan context only. |

### CR009 - VERIFIED

Key references checked: PDG 2025 quark-lepton compositeness review, arXiv:2006.12946, arXiv:2103.02708, HEPData ins1849964, arXiv:1112.4462, arXiv:1407.2410, arXiv:1412.6302, arXiv:1607.03669, arXiv:1707.02424, arXiv:1812.10443, arXiv:1706.03783, arXiv:1811.12260, arXiv:1609.08157.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| PDG/ATLAS contact-interaction endpoints: `LL+ 35.8`, `LL- 26.0`, `RR+ 35.5`, `RR- 26.5`, `LR+ 32.5`, `LR- 28.8 TeV`. | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-searches-quark-lep-compositeness.pdf | Y | Live PDG PDF returned 200; local PDG snapshot contains all six ATLAS endpoint values. |
| PDG/CMS endpoint checks: `LL destructive 23.9 TeV` and `RR constructive 36.4 TeV`. | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-searches-quark-lep-compositeness.pdf | Y | Local PDG snapshot contains the CMS endpoint range values. |
| CMS HEPData contact-interaction limit tables. | https://www.hepdata.net/record/ins1849964 | Y | Live HEPData record metadata matches CMS 2021; local extract supports the structured-table reference. |
| ATLAS/CMS full-Run-2 dataset metadata: `139` and `140 fb^-1` at 13 TeV. | https://arxiv.org/abs/2006.12946 and https://arxiv.org/abs/2103.02708 | N/A | Dataset metadata are policy carve-outs. |
| RS/EFT translation `4*pi/Lambda^2 ~ abs(g_q g_l)/M_V^2`, `O(10 TeV)`, and SMEFT theory references. | https://arxiv.org/abs/1706.03783 | N/A | EFT translation and theory context, not a directly measured catalog value. |
| Local quark-scan context `47.26 TeV` and `127.13 TeV`. | docs/quark_scan_methodology_note.tex | N/A | Local methodology comparison only. |

### CR010 - VERIFIED

Key references checked: PDG 2025 T-prime and B-prime listings, arXiv:1808.02343, arXiv:2209.07327, arXiv:2401.17165, arXiv:1706.03408, arXiv:0907.3155, arXiv:1301.4454, arXiv:2405.17605.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| `PDG2025:CR010:TB_doublet_T_mass_lower_limit`: `m_T > 1.37 TeV`. | https://pdg.lbl.gov/2025/listings/rpp2025-list-t-prime-quark.pdf | Y | Live PDG PDF returned 200; local combined PDG snapshot contains the doublet T row. |
| `PDG2025:CR010:TB_doublet_B_mass_lower_limit`: `m_B > 1.37 TeV`. | https://pdg.lbl.gov/2025/listings/rpp2025-list-b-prime-quark.pdf | Y | Live PDG PDF returned 200; local snapshot contains the doublet B row. |
| Simplified endpoints: `m_T > 1.70 TeV`; `m_B > 1.57`, `1.56`, and `1.54 TeV` for Wb/Hb/Wt/Zb benchmark corners. | PDG T-prime/B-prime URLs above | Y | Local PDG snapshot contains all four simplified endpoint values. |
| `CMS2023:CR010:T_all_third_generation_decays_uniform_lower_limit`: `m_T > 1.48 TeV`. | https://arxiv.org/abs/2209.07327 | Y | Live arXiv metadata matches CMS; local snapshot contains the 1.48 TeV envelope. |
| ATLAS 2018 combination dataset and collision energy: `36.1 fb^-1`, 13 TeV. | https://arxiv.org/abs/1808.02343 | N/A | Dataset metadata are policy carve-outs. |
| Branching-generalization and top-partner theory references; CMS review/projection material. | https://arxiv.org/abs/1301.4454 | N/A | Theory/review/projection context only. |

### CR011 - VERIFIED

Key references checked: arXiv:2503.11317, arXiv:2009.09429, arXiv:2312.00420, arXiv:2004.10612, arXiv:1611.02428, arXiv:2603.18630, arXiv:1109.1570, PDG 2025 WZ quartic-couplings review.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| `ATLAS2025:CR011:fiducial_sigma_WLWL_same_sign_upper_limit`: `sigma_fid < 0.45 fb` at 95% CL; expected limit `0.70 fb`; evidence with at least one longitudinal W at `3.3 sigma`. | https://arxiv.org/abs/2503.11317 | Y | Live arXiv metadata matches ATLAS; local snapshot contains the observed/expected limits and evidence statement. |
| `CMS2020:CR011:fiducial_sigma_WLWL_same_sign_upper_limit`: `sigma < 1.17 fb` at 95% CL. | https://arxiv.org/abs/2009.09429 | Y | Live arXiv metadata matches CMS; local snapshot contains the 1.17 fb limit. |
| ATLAS/CMS dataset metadata: `140 fb^-1`, 13 TeV and related VBS/aQGC datasets. | multiple arXiv URLs above | N/A | Dataset metadata are policy carve-outs. |
| PDG QGC/aQGC review and ATLAS 2026 coefficient-combination context. | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-wz-quartic-couplings.pdf | N/A | Framework/EFT context only; no CR011 coefficient value promoted. |
| Local `47.26/127.13 TeV` scan comparison and Contino composite-Higgs resonance baseline. | https://arxiv.org/abs/1109.1570 | N/A | Theory/local-scan context only. |

### CR012 - VERIFIED

Key references checked: PDG 2025 heavy-boson listing, CMS-B2G-20-009 public page and Table 2 PDF, HEPData ins2159368, arXiv:2210.00043, arXiv:2109.06055, arXiv:1906.05977, arXiv:1708.05379, arXiv:1705.09171, arXiv:1708.04445, arXiv:1506.00962, arXiv:1402.4431.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| `PDG2025:CR012:HVTB_Wprime_WZ_mass_lower`: `m(W') > 4.4 TeV` for HVT model B. | https://pdg.lbl.gov/2025/listings/rpp2025-list-heavy-bosons.pdf | Y | Live PDG PDF returned 200; local PDG snapshot contains the 4.4 TeV Wprime WZ row. |
| `CMS2023:CR012:HVTB_Vprime_VV_mass_lower`: mass-degenerate `V' -> VV` lower limit `4.5 TeV`. | https://cms-results.web.cern.ch/cms-results/public-results/publications/B2G-20-009/CMS-B2G-20-009_Table_002.pdf | Y | Live CMS Table 2 PDF returned 200; local CMS snapshot contains the 4.5 TeV value. |
| `CMS2023:CR012:HVTB_Vprime_VV_VH_mass_lower`: `4.8 TeV` when `VV+VH` categories are combined. | https://cms-results.web.cern.ch/cms-results/public-results/publications/B2G-20-009/CMS-B2G-20-009_Table_002.pdf | Y | Same Table 2/local snapshot contains the combined 4.8 TeV value. |
| `CMS2023:CR012:HVTB_Zprime_WW_excluded_intervals`: observed `Z' -> WW` exclusions in `1.3-3.1 TeV` and `3.3-3.5 TeV`. | https://cms-results.web.cern.ch/cms-results/public-results/publications/B2G-20-009/CMS-B2G-20-009_Table_002.pdf | Y | Same Table 2/local snapshot contains both disjoint intervals. |
| CMS-B2G-20-009 dataset/search range: `138 fb^-1`, 13 TeV, `1.3-6 TeV`; HEPData upper-limit tables. | https://arxiv.org/abs/2210.00043 and https://www.hepdata.net/record/ins2159368 | N/A | Dataset and structured-table metadata are policy carve-outs. |
| HVT framework `g_V=3`, local `47.26/127.13 TeV` scan context, and generic theory assumptions. | https://arxiv.org/abs/1402.4431 | N/A | Theory/model-normalization and local-scan context only. |

### CR013 - VERIFIED

Key references checked: PDG 2025 extra-dimensions listing/review, CMS EXO-22-024 public page, HEPData ins2787227, arXiv:2405.09320, arXiv:2102.13405, arXiv:1809.00327, arXiv:1504.05511, arXiv:1112.0688, arXiv:1210.8389, arXiv:hep-ph/9905221, arXiv:hep-ph/9909255.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| `PDG2025:CR013:CMS2024_RSG_diphoton_kMPl_0p1`: `m_G* > 4.8 TeV` for `k/Mbar_Pl = 0.1`. | https://pdg.lbl.gov/2025/listings/rpp2025-list-extra-dimensions.pdf | Y | Live PDG PDF returned 200; local PDG snapshot contains the 4.8 TeV CMS diphoton row. |
| `CMS2024:CR013:RSG_diphoton_exclusion_range_ktilde_0p01_0p2`: KK gravitons below `2.2-5.6 TeV` excluded for `0.01 < ktilde < 0.2`. | https://cms-results.web.cern.ch/cms-results/public-results/publications/EXO-22-024/ | Y | Live CMS public page returned 200; local public-page and arXiv snapshots contain the range. |
| `ATLAS2021:CR013` RS-graviton limits: `4.5 TeV`, `3.9 TeV`, and `2.2 TeV` for `k/Mbar_Pl = 0.1`, `0.05`, and `0.01`. | https://arxiv.org/abs/2102.13405 | Y | Live arXiv metadata matches ATLAS; local snapshot contains all three benchmark limits. |
| `CMS2018:CR013` exclusion range `2.3-4.6 TeV` for `ktilde = 0.01-0.2`. | https://arxiv.org/abs/1809.00327 | Y | Live arXiv metadata matches CMS; local snapshot contains the range. |
| `ATLAS2015:CR013:RSG_diphoton_kMPl_0p1`: `m_G* > 2.66 TeV`. | https://arxiv.org/abs/1504.05511 | Y | Live arXiv metadata matches ATLAS; local snapshot contains the 2.66 TeV value. |
| CMS/ATLAS dataset metadata: `138` and `139 fb^-1`, 13 TeV; HEPData limit tables. | CMS/ATLAS URLs above | N/A | Dataset and structured-table metadata are policy carve-outs. |
| RS model-dependence review and original RS/DHR theory references; local `47.26 TeV` comparison. | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-extra-dimensions.pdf | N/A | Theory/model-dependence and local-scan context only. |

### CR014 - VERIFIED

Key references checked: PDG 2025 top-quark review, arXiv:2303.15061, arXiv:2305.13439, CMS-B2G-25-005 public page, arXiv:2604.14058, arXiv:2304.01678, arXiv:1409.7339, arXiv:1702.06164, arXiv:1710.10614, arXiv:1811.02305, arXiv:1908.06463, arXiv:2106.11683, arXiv:2303.03864, arXiv:1602.01934, arXiv:1803.00939, arXiv:2010.05915, arXiv:2104.09512, arXiv:2507.05334.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| `PDG2025:CR014:ATLAS_4top_cross_section`: `22.5^{+6.6}_{-5.5} fb`. | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-top-quark.pdf | Y | Live PDG PDF returned 200; local PDG and ATLAS snapshots contain the ATLAS four-top cross section. |
| `PDG2025:CR014:CMS_4top_cross_section`: `17.7^{+4.4}_{-4.0} fb`. | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-top-quark.pdf | Y | Local PDG snapshot contains the CMS summary value; CMS arXiv snapshot contains the same central value with split uncertainties. |
| `CMSB2G25005:CR014:top_philic_vector_50pct_width_excluded_up_to`: vector mediator excluded up to `850 GeV` observed, `1000 GeV` expected for `50%` width. | https://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/B2G-25-005/index.html | Y | Live CMS public page returned 200; local snapshot contains the 850/1000 GeV statement. |
| CMS-B2G-25-005/arXiv:2604.14058 search range and datasets: `500 GeV-4 TeV`, `4-50%` widths, `138 fb^-1` at 13 TeV and `35 fb^-1` at 13.6 TeV. | https://arxiv.org/abs/2604.14058 | Y | Live arXiv metadata matches CMS note; local arXiv and public-result snapshots contain these values. Dataset components are metadata but support the result statement. |
| `ATLAS2024:CR014:top_philic_Zprime_cross_section_limit_range`: observed 95% CL upper-limit range `21-119 fb` and expected `14-86 fb`. | https://arxiv.org/abs/2304.01678 | Y | Live arXiv metadata matches ATLAS; local snapshot contains the limit ranges. |
| Post-2010 four-top search-history arXiv references and observation/evidence sequence. | multiple arXiv URLs above | Y | Live arXiv metadata matches manifest entries; no value-bearing mismatch found. |
| Local `47.26 TeV` quark-scan comparison and four-top EFT/top-philic theory references. | docs/quark_scan_methodology_note.tex and multiple arXiv URLs above | N/A | Local-scan, EFT, and theory-context material only. |
