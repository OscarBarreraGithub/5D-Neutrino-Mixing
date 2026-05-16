# Fact-check report: top_higgs_ew

Date: 2026-05-16T17:45:00-04:00
Agent: factcheck-codex-top_higgs_ew
Family: top_higgs_ew

## Summary

I re-fetched the cited PDG/PDGLive pages and PDFs, arXiv abstract/PDF pages, CMS/Gfitter/OSTI pages, and the local snapshots under `flavor_catalog/references/<process_id>/`. For arXiv sources, I checked title and author metadata against the manifests before checking the numerical claims in the abstract/PDF text or the local snapshot.

Policy carve-outs were applied to pure theory-normalization context such as CFW KK-gluon reference scales and dimensional EFT estimates. Those rows are marked `N/A` unless the process text claimed a concrete quoted theory prediction.

| Process | Verdict | Mismatches | Note |
|---|---:|---:|---|
| T001 | VERIFIED | 0 | Top FCNC `t -> Zc` limits, SM estimate, and HL-LHC projection matched. |
| T002 | VERIFIED | 0 | Top FCNC `t -> Zu` limits and SM estimate matched; RS scale context skipped by policy. |
| T005 | VERIFIED | 0 | `t -> cg` PDG/ATLAS/CMS limits and SM estimate matched. |
| T006 | VERIFIED | 0 | `t -> ug` PDG/ATLAS/CMS limits and SM estimate matched. |
| T007 | VERIFIED | 0 | `t -> Hc` PDG/ATLAS/CMS limits matched. |
| T010 | VERIFIED | 0 | `Zbb` PDG/LEP-SLC/FCC-ee values matched. |
| T015 | VERIFIED | 0 | `Z -> e mu` CMS/PDG/ATLAS limits matched. |
| T016 | VERIFIED | 0 | `Z -> e tau` PDG/ATLAS/CMS limits matched. |
| T017 | VERIFIED | 0 | `Z -> mu tau` PDG/ATLAS/CMS limits matched. |
| T018 | VERIFIED | 0 | `h -> mu tau` PDG/CMS/ATLAS/CMS-2015 values matched. |
| T019 | VERIFIED | 0 | `h -> e tau` PDG/CMS/ATLAS values matched. |
| T020 | MISMATCH | 1 | ATLAS 2020 `h -> e mu` limit is quoted as `< 6.1 x 10^-5`, but fetched ATLAS/PDG sources give `< 6.2 x 10^-5`. |
| EW001 | VERIFIED | 0 | PDG/Gfitter/HEPfit oblique parameter values matched. |
| EW002 | VERIFIED | 0 | PDG/FLAG/Hardy-Towner/CKM anomaly values matched; CFW scale context skipped by policy. |
| EW003 | VERIFIED | 0 | PDG/HFLAV/FLAG CKM values matched. |

## Fetch exceptions

| URL | Status | Impact |
|---|---:|---|
| https://doi.org/10.1103/PhysRevD.46.381 | HTTP 403 | Publisher DOI for the Peskin-Takeuchi theory reference. No measured catalog value depends on the blocked page; the process values were checked against PDG/Gfitter/HEPfit. |
| https://journals.aps.org/prd/abstract/10.1103/PhysRevD.110.075005 | HTTP 403 | Publisher page for an EW003 BSM-context reference. No measured catalog value depends on the blocked page; CKM values were checked against PDG/HFLAV/FLAG. |

## T001 - VERIFIED

Key references checked: arXiv:2301.11605, arXiv:1702.01404, arXiv:hep-ph/0409342, and arXiv:2601.14966 title/author metadata matched the manifest records.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| PDG2025:T001:tZc_left = `< 1.3e-4` | https://pdgweb.lbl.gov/2025/listings/rpp2025-list-t-quark.pdf | Y | PDG top listing gives the ATLAS `tZc` LH value as `< 0.13` in units of `10^-3`. |
| PDG2025:T001:tZc_right = `< 1.2e-4` | https://pdgweb.lbl.gov/2025/listings/rpp2025-list-t-quark.pdf | Y | PDG top listing gives the ATLAS `tZc` RH value as `< 0.12` in units of `10^-3`. |
| ATLAS2023:T001:dataset = `139 fb^-1 at sqrt(s) = 13 TeV` | https://arxiv.org/abs/2301.11605 | Y | ATLAS paper states full Run-2 139 fb^-1 at 13 TeV. |
| ATLAS2023:T001:tZc_left = `< 13e-5` | https://arxiv.org/abs/2301.11605 | Y | Paper gives `B(t -> Zc) < 13 x 10^-5` for left-handed coupling. |
| ATLAS2023:T001:tZc_right = `< 12e-5` | https://arxiv.org/abs/2301.11605 | Y | Paper gives `B(t -> Zc) < 12 x 10^-5` for right-handed coupling. |
| CMS2017:T001:dataset = `19.7 fb^-1 at sqrt(s) = 8 TeV` | https://arxiv.org/abs/1702.01404 | Y | CMS paper states 19.7 fb^-1 at 8 TeV. |
| CMS2017:T001:tZc = `< 0.049%` | https://arxiv.org/abs/1702.01404 | Y | CMS paper quotes `B(t -> Zc) < 0.049%`. |
| AguilarSaavedra2004:T001:SM = `1e-14` | https://arxiv.org/abs/hep-ph/0409342 | Y | Theory-prediction table gives SM `t -> cZ` at order `10^-14`. |
| PDG2025:T001:SM = order `1e-14` | https://pdgweb.lbl.gov/2025/listings/rpp2025-list-t-quark.pdf | Y | PDG top review context is consistent with the Aguilar-Saavedra SM estimate. |
| AirenFranceschini2026:T001:HL_LHC_projection = approximately `1e-6` | https://arxiv.org/abs/2601.14966 | Y | Projection paper exists and quotes HL-LHC sensitivity around `10^-6` for `t -> cZ`. |

## T002 - VERIFIED

Key references checked: arXiv:2301.11605, arXiv:1702.01404, arXiv:hep-ph/0409342, arXiv:0807.4937, and arXiv:0804.1954 title/author metadata matched the manifest records.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| PDG2025:T002:tZu_left = `< 6.2e-5` | https://pdgweb.lbl.gov/2025/listings/rpp2025-list-t-quark.pdf | Y | PDG top listing gives `< 0.062` in units of `10^-3`. |
| PDG2025:T002:tZu_right = `< 6.6e-5` | https://pdgweb.lbl.gov/2025/listings/rpp2025-list-t-quark.pdf | Y | PDG top listing gives `< 0.066` in units of `10^-3`. |
| ATLAS2023:T002:dataset = `139 fb^-1 at sqrt(s) = 13 TeV` | https://arxiv.org/abs/2301.11605 | Y | ATLAS paper states full Run-2 139 fb^-1 at 13 TeV. |
| ATLAS2023:T002:tZu_left = `< 6.2e-5` | https://arxiv.org/abs/2301.11605 | Y | Paper gives `B(t -> Zu) < 6.2 x 10^-5` for left-handed coupling. |
| ATLAS2023:T002:tZu_right = `< 6.6e-5` | https://arxiv.org/abs/2301.11605 | Y | Paper gives `B(t -> Zu) < 6.6 x 10^-5` for right-handed coupling. |
| CMS2017:T002:dataset = `19.7 fb^-1 at sqrt(s) = 8 TeV` | https://arxiv.org/abs/1702.01404 | Y | CMS paper states 19.7 fb^-1 at 8 TeV. |
| CMS2017:T002:tZu = `< 0.022%` | https://arxiv.org/abs/1702.01404 | Y | CMS paper quotes `B(t -> Zu) < 0.022%`. |
| AguilarSaavedra2004:T002:SM = `8e-17` | https://arxiv.org/abs/hep-ph/0409342 | Y | Theory-prediction table gives SM `t -> uZ` near `8 x 10^-17`. |
| Casagrande2008:T002:RS_suppression = approximately two orders below `t -> cZ` | https://arxiv.org/abs/0807.4937 | N/A | Theory scaling context, not a measured observable. |
| CFW2008:T002:KK_gluon_scale = about `21 TeV` | https://arxiv.org/abs/0804.1954 | N/A | Theory-normalization context skipped under the signoff policy. |

## T005 - VERIFIED

Key references checked: PDGLive Q007TUG, arXiv:2112.01302, arXiv:1610.03545, arXiv:hep-ph/0409342, and arXiv:0804.1954.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| PDG2025:T005:ATLAS2022:dataset = `139 fb^-1 at sqrt(s) = 13 TeV` | https://pdglive.lbl.gov/view/Q007TUG | Y | PDG/ATLAS entry identifies the 139 fb^-1 13 TeV dataset. |
| PDG2025:T005:ATLAS2022:cg_to_t_cross_section = `< 4.7 pb` | https://pdglive.lbl.gov/view/Q007TUG | Y | PDG top-gluon listing gives `sigma(cg -> t) < 4.7 pb`. |
| PDG2025:T005:ATLAS2022:w_leptonic_sum = `0.325` | https://arxiv.org/abs/2112.01302 | Y | ATLAS paper/snapshot contains the W leptonic-sum branching fraction used in the conversion. |
| PDG2025:T005:ATLAS2022:CuGct = `< 0.14 TeV^-2` | https://pdglive.lbl.gov/view/Q007TUG | Y | PDG lists `C_ct/Lambda^2 < 0.14 TeV^-2`. |
| PDG2025:T005:ATLAS2022:t_cg = `< 3.7e-4` | https://pdglive.lbl.gov/view/Q007TUG | Y | PDG lists `B(t -> cg) < 3.7 x 10^-4`. |
| CMS2017:T005:dataset = `5.0 fb^-1 at 7 TeV and 19.7 fb^-1 at 8 TeV` | https://arxiv.org/abs/1610.03545 | Y | CMS paper states both datasets. |
| CMS2017:T005:kappa_tcg = `< 1.8e-2 TeV^-1` | https://arxiv.org/abs/1610.03545 | Y | CMS paper quotes the `kappa_tcg` limit. |
| CMS2017:T005:t_cg = `< 4.1e-4` | https://arxiv.org/abs/1610.03545 | Y | CMS paper quotes `B(t -> cg) < 4.1 x 10^-4`. |
| AguilarSaavedra2004:T005:SM = `4.6e-12` | https://arxiv.org/abs/hep-ph/0409342 | Y | SM `t -> cg` estimate appears in the fetched theory paper. |
| CFW2008:T005:ordinary_RS_KK_gluon_context = about `21 TeV` | https://arxiv.org/abs/0804.1954 | N/A | Theory-normalization context skipped under the signoff policy. |
| CFW2008:T005:pNGB_Higgs_KK_gluon_context = about `33 TeV` | https://arxiv.org/abs/0804.1954 | N/A | Theory-normalization context skipped under the signoff policy. |

## T006 - VERIFIED

Key references checked: PDGLive Q007TUG, arXiv:2112.01302, arXiv:1610.03545, and arXiv:hep-ph/0409342.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| PDG2025:T006:ATLAS2022:dataset = `139 fb^-1 at sqrt(s) = 13 TeV` | https://pdglive.lbl.gov/view/Q007TUG | Y | PDG/ATLAS entry identifies the 139 fb^-1 13 TeV dataset. |
| PDG2025:T006:ATLAS2022:ug_to_t_cross_section = `< 3.0 pb` | https://pdglive.lbl.gov/view/Q007TUG | Y | PDG top-gluon listing gives `sigma(ug -> t) < 3.0 pb`. |
| PDG2025:T006:ATLAS2022:w_leptonic_sum = `0.325` | https://arxiv.org/abs/2112.01302 | Y | ATLAS paper/snapshot contains the W leptonic-sum branching fraction used in the conversion. |
| PDG2025:T006:ATLAS2022:CuGut = `< 0.057 TeV^-2` | https://pdglive.lbl.gov/view/Q007TUG | Y | PDG lists `C_ut/Lambda^2 < 0.057 TeV^-2`. |
| PDG2025:T006:ATLAS2022:t_ug = `< 6.1e-5` | https://pdglive.lbl.gov/view/Q007TUG | Y | PDG lists `B(t -> ug) < 6.1 x 10^-5`. |
| CMS2017:T006:dataset = `5.0 fb^-1 at 7 TeV and 19.7 fb^-1 at 8 TeV` | https://arxiv.org/abs/1610.03545 | Y | CMS paper states both datasets. |
| CMS2017:T006:kappa_tug = `< 4.1e-3 TeV^-1` | https://arxiv.org/abs/1610.03545 | Y | CMS paper quotes the `kappa_tug` limit. |
| CMS2017:T006:t_ug = `< 2.0e-5` | https://arxiv.org/abs/1610.03545 | Y | CMS paper quotes `B(t -> ug) < 2.0 x 10^-5`. |
| AguilarSaavedra2004:T006:t_ug_SM = approximately `3.6e-14` | https://arxiv.org/abs/hep-ph/0409342 | Y | Fetched theory paper gives the SM `t -> ug` branching ratio at `3.7 x 10^-14`, consistent with the catalog value. |

## T007 - VERIFIED

Key references checked: PDGLive Q007R01, arXiv:2404.02123, and arXiv:2407.15172 title/author metadata matched the manifest records.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| PDG2026:T007:tHc_combined = `< 3.4e-4` | https://pdgprod.lbl.gov/pdgprod/pdgLive/DataBlock.action?node=Q007R01&_eventName=showBR | Y | PDG datablock has table unit `10^-4` and `Hc < 3.4`. |
| ATLAS2024:T007:dataset = `140 fb^-1 at sqrt(s) = 13 TeV` | https://arxiv.org/abs/2404.02123 | Y | ATLAS paper states 140 fb^-1 at 13 TeV. |
| ATLAS2024:T007:tHc_multilepton = `< 3.3e-4` | https://arxiv.org/abs/2404.02123 | Y | ATLAS multilepton result quotes `< 3.3 x 10^-4`. |
| ATLAS2024:T007:tHc_combined = `< 3.4e-4` | https://arxiv.org/abs/2404.02123 | Y | ATLAS combined value quotes `< 3.4 x 10^-4`. |
| CMS2025:T007:dataset = `138 fb^-1 at sqrt(s) = 13 TeV` | https://arxiv.org/abs/2407.15172 | Y | CMS paper states 138 fb^-1 collected in 2016-2018. |
| CMS2025:T007:tHc_samesign = `< 0.043%` | https://arxiv.org/abs/2407.15172 | Y | CMS same-sign channel limit appears in the paper. |
| CMS2025:T007:tHc_combined = `< 0.037%` | https://arxiv.org/abs/2407.15172 | Y | CMS combined limit appears in the paper. |

## T010 - VERIFIED

Key references checked: PDG Z-boson listing, arXiv:hep-ex/0509008, and arXiv:2502.17281.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| `R_b^0 = 0.21629 +/- 0.00066` | https://pdg.lbl.gov/2025/listings/rpp2025-list-z-boson.pdf | Y | PDG Z listing gives `0.21629 +/- 0.00066`. |
| `A_FB^{0,b} = 0.0992 +/- 0.0016` | https://pdg.lbl.gov/2025/listings/rpp2025-list-z-boson.pdf | Y | PDG Z listing gives `9.92 +/- 0.16%`. |
| `A_b = 0.923 +/- 0.020` | https://pdg.lbl.gov/2025/listings/rpp2025-list-z-boson.pdf | Y | PDG Z listing gives `0.923 +/- 0.020`. |
| LEP/SLC final-combination pull for `A_FB^{0,b}` = `2.8` | https://arxiv.org/abs/hep-ex/0509008 | Y | LEP/SLC final combination gives a 2.8 sigma `A_FB^{0,b}` pull. |
| FCC-ee projected relative uncertainty for `R_b` and `A_FB^b` = `0.01%` | https://arxiv.org/abs/2502.17281 | Y | Projection source exists and quotes order `0.01%` total relative uncertainty. |

## T015 - VERIFIED

Key references checked: arXiv:2508.07512, PDG Z-boson listing, arXiv:2204.10783, arXiv:1408.5774, and arXiv:2107.10273.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| CMS2025:T015:zemu_limit = `< 1.9e-7` | https://arxiv.org/abs/2508.07512 | Y | CMS 138 fb^-1 paper quotes `B(Z -> e mu) < 1.9 x 10^-7`. |
| CMS2025:T015:dataset = `138 fb^-1 at sqrt(s) = 13 TeV` | https://arxiv.org/abs/2508.07512 | Y | CMS paper states 138 fb^-1 at 13 TeV. |
| PDG2025:T015:zemu_limit = `< 2.62e-7` | https://pdgweb.lbl.gov/2025/listings/rpp2025-list-z-boson.pdf | Y | PDG Z listing quotes `< 2.62 x 10^-7`. |
| ATLAS2023:T015:zemu_limit = `< 2.62e-7` | https://arxiv.org/abs/2204.10783 | Y | ATLAS Run-2 paper quotes `< 2.62 x 10^-7`. |
| ATLAS2023:T015:dataset = `139 fb^-1 at sqrt(s) = 13 TeV` | https://arxiv.org/abs/2204.10783 | Y | ATLAS paper states 139 fb^-1 at 13 TeV. |
| ATLAS2014:T015:zemu_limit = `< 7.5e-7` | https://arxiv.org/abs/1408.5774 | Y | ATLAS 8 TeV paper quotes `< 7.5 x 10^-7`. |
| ATLAS2014:T015:dataset = `20.3 fb^-1 at sqrt(s) = 8 TeV` | https://arxiv.org/abs/1408.5774 | Y | ATLAS paper states 20.3 fb^-1 at 8 TeV. |
| CalibbiMarcanoRoy2021:T015:teraZ = `O(1e12) Z decays` | https://arxiv.org/abs/2107.10273 | Y | Prospect paper discusses Tera-Z samples of order `10^12` to `3 x 10^12` Z decays. |

## T016 - VERIFIED

Key references checked: PDG Z-boson listing, arXiv:2105.12491, arXiv:2010.02566, arXiv:2508.07512, and arXiv:2107.10273.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| PDG2025:T016:zetau_limit = `< 5.0e-6` | https://pdgweb.lbl.gov/2025/listings/rpp2025-list-z-boson.pdf | Y | PDG Z listing quotes `< 5.0 x 10^-6`. |
| ATLAS2021:T016:zetau_combined_limit = `< 5.0e-6` | https://arxiv.org/abs/2105.12491 | Y | ATLAS combined result quotes `< 5.0 x 10^-6`. |
| ATLAS2021:T016:zetau_leptonic_limit = `< 7.0e-6` | https://arxiv.org/abs/2105.12491 | Y | ATLAS leptonic-only result quotes `< 7.0 x 10^-6`. |
| ATLAS2021:T016:dataset = `139 fb^-1 at sqrt(s) = 13 TeV` | https://arxiv.org/abs/2105.12491 | Y | ATLAS paper states 139 fb^-1 at 13 TeV. |
| ATLAS2020:T016:zetau_hadronic_limit = `< 8.1e-6` | https://arxiv.org/abs/2010.02566 | Y | ATLAS hadronic-tau paper quotes `< 8.1 x 10^-6`. |
| ATLAS2020:T016:dataset = `139 fb^-1 at 13 TeV and 20.3 fb^-1 at 8 TeV` | https://arxiv.org/abs/2010.02566 | Y | Paper states the Run-2 and Run-1 datasets. |
| CMS2025:T016:zetau_limit = `< 13.8e-6` | https://arxiv.org/abs/2508.07512 | Y | CMS paper quotes `B(Z -> e tau) < 13.8 x 10^-6`. |
| CMS2025:T016:dataset = `138 fb^-1 at sqrt(s) = 13 TeV` | https://arxiv.org/abs/2508.07512 | Y | CMS paper states 138 fb^-1 at 13 TeV. |
| CalibbiMarcanoRoy2021:T016:teraZ = `O(1e12) Z decays` | https://arxiv.org/abs/2107.10273 | Y | Prospect paper discusses Tera-Z samples of order `10^12` to `3 x 10^12` Z decays. |

## T017 - VERIFIED

Key references checked: PDG Z-boson listing, arXiv:2105.12491, and arXiv:2508.07512.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| PDG2025:T017:zmutau_limit = `< 6.5e-6` | https://pdgweb.lbl.gov/2025/listings/rpp2025-list-z-boson.pdf | Y | PDG Z listing quotes `< 6.5 x 10^-6`. |
| ATLAS2021:T017:zmutau_combined_limit = `< 6.5e-6` | https://arxiv.org/abs/2105.12491 | Y | ATLAS combined result quotes `< 6.5 x 10^-6`. |
| ATLAS2021:T017:zmutau_leptonic_tau_limit = `< 7.2e-6` | https://arxiv.org/abs/2105.12491 | Y | ATLAS leptonic-tau result quotes `< 7.2 x 10^-6`. |
| ATLAS2021:T017:dataset = `139 fb^-1 at sqrt(s) = 13 TeV` | https://arxiv.org/abs/2105.12491 | Y | ATLAS paper states 139 fb^-1 at 13 TeV. |
| CMS2025:T017:zmutau_limit = `< 12.0e-6` | https://arxiv.org/abs/2508.07512 | Y | CMS paper quotes `B(Z -> mu tau) < 12.0 x 10^-6`. |
| CMS2025:T017:dataset = `138 fb^-1 at sqrt(s) = 13 TeV` | https://arxiv.org/abs/2508.07512 | Y | CMS paper states 138 fb^-1 at 13 TeV. |

## T018 - VERIFIED

Key references checked: PDG Higgs review, arXiv:2105.03007, arXiv:2302.05225, arXiv:1502.07400, arXiv:1209.1397, and arXiv:0804.1954.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| PDG2025:T018:cms_run2 = `< 0.15%` | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-higgs-boson.pdf | Y | PDG Higgs review quotes CMS full Run-2 `h -> mu tau` limit `< 0.15%`. |
| PDG2025:T018:atlas_run2 = `< 0.18%` | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-higgs-boson.pdf | Y | PDG Higgs review quotes ATLAS full Run-2 limit `< 0.18%`. |
| CMS2021:T018:mutau_limit = `< 0.15%` | https://arxiv.org/abs/2105.03007 | Y | CMS paper quotes `< 0.15%`. |
| CMS2021:T018:dataset = `137 fb^-1 at sqrt(s) = 13 TeV` | https://arxiv.org/abs/2105.03007 | Y | CMS paper states 137 fb^-1 at 13 TeV. |
| ATLAS2023:T018:mutau_limit = `< 0.18%` | https://arxiv.org/abs/2302.05225 | Y | ATLAS paper quotes `< 0.18%`. |
| ATLAS2023:T018:dataset = `138 fb^-1 at sqrt(s) = 13 TeV` | https://arxiv.org/abs/2302.05225 | Y | ATLAS paper states 138 fb^-1 at 13 TeV. |
| ATLAS2023:T018:symmetry_difference = `0.25%` | https://arxiv.org/abs/2302.05225 | Y | ATLAS paper quotes the channel-asymmetry difference near `0.25 +/- 0.10%`. |
| CMS2015:T018:hint = `2.4 sigma` | https://arxiv.org/abs/1502.07400 | Y | CMS 8 TeV paper quotes the 2.4 sigma excess. |
| CMS2015:T018:limit = `< 1.51%` | https://arxiv.org/abs/1502.07400 | Y | CMS paper quotes the 95% CL limit `< 1.51%`. |
| CMS2015:T018:dataset = `19.7 fb^-1 at sqrt(s) = 8 TeV` | https://arxiv.org/abs/1502.07400 | Y | CMS paper states 19.7 fb^-1 at 8 TeV. |
| CMS2015:T018:yukawa_limit = `< 3.6e-3` | https://arxiv.org/abs/1502.07400 | Y | CMS paper quotes the LFV Yukawa constraint. |
| HarnikKoppZupan2012:order10 = order `10%` | https://arxiv.org/abs/1209.1397 | Y | Quoted theory prediction/context appears in the fetched paper. |
| CsakiFalkowskiWeiler2008:kk_gluon_rs = about `21 TeV` | https://arxiv.org/abs/0804.1954 | N/A | Theory-normalization context skipped under the signoff policy. |

## T019 - VERIFIED

Key references checked: PDG Higgs review, arXiv:2302.05225, arXiv:2105.03007, arXiv:1907.06131, and arXiv:1209.1397.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| PDG2025:T019:atlas_run2 = `< 0.20%` | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-higgs-boson.pdf | Y | PDG Higgs review quotes ATLAS full Run-2 `h -> e tau` limit `< 0.20%`. |
| PDG2025:T019:cms_run2 = `< 0.22%` | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-higgs-boson.pdf | Y | PDG Higgs review quotes CMS full Run-2 limit `< 0.22%`. |
| ATLAS2023:T019:etau_limit = `< 0.20%` | https://arxiv.org/abs/2302.05225 | Y | ATLAS paper quotes `< 0.20%`. |
| ATLAS2023:T019:dataset = `138 fb^-1 at sqrt(s) = 13 TeV` | https://arxiv.org/abs/2302.05225 | Y | ATLAS paper states 138 fb^-1 at 13 TeV. |
| CMS2021:T019:etau_limit = `< 0.22%` | https://arxiv.org/abs/2105.03007 | Y | CMS paper quotes `< 0.22%`. |
| CMS2021:T019:dataset = `137 fb^-1 at sqrt(s) = 13 TeV` | https://arxiv.org/abs/2105.03007 | Y | CMS paper states 137 fb^-1 at 13 TeV. |
| ATLAS2019:T019:etau_limit = `< 0.47%` | https://arxiv.org/abs/1907.06131 | Y | ATLAS partial Run-2 paper quotes `< 0.47%`. |
| ATLAS2019:T019:dataset = `36.1 fb^-1 at sqrt(s) = 13 TeV` | https://arxiv.org/abs/1907.06131 | Y | ATLAS paper states 36.1 fb^-1 at 13 TeV. |
| HarnikKoppZupan2012:T019:LFV_Higgs_EFT_allowance = order `10%` | https://arxiv.org/abs/1209.1397 | Y | Quoted theory prediction/context appears in the fetched paper. |

## T020 - MISMATCH

Key references checked: PDG Higgs review, CMS HIG-22-002 publication page, arXiv:1909.10235, arXiv:1209.1397, arXiv:0804.1954, and arXiv:0805.4652.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| PDG2025:T020:higgs_mass_hypothesis = `125 GeV` | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-higgs-boson.pdf | Y | PDG Higgs review states the full Run-2 searches at the 125 GeV Higgs hypothesis. |
| PDG2025:T020:cms_run2 = `< 4.4 x 10^-5` | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-higgs-boson.pdf | Y | PDG Higgs review quotes CMS `< 4.4 x 10^-5`, expected `< 4.7 x 10^-5`. |
| PDG2025:T020:atlas_run2 = `< 6.2 x 10^-5` | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-higgs-boson.pdf | Y | PDG Higgs review quotes ATLAS `< 6.2 x 10^-5`, expected `< 5.9 x 10^-5`. |
| CMS2023:T020:emu_limit = `< 4.4 x 10^-5` | https://cms-results.web.cern.ch/cms-results/public-results/publications/HIG-22-002/ | Y | CMS public page quotes the observed `H -> e mu` limit `< 4.4 x 10^-5`. |
| CMS2023:T020:dataset = `138 fb^-1 at sqrt(s) = 13 TeV` | https://cms-results.web.cern.ch/cms-results/public-results/publications/HIG-22-002/ | Y | CMS public page states 138 fb^-1 at 13 TeV. |
| CMS2023:T020:mass_scan = `110-160 GeV`; largest excess near `146 GeV` | https://cms-results.web.cern.ch/cms-results/public-results/publications/HIG-22-002/ | Y | CMS public page states the 110-160 GeV scan and the excess near 146 GeV. |
| ATLAS2020:T020:emu_limit = `< 6.1 x 10^-5` | https://arxiv.org/abs/1909.10235 | N | Fetched ATLAS paper gives observed `< 6.2 x 10^-5` and expected `< 5.9 x 10^-5`; PDG agrees with `< 6.2 x 10^-5`. |
| ATLAS2020:T020:dataset = `139 fb^-1 at sqrt(s) = 13 TeV` | https://arxiv.org/abs/1909.10235 | Y | ATLAS paper states full Run-2 139 fb^-1 at 13 TeV. |

The T020 mismatch is numerical and source-content based. The YAML entry `ATLAS2020:T020:emu_limit` claims `< 6.1 x 10^-5`, and the local source snapshot also contains the older/incorrect `6.1/5.8` wording. Re-fetching arXiv:1909.10235 and the PDG 2025 Higgs review gives the observed ATLAS 95% CL limit as `< 6.2 x 10^-5` with expected `< 5.9 x 10^-5`.

## EW001 - VERIFIED

Key references checked: PDG Standard Model review, Gfitter public oblique-parameter page, arXiv:2204.04204, and the Peskin-Takeuchi DOI metadata where accessible. The direct APS DOI returned HTTP 403, but it is a theory-reference citation and not the source of the catalog's measured numbers.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| PDG 2025 global oblique fit, `S` with `U=0` = `0.026 +/- 0.075` | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-standard-model.pdf | Y | PDG Table 10.8 gives `S = 0.026 +/- 0.075`. |
| PDG 2025 global oblique fit, `T` with `U=0` = `0.047 +/- 0.066` | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-standard-model.pdf | Y | PDG Table 10.8 gives `T = 0.047 +/- 0.066`. |
| PDG 2025 global oblique fit, fixed `U = 0` | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-standard-model.pdf | Y | PDG table labels this fit with `U = 0`. |
| PDG 2025 `rho(S,T) = 0.90` with `U=0` | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-standard-model.pdf | Y | PDG table gives correlation `0.90`. |
| PDG 2025 floating `S = 0.021 +/- 0.096` | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-standard-model.pdf | Y | PDG table gives floating-STU `S = 0.021 +/- 0.096`. |
| PDG 2025 floating `T = 0.04 +/- 0.12` | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-standard-model.pdf | Y | PDG table gives floating-STU `T = 0.04 +/- 0.12`. |
| PDG 2025 floating `U = 0.008 +/- 0.092` | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-standard-model.pdf | Y | PDG table gives floating-STU `U = 0.008 +/- 0.092`. |
| PDG 2025 warped-extra-dimension `S` coefficient = `30` | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-standard-model.pdf | Y | PDG review contains the quoted warped-context coefficient. |
| PDG 2025 one-sided `S` upper value = `0.18` | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-standard-model.pdf | Y | PDG warped-context passage gives the one-sided `S` value used for the estimate. |
| PDG 2025 warped `M_KK` lower-bound context = `3.2 TeV` | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-standard-model.pdf | Y | PDG warped-context passage gives `M_KK` around `3.2 TeV`. |
| Gfitter public oblique fit `S = 0.05 +/- 0.11` | https://project-gfitter.web.cern.ch/Oblique_Parameters/ | Y | Gfitter page quotes `S = 0.05 +/- 0.11`. |
| Gfitter public oblique fit `T = 0.09 +/- 0.13` | https://project-gfitter.web.cern.ch/Oblique_Parameters/ | Y | Gfitter page quotes `T = 0.09 +/- 0.13`. |
| Gfitter public oblique fit `U = 0.01 +/- 0.11` | https://project-gfitter.web.cern.ch/Oblique_Parameters/ | Y | Gfitter page quotes `U = 0.01 +/- 0.11`. |
| HEPfit/de Blas 2022 CDF `M_W = 80.4335 +/- 0.0094 GeV` | https://arxiv.org/abs/2204.04204 | Y | Paper quotes the CDF W-mass value. |
| HEPfit/de Blas 2022 ST fit `S = 0.100 +/- 0.073` | https://arxiv.org/abs/2204.04204 | Y | Paper quotes the standard-average ST fit value. |
| HEPfit/de Blas 2022 ST fit `T = 0.202 +/- 0.056` | https://arxiv.org/abs/2204.04204 | Y | Paper quotes the standard-average ST fit value. |
| HEPfit/de Blas 2022 STU fit `U = 0.134 +/- 0.087` | https://arxiv.org/abs/2204.04204 | Y | Paper quotes the standard-average STU fit value. |

## EW002 - VERIFIED

Key references checked: PDG Vud/Vus review, arXiv:2411.04268, OSTI Hardy-Towner page, arXiv:2111.05338, arXiv:2212.06862, and arXiv:0804.1954.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| PDG2025:EW002:Vud_superallowed = `0.97367 +/- 0.00032` | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-vud-vus.pdf | Y | PDG review gives `Vud = 0.97367(32)`. |
| PDG2025:EW002:Vus_fplus_product = `0.21656 +/- 0.00035` | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-vud-vus.pdf | Y | PDG review gives `f_+(0)\|Vus\| = 0.21656(35)`. |
| PDG2025:EW002:Vus_Kl3 = `0.22330 +/- 0.00053` | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-vud-vus.pdf | Y | PDG review gives `Vus = 0.22330` with total uncertainty `0.00053`. |
| PDG2025:EW002:Vus_kaon_average = `0.22431 +/- 0.00085` | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-vud-vus.pdf | Y | PDG review gives kaon average `0.22431(85)`. |
| PDG2025:EW002:first_row_sum = `0.9983` | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-vud-vus.pdf | Y | PDG review quotes the first-row sum `0.9983` with listed uncertainties. |
| FLAG2024:EW002:Vus_fplus_product = `0.21654` | https://arxiv.org/abs/2411.04268 | Y | FLAG 2024 quotes `\|Vus\| f_+(0) = 0.21654(41)`. |
| FLAG2024:EW002:fplus_Nf211 = `0.9698` | https://arxiv.org/abs/2411.04268 | Y | FLAG 2024 quotes `f_+(0) = 0.9698(17)`. |
| FLAG2024:EW002:Vus_from_fplus_Nf211 = `0.22328` | https://arxiv.org/abs/2411.04268 | Y | FLAG 2024 quotes `Vus = 0.22328(58)`. |
| FLAG2024:EW002:first_row_sum_fplus = `0.99802` | https://arxiv.org/abs/2411.04268 | Y | FLAG 2024 quotes first-row sum `0.99802(66)`. |
| FLAG2024:EW002:first_row_sum_fkfpi = `0.99888` | https://arxiv.org/abs/2411.04268 | Y | FLAG 2024 quotes first-row sum `0.99888(67)`. |
| HardyTowner2020:EW002:Vud = `0.97373 +/- 0.00031` | https://www.osti.gov/pages/biblio/1780888 | Y | OSTI full-text page quotes `Vud = 0.97373 +/- 0.00031`. |
| BCBCI2022:EW002:Delta_CKM_global = `-0.00195` | https://arxiv.org/abs/2111.05338 | Y | Bryman et al. paper quotes `Delta_CKM = (-19.5 +/- 5.3) x 10^-4`. |
| CKKM2022:EW002:Caa_significance_context = about `3 sigma` | https://arxiv.org/abs/2212.06862 | Y | Paper abstract describes the Cabibbo-angle anomaly at about 3 sigma. |
| CFW2008:EW002:rs_kk_gluon_bound = `21 TeV` | https://arxiv.org/abs/0804.1954 | N/A | Theory-normalization context skipped under the signoff policy. |
| CFW2008:EW002:composite_higgs_kk_bound = `33 TeV` | https://arxiv.org/abs/0804.1954 | N/A | Theory-normalization context skipped under the signoff policy. |

## EW003 - VERIFIED

Key references checked: PDG 2024 Vcb/Vub review, arXiv:2411.18639, arXiv:2411.04268, arXiv:2404.02123, arXiv:2407.15172, and the BSM context sources where accessible.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| PDG 2024 `\|V_cb\|` inclusive = `42.2 +/- 0.5` | https://pdg.lbl.gov/2024/reviews/rpp2024-rev-vcb-vub.pdf | Y | PDG review quotes inclusive `\|V_cb\| = 42.2 +/- 0.5` in units of `10^-3`. |
| PDG 2024 `\|V_cb\|` exclusive = `39.8 +/- 0.6` | https://pdg.lbl.gov/2024/reviews/rpp2024-rev-vcb-vub.pdf | Y | PDG review quotes exclusive `\|V_cb\| = 39.8 +/- 0.6` in units of `10^-3`. |
| PDG 2024 `\|V_cb\|` scaled average = `41.1 +/- 1.2` | https://pdg.lbl.gov/2024/reviews/rpp2024-rev-vcb-vub.pdf | Y | PDG review quotes average `41.1 +/- 1.2` in units of `10^-3`. |
| PDG 2024 `\|V_ub\|` inclusive = `4.13` | https://pdg.lbl.gov/2024/reviews/rpp2024-rev-vcb-vub.pdf | Y | PDG review quotes inclusive `\|V_ub\| = 4.13` with listed uncertainties. |
| PDG 2024 `\|V_ub\|` exclusive = `3.70` | https://pdg.lbl.gov/2024/reviews/rpp2024-rev-vcb-vub.pdf | Y | PDG review quotes exclusive `\|V_ub\| = 3.70 +/- 0.10 +/- 0.12`. |
| PDG 2024 `\|V_ub\|` scaled average = `3.82 +/- 0.20` | https://pdg.lbl.gov/2024/reviews/rpp2024-rev-vcb-vub.pdf | Y | PDG review quotes average `3.82 +/- 0.20`. |
| HFLAV 2023 preferred inclusive `\|V_cb\| = 41.97 +/- 0.48` | https://arxiv.org/abs/2411.18639 | Y | HFLAV report quotes inclusive `\|V_cb\| = 41.97 +/- 0.48`. |
| HFLAV 2023 combined exclusive `\|V_cb\| = 39.77 +/- 0.46` | https://arxiv.org/abs/2411.18639 | Y | HFLAV report quotes exclusive `\|V_cb\| = 39.77 +/- 0.46`. |
| HFLAV 2023 combined exclusive `\|V_ub\| = 3.43 +/- 0.12` | https://arxiv.org/abs/2411.18639 | Y | HFLAV report quotes exclusive `\|V_ub\| = 3.43 +/- 0.12`. |
| FLAG 2024 `B -> D* l nu` `\|V_cb\| = 39.23 +/- 0.65` | https://arxiv.org/abs/2411.04268 | Y | FLAG 2024 quotes the `B -> D* l nu` value. |
| FLAG 2024 `B -> pi l nu` `\|V_ub\| = 3.61 +/- 0.16` | https://arxiv.org/abs/2411.04268 | Y | FLAG 2024 quotes the `B -> pi l nu` value. |
| PDG 2024 inclusive-exclusive `\|V_cb\|` marginal consistency = `3.0 sigma` | https://pdg.lbl.gov/2024/reviews/rpp2024-rev-vcb-vub.pdf | Y | PDG review states the inclusive/exclusive `\|V_cb\|` tension is about 3 sigma. |

## Wave-7 follow-up addendum: T003, T004, T008, T012

Date: 2026-05-16T19:27:00-04:00
Agent: factcheck-codex-top_higgs_ew
Batch label: wave7_top_higgs_ew

I re-fetched the cited CMS/ATLAS public pages, arXiv abstract/PDF pages, PDG
pdgLive pages, the PDG 2025 Z-boson PDF, and the local snapshots under
`flavor_catalog/references/<process_id>/`. ArXiv author/title metadata matched
the manifests for the key references. Theory-normalization-only rows are marked
`N/A` under the L001 signoff policy.

| Process | Verdict | Mismatches | Unresolvable | Note |
|---|---:|---:|---:|---|
| T003 | VERIFIED | 0 | 0 | CMS/ATLAS/PDG `tqgamma` limits and the quoted SM `t -> c gamma` prediction matched. |
| T004 | VERIFIED | 0 | 0 | PDG/ATLAS/CMS `t -> u gamma` limits, datasets, and projection context matched. |
| T008 | VERIFIED | 0 | 0 | PDG/CMS/ATLAS `t -> H u` limits and datasets matched. |
| T012 | VERIFIED | 0 | 0 | PDG 2025 and LEP/SLC `Z -> c cbar` pole observables matched. |

### T003 - VERIFIED

Key references checked: CMS-TOP-21-013 public page and arXiv:2312.08229,
arXiv:2205.02537, arXiv:hep-ph/0409342, and arXiv:0804.1954. Titles and
authors match the manifest records.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| CMS2024:T003:dataset = `138 fb^-1 at sqrt(s) = 13 TeV` | https://cms-results.web.cern.ch/cms-results/public-results/publications/TOP-21-013/index.html | Y | CMS public page/arXiv abstract states 138 fb^-1 at 13 TeV. |
| CMS2024:T003:tcgamma = `< 1.51e-5`, expected `< 1.54e-5` | https://cms-results.web.cern.ch/cms-results/public-results/publications/TOP-21-013/index.html | Y | CMS abstract gives observed(expected) `B(t -> c gamma) < 1.51 (1.54) x 10^-5` at 95% CL. |
| CMS2024:T003:tugamma = `< 0.95e-5`, expected `< 1.20e-5` | https://cms-results.web.cern.ch/cms-results/public-results/publications/TOP-21-013/index.html | Y | Same CMS source gives `B(t -> u gamma) < 0.95 (1.20) x 10^-5`. |
| ATLAS2023:T003:dataset = `139 fb^-1 at sqrt(s) = 13 TeV` | https://arxiv.org/abs/2205.02537 | Y | ATLAS abstract states 139 fb^-1 at 13 TeV. |
| ATLAS2023:T003:tcgamma_left = `< 4.2e-5` | https://arxiv.org/abs/2205.02537 | Y | ATLAS abstract gives the left-handed `t -> c gamma` limit as `4.2 x 10^-5`. |
| ATLAS2023:T003:tcgamma_right = `< 4.5e-5` | https://arxiv.org/abs/2205.02537 | Y | ATLAS abstract gives the right-handed `t -> c gamma` limit as `4.5 x 10^-5`. |
| PDG2026:T003:gammaq_summary = `< 9.5e-6` | https://pdgprod.lbl.gov/pdgprod/pdgLive/Particle.action?home=&node=Q007 | Y | pdgLive top listing gives `t -> gamma q (q = u,c) < 9.5 x 10^-6` at 95% CL; footnote confirms the ratio is to `Gamma(t -> Wb)`. |
| AguilarSaavedra2004:T003:SM = `4.6e-14` | https://arxiv.org/abs/hep-ph/0409342 | Y | Fetched PDF gives `Br(t -> c gamma) = 4.6 x 10^-14` with the quoted uncertainty components. |
| AguilarSaavedra2004:T003:dipole_normalization = `0.428 * lambda_qt^2` | https://arxiv.org/abs/hep-ph/0409342 | N/A | Formula appears in the source, but it is a theory normalization convention rather than a measured observable. |
| CsakiFalkowskiWeiler2008WarpedFlavor metadata | https://arxiv.org/abs/0804.1954 | N/A | Key-reference metadata checked; no measured T003 value is taken from this context source. |

### T004 - VERIFIED

Key references checked: PDGLive Q007, arXiv:2205.02537, CMS-TOP-21-013 /
arXiv:2312.08229, arXiv:0804.1954, and arXiv:1709.03975. Titles and authors
match the manifest records.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| PDG2026:T004:tgammaq_combined = `< 9.5e-6` | https://pdgprod.lbl.gov/pdgprod/pdgLive/Particle.action?home=&node=Q007 | Y | pdgLive top listing gives the combined `t -> gamma q (q = u,c)` row as `< 9.5 x 10^-6` at 95% CL. |
| ATLAS2023:T004:tugamma_left = `< 0.85e-5` | https://arxiv.org/abs/2205.02537 | Y | ATLAS abstract gives the left-handed `t -> u gamma` limit as `0.85 x 10^-5`. |
| ATLAS2023:T004:tugamma_right = `< 1.2e-5` | https://arxiv.org/abs/2205.02537 | Y | ATLAS abstract gives the right-handed `t -> u gamma` limit as `1.2 x 10^-5`. |
| CMS2024:T004:tugamma_observed = `< 0.95e-5`, expected `< 1.20e-5` | https://cms-results.web.cern.ch/cms-results/public-results/publications/TOP-21-013/ | Y | CMS public page/arXiv abstract gives observed(expected) `B(t -> u gamma) < 0.95 (1.20) x 10^-5`. |
| ATLAS2023:T004:dataset = `139 fb^-1 at sqrt(s) = 13 TeV` | https://arxiv.org/abs/2205.02537 | Y | ATLAS abstract states 139 fb^-1 at 13 TeV. |
| CMS2024:T004:dataset = `138 fb^-1 at sqrt(s) = 13 TeV` | https://cms-results.web.cern.ch/cms-results/public-results/publications/TOP-21-013/ | Y | CMS public page/arXiv abstract states 138 fb^-1 at 13 TeV. |
| AguilarSaavedra2017:T004:SM_and_projection_context | https://arxiv.org/abs/1709.03975 | Y | Source exists and quotes SM top-FCNC rates around `10^-14` for charm photon/Z modes, up-quark modes one order lower, and `10^-5`/`10^-6` HL-LHC/FCC-hh projection levels. |
| CFW2008:T004:rs_flavor_context | https://arxiv.org/abs/0804.1954 | N/A | Key-reference metadata checked; no measured T004 value is taken from this context source. |

### T008 - VERIFIED

Key references checked: PDGLive Q007R00, arXiv:2111.02219, arXiv:2404.02123,
arXiv:2407.15172, and arXiv:0804.1954. Titles and authors match the manifest
records.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| PDG2026:T008:tHu_headline = `< 1.9e-4` | https://pdgprod.lbl.gov/pdgprod/pdgLive/DataBlock.action?node=Q007R00&_eventName=showBR | Y | pdgLive datablock has units `10^-4` and the CMS `t -> H u` row `< 1.9` at 95% CL. |
| CMS2022:T008:dataset = `137 fb^-1 at sqrt(s) = 13 TeV` | https://arxiv.org/abs/2111.02219 | Y | CMS diphoton abstract states 137 fb^-1 at 13 TeV. |
| CMS2022:T008:tHu_diphoton = `< 0.019%`, expected `< 0.031%` | https://arxiv.org/abs/2111.02219 | Y | CMS abstract gives observed(expected) `B(t -> Hu) < 0.019 (0.031)%`. |
| ATLAS2024:T008:dataset = `140 fb^-1 at sqrt(s) = 13 TeV` | https://arxiv.org/abs/2404.02123 | Y | ATLAS multilepton abstract states 140 fb^-1 at 13 TeV. |
| ATLAS2024:T008:tHu_multilepton = `< 2.8e-4`, expected `< 3.0e-4` | https://arxiv.org/abs/2404.02123 | Y | ATLAS abstract gives observed(expected) `B(t -> Hu) < 2.8 (3.0) x 10^-4`. |
| ATLAS2024:T008:tHu_combined = `< 2.6e-4`, expected `< 1.8e-4` | https://arxiv.org/abs/2404.02123 | Y | ATLAS abstract gives combined observed(expected) `B(t -> Hu) < 2.6 (1.8) x 10^-4`. |
| CMS2025:T008:dataset = `138 fb^-1 at sqrt(s) = 13 TeV, collected in 2016-2018` | https://arxiv.org/abs/2407.15172 | Y | CMS same-sign dilepton abstract states 2016-2018 data, 13 TeV, 138 fb^-1. |
| CMS2025:T008:tHu_samesign = `< 0.072%`, expected `< 0.059%` | https://arxiv.org/abs/2407.15172 | Y | CMS abstract gives observed(expected) `B(t -> Hu) < 0.072% (0.059%)`. |
| CMS2025:T008:tHu_combined = `< 0.019%`, expected `< 0.027%` | https://arxiv.org/abs/2407.15172 | Y | CMS abstract gives combined observed(expected) `B(t -> Hu) < 0.019% (0.027%)`. |
| CsakiFalkowskiWeiler2008WarpedFlavor metadata | https://arxiv.org/abs/0804.1954 | N/A | Key-reference metadata checked; no measured T008 value is taken from this context source. |

### T012 - VERIFIED

Key references checked: PDG 2025 Z-boson PDF, arXiv:hep-ex/0509008,
arXiv:0804.1954, arXiv:0807.4937, arXiv:1401.2447, and arXiv:2107.00616.
Titles and authors match the manifest records.

| Claim | source_url | Verified? | Note |
|---|---|---:|---|
| `R_c^0 = 0.1721 +/- 0.0030` | https://pdg.lbl.gov/2025/listings/rpp2025-list-z-boson.pdf | Y | PDG 2025 Z listing gives `R_c = Gamma(cc)/Gamma(hadrons) = 0.1721 +/- 0.0030 OUR FIT`. |
| `A_FB^{0,c} = 0.0707 +/- 0.0035` | https://pdg.lbl.gov/2025/listings/rpp2025-list-z-boson.pdf | Y | PDG listing gives `7.07 +/- 0.35%`; catalog converts this to `0.0707 +/- 0.0035`. |
| `A_c = 0.670 +/- 0.027` | https://pdg.lbl.gov/2025/listings/rpp2025-list-z-boson.pdf | Y | PDG listing gives `0.670 +/- 0.027 OUR FIT`. |
| LEP/SLC final combination repeats `R_c^0 = 0.1721 +/- 0.0030` | https://arxiv.org/abs/hep-ex/0509008 | Y | Fetched LEP/SLC PDF contains the same `R_c` value. |
| LEP/SLC final combination repeats `A_FB^{0,c} = 0.0707 +/- 0.0035` | https://arxiv.org/abs/hep-ex/0509008 | Y | Fetched LEP/SLC PDF contains the same charm forward-backward asymmetry. |
| LEP/SLC final combination repeats `A_c = 0.670 +/- 0.027` | https://arxiv.org/abs/hep-ex/0509008 | Y | Fetched LEP/SLC PDF contains the same `A_c` value. |
| cfw_2008_rs_flavor metadata | https://arxiv.org/abs/0804.1954 | N/A | Context reference checked; no measured T012 value is taken from this paper. |
| casagrande_2008_rs_ewpt metadata | https://arxiv.org/abs/0807.4937 | N/A | Context reference checked; no measured T012 value is taken from this paper. |
| freitas_2014_z_widths metadata | https://arxiv.org/abs/1401.2447 | N/A | Theory-update reference checked; no measured T012 value is taken from this paper. |
| alcaraz_2021_z_lineshape_fcc metadata | https://arxiv.org/abs/2107.00616 | N/A | Projection/context reference checked; no current measured T012 value is taken from this paper. |

No mismatches or unresolvable URLs were found for the wave-7 follow-up batch.
