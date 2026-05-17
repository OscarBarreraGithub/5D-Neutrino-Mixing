# Flavor Catalog — Fact-Check Status (v0 → v0.1)

**Goal**: independently verify that each cited paper EXISTS *and* actually
CONTAINS the numerical value(s) the catalog entry claims to take from it.
Not just URL HTTP-200; the paper content must match.

**Method**: per-family codex (gpt-5.5 xhigh) fact-checkers with WebFetch.
Each agent reads its family's `.yaml` source_manifest, fetches the cited
arXiv / PDG / HFLAV / FLAG / experiment-collaboration pages, and verifies the
claimed numerical values appear at the cited source.

**Prior review credit**: the original PKA/WA/CA cycles + Opus round-1/round-2
sign-offs verified internal consistency (every TeX claim resolves to a manifest
entry; sha256 of snapshot matches; status_history transitions clean). They did
NOT independently re-fetch every source. This pass closes that gap.

**Status legend**:
- `PENDING` — not yet fact-checked.
- `VERIFIED` — all cited sources fetched; all claimed values match.
- `PARTIAL` — most sources verified; some unresolvable (e.g. paywalled or 404).
- `MISMATCH` — at least one claimed value does not match the cited source.
- `FAILED` — multiple claims unverifiable.

Each per-family fact-checker writes a detailed report at
`flavor_catalog/audits/factcheck_<family>.md`. Then an aggregator compiles
the rows below.

## Per-process status table

| process_id | family | status | date | fact-checker | mismatches/notes |
|---|---|---|---|---|---|
| K001 | kaon | VERIFIED | 2026-05-16 | factcheck-codex-kaon | PDG epsilon_K, BGS SM prediction, and FLAG bag values match. |
| K002 | kaon | VERIFIED | 2026-05-16 | factcheck-codex-kaon | PDG Delta m_K fits, lattice inputs, and KTeV values match. |
| K003 | kaon | VERIFIED | 2026-05-16 | factcheck-codex-kaon | PDG, KTeV, NA48, RBC/UKQCD, and Aebischer-Buras epsilon-prime values match. |
| K004 | kaon | VERIFIED | 2026-05-16 | factcheck-codex-kaon | NA62 K+ -> pi+ nu nubar measurements and SM predictions match. |
| K005 | kaon | VERIFIED | 2026-05-16 | factcheck-codex-kaon | PDG/KOTO KL -> pi0 nu nubar limit, NA62 context, and SM values match. |
| K006 | kaon | VERIFIED | 2026-05-16 | factcheck-codex-kaon | PDG KL -> mu mu value, BNL E871 inputs, and supporting theory values match. |
| K008 | kaon | VERIFIED | 2026-05-16 | factcheck-codex-kaon | KL -> pi0 e e PDG/KTeV/NA48/ISU values match; RS context skipped by policy. |
| K009 | kaon | VERIFIED | 2026-05-16 | factcheck-codex-kaon | KL -> pi0 mu mu PDG/KTeV/NA48/ISU values match; RS context skipped by policy. |
| K010 | kaon | VERIFIED | 2026-05-16 | factcheck-codex-kaon | KS -> pi0 e e PDG/NA48 branching values and event counts match. |
| K012 | kaon | VERIFIED | 2026-05-16 | factcheck-codex-kaon | KS -> mu mu PDG/LHCb limits and theory-context values match. |
| K013 | kaon | VERIFIED | 2026-05-16 | factcheck-codex-kaon | KL -> pi0 gamma gamma PDG average, KTeV/NA48 inputs, and NA48/2 values match. |
| K017 | kaon | VERIFIED | 2026-05-16 | factcheck-codex-kaon | R_K PDG API, NA62, KLOE, and SM prediction values match. |
| K018 | kaon | VERIFIED | 2026-05-16 | factcheck-codex-kaon | K_l3 PDG, FlaviaNet, and FLAG V_us f_+(0) values match. |
| C001 | charm | VERIFIED | 2026-05-16 | factcheck-codex-charm | HFLAV/PDG x_D, y_D, Delta m_D, Delta y, and Delta Gamma/Gamma values match. |
| C002 | charm | VERIFIED | 2026-05-16 | factcheck-codex-charm | HFLAV/PDG q/p, phi_D, no-indirect-CPV test, and metadata checks match. |
| C003 | charm | VERIFIED | 2026-05-16 | factcheck-codex-charm | HFLAV direct/indirect CPV fit and LHCb Delta A_CP values match. |
| C004 | charm | VERIFIED | 2026-05-16 | factcheck-codex-charm | D0 -> mu mu PDG/CMS/LHCb limits and long-distance theory inputs match. |
| C005 | charm | VERIFIED | 2026-05-16 | factcheck-codex-charm | D0 -> e e PDG, Belle, and BABAR limits match. |
| C006 | charm | VERIFIED | 2026-05-16 | factcheck-codex-charm | D0 -> e mu PDG, LHCb, Belle, and BABAR LFV limits match. |
| C007 | charm | VERIFIED | 2026-05-16 | factcheck-codex-charm | D+ -> pi+ mu mu PDG/LHCb values, datasets, and search scope match. |
| C008 | charm | VERIFIED | 2026-05-16 | factcheck-codex-charm | D+ -> pi+ e mu PDG/LHCb/BABAR charge-mode limits match. |
| B001 | beauty | VERIFIED | 2026-05-16 | factcheck-codex-beauty | Delta m_d HFLAV, companion x_d/chi_d, and Belle II values match. |
| B002 | beauty | VERIFIED | 2026-05-16 | factcheck-codex-beauty | sin(2 beta) HFLAV, LHCb, Belle II, and penguin-shift values match. |
| B003 | beauty | VERIFIED | 2026-05-16 | factcheck-codex-beauty | Delta m_s HFLAV/PDG/LHCb values and companion quantities match. |
| B004 | beauty | VERIFIED | 2026-05-16 | factcheck-codex-beauty | phi_s, Delta Gamma_s, and LHCb Run-2 values match. |
| B005 | beauty | VERIFIED | 2026-05-16 | factcheck-codex-beauty | Bs -> mu mu PDG/HFLAV average, SM prediction, and experiment inputs match. |
| B006 | beauty | VERIFIED | 2026-05-16 | factcheck-codex-beauty | B0 -> mu mu PDG/HFLAV limits, SM prediction, and supporting inputs match. |
| B009 | beauty | VERIFIED | 2026-05-16 | factcheck-codex-beauty | B+ -> tau nu HFLAV/PDG/UTfit and B-factory inputs match. |
| B011 | beauty | VERIFIED | 2026-05-16 | factcheck-codex-beauty | B -> Xs gamma HFLAV/PDG averages, SM predictions, and projection context match. |
| B015 | beauty | VERIFIED | 2026-05-16 | factcheck-codex-beauty | B -> Xs l l HFLAV, BaBar/Belle inputs, and theory bin predictions match. |
| B016 | beauty | VERIFIED | 2026-05-16 | factcheck-codex-beauty | B+ -> K+ l l and B0 -> K0 l l HFLAV averages and key references match. |
| B017 | beauty | VERIFIED | 2026-05-16 | factcheck-codex-beauty | Inclusive/exclusive b -> s l l averages and LFU-ratio context match. |
| B018 | beauty | VERIFIED | 2026-05-16 | factcheck-codex-beauty | R_K HFLAV bin averages, LHCb 2021 value, and SM context match. |
| B019 | beauty | VERIFIED | 2026-05-16 | factcheck-codex-beauty | R_K* HFLAV bin averages, LHCb values, and SM theory context match. |
| B021 | beauty | VERIFIED | 2026-05-16 | factcheck-codex-beauty | Lambda_b -> Lambda mu mu PDG/LHCb/CDF branching and angular values match. |
| B022 | beauty | VERIFIED | 2026-05-16 | factcheck-codex-beauty | B+ -> K+ nu nubar HFLAV/PDG/Belle II, BaBar limit, and SM prediction match. |
| B023 | beauty | VERIFIED | 2026-05-16 | factcheck-codex-beauty | B -> K* nu nubar PDG/HFLAV limits and SM prediction match. |
| B025 | beauty | VERIFIED | 2026-05-16 | factcheck-codex-beauty | R(D) HFLAV CKM2025 average, SM comparator, tension, and inputs match. |
| B026 | beauty | VERIFIED | 2026-05-16 | factcheck-codex-beauty | R(D*) HFLAV CKM2025 average, SM comparator, tension, and inputs match. |
| B032 | beauty | VERIFIED | 2026-05-16 | factcheck-codex-beauty | B -> K pi branching fractions, direct CP values, PDG CP parameters, and sum-rule inputs match. |
| B033 | beauty | VERIFIED | 2026-05-16 | factcheck-codex-beauty | S(phi K0), C(phi K0), sin(2 beta) comparator, and Belle II inputs match. |
| B034 | beauty | VERIFIED | 2026-05-16 | factcheck-codex-beauty | Bs -> phi phi CP parameters, branching fraction, and earlier LHCb/CDF values match. |
| T001 | top_higgs_ew | VERIFIED | 2026-05-16 | factcheck-codex-top_higgs_ew | t -> Zc PDG/ATLAS/CMS limits, SM estimate, and HL-LHC projection match. |
| T002 | top_higgs_ew | VERIFIED | 2026-05-16 | factcheck-codex-top_higgs_ew | t -> Zu PDG/ATLAS/CMS limits and SM estimate match; RS context skipped by policy. |
| T005 | top_higgs_ew | VERIFIED | 2026-05-16 | factcheck-codex-top_higgs_ew | t -> cg PDG/ATLAS/CMS limits and SM estimate match. |
| T006 | top_higgs_ew | VERIFIED | 2026-05-16 | factcheck-codex-top_higgs_ew | t -> ug PDG/ATLAS/CMS limits and SM estimate match. |
| T007 | top_higgs_ew | VERIFIED | 2026-05-16 | factcheck-codex-top_higgs_ew | t -> Hc PDG/ATLAS/CMS limits and datasets match. |
| T010 | top_higgs_ew | VERIFIED | 2026-05-16 | factcheck-codex-top_higgs_ew | Zbb PDG, LEP/SLC, and FCC-ee values match. |
| T015 | top_higgs_ew | VERIFIED | 2026-05-16 | factcheck-codex-top_higgs_ew | Z -> e mu CMS/PDG/ATLAS limits and Tera-Z prospect match. |
| T016 | top_higgs_ew | VERIFIED | 2026-05-16 | factcheck-codex-top_higgs_ew | Z -> e tau PDG/ATLAS/CMS limits and Tera-Z prospect match. |
| T017 | top_higgs_ew | VERIFIED | 2026-05-16 | factcheck-codex-top_higgs_ew | Z -> mu tau PDG/ATLAS/CMS limits and datasets match. |
| T018 | top_higgs_ew | VERIFIED | 2026-05-16 | factcheck-codex-top_higgs_ew | h -> mu tau PDG/CMS/ATLAS/CMS-2015 values and theory context match. |
| T019 | top_higgs_ew | VERIFIED | 2026-05-16 | factcheck-codex-top_higgs_ew | h -> e tau PDG/CMS/ATLAS values and theory context match. |
| T020 | top_higgs_ew | VERIFIED | 2026-05-16 | factcheck-codex-top_higgs_ew + factcheck-followup | Family report found an ATLAS h -> e mu 6.1e-5 mismatch; follow-up fixed the catalog to 6.2e-5 observed and 5.9e-5 expected, so current status is VERIFIED. |
| EW001 | top_higgs_ew | VERIFIED | 2026-05-16 | factcheck-codex-top_higgs_ew | PDG/Gfitter/HEPfit oblique S, T, U values match; blocked DOI was theory-only. |
| EW002 | top_higgs_ew | VERIFIED | 2026-05-16 | factcheck-codex-top_higgs_ew | PDG/FLAG/Hardy-Towner/CKM-anomaly values match; RS context skipped by policy. |
| EW003 | top_higgs_ew | VERIFIED | 2026-05-16 | factcheck-codex-top_higgs_ew | PDG/HFLAV/FLAG V_cb and V_ub values match; blocked APS page was BSM context only. |
| L001 | charged_lepton | VERIFIED | 2026-05-16 | factcheck-codex-charged_lepton | mu -> e gamma MEG II, PDG, MEG II 2024, and MEG 2016 values match. |
| L002 | charged_lepton | VERIFIED | 2026-05-16 | factcheck-codex-charged_lepton | mu -> 3e PDG/SINDRUM limit and Mu3e sensitivity values match. |
| L003 | charged_lepton | VERIFIED | 2026-05-16 | factcheck-codex-charged_lepton | mu-e conversion Al projections and Au benchmark match; no direct published Al limit claimed. |
| L004 | charged_lepton | VERIFIED | 2026-05-16 | factcheck-codex-charged_lepton | mu-e conversion Au PDG/SINDRUM II and COMET/Mu2e context values match. |
| L005 | charged_lepton | VERIFIED | 2026-05-16 | factcheck-codex-charged_lepton | mu-e conversion Ti PDG/SINDRUM II and COMET/Mu2e context values match. |
| L006 | charged_lepton | VERIFIED | 2026-05-16 | factcheck-codex-charged_lepton | Muonium-antimuonium PDG/Willmann/MACE values match; 8.2e-11 arXiv vs 8.3e-11 PDG nuance documented. |
| L007 | charged_lepton | VERIFIED | 2026-05-16 | factcheck-codex-charged_lepton | tau -> mu gamma PDG, Belle, BaBar, and Belle II projection values match. |
| L008 | charged_lepton | VERIFIED | 2026-05-16 | factcheck-codex-charged_lepton | tau -> e gamma PDG, BaBar, Belle, and Belle II projection values match. |
| L009 | charged_lepton | VERIFIED | 2026-05-16 | factcheck-codex-charged_lepton | tau -> 3mu PDG, Belle II, LHCb, Belle, and BaBar values match. |
| L010 | charged_lepton | VERIFIED | 2026-05-16 | factcheck-codex-charged_lepton | tau -> 3e PDG, Belle, BaBar, and Belle II prospect values match. |
| L023 | charged_lepton | VERIFIED | 2026-05-16 | factcheck-codex-charged_lepton | Neutrino trident CHARM-II, CCFR, NuTeV, Altmannshofer, DUNE, and Belle II values match. |
| E001 | edm_neutrino | VERIFIED | 2026-05-16 | factcheck-codex-edm_neutrino | Electron EDM PDG, Roussy 2023, ACME 2018, and context metadata match. |
| E002 | edm_neutrino | VERIFIED | 2026-05-16 | factcheck-codex-edm_neutrino | Muon EDM PDG, Bennett 2009, PSI proposal, and Renga 2024 projection values match. |
| E004 | edm_neutrino | VERIFIED | 2026-05-16 | factcheck-codex-edm_neutrino | Neutron EDM PDG, Abel 2020, Pendlebury, and key theory/lattice metadata match. |
| E006 | edm_neutrino | VERIFIED | 2026-05-16 | factcheck-codex-edm_neutrino | Hg-199 EDM Graner, PDG neutron cross-reference rows, and Sahoo interpretation values match. |
| E007 | edm_neutrino | VERIFIED | 2026-05-16 | factcheck-codex-edm_neutrino | Ra-225 and Xe-129 EDM limits plus Argonne/DFG/KU Leuven projection context match. |
| E008 | edm_neutrino | VERIFIED | 2026-05-16 | factcheck-codex-edm_neutrino | qCEDM neutron/Hg anchors, derived bounds, and quoted response formulae match. |
| E009 | edm_neutrino | PARTIAL | 2026-05-16 | factcheck-codex-edm_neutrino | All value-bearing neutron EDM and Weinberg-operator claims verified; INSPIRE Weinberg 1989 key-reference URL rendered JS-only, cross-checked via APS and local snapshot. |
| CR001 | collider_rs | VERIFIED | 2026-05-17 | factcheck-codex-collider_rs-w9 | CMS 2026 KK-gluon to ttbar 5.5 TeV exclusion, CMS 2019, PDG Live, and ATLAS historical limits match; scan/theory context skipped by policy. |
| CR002 | collider_rs | VERIFIED | 2026-05-17 | factcheck-codex-collider_rs-w9 | PDG/ATLAS/CMS T5/3 pair-production mass limits match; projection and scan-context items skipped by policy. |
| CR003 | collider_rs | VERIFIED | 2026-05-17 | factcheck-codex-collider_rs-w9 | PDG Live and ATLAS/CMS charge-2/3 top-partner mass-limit rows match. |
| CR004 | collider_rs | VERIFIED | 2026-05-17 | factcheck-codex-collider_rs-w9 | PDG 2025 bottom-partner bH, bZ, and tW mass limits match. |
| CR005 | collider_rs | VERIFIED | 2026-05-17 | factcheck-codex-collider_rs-w9 | CMS/ATLAS/PDG high-mass dilepton mass and fiducial cross-section limits match. |
| CR006 | collider_rs | VERIFIED | 2026-05-17 | factcheck-codex-collider_rs-w9 | PDG/CMS charged-resonance lnu and tb mass limits match. |
| CR007 | collider_rs | VERIFIED | 2026-05-17 | factcheck-codex-collider_rs-w9 | PDG RS/bulk-graviton and ATLAS bulk-graviton limits match. |
| CR008 | collider_rs | VERIFIED | 2026-05-17 | factcheck-codex-collider_rs-w9 | ATLAS singlet-T and CMS all-branching-envelope mass limits match. |
| CR009 | collider_rs | VERIFIED | 2026-05-17 | factcheck-codex-collider_rs-w9 | PDG ATLAS/CMS compositeness-scale endpoints match; EFT translation skipped by policy. |
| CR010 | collider_rs | VERIFIED | 2026-05-17 | factcheck-codex-collider_rs-w9 | PDG VLQ doublet and simplified T/B endpoint limits match. |
| CR011 | collider_rs | VERIFIED | 2026-05-17 | factcheck-codex-collider_rs-w9 | ATLAS/CMS longitudinal same-sign WW fiducial cross-section limits match; aQGC/EFT context skipped by policy. |
| CR012 | collider_rs | VERIFIED | 2026-05-17 | factcheck-codex-collider_rs-w9 | PDG/CMS HVT diboson mass limits and CMS Table 2 intervals match. |
| CR013 | collider_rs | VERIFIED | 2026-05-17 | factcheck-codex-collider_rs-w9 | PDG/CMS/ATLAS diphoton RS-graviton mass limits match. |
| CR014 | collider_rs | VERIFIED | 2026-05-17 | factcheck-codex-collider_rs-w9 | PDG four-top cross sections, CMS top-philic vector limit, and ATLAS cross-section-limit range match. |

## Summary

Status totals across 89 process entries:

| status | count |
|---|---:|
| VERIFIED | 88 |
| PARTIAL | 1 |
| MISMATCH | 0 |
| FAILED | 0 |

Family totals:

| family | VERIFIED | PARTIAL | MISMATCH | FAILED | total |
|---|---:|---:|---:|---:|---:|
| kaon | 13 | 0 | 0 | 0 | 13 |
| charm | 8 | 0 | 0 | 0 | 8 |
| beauty | 21 | 0 | 0 | 0 | 21 |
| top_higgs_ew | 15 | 0 | 0 | 0 | 15 |
| charged_lepton | 11 | 0 | 0 | 0 | 11 |
| edm_neutrino | 6 | 1 | 0 | 0 | 7 |
| collider_rs | 14 | 0 | 0 | 0 | 14 |

T020 was originally reported as `MISMATCH` in
`factcheck_top_higgs_ew.md`: the family audit found the ATLAS `h -> e mu`
limit had been quoted as `< 6.1 x 10^-5` while the live ATLAS/PDG sources give
`< 6.2 x 10^-5` observed and `< 5.9 x 10^-5` expected. A follow-up fix updated
the catalog sidecar/prose and set the process fact-check verdict to
`VERIFIED`, so T020 is counted as `VERIFIED` here.

Remaining issues:

- `E009` remains `PARTIAL` because the INSPIRE URL for the Weinberg 1989
  key reference rendered as a JavaScript-only page under WebFetch. The same
  metadata was cross-checked against APS and the local snapshot, and no
  value-bearing catalog claim depends on the unresolved INSPIRE page.
- `EW001` and `EW003` had blocked publisher pages for theory/context references
  (`HTTP 403`), but their measured/catalog numerical values were verified
  against PDG, Gfitter, HEPfit, HFLAV, and FLAG sources.
- `B009` required a TLS/certificate workaround for the UTfit page; the content
  was fetched and the quoted value was visible.
- `L006` retains the documented convention/source nuance where the Willmann
  arXiv abstract says `8.2e-11` while PDG and the catalog use the canonical
  `8.3e-11` value.

Recommendation: ready to tag `flavor-catalog-v0.1`. There are no remaining
value-bearing mismatches or failed rows. The only open status is the E009
metadata-access caveat, which is documented, cross-checked through an alternate
source, and non-blocking for the numerical catalog.
