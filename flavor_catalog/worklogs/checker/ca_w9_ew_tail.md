# CA Worklog: ca_w9_ew_tail
**Date**: 2026-05-17
**Family**: collider_rs
**Batch ID**: Wave-9 EW-precision-tail, PRIMARY collider_rs
**Cycle**: 1
**Checker agent**: CA-CA-w9-ew-tail
**Process IDs**: CR009 CR011 CR012 CR013

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| CR009 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| CR011 | PASS | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| CR012 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| CR013 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| CR009 | llqq contact-interaction scale Lambda, combined dielectron+dimuon LL constructive | Lambda_LL^+ > 35.8 TeV | PDG 2025 compositeness review quoting ATLAS arXiv:2006.12946 | 4bed92264b2c |
| CR009 | llqq contact-interaction scale Lambda, combined dielectron+dimuon LL destructive | Lambda_LL^- > 26.0 TeV | PDG 2025 compositeness review quoting ATLAS arXiv:2006.12946 | 4bed92264b2c |
| CR009 | llqq contact-interaction scale Lambda, combined dielectron+dimuon RR constructive | Lambda_RR^+ > 35.5 TeV | PDG 2025 compositeness review quoting ATLAS arXiv:2006.12946 | 4bed92264b2c |
| CR009 | llqq contact-interaction scale Lambda, combined dielectron+dimuon RR destructive | Lambda_RR^- > 26.5 TeV | PDG 2025 compositeness review quoting ATLAS arXiv:2006.12946 | 4bed92264b2c |
| CR009 | llqq contact-interaction scale Lambda, combined dielectron+dimuon LR constructive | Lambda_LR^+ > 32.5 TeV | PDG 2025 compositeness review quoting ATLAS arXiv:2006.12946 | 4bed92264b2c |
| CR009 | llqq contact-interaction scale Lambda, combined dielectron+dimuon LR destructive | Lambda_LR^- > 28.8 TeV | PDG 2025 compositeness review quoting ATLAS arXiv:2006.12946 | 4bed92264b2c |
| CR009 | llqq contact-interaction scale Lambda, CMS lower endpoint, LL destructive | Lambda_LL > 23.9 TeV | PDG 2025 compositeness review quoting CMS arXiv:2103.02708 | 4bed92264b2c |
| CR009 | llqq contact-interaction scale Lambda, CMS upper endpoint, RR constructive | Lambda_RR > 36.4 TeV | PDG 2025 compositeness review quoting CMS arXiv:2103.02708 | 4bed92264b2c |
| CR011 | Fiducial cross-section upper limit for pp -> jj W_L^+/- W_L^+/- | sigma_fid(pp -> jj W_L^+/- W_L^+/-) < 0.45 fb (95% CL) | ATLAS Collaboration, Phys. Rev. Lett. 135 (2025) 111802 | 449b1c2c9f28 |
| CR011 | Production cross-section upper limit for longitudinally polarized same-sign W^+/- W^+/- pairs | sigma(pp -> jj W_L^+/- W_L^+/-) < 1.17 fb (95% CL) | CMS Collaboration, Phys. Lett. B 812 (2020) 136018 | 6ceb6660bd95 |
| CR012 | HVT model B W' -> WZ observed resonance-mass lower limit | m(W') > 4.4 TeV | PDG 2025 heavy-boson listing, document TUMASYAN 23AP; CMS-B2G-20-009 Table 2 | 987f6211a9b5 |
| CR012 | HVT model B mass-degenerate V' -> VV observed resonance-mass lower limit | m(V') > 4.5 TeV for the mass-degenerate VV category | CMS-B2G-20-009 Table 2 | 5560f5834e4b |
| CR012 | HVT model B mass-degenerate V' -> VV+VH observed resonance-mass lower limit | m(V') > 4.8 TeV when VV and VH channels are combined | CMS-B2G-20-009 Table 2; PDG 2025 comment for TUMASYAN 23AP | 5560f5834e4b |
| CR012 | HVT model B Z' -> WW observed excluded resonance-mass intervals | m(Z') excluded in 1.3-3.1 TeV and 3.3-3.5 TeV | CMS-B2G-20-009 Table 2 | 5560f5834e4b |
| CR013 | Lower limit on lightest RS KK graviton mass in pp -> G* -> gamma gamma | m_G* > 4.8 TeV | PDG 2025 Extra Dimensions listing, document HAYRAPETYAN 24AJ | ce4823040682 |
| CR013 | RS graviton mass-exclusion range in pp -> G* -> gamma gamma | m_G* below 2.2-5.6 TeV excluded for 0.01 < ktilde < 0.2 | CMS EXO-22-024 public result and arXiv:2405.09320 | f3ee8d2e3e10 |
| CR013 | Lower limit on RS1 graviton mass in pp -> G* -> gamma gamma | m_G* > 4.5 TeV for k/M_Pl = 0.1 | ATLAS, Phys. Lett. B 822 (2021) 136651, arXiv:2102.13405 | eb98469106b1 |
| CR013 | Lower limit on RS1 graviton mass in pp -> G* -> gamma gamma | m_G* > 3.9 TeV for k/M_Pl = 0.05 | ATLAS, Phys. Lett. B 822 (2021) 136651, arXiv:2102.13405 | eb98469106b1 |
| CR013 | Lower limit on RS1 graviton mass in pp -> G* -> gamma gamma | m_G* > 2.2 TeV for k/M_Pl = 0.01 | ATLAS, Phys. Lett. B 822 (2021) 136651, arXiv:2102.13405 | eb98469106b1 |
| CR013 | RS graviton mass-exclusion range in pp -> G* -> gamma gamma | m_G* below 2.3-4.6 TeV excluded for ktilde = 0.01-0.2 | CMS, Phys. Rev. D 98 (2018) 092001, arXiv:1809.00327 | 6cbc675c212d |
| CR013 | Lower limit on lightest RS graviton mass in pp -> G* -> gamma gamma | m_G* > 2.66 TeV for k/M_Pl = 0.1 | ATLAS, Phys. Rev. D 92 (2015) 032004, arXiv:1504.05511 | fa8a6d18d563 |

## Per-process issues
- CR009: CHK-1 fails. The TeX post-2010 history quotes historical measured contact-interaction limit ranges, including \(O(10~\mathrm{TeV})\), \(15.4\)--\(26.3~\mathrm{TeV}\), \(16.7\)--\(25.2~\mathrm{TeV}\), \(24\)--\(40~\mathrm{TeV}\), and \(20\)--\(32~\mathrm{TeV}\), outside `pdg_or_equivalent.values` (`CR009.tex:71`, `CR009.tex:74`, `CR009.tex:80`-`83`). These are measured experimental exclusion limits, not luminosity/theory/EFT-translation carve-outs. Rework should either remove the historical numeric ranges from prose/sidecar `why_matters` text or promote them as structured `pdg_or_equivalent.values` entries with full strict metadata.
- CR011: CHK-2 fails. The TeX `Key references` section lists snapshot-name stems such as `atlas_2025_arxiv2503_11317` and `cms_2020_arxiv2009_09429`, but `source_manifest.yaml` keys are `ATLAS2025_LongitudinalWW`, `CMS2020_PolarizedSameSignWW`, `ATLAS2024_SameSignWWjj`, `ATLAS2020_ZZjj`, `ATLAS2017_WWjjAQGC`, `ATLAS2026_AQGCCombination`, `PDG2025_WZQuarticCouplings`, and `Contino2011_CompositeHiggsResonances`. Rework should make the key-reference list resolve to manifest keys.
- CR012: None.
- CR013: None.

## Evidence notes
- CHK-1: Applied the L001/B001/B021 carve-out. Dataset luminosities, \(\sqrt{s}\), expected limits carried inside an observed-limit block, HEPData table availability, theory frameworks, EFT translations, and the local 47.26/127.13 TeV methodology-note comparison were not required in `pdg_or_equivalent.values`. CR012 and CR013 measured mass-exclusion claims are in `pdg_or_equivalent.values` with year, value, uncertainty/null, units, source URL, access date, snapshot path, and sha256. CR011 measured cross-section upper limits are likewise promoted. CR009 current ATLAS/CMS promoted limits pass, but the historical measured contact-limit ranges in the TeX do not.
- CHK-2: `source_manifest.yaml` validation found all manifest entries for CR009, CR011, CR012, and CR013 have non-empty `snapshot_path` values tracked by `git ls-files`; `sha256sum -c sha256sums.txt` passed in all four reference directories. CR009, CR012, and CR013 TeX key-reference tokens exactly match manifest keys. CR011 fails because the TeX key-reference tokens are snapshot stems rather than manifest keys.
- CHK-3: `find flavor_catalog/references/CR009 flavor_catalog/references/CR011 flavor_catalog/references/CR012 flavor_catalog/references/CR013 -type f -iname '*.pdf'` returned no files. The directories contain text/XML/HTML/JSON snapshots, manifests, and sha256 lists only.
- CHK-4: Before CA mutation, all four sidecars showed `WRITER-INITIATED -> WRITER-DONE` with ISO 8601 `timestamp` and `at` fields. CR009 and CR011 now receive `WRITER-REWORK`; CR012 and CR013 receive `CHECKER-DONE`.
- CHK-5: The `NO` code-coverage claims are consistent. Direct CR009 DY/contact and CR011 VBS/aQGC greps returned no implementation hits. The CR012 collider/recast grep returned only `tests/test_alpha_s.py:89`, a CMS/RunDec alpha-s example. The CR013 direct-search grep returned only false-positive `m_GeV` field names in `quarkConstraints/paper_0710_1869/verifier.py`, not diphoton/RS-graviton collider logic. The cited adjacent evidence lines exist: `quarkConstraints/scan.py:359`, `:376`, `:377`, `:379`; `scanParams/scan.py:399`, `:429`, `:523`, `:528`; `quarkConstraints/deltaf2.py:1`, `:320`, `:459`, `:557`; `solvers/bessel.py:5`-`:6`; and `tests/test_alpha_s.py:88`-`:89`.
- CHK-6: `HIGH` is consistent for all four. CR009 needs a binned high-mass dilepton likelihood or equivalent SMEFT/collider reinterpretation; CR011 needs polarization-sensitive VBS simulation/recast machinery; CR012 needs resonance production, widths, branching fractions, and collider likelihood/recast support; CR013 needs a spin-2 spectrum/branching model and diphoton upper-limit reinterpretation.
- CHK-7: No entry contradicts the rc1.1/methodology-note low-energy flavor result. The entries repeatedly state that LHC EW-tail/diphoton/diboson reach is weaker than the default anarchic-flavor \(M_{KK}^{\min}(p50,g_{s*}=3)=47.26~\mathrm{TeV}\) comparison, and the methodology-note lines checked retain 47.26 TeV and 127.13 TeV.
- CHK-8: No unresolved `\cite`, `\ref`, `\textbf{CHECK}`, `TODO`, `FIXME`, missing-reference marker, or undefined-reference marker was found in the four TeX files. Lightweight balance checks gave matched `\(`/`\)`, matched `\[`/`\]`, and brace delta 0 for each TeX file.
