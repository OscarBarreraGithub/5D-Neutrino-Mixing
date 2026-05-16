# CA Worklog: ca_w6_edm_v2
**Date**: 2026-05-16
**Family**: edm_neutrino
**Process IDs**: E009

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| E009 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| E009 | \(|d_n|\) direct anchor | \(<1.8\times10^{-26}\ e\,\mathrm{cm}\), 90% CL | PDG Live 2026 S017EDM | ea22e4f3da86 |
| E009 | \(d_n\) primary measurement | \((0.0\pm1.1_{\rm stat}\pm0.2_{\rm sys})\times10^{-26}\ e\,\mathrm{cm}\) | Abel et al. 2020, arXiv:2001.11966 | a78b99dfd396 |
| E009 | Pospelov-Ritz Weinberg normalization | \(|d_n(w)|\simeq e\,22~\mathrm{MeV}\,w(1~\mathrm{GeV})\) | Pospelov and Ritz 2005, arXiv:hep-ph/0504231 | e07e62e5a797 |
| E009 | Derived \(w(1~\mathrm{GeV})\) benchmark | \(<4.1\times10^{-11}\ \mathrm{GeV}^{-2}\), central one-source translation | PDG Live 2026 S017EDM plus Pospelov-Ritz 2005 | ea22e4f3da86 / e07e62e5a797 |
| E009 | Haisch-Hala \(O_6\) normalization | \((d_n/e)_{O_6}=74(1\pm0.5)~\mathrm{MeV}\) | Haisch and Hala 2019, arXiv:1909.08955 | ae1c67f2535e |
| E009 | Derived \(C_6\) benchmark | \(<1.2\times10^{-11}\ \mathrm{GeV}^{-2}\), central one-source translation | PDG Live 2026 S017EDM plus Haisch-Hala 2019 | ea22e4f3da86 / ae1c67f2535e |

## Issues (if any)
- None.

## Evidence notes
- CHK-1 PASS: `rg -n "[0-9]" flavor_catalog/processes/edm_neutrino/E009.tex` finds the neutron EDM limit and Abel measurement in lines 22-28, plus the Pospelov-Ritz and Haisch-Hala benchmark numbers in lines 53-60. After WA-v2, all value-bearing numerical claims have corresponding `pdg_or_equivalent` entries with year, value, uncertainty or CL, source URL, access date, and sha256 metadata. Process IDs, source years, file-line citations, and \(\Delta F=2\) are not treated as measured/equivalent observable values.
- CHK-2 PASS: every process-local source key cited in the TeX matches an entry in `flavor_catalog/references/E009/source_manifest.yaml`: `PDG2026:NeutronEDM`, `Abel:NeutronEDM2020`, `Weinberg:ThreeGluon1989`, `PospelovRitz:EDMReview2005`, `ChuppRamseyMusolf:EDMGlobal2014`, `KoenigNeubertStraub:CompositeDipoles2014`, `HaischHala:WeinbergSumRules2019`, `Bhattacharya:WeinbergLattice2022`, and `CsakiFalkowskiWeiler:RSFlavor2008`. Each manifest entry has a non-empty `snapshot_path`, and `git ls-files --error-unmatch` confirms each snapshot is tracked under `flavor_catalog/references/E009/`.
- CHK-3 PASS: `ls -la flavor_catalog/references/E009` and `git ls-files flavor_catalog/references/E009` show only `.txt` snapshots plus `source_manifest.yaml`; `find flavor_catalog/references/E009 -maxdepth 1 -type f -iname '*.pdf' -print` returned no files.
- CHK-4 PASS: `E009.yaml` contains ordered `WRITER-INITIATED`, `PKA-DONE`, `WRITER-DONE`, `WRITER-REWORK`, and cycle-2 `WRITER-DONE` entries with ISO 8601 timestamps. Before this CA update, the most recent entry was `WRITER-DONE` from `wa_w6_edm_v2`.
- CHK-5 PASS: the prompt's literal `rg -l -E "<process keyword>" ...` form errors in this ripgrep because `-E` is parsed as an encoding flag. The equivalent focused command `rg -l -e "Weinberg|three.?gluon|neutron.?EDM|CP.?odd.*gluon|electric dipole moment" quarkConstraints/ qcd/ flavorConstraints/ neutrinos/ yukawa/ warpConfig/ solvers/ scanParams/ tests/ 2>&1` returned only `qcd/alphaS.ipynb`, whose matching line is a prose note about a two-loop neutron EDM calculation, not a live constraint implementation. A non-notebook implementation grep for `dipole|muToEGamma|check_mu_to_e_gamma|EDM` found only the unrelated LFV dipole path; the TeX-cited nearby lines `flavorConstraints/muToEGamma.py:3`, `:21`, and `:81` exist.
- CHK-6 PASS: `HIGH` is consistent with the rubric. E009 needs new CP-odd colored operators, threshold matching, RG/operator mixing, and a hadronic matrix-element convention; it is not an existing \(\Delta F=2\) SLL/SLR/VLL/VRR/LR1/LR2 use case.
- CHK-7 PASS: E009 is a companion-catalog EDM/EFT entry. The CFW 2008/rc1.1 material is used only as RS-flavor context, and no load-bearing rc1.1 paper number is silently revised.
- CHK-8 PASS: `rg -n "\\\\(cite|ref|label)|TODO|CHECK|\\?\\?|undefined|Missing|textbf\\{CHECK\\}" flavor_catalog/processes/edm_neutrino/E009.tex` returned no unresolved-reference or CHECK markers. Math-mode usage in the displayed values and operator notation is cosmetically acceptable.
