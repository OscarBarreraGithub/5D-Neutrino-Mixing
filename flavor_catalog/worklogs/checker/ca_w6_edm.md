# CA Worklog: ca_w6_edm
**Date**: 2026-05-16
**Family**: edm_neutrino
**Process IDs**: E009

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| E009 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| E009 | \(|d_n|\) direct anchor | \(<1.8\times10^{-26}\ e\,\mathrm{cm}\), 90% CL | PDG Live 2026 S017EDM | ea22e4f3da86 |
| E009 | \(d_n\) primary measurement | \((0.0\pm1.1_{\rm stat}\pm0.2_{\rm sys})\times10^{-26}\ e\,\mathrm{cm}\) | Abel et al. 2020, arXiv:2001.11966 | a78b99dfd396 |
| E009 | Pospelov-Ritz Weinberg normalization | \(|d_n(w)|\simeq e\,22~\mathrm{MeV}\,w(1~\mathrm{GeV})\) | Pospelov and Ritz 2005, arXiv:hep-ph/0504231 | e07e62e5a797 |
| E009 | Derived \(w(1~\mathrm{GeV})\) benchmark | \(<4.1\times10^{-11}\ \mathrm{GeV}^{-2}\), central one-source translation | PDG 2026 neutron anchor plus Pospelov-Ritz 2005 | ea22e4f3da86 / e07e62e5a797 |
| E009 | Haisch-Hala \(O_6\) normalization | \((d_n/e)_{O_6}=74(1\pm0.5)~\mathrm{MeV}\) | Haisch and Hala 2019, arXiv:1909.08955 | ae1c67f2535e |
| E009 | Derived \(C_6\) benchmark | \(<1.2\times10^{-11}\ \mathrm{GeV}^{-2}\), central one-source translation | PDG 2026 neutron anchor plus Haisch-Hala 2019 | ea22e4f3da86 / ae1c67f2535e |

## Check evidence
- CHK-1 failed under the strict sidecar-placement rule. The measured neutron EDM limit and Abel measurement in `E009.tex` lines 22-28 are in `pdg_or_equivalent` with year, value, uncertainty/CL, source URL, access date, snapshot path, and sha256. However, `E009.tex` lines 53-60 also quote the Pospelov-Ritz \(22~\mathrm{MeV}\), \(1~\mathrm{GeV}\), \(4.1\times10^{-11}\ \mathrm{GeV}^{-2}\), Haisch-Hala \(74(1\pm0.5)~\mathrm{MeV}\), and \(1.2\times10^{-11}\ \mathrm{GeV}^{-2}\) benchmark numbers. Those values trace to `paper_era_reference` or `auxiliary_theory_inputs`, not to `pdg_or_equivalent` entries as required by CHK-1.
- CHK-2 passed. Every process-local source key cited in the TeX resolves in `flavor_catalog/references/E009/source_manifest.yaml`, and each manifest entry has a non-empty `snapshot_path`. `git ls-files flavor_catalog/references/E009` confirms all snapshots are tracked under the process reference directory.
- CHK-3 passed. `ls flavor_catalog/references/E009/` shows only `.txt` snapshots plus `source_manifest.yaml`; `find flavor_catalog/references/E009 -iname '*.pdf'` returned no PDF files.
- CHK-4 passed for the pre-check state. `E009.yaml` contained ordered `WRITER-INITIATED`, `PKA-DONE`, and terminal `WRITER-DONE` entries with ISO 8601 timestamps before this CA update.
- CHK-5 passed. The prompt's literal `rg -l -E "<process keyword>" ...` form fails in this ripgrep because `-E` is an encoding flag, so I reran the intended search with `rg -l -e "Weinberg|three.?gluon|neutron.?EDM|CP.?odd.*gluon"` over the required implementation/test directories. It returned no files. A broader non-notebook `EDM|electric dipole|dipole|muToEGamma` search found only unrelated `mu->e gamma` dipole code; the TeX-cited nearby lines `flavorConstraints/muToEGamma.py:3`, `:21`, and `:81` exist.
- CHK-6 passed. `HIGH` is consistent because E009 needs new CP-odd colored operators, threshold/RG mixing, and a hadronic matrix-element convention; it is not covered by the existing Delta F=2 basis or the muon LFV dipole path.
- CHK-7 passed. E009 is a companion-catalog EDM entry and does not revise a load-bearing rc1.1 paper number.
- CHK-8 passed. No `\textbf{CHECK}`, `\cite`, `\ref`, `TODO`, `FIXME`, `??`, or simple brace/math delimiter imbalance was found; `chktex`/`lacheck` are not installed in this environment.

## Issues (if any)
- E009: CHK-1 rework required. Either move the TeX lines 53-60 benchmark numeric claims into value-bearing `pdg_or_equivalent` entries with the required metadata, or soften/remove those numerical benchmark claims from the TeX in the next WA cycle. I did not edit the TeX per CA hard rules.
