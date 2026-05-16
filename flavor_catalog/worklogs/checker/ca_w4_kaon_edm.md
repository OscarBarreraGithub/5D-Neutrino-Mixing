# CA Worklog: ca_w4_kaon_edm
**Date**: 2026-05-16
**Family**: kaon_edm
**Process IDs**: K017 E004 E006 E008

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| K017 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| E004 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| E006 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| E008 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| K017 | \(R_K\) PDG average | \((2.488\pm0.009)\times10^{-5}\) | PDG API S010R28 / PDG2025:K017RK | c00ed1ff90dd |
| K017 | \(R_K\) NA62 input | \((2.488\pm0.007_{\rm stat}\pm0.007_{\rm syst})\times10^{-5}\) | NA62 2013 | b6389e671653 |
| K017 | \(R_K\) KLOE input | \((2.493\pm0.025_{\rm stat}\pm0.019_{\rm syst})\times10^{-5}\) | KLOE 2009 | c473017a6d6d |
| K017 | \(R_K^{\rm SM}\) | \((2.477\pm0.001)\times10^{-5}\) | Cirigliano-Rosell 2007 | 9178029253ec |
| E004 | \(|d_n|\) direct limit | \(<1.8\times10^{-26}\ e\,\mathrm{cm}\) at 90% CL | PDG Live 2026 S017EDM / Abel 2020 | 211affe27bf5 |
| E004 | \(d_n\) primary measurement | \((0.0\pm1.1_{\rm stat}\pm0.2_{\rm sys})\times10^{-26}\ e\,\mathrm{cm}\) | Abel et al. 2020 | b413e25a28af |
| E004 | \(|d_n|\) superseded direct limit | \(<3.0\times10^{-26}\ e\,\mathrm{cm}\) at 90% CL | PDG Live 2026 S017EDM / Pendlebury 2015 | 211affe27bf5 |
| E006 | \(|d_{\mathrm{Hg}}|\) direct limit | \(<7.4\times10^{-30}\ e\,\mathrm{cm}\) at 95% CL | Graner et al. 2016 | 44522bf9473b |
| E006 | PDG Hg-derived neutron rows | Graner \(<0.16\), Sahoo \(<0.22\) in \(10^{-25}\ e\,\mathrm{cm}\), both 95% CL | PDG Live 2026 S017EDM cross-reference | 15f7baddbe77 |
| E006 | Sahoo Hg-derived translation limits | \(d_n<2.2\times10^{-26}\), \(d_p<2.1\times10^{-25}\), \(|\bar\theta|<1.1\times10^{-10}\), \(|\tilde d_u-\tilde d_d|<5.5\times10^{-27}\) | Sahoo 2017 | 88e10fc1c214 |
| E006 | Hg post-2008 improvement | factor of 4 | Graner et al. 2016 snapshot | 44522bf9473b |
| E008 | \(|d_n|\) anchor | \(<1.8\times10^{-26}\ e\,\mathrm{cm}\) at 90% CL | PDG Live 2026 S017EDM | 2478ea90f896 |
| E008 | \(d_{\rm Hg}\) anchor | \((2.20\pm2.75_{\rm stat}\pm1.48_{\rm syst})\times10^{-30}\ e\,\mathrm{cm}\); \(|d_{\rm Hg}|<7.4\times10^{-30}\ e\,\mathrm{cm}\) at 95% CL | Graner et al. 2016 | 8c99b7a6bd9c |
| E008 | neutron qCEDM benchmark | \(|\tilde d_d+0.5\tilde d_u|<1.6\times10^{-26}\ \mathrm{cm}\), with \((1\pm0.5)\) normalization | Pospelov-Ritz 2000 plus current neutron anchor | 9ca4f0955eb5 |
| E008 | Hg qCEDM benchmark | \(|\tilde d_u-\tilde d_d|<1.1\times10^{-27}\ \mathrm{cm}\) from \(7\times10^{-3}e(\tilde d_u-\tilde d_d)\) | Olive-Pospelov-Ritz-Santoso 2005 plus Graner anchor | 6acee076efae |
| E008 | Hg post-2008 improvement | factor of 4 | Graner et al. 2016 | 8c99b7a6bd9c |

## Check Evidence
- CHK-1: Grepped each TeX with `rg -n "[0-9]"` and compared physical numerical claims to sidecar `pdg_or_equivalent` blocks and local snapshots. E004 and E008 pass. K017 and E006 fail the strict placement/structured-metadata rule listed below, although the snapshot values themselves match.
- CHK-2: All TeX source keys resolve to the corresponding `flavor_catalog/references/<process_id>/source_manifest.yaml`. Each manifest entry has a non-empty `snapshot_path`; `test -s` plus `git ls-files --error-unmatch` confirmed every snapshot is tracked under its process reference directory.
- CHK-3: `ls flavor_catalog/references/{K017,E004,E006,E008}/` showed only `.txt` snapshots and `source_manifest.yaml`; `find ... -iname '*.pdf'` returned no PDF files.
- CHK-4: Before this CA update, all four sidecars contained `WRITER-INITIATED` followed by `WRITER-DONE` with ISO 8601 timestamps, and the most recent pre-check state was `WRITER-DONE`.
- CHK-5: The installed `rg` treats `-E` as an encoding flag, so the semantic command was run with `-e`. Focused non-notebook greps found no K017 \(R_K/K_{\ell2}\), E004 neutron-EDM/qCEDM/Weinberg, E006 Hg/Schiff/diamagnetic, or E008 qCEDM implementation. Generic dipole greps found only `flavorConstraints/muToEGamma.py:3`, `:21`, and `:81`, matching the TeX near-miss notes; K017's cited `quarkConstraints/deltaf2.py:209` and `quarkConstraints/modern/phenomenology.py:23` also exist.
- CHK-6: The implementation-difficulty labels are consistent. K017 is MEDIUM because it needs a new charged-current leptonic observable but no new lattice matrix elements. E004, E006, and E008 are HIGH because EDM/qCEDM constraints need new CP-odd matching, RG/threshold treatment, and hadronic/nuclear or atomic inputs.
- CHK-7: No load-bearing rc1.1 number is contradicted. Searches outside `flavor_catalog/` found only existing CFW convention/audit text for 21/33 TeV and plan rows; these entries are companion-catalog EDM/LFU additions and do not silently revise the paper.
- CHK-8: No `\textbf{CHECK}`, TODO/FIXME, `\cite`, `\ref`, unresolved-marker hits, or simple brace imbalance was found in the four TeX files.

## Issues (if any)
- K017: CHK-1 fails under the strict instruction that every numerical claim in the TeX trace to a `pdg_or_equivalent` entry with year, value, uncertainty/CL where applicable, source URL, access date, snapshot path, and sha256. The PDG average, NA62 input, KLOE input, and Cirigliano-Rosell SM prediction blocks contain values, uncertainties, URLs, access dates, snapshot paths, and hashes, but lack explicit `year` fields.
- E006: CHK-1 fails under the same strict placement rule. The direct Hg limit, PDG cross-reference rows, and Sahoo translation limits are represented in `pdg_or_equivalent`, but the TeX post-2008 claim that Graner 2016 improved the previous Hg bound by a factor of \(4\) appears only in the source snapshot and `post_2008_context`, not in a value-bearing `pdg_or_equivalent` entry with full strict metadata.
