# CA Worklog: ca_w5b_top_higgs_ew_v2
**Date**: 2026-05-16
**Family**: top_higgs_ew
**Process IDs**: T006

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| T006 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| T006 | ATLAS dataset for \(u g \to t\) | \(139~\mathrm{fb}^{-1}\) at \(13~\mathrm{TeV}\) | PDG 2025 Q007TUG / ATLAS 2022 | 137350422325 |
| T006 | \(\sigma(ug\to t)\mathcal{B}(t\to bW)\mathcal{B}(W\to\ell\nu)\) | \(<3.0~\mathrm{pb}\), 95% CL | PDG 2025 Q007TUG / ATLAS 2022 | 137350422325 |
| T006 | \(\mathcal{B}(W\to\ell\nu)\), ATLAS lepton sum | 0.325 | ATLAS 2022 arXiv:2112.01302 | ed435aaa49e3 |
| T006 | \(|C^{ut}_{uG}|/\Lambda^2\) | \(<0.057~\mathrm{TeV}^{-2}\), 95% CL | PDG 2025 Q007TUG / ATLAS 2022 | 137350422325 |
| T006 | \(\mathcal{B}(t\to ug)\), ATLAS EFT interpretation | \(<6.1\times10^{-5}\), 95% CL | PDG 2025 Q007TUG / ATLAS 2022 | 137350422325 |
| T006 | CMS dataset for \(t\)-channel single-top FCNC search | \(5.0~\mathrm{fb}^{-1}\) at \(7~\mathrm{TeV}\) and \(19.7~\mathrm{fb}^{-1}\) at \(8~\mathrm{TeV}\) | CMS 2017 arXiv:1610.03545 | ce1b34051294 |
| T006 | \(|\kappa_{tug}|/\Lambda\) | \(<4.1\times10^{-3}~\mathrm{TeV}^{-1}\), 95% CL | CMS 2017 arXiv:1610.03545 | ce1b34051294 |
| T006 | \(\mathcal{B}(t\to ug)\), CMS interpretation | \(<2.0\times10^{-5}\), 95% CL | CMS 2017 arXiv:1610.03545 | ce1b34051294 |
| T006 | SM context estimate for \(\mathcal{B}(t\to ug)\) | \(\simeq3.6\times10^{-14}\) | Aguilar-Saavedra 2004 arXiv:hep-ph/0409342 | ed810970e049 |

## Evidence notes
- CHK-1: Grepped numeric claims in `T006.tex`; each load-bearing value cited in the PDG/equivalent section is present in `T006.yaml` under `pdg_or_equivalent.values` with year, value, uncertainty field, source URL, access date, and sha256. The WA-v2 fix correctly promoted `AguilarSaavedra2004:T006:t_ug_SM` into `pdg_or_equivalent.values`.
- CHK-2: The TeX key references `PDG2025TopFCNCTUG`, `ATLAS2022TopFCNCTUG`, `CMS2017TopFCNCTUG`, `AguilarSaavedra2004TopFCNC`, and `CsakiFalkowskiWeiler2008WarpedFlavor` all resolve to entries in `flavor_catalog/references/T006/source_manifest.yaml`. Each entry has a non-empty `snapshot_path` under `flavor_catalog/references/T006/`, and `git ls-files flavor_catalog/references/T006` shows the snapshots are tracked.
- CHK-3: `ls flavor_catalog/references/T006/` shows only `.txt` snapshots plus `source_manifest.yaml`; no PDF files are tracked.
- CHK-4: `status_history` contains `WRITER-INITIATED` then `WRITER-DONE`, then the previous cycle `WRITER-REWORK`, then WA-v2 `WRITER-DONE`, all with ISO 8601 timestamps. This CA appended `CHECKER-DONE` as the most recent entry.
- CHK-5: The user-specified `rg -l -E "<process keyword>" ...` form is incompatible with installed ripgrep because `-E` is parsed as encoding. The equivalent `rg -l -e "t.?->.?u.?g|t.?->.?g.?u|tug|utg|u.?g.?->.?t|tqg|uG" ... -g "!*.ipynb"` returned no implementation files. A broader FCNC grep returned only generic prose/comment hits at `quarkConstraints/PAPER_0710_1869.md:35` and `tests/test_paper_couplings.py:258`. `quarkConstraints/modern/phenomenology.py:23` enumerates only `epsilon_K`, `K`, `B_d`, `B_s`, and `D0`; line 167 states the surface has no numeric EFT backend.
- CHK-6: `HIGH` is consistent with the rubric: a live constraint would need a top chromomagnetic FCNC operator convention, RS-to-SMEFT matching, \(t\to ug\) width and \(ug\to t\) production normalization, collider/PDF reinterpretation, and QCD-running choices outside the existing \(\Delta F=2\) operator basis.
- CHK-7: No T006 load-bearing values are in rc1.1. The only CFW context is explicitly paper-era baseline context and does not revise the rc1.1 quark-scan numbers.
- CHK-8: No `\textbf{CHECK}`, TODO, `\cite`, or `\ref` tokens were found in `T006.tex`; math mode and escaped underscores are acceptable on inspection.

## Issues (if any)
- None.
