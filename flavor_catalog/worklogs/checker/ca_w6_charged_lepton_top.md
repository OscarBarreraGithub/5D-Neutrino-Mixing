# CA Worklog: ca_w6_charged_lepton_top
**Date**: 2026-05-16
**Family**: charged_lepton_top
**Process IDs**: L023 T020

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| L023 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| T020 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| L023 | CHARM-II \(\sigma_{\rm exp}/\sigma_{\rm SM}\) | \(1.58\pm0.64\) | DUNE trident compilation, arXiv:1902.06765 | 19acfb9bc4c3 |
| L023 | CCFR \(\sigma_{\rm exp}/\sigma_{\rm SM}\) | \(0.82\pm0.28\) | CCFR PRL/PubMed snapshot plus DUNE compilation | 903c8a531035 |
| L023 | CCFR corrected trident events vs SM | \(37.0\pm12.4\) vs \(45.3\pm2.3\) | CCFR PRL/PubMed snapshot | 903c8a531035 |
| L023 | NuTeV \(\sigma_{\rm exp}/\sigma_{\rm SM}\) | \(0.72^{+1.73}_{-0.72}\) | DUNE compilation and NuTeV arXiv snapshot | 19acfb9bc4c3 |
| L023 | Belle-II auxiliary reach inputs | \(\sqrt{s}=10.58\,\mathrm{GeV}\), \(50\,\mathrm{ab}^{-1}\) | Kaneta-Shimomura arXiv:1701.00156 | 9f1614fe9aa0 |
| L023 | DUNE projection wording | about \(25\%\) after about 3 years per beam mode; \(40\%\) baseline and \(25\%\) improved contour | DUNE trident compilation, arXiv:1902.06765 | 19acfb9bc4c3 |
| T020 | PDG/CMS \(\mathcal{B}(h\to e\mu)\) | \(<4.4\times10^{-5}\), expected \(<4.7\times10^{-5}\), 95% CL | PDG 2025 Higgs LFV review | 8178cb0c1055 |
| T020 | PDG/ATLAS \(\mathcal{B}(h\to e\mu)\) | \(<6.2\times10^{-5}\), expected \(<5.9\times10^{-5}\), 95% CL | PDG 2025 Higgs LFV review | 8178cb0c1055 |
| T020 | ATLAS primary \(\mathcal{B}(H\to e\mu)\) | \(<6.1\times10^{-5}\), expected \(<5.8\times10^{-5}\); \(139\,\mathrm{fb}^{-1}\) at 13 TeV | ATLAS arXiv:1909.10235 snapshot | e1e4c131be4c |
| T020 | CMS direct search dataset and scan | \(138\,\mathrm{fb}^{-1}\) at 13 TeV; 110--160 GeV scan; largest excess near 146 GeV with 3.8/2.8 sigma | CMS HIG-22-002 public page | 072412def18c |
| T020 | Higgs mass hypothesis in direct searches | \(125\,\mathrm{GeV}\) | PDG/CMS/ATLAS snapshots | 8178cb0c1055 |

## Issues (if any)
- L023: CHK-1 fails under the strict batch rule. `L023.tex` lines 50, 53, and 57--59 contain numerical claims for Altmannshofer 2014 95% CL contours, Belle-II \(\sqrt{s}=10.58\,\mathrm{GeV}\) and \(50\,\mathrm{ab}^{-1}\), and DUNE \(25\%\)/\(40\%\)/\(25\%\) projection wording. These values are traceable in local snapshots and YAML `auxiliary_context`, but not in YAML `pdg_or_equivalent` entries with the required year/value/source URL/access-date/sha256 metadata. Rework should duplicate those retained TeX numbers into `pdg_or_equivalent` entries or remove the numerics from the TeX prose.
- T020: CHK-1 fails under the strict batch rule. `T020.tex` line 14 quotes a \(125\,\mathrm{GeV}\) Higgs mass hypothesis. The number appears in local snapshots and `pdg_or_equivalent.value_summary`, but not as a `pdg_or_equivalent.values[]` entry with value metadata. Rework should add a per-value entry for the 125 GeV mass hypothesis or remove the numeric claim from the TeX prose.

## Checklist evidence
- CHK-2: All process-local source keys cited in `L023.tex` and `T020.tex` appear in the corresponding `source_manifest.yaml`. Each manifest entry has a non-empty `snapshot_path`, and `git ls-files` confirms every snapshot path is tracked. Value-ID anchors such as `PDG2025:T020:cms_run2` and `CMS2023:T020:mass_scan` resolve to YAML `pdg_or_equivalent.values[]`, not manifest source keys.
- CHK-3: `ls flavor_catalog/references/L023/` and `ls flavor_catalog/references/T020/` showed only text snapshots, `source_manifest.yaml`, and for T020 `sha256sums.txt`. `find ... -name '*.pdf'` returned no PDFs.
- CHK-4: Both sidecars had ordered ISO-8601 `WRITER-INITIATED` and terminal `WRITER-DONE` entries before this CA update. L023 also retained its intervening `PKA-DONE` entry.
- CHK-5: The literal requested `rg -l -E "<process keyword>" ...` form errors in this ripgrep because `-E` is parsed as an encoding flag. Equivalent `rg -l -e` focused searches over the required directories returned no L023 trident implementation and no T020 \(h\to e\mu\)/`Y_e_mu`/`Y_mu_e` implementation. T020's cited nearby lines exist at `flavorConstraints/muToEGamma.py:1`, `flavorConstraints/muToEGamma.py:75`, `scanParams/scan.py:523`, `scanParams/scan.py:524`, `yukawa/charged_lepton.py:1`, and `yukawa/charged_lepton.py:36`.
- CHK-6: HIGH is consistent for both processes. L023 requires a new neutrino-trident mode calculation with lepton-current operators, target/form-factor handling, flux integration, and experiment-specific rates. T020 requires new charged-lepton Higgs-Yukawa misalignment matching plus Higgs production/width interpretation. Neither is a low-risk update to the existing \(\Delta F=2\) basis or the existing \(\mu\to e\gamma\) path.
- CHK-7: Searches outside `flavor_catalog/` found only plan rows for L023 and T020, not rc1.1 load-bearing values to contradict.
- CHK-8: `rg` found no `CHECK`, `TODO`, `FIXME`, `??`, `\cite`, or `\ref` markers in either TeX file, and no obvious math-mode or bracing defects were observed.
