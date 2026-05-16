# CA Worklog: ca_w5b_top_higgs_ew
**Date**: 2026-05-16
**Family**: top_higgs_ew
**Process IDs**: T006, T016, T017

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| T006 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| T016 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| T017 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| T006 | ATLAS dataset for ug -> t | 139 fb^-1 at sqrt(s) = 13 TeV | PDG 2025 / ATLAS 2022 | 137350422325 |
| T006 | sigma(ug -> t) B(t -> bW) B(W -> l nu) | < 3.0 pb, 95% CL | PDG 2025 / ATLAS 2022 | 137350422325 |
| T006 | B(W -> l nu) leptonic sum | 0.325 | ATLAS 2022 | ed435aaa49e3 |
| T006 | \|C_uG^{ut}\| / Lambda^2 | < 0.057 TeV^-2, 95% CL | PDG 2025 / ATLAS 2022 | 137350422325 |
| T006 | B(t -> u g), ATLAS EFT interpretation | < 6.1e-5, 95% CL | PDG 2025 / ATLAS 2022 | 137350422325 |
| T006 | CMS dataset | 5.0 fb^-1 at 7 TeV and 19.7 fb^-1 at 8 TeV | CMS 2017 | ce1b34051294 |
| T006 | \|kappa_tug\| / Lambda | < 4.1e-3 TeV^-1, 95% CL | CMS 2017 | ce1b34051294 |
| T006 | B(t -> u g), CMS interpretation | < 2.0e-5, 95% CL | CMS 2017 | ce1b34051294 |
| T006 | SM B(t -> u g) context | approximately 3.6e-14 | Aguilar-Saavedra 2004 | ed810970e049 |
| T016 | B(Z -> e tau), PDG/ATLAS combined | < 5.0e-6, 95% CL | PDG 2025 / ATLAS 2021 | d6c641885851 |
| T016 | B(Z -> e tau), ATLAS leptonic tau only | < 7.0e-6, 95% CL | ATLAS 2021 | 0869c301a050 |
| T016 | CMS dataset for Z -> e tau | 138 fb^-1 at sqrt(s) = 13 TeV | CMS 2025 | 9370bd0fcb6e |
| T016 | B(Z -> e tau), CMS observed expected | < 13.8e-6; expected < 11.4e-6, 95% CL | CMS 2025 | 9370bd0fcb6e |
| T016 | Tera-Z scale for LFV Z studies | O(1e12) Z decays | Calibbi-Marcano-Roy 2021 | b6bc6608c1cd |
| T017 | B(Z -> mu tau), PDG/ATLAS combined | < 6.5e-6, 95% CL | PDG 2025 / ATLAS 2021 | 38a914a87c02 |
| T017 | ATLAS dataset for Z -> mu tau | 139 fb^-1 at sqrt(s) = 13 TeV | ATLAS 2021 | cd52f80f1ce5 |
| T017 | B(Z -> mu tau), ATLAS leptonic tau only | < 7.2e-6, 95% CL | ATLAS 2021 | cd52f80f1ce5 |
| T017 | CMS dataset for Z -> mu tau | 138 fb^-1 at sqrt(s) = 13 TeV | CMS 2025 | 0dee8d38e4f2 |
| T017 | B(Z -> mu tau), CMS observed expected | < 12.0e-6; expected < 5.3e-6, 95% CL | CMS 2025 | 0dee8d38e4f2 |
| T017 | Tera-Z scale for LFV Z studies | O(10^12) Z decays | Calibbi-Marcano-Roy 2021 | d9a926db4f51 |

## Evidence notes
- CHK-1: Grepped all three TeX files for numerical claims. T016 and T017 physics numbers map to `pdg_or_equivalent.values` with year, value, uncertainty/null uncertainty, source URL, access date, snapshot path, and sha256. T006's measured ATLAS/CMS values satisfy this, but the TeX also quotes the SM context estimate `B(t -> u g) ~= 3.6e-14`; that value is verified in `theory_context`, not under `pdg_or_equivalent.values`, so T006 fails the strict placement rule.
- CHK-2: Source-key references in the TeX resolve to `flavor_catalog/references/<process_id>/source_manifest.yaml`. Each manifest entry has a non-empty process-local `snapshot_path`, and `git ls-files` confirms the referenced `.txt` snapshots are tracked.
- CHK-3: `ls flavor_catalog/references/T006`, `T016`, and `T017` shows only `.txt` snapshots, optional `sha256sums.txt`, and `source_manifest.yaml`; `find ... -name '*.pdf'` returned no PDFs.
- CHK-4: Before CA edits, all three sidecars contained `WRITER-INITIATED` followed by `WRITER-DONE` with ISO 8601 timestamps, and the most recent pre-check state was `WRITER-DONE`.
- CHK-5: The prompt's literal `rg -l -E` form fails in this environment because this ripgrep uses `-E` as the encoding flag. The semantic searches were rerun with `rg -l -e`. Non-notebook process-specific searches for `tug`/`utg`/`tqg`/`C_uG`, `Z -> e tau`/`Zetau`, and `Z -> mu tau`/`Zmutau` found no observable implementations. The generic evidence cited in the TeX exists at `quarkConstraints/PAPER_0710_1869.md:35`, `tests/test_paper_couplings.py:258`, `flavorConstraints/muToEGamma.py:75`, `scanParams/scan.py:524`, and `quarkConstraints/modern/phenomenology.py:23` and `:167`.
- CHK-6: `HIGH` is consistent for all three processes. T006 needs a new top chromomagnetic FCNC convention, RS-to-SMEFT matching, production/width normalization, and recast/RG choices. T016 and T017 need LFV electroweak `Z l_i l_j` operator conventions, lepton-sector matching, and a combination policy with tau LFV.
- CHK-7: No contradictions with rc1.1 load-bearing numbers were found. Searches outside `flavor_catalog/` found only plan rows or unrelated examples, and the checked values are companion catalog additions.
- CHK-8: No `CHECK`, `TODO`, `\ref{}`, `\cite{}`, `\label{}`, or unresolved-reference markers were found in the three TeX files; `git diff --check` on the relevant TeX paths was clean.

## Issues (if any)
- T006: CHK-1 fails under the strict CA placement rule. Move or duplicate the SM `B(t -> u g) ~= 3.6e-14` context value into `pdg_or_equivalent.values` with explicit `year`, `value`, `uncertainty`, `source_url`, `access_date`, `snapshot_path`, and `sha256`, or remove that numerical claim from the TeX in the next WA pass.
