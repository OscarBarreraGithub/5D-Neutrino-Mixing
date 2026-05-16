# CA Worklog: ca_wave1_top_higgs_ew
**Date**: 2026-05-16
**Family**: top_higgs_ew
**Process IDs**: T001 T010

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| T001 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| T010 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| T001 | B(t -> Z c), left-handed benchmark | < 1.3e-4 | PDG 2025 top-quark listing | a9f404618edc |
| T001 | B(t -> Z c), right-handed benchmark | < 1.2e-4 | PDG 2025 top-quark listing | a9f404618edc |
| T001 | B(t -> Z c), CMS comparison limit | < 0.049% | CMS 2017 / arXiv:1702.01404 | 4457f9b0ce55 |
| T001 | ATLAS dataset | 139 fb^-1 at sqrt(s) = 13 TeV | ATLAS 2023 / arXiv:2301.11605 | f7c27db54a6e |
| T001 | CMS dataset | 19.7 fb^-1 at sqrt(s) = 8 TeV | CMS 2017 / arXiv:1702.01404 | 4457f9b0ce55 |
| T001 | SM B(t -> c Z) estimate | about 1e-14 | Aguilar-Saavedra 2004 / PDG 2025 | 104a32f81cfd |
| T001 | HL-LHC projected equivalent B(t -> c Z) reach | about 1e-6 | Airen-Franceschini 2026 / arXiv:2601.14966 | e17cbee2269b |
| T010 | R_b^0 | 0.21629 +/- 0.00066 | PDG 2025 Z-boson listing | 0d0b9e2377ce |
| T010 | A_FB^{0,b} | 0.0992 +/- 0.0016; source listing 9.92 +/- 0.16% | PDG 2025 Z-boson listing | 0d0b9e2377ce |
| T010 | A_b | 0.923 +/- 0.020 | PDG 2025 Z-boson listing | 0d0b9e2377ce |
| T010 | A_FB^{0,b} LEP/SLC pull | 2.8 sigma | LEP/SLC final Z-resonance combination | fe6c50814880 |
| T010 | FCC-ee projected relative uncertainty for R_b and A_FB^b | order 0.01% | Roehrig et al. 2025 / arXiv:2502.17281 | 035f06739e5f |

## Issues (if any)
- T001: CHK-1 fails under the batch instruction that every numerical claim in the TeX trace to a `pdg_or_equivalent` entry with year, value, uncertainty, source URL, access date, and sha256. The headline PDG left/right limits satisfy this, but TeX lines 50-60 and 64 cite contextual numerical claims (19.7 fb^-1 at 8 TeV, 139 fb^-1 at 13 TeV, 0.049%, 10^-4 level, 10^-6 projection, and 10^-14 SM estimate) that are only in `auxiliary_values`; those auxiliary entries also lack explicit `year` fields.
- T010: CHK-1 fails under the same strict placement rule. The current PDG values for R_b^0, A_FB^{0,b}, and A_b are in `pdg_or_equivalent`, but TeX lines 28-29 and 50 cite the 2.8 sigma LEP/SLC pull and order-0.01% FCC-ee projection, which are in `additional_numerical_context` rather than `pdg_or_equivalent`.
- T001: CHK-2 PASS. Key references `PDG2025TopFCNCTZQ`, `ATLAS2023TopFCNCTZQ`, `CMS2017TZQFCNC`, `AguilarSaavedra2004TopFCNC`, `CsakiFalkowskiWeiler2008WarpedFlavor`, and `AirenFranceschini2026RSTopFCNC` resolve to source-manifest entries; each manifest entry has a non-empty tracked text snapshot. All source URLs returned HTTP 200, and computed snapshot sha256 values match the sidecar/source manifest.
- T010: CHK-2 PASS. Key references `pdg_2025_z_boson`, `lepslc_2006_z_resonance`, `cfw_2008_rs_flavor`, `casagrande_2008_rs_ewpt`, `freitas_2014_z_widths`, and `fcc_ee_2025_zbb_projection` resolve to source-manifest entries; each manifest entry has a non-empty tracked text snapshot. All source URLs returned HTTP 200, and computed snapshot sha256 values match the sidecar/source manifest.
- T001/T010: CHK-3 PASS. `ls flavor_catalog/references/T001/` and `ls flavor_catalog/references/T010/` show only text snapshots, `sha256sums.txt` for T001, and source manifests; `find ... -name '*.pdf'` returned no PDF files.
- T001/T010: CHK-4 PASS on pre-CA state. Each YAML sidecar had `WRITER-INITIATED` followed by `WRITER-DONE` with ISO 8601 timestamps, and the most recent pre-checker state was `WRITER-DONE`.
- T001: CHK-5 PASS. The installed `rg` treats `-E` as an encoding flag, so the semantic command was run with `-e`. Exact keywords `t.?->.?c.?Z|t.?->.?Z.?c|tcZ|tZc` have no non-notebook hits in the required directories; broader FCNC text only finds generic prose/comments in `quarkConstraints/PAPER_0710_1869.md:35` and `tests/test_paper_couplings.py:258`.
- T010: CHK-5 PASS. Exact non-notebook keywords `\bR_b\b|\bA_FB\b|\bAfb\b|\bA_b\b|\bZbb\b|Z.?->.?b` have no hits in the required directories. A broader Z search finds only generic `M_Z` support, including `qcd/running.py:3`, `qcd/constants.py:11`, and `quarkConstraints/qcd_running.py:100`.
- T001/T010: CHK-6 PASS. `HIGH` is consistent with the rubric: T001 needs new top electroweak FCNC matching/width interpretation, and T010 needs a new electroweak-pole likelihood plus RS matching for Z b_L and Z b_R coupling shifts.
- T001/T010: CHK-7 PASS. Searching outside `flavor_catalog/` found only plan rows or unrelated `PR_benchmark` text, not rc1.1 paper claims contradicted by these catalog entries.
- T001/T010: CHK-8 PASS. No `\textbf{CHECK}`, `\cite`, `\ref`, `TODO`, or unresolved-reference markers were found; `git diff --check` on the two TeX files was clean.
