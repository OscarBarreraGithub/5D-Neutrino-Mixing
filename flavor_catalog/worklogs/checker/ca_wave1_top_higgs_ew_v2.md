# CA Worklog: ca_wave1_top_higgs_ew_v2
**Date**: 2026-05-16
**Family**: top_higgs_ew
**Process IDs**: T001 T010

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| T001 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| T010 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| T001 | B(t -> Z c), left-handed benchmark | < 1.3e-4, 95% CL | PDG 2025 top-quark listing | a9f404618edc |
| T001 | B(t -> Z c), right-handed benchmark | < 1.2e-4, 95% CL | PDG 2025 top-quark listing | a9f404618edc |
| T001 | B(t -> Z c), CMS comparison limit | < 0.049%, 95% CL | CMS 2017 / arXiv:1702.01404 | 4457f9b0ce55 |
| T001 | CMS dataset | 19.7 fb^-1 at sqrt(s) = 8 TeV | CMS 2017 / arXiv:1702.01404 | 4457f9b0ce55 |
| T001 | ATLAS dataset | 139 fb^-1 at sqrt(s) = 13 TeV | ATLAS 2023 / arXiv:2301.11605 | f7c27db54a6e |
| T001 | ATLAS t -> Z c reach | 10^-4 level, represented by < 13e-5 and < 12e-5 benchmarks | ATLAS 2023 / arXiv:2301.11605 | f7c27db54a6e |
| T001 | SM B(t -> c Z) estimate | about 1e-14 | Aguilar-Saavedra 2004 and PDG 2025 | 104a32f81cfd |
| T001 | HL-LHC projected equivalent B(t -> c Z) reach | about 1e-6 | Airen-Franceschini 2026 / arXiv:2601.14966 | e17cbee2269b |
| T010 | R_b^0 | 0.21629 +/- 0.00066 | PDG 2025 Z-boson listing | 0d0b9e2377ce |
| T010 | A_FB^{0,b} | 0.0992 +/- 0.0016; source listing 9.92 +/- 0.16% | PDG 2025 Z-boson listing | 0d0b9e2377ce |
| T010 | A_b | 0.923 +/- 0.020 | PDG 2025 Z-boson listing | 0d0b9e2377ce |
| T010 | LEP/SLC final-combination pull for A_FB^{0,b} | 2.8 sigma | LEP/SLC final Z-resonance combination | fe6c50814880 |
| T010 | FCC-ee projected relative uncertainty for R_b and A_FB^b | order 0.01% | Roehrig et al. 2025 / arXiv:2502.17281 | 035f06739e5f |

## Issues (if any)
- None.

## Evidence notes
- T001 CHK-1: all load-bearing TeX numerical claims now trace to `pdg_or_equivalent.values` entries with year, value, uncertainty, source URL, access date, and sha256: PDG left/right limits, CMS comparison limit and dataset, ATLAS dataset and benchmarks, SM estimate, and HL-LHC projection.
- T010 CHK-1: the PDG Z-pole values, LEP/SLC `2.8 sigma` pull, and FCC-ee order-`0.01%` projection are all in `pdg_or_equivalent` with year, value, uncertainty, source URL, access date, snapshot path, and sha256.
- CHK-2: bibliographic keys and snapshot references used in the TeX resolve to `source_manifest.yaml` entries for both processes. Each manifest entry has a non-empty `snapshot_path`, and `git ls-files` confirms the referenced text snapshots are tracked under the process reference directory.
- CHK-3: `ls flavor_catalog/references/T001/`, `ls flavor_catalog/references/T010/`, and `find ... -name '*.pdf'` show no publisher PDFs; the reference directories contain text snapshots, manifests, and T001's `sha256sums.txt`.
- CHK-4: both sidecars contain `WRITER-INITIATED` followed by the original `WRITER-DONE`, the CA `WRITER-REWORK`, and the WA-v2 `WRITER-DONE`, all with ISO 8601 timestamps. The most recent pre-check state was `WRITER-DONE`.
- CHK-5: the requested `rg -l -E` form was run and fails because this ripgrep treats `-E` as an encoding flag, so the semantic searches were rerun with `-e`. For T001, exact `t -> cZ` / `tZc` / `tcZ` / `Ztc` non-notebook searches have no hits; broader FCNC hits are only generic prose/comments. For T010, exact `R_b` / `A_FB` / `A_b` / `Zbb` / `Z -> b` non-notebook search finds only `scanParams/THEORY_PRIORS.md` via the unrelated `PR_benchmark` substring, while broader Z searches find generic `M_Z` QCD support, not a Zbb likelihood.
- CHK-6: `HIGH` is consistent for both processes. T001 needs a new electroweak top-FCNC matching and width or production interpretation, and T010 needs a new electroweak-pole likelihood plus RS matching for `Z b_L` and `Z b_R` coupling shifts.
- CHK-7: searching outside `flavor_catalog/` found only plan rows or unrelated `PR_benchmark` text; no rc1.1 load-bearing paper number is contradicted.
- CHK-8: no `\textbf{CHECK}`, `\cite`, `\ref`, `TODO`, or unresolved-reference markers were found; `git diff --check` on the touched YAML and relevant TeX paths was clean.
