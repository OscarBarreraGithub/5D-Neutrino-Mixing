# CA Worklog: ca_w8_B_rare_leptonic
**Date**: 2026-05-17
**Family / batch**: beauty / CA-w8-B-rare-leptonic
**Cycle**: 1
**Process IDs**: B007 B008

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | CHK-W8 | overall |
|---|---|---|---|---|---|---|---|---|---|---|
| B007 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| B008 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| B007 | BR(B_s0 -> e+e-) current limit | < 9.4e-9, 90% CL | PDG live/API 2026 S086.20 | bec216f3aa33 |
| B007 | BR(B0 -> e+e-) current limit | < 2.5e-9, 90% CL | PDG live/API 2026 S042.6 | bec216f3aa33 |
| B007 | BR(B_s0 -> e+e-) LHCb companion limit | < 11.2e-9, 95% CL | LHCb 2020 / arXiv:2003.03999 | dd07d58ad222 |
| B007 | BR(B0 -> e+e-) LHCb companion limit | < 3.0e-9, 95% CL | LHCb 2020 / arXiv:2003.03999 | dd07d58ad222 |
| B007 | BR(B_s0 -> e+e-) historical CDF limit | < 2.8e-7, 90% CL | CDF 2009 / arXiv:0901.3803 | 151fba3738b |
| B007 | BR(B0 -> e+e-) historical CDF limit | < 8.3e-8, 90% CL | CDF 2009 / arXiv:0901.3803 | 151fba3738b |
| B007 | BR(B0 -> e+e-) historical BABAR limit | < 1.13e-7, 90% CL | BABAR 2008 / arXiv:0712.1516 | dbaf89484389 |
| B007 | BR(B0 -> e+e-) historical Belle row | < 1.9e-7, 90% CL | PDG live/API 2026 S042.6 | bec216f3aa33 |
| B007 | B_s0 lifetime-assumption limit shift | +/- 2.4% | LHCb 2020 / arXiv:2003.03999 | dd07d58ad222 |
| B008 | BR(B_s0 -> tau+tau-) current limit | < 6.8e-3, 95% CL | PDG live/API 2026 S086.130 | ef1c4ed4f001 |
| B008 | BR(B0 -> tau+tau-) current limit | < 2.1e-3, 95% CL | PDG live/API 2026 S042.336 | 4f20e368d910 |
| B008 | BR(B0 -> tau+tau-) historical BABAR limit | < 4.1e-3, 90% CL | BABAR 2006 / arXiv:hep-ex/0511015 | d0acf4944c8c |
| B008 | LHCb B_s0/B0 -> tau+tau- dataset metadata | 3 fb^-1 at 7 and 8 TeV | LHCb 2017 / arXiv:1703.02508 | b2b79597fa8a |
| B008 | SM BR(B_s -> tau+tau-) theory input | (7.73 +/- 0.49)e-7 | Bobeth et al. 2014 / arXiv:1311.0903 | dda74abee320 |
| B008 | SM BR(B_d -> tau+tau-) theory input | (2.22 +/- 0.19)e-8 | Bobeth et al. 2014 / arXiv:1311.0903 | dda74abee320 |

## Issues (if any)
- None.

## Evidence notes
- CHK-1: B007's measured experimental limits in TeX match `pdg_or_equivalent.values` entries for the PDG 90% CL limits and the LHCb 95% CL companion limits, with value or upper_limit, uncertainty null for limits, year, units, source URL, access date, snapshot path, and sha256. Historical CDF/BABAR/Belle measured limits are also structured in `pdg_or_equivalent.values`. The CFW KK reference-scale numbers, SM predictions, dataset metadata, and the LHCb lifetime-assumption shift were treated under the L001/B001/B021 carve-out and not required as measured-observable entries.
- CHK-1: B008's measured experimental limits in TeX match `pdg_or_equivalent.values` entries for the PDG/LHCb 95% CL limits, and the historical BABAR measured limit is also structured there. The LHCb 3 fb^-1, 7/8 TeV dataset descriptor and Bobeth SM predictions are kept in supporting/auxiliary fields with provenance, consistent with the carve-out.
- CHK-2: all TeX reference keys resolve to process-local manifest keys. B007 keys checked: `PDG2026_BsBdEe`, `LHCb2020_BsBdEe`, `BobethEtAl2013_BqllSM`, `CDF2009_BsBdEe`, `FleischerJaarsmaTetlalmatziXolocotzi2017_BqllNP`, `BaBar2008_BdEe`, and `CsakiFalkowskiWeiler2008_RSFlavor`. B008 keys checked: `PDG2026_BsTauTau`, `PDG2026_BdTauTau`, `AaijEtAl2017_BqTauTau`, `AubertEtAl2006_BdTauTau`, `LeesEtAl2017_BKTauTauRelated`, `BobethEtAl2014_BqTauTauSM`, `CapdevilaEtAl2018_BSTauTauNP`, `BordoneNavarro2023_BSTauTauNP`, and `CsakiFalkowskiWeiler2008_RSFlavor`. Every manifest entry has a non-empty `snapshot_path`, and `git ls-files` confirms the text snapshots are tracked.
- CHK-3: `find flavor_catalog/references/B007 flavor_catalog/references/B008 -type f` shows only text snapshots, `sha256sums.txt`, and `source_manifest.yaml`; `find ... -iname '*.pdf'` found no publisher PDFs.
- CHK-4: both sidecars show legal `DRAFT -> WRITER-INITIATED -> WRITER-DONE` histories with ISO 8601 timestamps before this CA transition.
- CHK-5: B007 `code_coverage.status: NO` is supported. Exact focused greps for `B_s0/B0 -> e+e-`, dielectron, and `b -> s,d e+e-` rare-leptonic patterns in implementation/test directories returned no hits. The cited nearby evidence lines exist and are neutral-B mixing: `quarkConstraints/deltaf2.py:225`, `:239`, `:903`, and `:922`.
- CHK-5: B008 `code_coverage.status: NO` is supported. Exact focused greps for `B_s0/B0 -> tau+tau-`, ditau, tautau, and `b -> s,d tau+tau-` rare-leptonic patterns in implementation/test directories returned no hits. The cited nearby evidence lines exist: B mixing at `quarkConstraints/deltaf2.py:903` and `:922`, bridge B mixing at `quarkConstraints/modern/phenomenology.py:646` and `:657`, unrelated `muToEGamma.py:75`, and `scanParams/README.md:161`.
- CHK-6: `MEDIUM` is consistent for both. Each needs a new Delta B = 1 pure-leptonic branching-ratio evaluator and standard Wilson/input plumbing, but not a new RG calculation, lattice calculation, or exclusive form-factor/angular likelihood.
- CHK-7: searches over `flavor_catalog/`, `docs/`, the plan rows, and methodology-note/rc1.1 context found no load-bearing B007/B008 numerical contradiction. Existing references outside the process files list these as previously NOT-SEEN/deferred targets rather than supplying conflicting values.
- CHK-8: no `\cite`, `\ref`, `\textbf{CHECK}`, TODO/FIXME, unresolved markers, or `??` markers were found in the B007/B008 TeX files. A simple delimiter pass found balanced braces and zero raw dollar delimiters; display math uses paired `\[` / `\]`.
- CHK-W8: both files live under `flavor_catalog/processes/secondary/beauty/`. Both sidecars carry `priority_tier: SECONDARY`, non-empty `priority_rationale`, and `promoted_in_wave: 8`. Both TeX files carry a SECONDARY-tier note pointing at `flavor_catalog/PRIORITY_TIERS.md`.
