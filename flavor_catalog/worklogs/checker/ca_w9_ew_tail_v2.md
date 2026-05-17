# CA Worklog: ca_w9_ew_tail_v2
**Date**: 2026-05-17
**Family**: collider_rs
**Batch**: ew_tail
**Cycle**: 2
**Checker agent**: CA-v2-w9-ew-tail
**Process IDs**: CR009 CR011

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| CR009 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| CR011 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

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

## Per-process issues
- CR009: None.
- CR011: None.

## Evidence notes
- CHK-1: PASS. CR009 cycle-1 historical measured contact-range numerals are absent from `CR009.tex` and `CR009.yaml` prose after the WA-v2 fix; the remaining current ATLAS/CMS observed limits are promoted under `pdg_or_equivalent.values` with value, year, CL/uncertainty-null, units, source URL, access date, snapshot path, and sha256. Applied the L001/T003 carve-out to theory-normalization and schematic RS-matching numerals such as \(4\pi/\Lambda^2\), \(O(10~\mathrm{TeV})\), and the local 47.26/127.13 TeV comparison. CR011 measured 0.45 fb and 1.17 fb limits are promoted; dataset luminosities, expected limits, significance values, and EFT/aQGC context remain supporting metadata.
- CHK-2: PASS. CR009 and CR011 TeX Key references resolve to canonical keys in their process-local `source_manifest.yaml` files. The CR011 cycle-1 snapshot-stem mismatch is fixed: the Key references now use `ATLAS2025_LongitudinalWW`, `CMS2020_PolarizedSameSignWW`, `ATLAS2024_SameSignWWjj`, `ATLAS2020_ZZjj`, `ATLAS2017_WWjjAQGC`, `ATLAS2026_AQGCCombination`, `PDG2025_WZQuarticCouplings`, and `Contino2011_CompositeHiggsResonances`. `git ls-files --error-unmatch` confirmed all manifest snapshot paths are tracked, and `sha256sum -c sha256sums.txt` passed in both reference directories.
- CHK-3: PASS. `find flavor_catalog/references/CR009 flavor_catalog/references/CR011 -type f` with PDF/office-document filters returned no files. Both directories contain text snapshots, `source_manifest.yaml`, and `sha256sums.txt` only.
- CHK-4: PASS. Before this checker update, both YAML sidecars had the legal cycle-2 path ending in `WRITER-DONE` from WA-v2 with ISO 8601 `at` values. This check appends `CHECKER-DONE`, sets `checker_passed_at`, and updates `last_updated_at`.
- CHK-5: PASS. Targeted implementation searches across `quarkConstraints`, `qcd`, `flavorConstraints`, `neutrinos`, `yukawa`, `warpConfig`, `solvers`, `scanParams`, and `tests` found no Drell-Yan/contact-interaction, VBS/polarization/aQGC, or collider-reinterpretation implementation. The broad ATLAS/CMS/collider grep returned only the unrelated CMS/RunDec alpha-s example at `tests/test_alpha_s.py:88`-`89`. The cited adjacent low-energy flavor and KK-scale bookkeeping lines exist in `quarkConstraints/scan.py`, `quarkConstraints/deltaf2.py`, and `scanParams/scan.py`.
- CHK-6: PASS. `HIGH` remains consistent. CR009 needs a binned high-mass dilepton likelihood or SMEFT/collider reinterpretation with RS neutral-vector matching and EFT-validity logic. CR011 needs polarization-sensitive VBS simulation/recast machinery, model-dependent resonance spectrum assumptions, acceptance, and likelihood support.
- CHK-7: PASS. Neither entry contradicts the rc1.1/methodology-note low-energy flavor result. `docs/quark_scan_methodology_note.tex` still carries the 47.26 TeV p50 and 127.13 TeV p95 \(\gs=3\) crossings, and both entries frame the collider-tail information as cross-check or non-anarchic-variant leverage rather than replacing the low-energy flavor bound.
- CHK-8: PASS. Fixed-string searches found no `CHECK`, `TODO`, `FIXME`, `MISSING`, unresolved `\cite{}`, unresolved `\ref{}`, or `??` markers in `CR009.tex` or `CR011.tex`. Lightweight delimiter counts are balanced: CR009 has `\(`/`\)` = 39/39, `\[`/`\]` = 1/1, braces = 82/82, dollar count = 0; CR011 has `\(`/`\)` = 36/36, `\[`/`\]` = 1/1, braces = 59/59, dollar count = 0. `git diff --check` was clean before checker-status edits.
