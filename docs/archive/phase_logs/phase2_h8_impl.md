# Phase 2 Hole #8 Implementation: Zero-Pass Finite Statistics

Date: 2026-05-15
Branch: `scan/zero-pass-statistics`
Base branch: `paper/quark-scan-2026q2`

## Inventory

| Run | Run dir | N | k | Wilson 95% `p_UL` |
|---|---|---:|---:|---:|
| Run 3 moreUV | `scan_outputs/rs_anarchy_run3_moreUV_20260515T085324` | 1,600,000 | 0 | 2.304e-6 |
| Run 3 moreIR | `scan_outputs/rs_anarchy_run3_moreIR_20260515T085324` | 1,600,000 | 0 | 2.304e-6 |
| Run C CFW-like gate | `scan_outputs/rs_anarchy_runC_20260515T085323` | 4,000,000 | 0 | 9.216e-7 |

Full provenance, gates, seeds, tile counts, c-patterns, Y-priors, and
`tile_summary.json` SHA-256 hashes are recorded in
`docs/audits/zero_pass_inventory.md`.

## Methodology Note Locations Updated

| Location | Change |
|---|---|
| `docs/quark_scan_methodology_note.tex:487-493` | Replaced the absolute 2-sigma zero-pass wording with finite-ensemble caveat language. |
| `docs/quark_scan_methodology_note.tex:739-743` | Updated Run 3 moreUV/moreIR figure caption to quote N and Wilson 95% upper limits. |
| `docs/quark_scan_methodology_note.tex:757-769` | Replaced "exactly zero" wording with explicit `N=1,600,000`, `p_UL=2.304e-6`, and finite-ensemble interpretation. |
| `docs/quark_scan_methodology_note.tex:853-858` | Added Run C factor-1.5/2.5 zero-pass finite-statistics bound. |
| `docs/phase_logs/invalidation_gate_rerun.md:93-95` | Replaced zero-pass shorthand with N and `p_UL` values. |
| `docs/phase_logs/invalidation_gate_signoff.md:76-82` | Replaced absolute zero-pass summary with finite-ensemble bounds. |

## Code And Figures

- Added `quarkConstraints/finite_stats.py::wilson_upper_limit`.
- Added `tests/test_finite_stats.py`.
- Updated `scripts/rs_anarchy_mkk_min_hist_by_cvals.py` so empty-run legend
  entries show `N` and Wilson 95% `p_UL`; regenerated
  `results/figures/quark/rs_anarchy_mkk_min_hist_by_cvals.{pdf,png}`.
- Updated `scripts/rs_anarchy_cfw_comparison.py` so Run C-style zero stored
  PDG-pass fallback plots carry the finite-ensemble upper-limit note.

## Verification

| Check | Result |
|---|---|
| `pytest -q tests/test_finite_stats.py` | 3 passed in 0.05 s |
| `pytest -q tests/test_finite_stats.py tests/test_cfw_comparison.py` | 5 passed in 76.05 s |
| Run 3 c-value plot regeneration | exited 0; moreUV/moreIR legend entries show `N=1.6M` and `p_UL=2.3e-6` |
| Run C CFW no-plot fallback check | exited 0; summary includes `zero_pass_note = stored gate: N=4.0M, p_UL=9.2e-7` |
| `pdflatex -interaction=nonstopmode quark_scan_methodology_note.tex` from `docs/` | rebuilt `docs/quark_scan_methodology_note.pdf` |

Hole #8 ready for peer review.
