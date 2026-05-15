# Invalidation-Gate RS-Anarchy Rerun

Date: 2026-05-15
Branch: `invalidation/gate-rerun`
Base branch: `paper/quark-scan-2026q2`

## SLURM submissions

All nine jobs were submitted on `serial_requeue` with the existing dispatcher
scripts and completed with exit code `0:0`.

| Run | Job ID | Output directory | Final state | Elapsed |
|---|---:|---|---|---:|
| RUNA baseline | 13024513 | `scan_outputs/rs_anarchy_runA_20260515T085316` | COMPLETED | 01:05:53 |
| Run 3 baseline | 13024514 | `scan_outputs/rs_anarchy_run3_baseline_20260515T085324` | COMPLETED | 00:14:47 |
| Run 3 qtop_shifted | 13024515 | `scan_outputs/rs_anarchy_run3_qtop_shifted_20260515T085324` | COMPLETED | 00:15:25 |
| Run 3 moreUV | 13024516 | `scan_outputs/rs_anarchy_run3_moreUV_20260515T085324` | COMPLETED | 00:14:03 |
| Run 3 moreIR | 13024517 | `scan_outputs/rs_anarchy_run3_moreIR_20260515T085324` | COMPLETED | 00:14:00 |
| Run B narrow_uniform | 13024518 | `scan_outputs/rs_anarchy_runB_narrow_uniform_20260515T085324` | COMPLETED | 00:16:21 |
| Run B wide_uniform | 13024519 | `scan_outputs/rs_anarchy_runB_wide_uniform_20260515T085325` | COMPLETED | 00:13:05 |
| Run B gaussian_3sigma | 13024520 | `scan_outputs/rs_anarchy_runB_gaussian_3sigma_20260515T085324` | COMPLETED | 00:15:19 |
| Run C CFW-like gate | 13024521 | `scan_outputs/rs_anarchy_runC_20260515T085323` | COMPLETED | 00:37:23 |

No job remained pending after 20 minutes and no manual requeue was needed.

## RUNA halt trigger

The M_KK = 3 TeV PDG-pass count is unchanged, as expected for a
Delta F=2-only invalidation:

| Quantity | Pre-audit RUNA | Post-audit RUNA | Change |
|---|---:|---:|---:|
| `n_pdg_pass` at 3 TeV | 164,080 | 164,080 | 0.000% |
| `pdg_pass_fraction` at 3 TeV | 0.16408 | 0.16408 | 0.000% |

The 3 TeV `max_ratio_percentiles_pdg` are much larger post-audit:

| Percentile | Old | New | New/old |
|---|---:|---:|---:|
| p05 | 0.1088 | 1.2319 | 11.33 |
| p25 | 0.3894 | 6.8108 | 17.49 |
| p50 | 0.9677 | 16.9439 | 17.51 |
| p75 | 2.1995 | 38.0196 | 17.29 |
| p95 | 6.8278 | 119.1491 | 17.45 |

The signed-off 22.49x cumulative factor is a benchmark epsilon_K ratio,
not a constant multiplier for every anarchic draw. A direct old/new
streaming comparison of the 3 TeV RUNA PDG-passing draws gives median
per-draw epsilon_K shift 18.55x and median max-ratio shift 18.11x. The
draw-dependent operator mixture and max-over-systems observable explain
the difference from the benchmark value.

## Before/after headline numbers

The headline values are computed from pooled PDG-passing draws via
`M_KK_min = M_KK_tile * sqrt(max_ratio)`.

| Quantity | Pre-audit RUNA | Post-audit RUNA | New/old |
|---|---:|---:|---:|
| p50, perturbative `g_s` | 3.70 TeV | 16.54 TeV | 4.47 |
| p50, `g_s*=3` | 10.58 TeV | 47.26 TeV | 4.47 |
| p95, perturbative `g_s` | 9.98 TeV | 44.50 TeV | 4.46 |
| p95, `g_s*=3` | 28.50 TeV | 127.13 TeV | 4.46 |

Finite-Monte-Carlo Wilson-score 95% intervals for RUNA:

| Crossing | Perturbative `g_s` | `g_s*=3` |
|---|---:|---:|
| p50 | [16.52, 16.56] TeV | [47.19, 47.32] TeV |
| p95 | [44.41, 44.58] TeV | [126.88, 127.38] TeV |

BGS NP-budget band on the p50 headline:

| Convention | Central | Budget band |
|---|---:|---:|
| Perturbative `g_s` | 16.54 TeV | -8.74 / +24.28 TeV |
| `g_s*=3` | 47.26 TeV | -24.98 / +69.37 TeV |

Thus the methodology-note headline quote is
`M_KK^min p50 = 47.26^{+69.37}_{-24.98} TeV` in the `g_s*=3`
convention, conditional on the central BGS + LO-only conventions.

## Follow-up crossing comparison

| Run | Old p50 `g_s*=3` | New p50 `g_s*=3` | Old p95 `g_s*=3` | New p95 `g_s*=3` |
|---|---:|---:|---:|---:|
| RUNA | 10.58 | 47.26 | 28.50 | 127.13 |
| Run 3 baseline | 10.57 | 47.26 | 28.44 | 126.76 |
| Run 3 qtop_shifted | 10.22 | 45.73 | 27.00 | 120.39 |
| Run B narrow_uniform | 11.03 | 49.79 | 29.12 | 132.35 |
| Run B wide_uniform | 11.91 | 52.76 | 32.27 | 139.34 |
| Run B gaussian_3sigma | 10.89 | 48.32 | 30.97 | 139.55 |
| Run 3 moreUV | zero PDG passes | zero PDG passes | zero PDG passes | zero PDG passes |
| Run 3 moreIR | zero PDG passes | zero PDG passes | zero PDG passes | zero PDG passes |
| Run C | zero PDG passes | zero PDG passes | zero PDG passes | zero PDG passes |

`scan_outputs/followup_crossings_summary.json` was regenerated with the
post-audit run paths and retains the pre-audit values under
`pre_audit_reference`.

## Per-system changes

At RUNA's 3 TeV tile, per-system pass fractions among PDG-passing draws
changed as follows:

| System | Old pass fraction | New pass fraction | Change |
|---|---:|---:|---:|
| epsilon_K | 51.224% | 4.089% | -47.135 pp |
| Delta_m_K | 100.000% | 99.999% | -0.001 pp |
| Delta_m_Bd | 99.938% | 99.993% | +0.055 pp |
| Delta_m_Bs | 99.923% | 100.000% | +0.077 pp |
| Delta_m_D0 | 99.786% | 98.932% | -0.854 pp |

The post-audit RUNA binding-system fractions over all pooled PDG-passing
draws are:

| Binding system | Fraction |
|---|---:|
| epsilon_K | 99.5001% |
| Delta_m_K | 0.0005% |
| Delta_m_Bd | 0.0152% |
| Delta_m_Bs | 0.0453% |
| Delta_m_D0 | 0.4389% |

epsilon_K still binds, and more strongly than before.

## Plot and document regeneration

The pre-audit figure snapshot is stored at
`results/figures/quark_pre_audit_constants/`. The canonical figures under
`results/figures/quark/` were regenerated from the new run directories,
including the RUNA histogram, `g_s*=3` histogram, summary plots, gate
sensitivity, c-value overlay, Y-prior overlays, PDG-tightness plot, and
CFW comparison. `docs/quark_scan_methodology_note.pdf` was rebuilt with
`pdflatex` and now has 18 pages.

Invalidation gate cleared; final manifest sign-off can proceed.
