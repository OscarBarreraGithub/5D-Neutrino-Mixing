# W9 — Custodial-ON 1M scan + comparison-ready save PLAN (codex author)

You are authoring an IMPLEMENTATION + RUN PLAN only. Do NOT write production code or launch jobs
in this run. Output the plan to `.orchestration/runs/W9-CUSTODIAL-SCAN/plan_codex.md`, end `PLAN-READY`.

## Goal
Produce a custodial-ON quark-only 1M scan that is APPLES-TO-APPLES with the existing non-custodial
1M scan, and save BOTH in a structure + documented schema so a future UI can plot
**no-custodial vs custodial** M_KK survival side by side.

## Ground truth to read
- Non-custodial baseline: `scan_outputs/wq_quarkonly_1M_20128400/` (its `scan_plan.json` has the
  EXACT grid: base_seed 202606040000, r_grid [0.05,0.1,0.25,0.5,1.0], mkk_tev
  [1,2,3,5,7,10,15,20,30,50], 20000 draws/mkk/r, seed formulas). Custodial run MUST reuse the
  SAME grid + SAME seed formulas so points are identical draw-for-draw.
- Harness: `scripts/run_full_catalog_scan.py` (has `--quark-only`; NO `--ew-model` flag yet).
- sbatch: `scripts/run_wq_quarkonly_1m_array.sbatch` (the 50-task array driver).
- Custodial branch: `ew_model="custodial_rs_plr"` in `quarkConstraints/rs_ew_couplings.py`;
  oblique custodial in `oblique_stu.py`. The point_builder/harness must be able to select it.
- Existing analysis: `scan_outputs/wq_quarkonly_1M_20128400/analysis/` (mirror its outputs).
- Will run AFTER W7 (μ→eγ) and W8 (custodial PR2) land — so plan for the custodial branch to
  include PR2's loop term if available, but the scan is QUARK-ONLY (μ→eγ is lepton sector, not in it).

## What the plan must contain
1. **Harness wiring (code change, dual-gated):** add an `--ew-model {minimal_rs,custodial_rs_plr}`
   flag that flows through point_builder so EW001/T010/T011 evaluate in the chosen mode. Default
   `minimal_rs` → existing behavior byte-identical (config hash unchanged). Tests for both modes.
2. **Run config:** a custodial sbatch (or a parameter to the existing one) using the SAME grid +
   seeds as the baseline, output root e.g. `scan_outputs/wq_quarkonly_1M_custodial_<jobid>/`.
   Resource estimate per tile + total; checkpoint/merge procedure (mirror the baseline's).
3. **Comparison-ready save (the user's key requirement):** a documented schema + a `comparison/`
   manifest that pairs the two runs: identical (r, M_KK, draw_seed) keys, per-(r,M_KK) survival
   fractions for BOTH runs, which constraint vetoes each, and a README so a UI can load both and
   plot survival-vs-M_KK (minimal vs custodial) per r. Specify exact file formats (JSON/parquet/csv).
4. **Analysis/plots:** mirror the baseline analysis on the custodial run + a side-by-side
   minimal-vs-custodial headline figure (survival vs M_KK per r; expect the floor to relax from
   ~25-30 TeV physical toward ~2-3 TeV under custodial).
5. **Isolation:** files touched (harness flag overlaps NOTHING in W7/W8 ideally — flag if it does);
   note that the actual SLURM run + analysis happen after W7/W8 commits.

## Constraints
Dual-gate: a second codex + Opus review this plan. Cite file:line, exact grid/seed reuse, file
formats, tests. Determinism, no fabricated physics, honest tags. End with `PLAN-READY`.
