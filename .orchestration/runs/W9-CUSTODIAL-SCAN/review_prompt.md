# W9 plan REVIEW — independent critique

Independently critique the plan at `.orchestration/runs/W9-CUSTODIAL-SCAN/plan_codex.md`.
Also read `.orchestration/runs/W9-CUSTODIAL-SCAN/plan_prompt.md` (the requirements) and verify
against the REAL repo. Do NOT rubber-stamp. Check specifically:

1. **Apples-to-apples reuse:** does the plan reuse the EXACT grid + seed formulas from
   `scan_outputs/wq_quarkonly_1M_20128400/scan_plan.json` so the custodial run is draw-for-draw
   identical to the baseline? Any seed/grid drift = BLOCKING.
2. **Harness `--ew-model` flag:** is the wiring through `scripts/run_full_catalog_scan.py` →
   point_builder → EW001/T010/T011 correct and does default `minimal_rs` stay byte-identical
   (config hash unchanged)? Are both-mode tests specified?
3. **Comparison-ready save:** is the `comparison/` schema concrete and sufficient for a UI to plot
   minimal-vs-custodial survival-vs-M_KK per r (paired identical-seed keys, per-(r,M_KK) survival
   for BOTH runs, which constraint vetoes, README)? File formats specified?
4. **Correctness/honesty:** determinism, honest rigorous|proxy|partial tags, no fabricated physics,
   resource estimate sane, checkpoint/merge mirrors baseline.
5. **Isolation:** files-touched list accurate; overlaps with W7 (lepton)/W8 (custodial PR2) flagged.

End your review with EXACTLY one line: `VERDICT: APPROVE` or `VERDICT: NEEDS-FIXES`,
followed by a numbered list of any required fixes (specific, file:line where possible).
