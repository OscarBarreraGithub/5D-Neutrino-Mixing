# WQ-1M FOLLOW-UP FIX (codex author, gpt-5.x xhigh)

Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. A real defect was found in the COMMITTED `scripts/run_wq_quarkonly_1m_array.sbatch` (commit 2c0db4a): each of the 50 array tasks computes its OWN `TIMESTAMP="$(date ...)"`, so the output root `scan_outputs/wq_quarkonly_1M_${TIMESTAMP}` DIFFERS per task — the 1M run scattered across 10 timestamped roots instead of one. (No data loss this time — seeds are disjoint so all 1M rows are unique; the orchestrator merged them. But a re-run would scatter again, and two tasks sharing the same second could collide on an output path and resume-skip each other → silent data loss.)

## FIX (minimal, surgical)
Make ALL array tasks share ONE output root. Replace the per-task `date`-derived root with a STABLE identifier common to the whole array — `SLURM_ARRAY_JOB_ID` is identical across all tasks of one array submission. E.g. `OUTPUT_ROOT="${WQ_OUTPUT_ROOT:-${REPO_ROOT}/scan_outputs/wq_quarkonly_1M_${SLURM_ARRAY_JOB_ID}}"`. Keep the `${WQ_OUTPUT_ROOT}` override. Keep the per-task `${WQ_R_LABEL}/shard-${WQ_DRAW_SHARD_PADDED}` subdir layout and ALL seed logic (base/stride/draws) UNCHANGED — seed-disjointness must be preserved exactly. Guard against an unset `SLURM_ARRAY_JOB_ID` (non-array invocation) with a sane fallback. Ensure the scan_plan.json / seed_disjointness.md sidecar writes still land once under the shared root (idempotent / first-writer-wins is fine).

Do NOT change the harness, the analysis script, or any physics. Touch ONLY the sbatch (and a tiny test if practical). Confirm `python -m pytest tests/ -q` stays green (this is an sbatch change; tests should be unaffected).

## OUTPUT (<=10 lines)
The exact diff (before/after lines); confirm seed logic UNCHANGED; confirm all tasks now resolve to one root; confirm the unset-SLURM_ARRAY_JOB_ID fallback; pytest counts. End with: WQ-1M-FIX-DONE.
