```diff
-#     python scripts/analyze_wq_quarkonly.py --input scan_outputs/wq_quarkonly_1M_<ts>
+#     python scripts/analyze_wq_quarkonly.py --input scan_outputs/wq_quarkonly_1M_<array_job_id>
-TIMESTAMP="${WQ_RUN_TIMESTAMP:-$(date -u +%Y%m%dT%H%M%S)}"
-OUTPUT_ROOT="${WQ_OUTPUT_ROOT:-${REPO_ROOT}/scan_outputs/wq_quarkonly_1M_${TIMESTAMP}}"
+WQ_OUTPUT_RUN_ID="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"
+OUTPUT_ROOT="${WQ_OUTPUT_ROOT:-${REPO_ROOT}/scan_outputs/wq_quarkonly_1M_${WQ_OUTPUT_RUN_ID}}"
```
Seed logic UNCHANGED: no edits to task plan eval, base/stride/draw vars, or scan invocation.
All array tasks now share `wq_quarkonly_1M_${SLURM_ARRAY_JOB_ID}`; task-0 sidecars still land under that shared root.
Unset fallback: `${SLURM_JOB_ID:-manual}`. Pytest: `1713 passed, 1 skipped in 799.79s`.
WQ-1M-FIX-DONE