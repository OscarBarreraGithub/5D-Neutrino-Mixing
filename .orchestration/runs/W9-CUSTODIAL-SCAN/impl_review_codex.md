Findings: none.

Verified against the W9 plan:
- Minimal payloads omit `ew_model`; pinned hashes remain `45e21a07585f7489` and `c6939cc65d71f86a`.
- Worker round-trip preserves `ew_model` for minimal and custodial configs.
- Harness build sites thread `cfg.ew_model` through spectrum `model_label` and builder `ew_model`, including universal-c sanity.
- Comparison builder normalizes `seed`, `params.M_KK`, and `quark_fit_r`, pairs on `(r, mkk_tev, draw_seed)`, writes README/schema, survival aggregates, and veto classes `rigorous|proxy|not_evaluated`.
- Custodial sbatch reuses `scripts/wq_quarkonly_1m_plan.py`; planner output exactly matches baseline `scan_plan.json`, and the command adds `--ew-model custodial_rs_plr --quark-only`.
- No fabricated production data found; builder streams raw JSONL inputs and tests use temp synthetic fixtures only.

Requested test command passed: `27 passed in 11.43s`.

VERDICT: APPROVE