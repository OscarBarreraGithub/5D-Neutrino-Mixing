# W2 PHASE 5 — SUB-STEP 5c: EW003 (codex, gpt-5.x xhigh). Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. Implement sub-step 5c of the DUAL-APPROVED plan `.orchestration/runs/W2-P5/plan.md` (step 13/17), using the 5a builder (committed 3f7762a: `rs_charged_current` epsilon) + the 5b charged_current adapter. First a SHORT plan, then implement. Codex + Opus dual-review (both must APPROVE). SCOPE = EW003 ONLY (last Phase-5 sub-step).

EW003 (|Vcb|/|Vub| inclusive-vs-exclusive):
- ADD charged-current diagnostics for `epsilon_cb` and `epsilon_ub` (from rs_charged_current) to EW003's result diagnostics.
- KEEP the main inclusive-vs-exclusive PULL DATA-LEVEL (the existing rigorous data tension) UNLESS a covariance/scheme input is supplied — do NOT convert it to a naive NP veto. Reason (per plan): a common vector rescaling `|1+epsilon|` CANCELS within each inclusive/exclusive determination of the same |V_ij|, so a universal CC shift does NOT move the incl-vs-excl RATIO; only a method-dependent (incl≠excl) shift would, which requires a covariance/scheme model not available. The correlated CKM/form-factor uncertainties must NOT be treated as independent quadrature.
- So EW003 stays PARTIAL/data-level: the SM-vs-data pull is rigorous; the RS-CC effect on the incl-vs-excl ratio is a documented diagnostic (≈cancels in v1), with a NEEDS-HUMAN note that a rigorous incl-vs-excl differential treatment needs a covariance/scheme input. Honest — do NOT claim a rigorous NP veto.

GRACEFUL DEGRADATION: if rs_charged_current present, surface the epsilon_cb/epsilon_ub diagnostics; if absent, behave as before (data-level pull), non-crashing.

TESTS: confirm (a) the incl-vs-excl PULL is unchanged vs the committed value (data-level preserved); (b) with a point, the epsilon_cb/epsilon_ub diagnostics appear and a UNIVERSAL CC shift cancels in the ratio (≈no change to the pull) — assert the cancellation; (c) absent-path non-crashing; (d) the NEEDS-HUMAN covariance note is present. Update tests accordingly — enumerate, no silent coverage loss. `python -m pytest tests/ -q` stays green.

CONSTRAINTS: physics ONLY via adapters; ConstraintResult numeric fields real finite floats; touch ONLY EW003 + its adapter usage/tests.

OUTPUT (<=12 lines): short plan; what changed in EW003; the data-level pull unchanged + epsilon diagnostics + universal-cancellation check + covariance NEEDS-HUMAN note; test-count change; pytest counts. End with: P5C-DONE.
