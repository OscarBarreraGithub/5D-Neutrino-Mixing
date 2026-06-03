# W2 PHASE 3 — SUB-STEP 3b: Z-pole rewire (codex, gpt-5.x xhigh). Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. Implement EXACTLY sub-step 3b of the DUAL-APPROVED plan `.orchestration/runs/W2-P3/plan.md` (steps 6-8) using the 3a builder (committed 7813e6c: `flavor_catalog_constraints/rs_ew_builder.build_from_rs_ew_inputs`, `quarkConstraints/rs_ew_couplings` with `z_delta_g_{L,R}_{u,d}`). First a SHORT plan, then implement. Codex + Opus dual-review (both must APPROVE).

REWIRE (Z-pole quark neutral currents):
- **T010, T011** (Z→bb̄: R_b / A_b^FB etc.): replace the proxy NP path (`zpole_evaluate_zbb_with_proxy(...)`) with the rigorous shift: `zpole_shifted_couplings(zpole_sm_couplings("b"), delta_g_left=z_delta_g_L_d[2,2], delta_g_right=z_delta_g_R_d[2,2])` then the existing `zpole_evaluate_quark`. Keep SM radiators/anchors/budgets/pulls UNCHANGED. **KEEP a PARTIAL/NEEDS-HUMAN diagnostic** on T010/T011 for the missing fermion-KK/custodial/BKT Zbb completion (Phase 6) — the gauge-NC piece is now rigorous but classic Zbb fermion-mixing is not.
- **T012** (Z→cc̄: R_c / A_c): replace `zpole_evaluate_zcc_with_proxy(...)` with the same `shifted_couplings` path using `z_delta_g_L_u[1,1]` / `z_delta_g_R_u[1,1]` (charm). Remove the proxy NEEDS-HUMAN flag for T012 minimal-RS gauge NC (now rigorous). (T013 has no standalone file; covered by T012.)

GRACEFUL DEGRADATION (critical — keep the whole suite green): the rewired constraints read the rigorous `rs_ew_couplings` from the ParameterPoint WHEN PRESENT (a point built via `build_from_rs_ew_inputs`). When ABSENT (old-style points from `build_from_quark_couplings`, which the existing tests + much of the suite use), the constraint must DEGRADE GRACEFULLY: return a non-vetoing result with `evaluated=False` + a clear "rs_ew_couplings not provided" note — NOT crash, NOT a silent fake pass. (The 100M scan will always provide rs_ew_couplings.)

TESTS (update the affected constraint tests):
- Add a RIGOROUS-PATH test: build a point via `build_from_rs_ew_inputs`, confirm T010/T011/T012 now consume z_delta_g (cross-check the shifted Z-bb̄/Z-cc̄ couplings independently of the adapter).
- Add/adjust the ABSENT-PATH test: an old-style point ⇒ non-vetoing evaluated=False (no crash).
- SM-limit: a universal-c / a(c)=a_ref point ⇒ rewired constraints recover their committed SM-only output (δg=0).
- Keep every other constraint's behavior unchanged; `python -m pytest tests/ -q` stays green (≈1652+).

CONSTRAINTS: reach physics ONLY via the zpole adapter; numeric ConstraintResult fields real finite floats (complex couplings in diagnostics); do not touch unrelated constraints/cores. One sub-step = T010/T011/T012 only.

OUTPUT (<=14 lines): short plan; the rewired files (file:line); the rigorous-path shifted Z-bb̄/Z-cc̄ coupling numbers + SM-limit recovery + the absent-path non-vetoing behavior; pytest counts. End with: P3B-DONE.
