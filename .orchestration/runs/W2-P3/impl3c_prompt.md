# W2 PHASE 3 — SUB-STEP 3c: FCNC-Z rewire (codex, gpt-5.x xhigh). Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. Implement EXACTLY sub-step 3c of the DUAL-APPROVED plan `.orchestration/runs/W2-P3/plan.md` (step 9), using the 3a builder (committed: `quarkConstraints/rs_ew_couplings` z_delta_g_{L,R}_d) and following the 3b rewire pattern (committed 7c0ecff). First a SHORT plan, then implement. Codex + Opus dual-review (both must APPROVE).

REWIRE (FCNC-Z, down-sector):
- **T014** (Z→b s̄+b̄ s, b d̄+b̄ d, s d̄+s̄ d): replace the proxy path `zpole_down_fcnc_branching_fraction_with_proxy(...)` with `zpole_down_fcnc_branching_fraction_from_couplings(delta_g_left=z_delta_g_L_d[i,j], delta_g_right=z_delta_g_R_d[i,j])` for each channel (bs: [1,2]/[2,1]; bd: [0,2]/[2,0]; sd: [0,1]/[1,0] — use the correct generation indices for the down mass-basis off-diagonal entries). Keep the SM-zero / direct-width policy, the budgets, and the B<2.9e-3 limits UNCHANGED. Remove the proxy NEEDS-HUMAN flag (down-sector FCNC-Z is now rigorous from the gauge-NC misalignment).

GRACEFUL DEGRADATION (same as 3b): use rigorous `rs_ew_couplings` when present; when ABSENT (old-style points), return a non-vetoing result `evaluated=False` + `missing_extra=rs_ew_couplings` (no crash, no fake pass). The 100M scan always provides it.

TESTS (update T014 tests):
- Rigorous-path: build a point via `build_from_rs_ew_inputs`, confirm the FCNC width/BR now comes from the OFF-DIAGONAL z_delta_g_d[i,j] (cross-check the off-diagonal coupling → width independently of the adapter).
- SM-limit: universal-c / a(c)=a_ref ⇒ off-diagonal z_delta_g_d=0 ⇒ FCNC BR=0 ⇒ non-vetoing pass (recover committed SM-zero behavior).
- Absent-path: old-style point ⇒ non-vetoing evaluated=False.
- Keep all other constraints unchanged; if proxy-only T014 tests are removed, replace with rigorous/absent/SM-limit equivalents (do NOT silently drop coverage — enumerate the change). `python -m pytest tests/ -q` stays green.

CONSTRAINTS: reach physics ONLY via the zpole adapter; numeric ConstraintResult fields real finite floats (complex couplings in diagnostics); touch ONLY T014 + its tests (+ a helper if needed). 

OUTPUT (<=14 lines): short plan; rewired file (file:line); the rigorous-path off-diagonal FCNC coupling → BR for one channel (e.g. bs) + SM-limit=0 + absent-path non-vetoing; the test-count change (enumerated); pytest counts. End with: P3C-DONE.
