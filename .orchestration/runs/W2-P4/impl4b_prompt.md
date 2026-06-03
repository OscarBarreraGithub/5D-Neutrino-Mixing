# W2 PHASE 4 — SUB-STEP 4b: Z-LFV rewire (codex, gpt-5.x xhigh). Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. Implement sub-step 4b of the DUAL-APPROVED plan `.orchestration/runs/W2-P4/plan.md` (step 12), using the 4a lepton builder (committed 7265c8d: `quarkConstraints/rs_ew_couplings` `z_delta_g_{L,R}_e`) and the 3b/3c/3d rewire+degradation pattern. First a SHORT plan, then implement. Codex + Opus dual-review (both must APPROVE). SCOPE = Z-LFV ONLY (T015/T016/T017).

REWIRE (lepton-flavor-violating Z decays):
- **T015** (Z→eμ), **T016** (Z→eτ), **T017** (Z→μτ): replace the proxy path with the rigorous off-diagonal charged-lepton coupling via `z_lfv_branching_fraction_from_couplings(delta_g_left=z_delta_g_L_e[i,j], delta_g_right=z_delta_g_R_e[i,j])` (eμ:[0,1], eτ:[0,2], μτ:[1,2]). Use the adapter-exported name (`quarkConstraints/zpole_lfv.py`).
- **IMPORTANT v1 PHYSICS**: with the diagonal charged-lepton fit (U_e=I, universal c_L), `z_delta_g_e` is DIAGONAL ⇒ off-diagonal=0 ⇒ tree-level Z-LFV NP = 0. This is CORRECT (loop/dipole-spurion LFV is deferred). So in v1 these constraints become RIGOROUS-tree but trivially NON-VETOING (BR_NP=0). KEEP a PARTIAL/NEEDS-HUMAN diagnostic noting loop-induced Z-LFV (from the lfv_dipole_spurion) is deferred to Phase 7 — do NOT claim full Z-LFV rigor, and do NOT fake a nonzero bound. Remove only the OLD overlap-proxy NEEDS-HUMAN, replacing it with the honest "tree-level rigorous (=0 for diagonal fit); loop-LFV deferred" status.
- SM (Z→ℓ_iℓ_j is zero in SM) / anchors / limits UNCHANGED.

GRACEFUL DEGRADATION (as prior sub-steps): rigorous `rs_ew_couplings` when present; ABSENT ⇒ non-vetoing `evaluated=False` + `missing_extra` (no crash/fake pass).

TESTS:
- Rigorous v1-path: point via `build_from_rs_ew_inputs` (diagonal fit) ⇒ off-diagonal z_delta_g_e=0 ⇒ BR_NP=0, non-vetoing pass (rigorous-tree).
- LFV-LIVE path: feed a NON-diagonal U_e / non-universal lepton-c toy point ⇒ off-diagonal z_delta_g_e≠0 ⇒ T015/T016/T017 produce NONZERO BR (cross-check independently of the adapter) ⇒ confirm the rewire actually bites when LFV exists, and scales as (m_Z/M_KK)^4.
- Absent-path: old-style point ⇒ non-vetoing evaluated=False.
- Replace proxy-only tests with these — enumerate, no silent coverage loss. `python -m pytest tests/ -q` stays green.

CONSTRAINTS: physics ONLY via the zpole_lfv adapter; ConstraintResult numeric fields real finite floats (complex couplings in diagnostics); touch ONLY T015/T016/T017 + their adapter/tests.

OUTPUT (<=14 lines): short plan; rewired files + IDs; v1-zero (rigorous-tree, non-vetoing) + LFV-live-toy nonzero BR + scaling + absent-path; the PARTIAL/loop-deferred flag; test-count change; pytest counts. End with: P4B-DONE.
