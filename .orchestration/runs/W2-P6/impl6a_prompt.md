# W2 PHASE 6 — SUB-STEP 6a: Zbb fermion-KK mixing (codex, gpt-5.x xhigh). Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. Implement EXACTLY sub-step 6a of the DUAL-APPROVED plan `.orchestration/runs/W2-P6/plan.md`, building on the committed builder (cbf3529). First a SHORT plan, then implement code+tests. Codex + Opus dual-review (both must APPROVE). SCOPE = Zbb fermion-mixing + T010/T011 ONLY (Higgs-LFV = 6b separate).

BUILD (per the approved plan, minimal-RS Casagrande ZMA — NOT a fake full tower):
- REMOVE the hard guard in `flavor_catalog_constraints/rs_ew_builder.py` (~lines 53-54: `if include_fermion_kk_mixing: raise ValueError("...deferred beyond Phase 3a")`) and implement the real path.
- New Zbb fermion-mixing helper/dataclass (e.g. in quarkConstraints/rs_ew_couplings.py or a sibling) computing the minimal non-custodial Casagrande pieces:
  `delta g_L^b,ferm = +(m_b^2 / (2 Lambda_IR^2)) * B_d`, `delta g_R^b,ferm = -(m_b^2 / (2 Lambda_IR^2)) * B_Q`,
  with `B` ∋ `1/(1-2c)*(1/F^2 - 1 + F^2/(3+2c))` + the light-generation Yukawa-ratio sums:
  g_L uses `sum_i |Y_d[2,i]|^2/|Y_d[2,2]|^2` with `c_d[i]`; g_R uses `sum_i |Y_d[i,2]|^2/|Y_d[2,2]|^2` with `c_Q[i]`.
  (Add a comment pinning the L↔R B_d/B_Q cross-assignment.) Use `Lambda_IR = spectrum.lambda_ir_gev` (NOT kk_ew_mass/x1); M_KK in Casagrande = geometric Lambda_IR. Inputs: `bulk_state.c_Q[2]/c_d[2]/F_Q/F_d/Y_d_bulk_basis`, `m_b = fit_result.masses_down[2]`.
- MERGE the fermion piece into `rs_ew_couplings.z_delta_g_L_d[2,2]` and `z_delta_g_R_d[2,2]` (add to the existing gauge-NC piece), atomically, with `fermion_kk_mixing_included=True` metadata. The `m_t^2/M_KK^2` top-partner/custodial Zb_L term is NEEDS-HUMAN: do NOT compute it; flag `custodial_toppartner_zbL_needs_human=True`.
- Gate via `include_fermion_kk_mixing=True` (opt-in); fail loud if bulk_state Y_d/c_Q/c_d or m_b absent.

REWIRE T010/T011: they already read `z_delta_g_L/R_d[2,2]` (now gauge+fermion); REMOVE the old PARTIAL/NEEDS-HUMAN-for-fermion-KK wording, replace with `minimal_rs_tree_complete=True` + `custodial_variant_needs_human=True` (custodial/top-partner Zb_L still flagged). SM anchors + zpole observable formulas UNCHANGED.

6a GATE TESTS (recompute independently):
- alignment / `m_b->0` ⇒ fermion piece = 0 (only gauge piece remains).
- `1/M_KK^2` scaling of the fermion piece (double Lambda_IR ⇒ ÷4).
- Casagrande SIGN: `delta g_L^b,ferm > 0`, `delta g_R^b,ferm < 0` (for IR-composite b); state the magnitude (~per-mille-ish × m_b^2/Lambda_IR^2).
- universal-c / SM-limit ⇒ T010/T011 recover committed SM (gauge δg=0 AND fermion piece consistent).
- T010/T011 diagnostics now show minimal_rs_tree_complete + custodial_variant_needs_human; the top-partner term is NOT silently added.
- determinism; finite; no double-count of the gauge piece.

CONSTRAINTS: physics ONLY via adapters; ConstraintResult numeric fields real finite floats; touch the builder/coupling module + T010/T011 + tests. `python -m pytest tests/ -q` stays green + new tests.

OUTPUT (<=16 lines): short plan; files; ACTUAL numbers (delta g_L/R^b,ferm sample + sign, m_b->0=0, scaling, T010/T011 SM-limit); confirm top-partner NEEDS-HUMAN not computed; pytest counts. End with: P6A-DONE.
