# W2 PHASE 6 — SUB-STEP 6b: Higgs-LFV (codex, gpt-5.x xhigh). Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. Implement EXACTLY sub-step 6b of the DUAL-APPROVED plan `.orchestration/runs/W2-P6/plan.md`, building on 6a (committed 85c06ea) + the 4a lepton builder. First a SHORT plan, then implement code+tests. Codex + Opus dual-review (both must APPROVE). SCOPE = Higgs-LFV (T018/T019/T020) ONLY. This is the LAST P6 sub-step.

BUILD (per the approved plan):
- New `quarkConstraints/rs_higgs_yukawas.py` with frozen `RSHiggsYukawaCouplings`; add `rs_higgs_yukawas` to `KNOWN_EXTRA_KEYS`; add `include_higgs_yukawas=True` kwarg to `build_rs_ew_extras` (populate from `lepton_mass_basis_couplings`).
- Higgs-LFV Yukawa matrix (lepton analog of Casagrande Δg_h): `Y_ij = -[(m_i/v)(delta_E)_{ij} + (delta_L)_{ij}(m_j/v)]` where `delta_L = x_e U_e_L† diag[B(c_E)] U_e_L x_e`, `delta_E = x_e U_e_R† diag[B(c_L)] U_e_R x_e`, `x_e = diag(m_e,m_mu,m_tau)/Lambda_IR`, `B(c)=1/(1-2c)(1/F²-1+F²/(3+2c))`. Lepton masses from `yukawa.constants.LEPTON_MASSES`; `Lambda_IR = spectrum.lambda_ir_gev`; v from the point. Expose `higgs_yukawa_matrix` (the attribute the adapter reads via _first_present_attr), units dimensionless, matching_assumption, includes_fermion_kk_mixing, diagnostics.
- **v1 HONESTY**: with the diagonal charged-lepton fit (U_e_L=U_e_R=I, Y_E_bar_matrix=diag, scalar/broadcast c_L), both delta_L and delta_E are DIAGONAL ⇒ Y_ij off-diagonal = EXACTLY 0 ⇒ T018/T019/T020 evaluate ZERO (rigorous-tree, non-vetoing). Nonzero tree Higgs-LFV requires non-diagonal Y_E/rotations (future anarchic lepton fit). So in v1 these are rigorous-but-zero; keep an honest note that generic RS Higgs-LFV needs a non-diagonal charged-lepton structure (deferred). Do NOT fake a nonzero bound.

REWIRE T018/T019/T020 (h→eμ/eτ/μτ): change the required extra to `rs_higgs_yukawas`; pass it to the higgs-LFV adapter (`higgs_lfv_branching_fraction_with_proxy`; core `h_lfv_branching_fraction_with_proxy`); KEEP the same SM=0 and BR formula `m_h(|Y_ij|²+|Y_ji|²)/(8π Γ_h)`. Graceful degradation: absent rs_higgs_yukawas ⇒ non-vetoing evaluated=False missing_extra.

6b GATE TESTS:
- Diagonal-v1: build via `build_from_rs_ew_inputs(...,include_higgs_yukawas=True)` ⇒ Y_h off-diagonal EXACTLY 0 ⇒ T018/T019/T020 predicted=0 (non-vetoing, rigorous-tree). 
- LFV-LIVE: a non-diagonal-U_e / anarchic-Y_E toy ⇒ nonzero Y_ij ⇒ nonzero BR matching the higgs_lfv core (cross-check independent of the adapter); v²/M_KK² scaling.
- Alignment/m_ℓ→0 sanity; absent-path non-vetoing; SM=0 unchanged; determinism/finite/immutable.
- Replace proxy-only tests — enumerate, no silent coverage loss.

CONSTRAINTS: physics ONLY via adapters; ConstraintResult numeric fields real finite floats; touch the new module + builder + T018/T019/T020 + tests. `python -m pytest tests/ -q` stays green.

OUTPUT (<=14 lines): short plan; new module + API; diagonal-v1 Y_h=0 + T018/19/20=0 + LFV-live-toy nonzero Y_ij/BR + scaling + absent; test-count change; pytest counts. End with: P6B-DONE.
