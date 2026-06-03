# W2 PHASE 3 — SUB-STEP 3a (codex, gpt-5.x xhigh). Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. Implement EXACTLY sub-step 3a of the DUAL-APPROVED plan `.orchestration/runs/W2-P3/plan.md` (and design `.orchestration/rs_ew_sector_design_CONSENSUS.md` §3,§4,§6). First a SHORT plan, then implement code+tests. Codex + Opus will dual-review (both must APPROVE). **3a does NOT rewire any constraint** — builder + quark Z-matrices + contacts + the 3a gate tests ONLY (rewiring is 3b/3c/3d).

BUILD (per the approved plan):
- Add `rs_ew_spectrum`, `rs_ew_couplings`, `rs_semileptonic_wilsons` to `KNOWN_EXTRA_KEYS` (fail-loud `make_point` preserved).
- `build_from_rs_ew_inputs(quark_fit_result, Lambda_IR, k, n_gauge_modes, ew_model="minimal_rs", ...)`: builds `RSEWSpectrum` (from `quarkConstraints/rs_ew_spectrum.py`), sets `kk_ew_mass_gev=spectrum.kk_ew_mass_gev`; reads c-values from `quark_fit_result.bulk_state.c_Q/.c_u/.c_d` and rotations `.U_L_u/.U_L_d/.U_R_u/.U_R_d`.
- Quark mass-basis Z matrices: `z_delta_g_A^f = s_Z * g_A^{f,SM} * (m_Z^2/M_KK^2) * [U^dagger diag(a(c_f)-a_ref) U]` (NOT raw overlap), with u_L:(U_L_u,c_Q), d_L:(U_L_d,c_Q), u_R:(U_R_u,c_u), d_R:(U_R_d,c_d); store dimensionless additive shifts in the zpole convention (NO extra g_Z); `z_total = g^SM·I + z_delta_g`. Must be Hermitian.
- Neutral contacts (GeV^-2): `(g_Z^2/m_Z^2)[(g_qA^SM δ_ij + z_delta_g_qA)(g_lB^SM δ_ab) − g_qA^SM g_lB^SM δ_ij δ_ab] + Σ_V g_Vq g_Vl/M_V^2`, using SM charged-lepton couplings g_l^SM (charged-lepton δg_l=0 in Phase 3; heavy-vector lepton g_Vl deferred to Phase 4 — flag `includes_heavy_neutral_lepton=False`). 
- `rs_semileptonic_wilsons.{b_to_s_ll,b_to_d_ll,s_to_d_ll,c_to_u_ll}` C9/C10/C9p/C10p built from the contacts. IMPORTANT: these are GeV^-2 / dimensionless additive Wilsons carrying their OWN 1/m_Z^2·1/M_V^2 — they must be consumable WITHOUT `_wilson_prefactor` and WITHOUT any second 1/M_KK^2 (3d will inject them additively at the core's wilson-value consumption point).

3a GATE TESTS (run BEFORE any rewiring; recompute independently, not call-through):
- Schema fail-loud (unknown extra raises); Hermiticity of z_total/z_delta_g matrices.
- Full-subtraction SM limit: a(c)=a_ref / identity rotations ⇒ all z_delta_g=0, all contacts/Wilsons=0.
- Per-species family-universal c ⇒ zero off-diagonal (FCNC) but nonzero diagonal shifts allowed.
- IR b_R sign+magnitude: small c_{bR} ⇒ z_delta_g_R^b < 0, |·|~1e-3 at M_KK~3 TeV, scaling 1/M_KK^2.
- Contact units (GeV^-2); no mutation of ParameterPoint.extras; determinism; finite.
- Independent manual recompute of one b→sμμ, one s→dee, one c→uμμ Wilson (from spectrum+rotations, not the new builder).

CONSTRAINTS: ADD a new builder module/functions + new tests; do NOT modify the 103 constraints or existing cores beyond adding the extras to KNOWN_EXTRA_KEYS. `python -m pytest tests/ -q` stays green + new tests pass. Numeric result fields real finite floats; complex matrices live in extras/diagnostics (not ConstraintResult numeric fields).

OUTPUT (<=16 lines): short plan; new module path + public API; the ACTUAL asserted numbers (a sample z_delta_g_R^b @3TeV, the SM-limit=0 check, a sample C9/C10, the independent-recompute match); pytest counts. End with: P3A-DONE.
