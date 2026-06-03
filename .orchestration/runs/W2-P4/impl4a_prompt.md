# W2 PHASE 4 — SUB-STEP 4a (codex, gpt-5.x xhigh). Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. Implement EXACTLY sub-step 4a of the DUAL-APPROVED plan `.orchestration/runs/W2-P4/plan.md`, building on the 3a quark builder (committed; `quarkConstraints/rs_ew_couplings.py`, `rs_semileptonic_wilsons.py`, `flavor_catalog_constraints/rs_ew_builder.py`). First a SHORT plan, then implement code+tests. Codex + Opus dual-review (both must APPROVE). **4a does NOT rewire any constraint** — lepton builder + lepton/neutrino Z-matrices + completed contacts + Wilson bundle extension + 4a gate tests ONLY.

BUILD (per the approved plan):
- Lepton inputs from `yukawa.compute_yukawas.compute_all_yukawas` → `YukawaResult` (read c_L,c_E,c_N,M_N from `.params`; Y_E_bar,Y_N_bar are (3,) vectors → use `np.diag`; Y_N_matrix; f_L/f_E/f_N/f_N_UV) + `neutrinos.neutrinoValues.get_pmns()`. Charged-lepton fit is DIAGONAL ⇒ set `U_e_L=U_e_R=I_3` (no second PMNS rotation); current basis IS the charged-lepton mass basis.
- `lepton_mass_basis_couplings` typed extra `RSLeptonMassBasisCouplings`: c_L[3],c_E[3],c_N[3], overlaps, `np.diag(Y_E_bar)`, Y_N_bar (vector), Y_N_matrix, Y_N_bar_matrix=2k·Y_N_matrix, pmns, U_e_L=U_e_R=I, U_nu basis metadata, lfv_dipole_spurion=Y_N_bar_matrix·Y_N_bar_matrix†, matching-status diagnostics.
- Lepton/neutrino Z matrices: `z_delta_g_{L,R}_e` (charged leptons, c_L/c_E via U_e=I) and `z_delta_g_L_nu` (active ν, LH, uses c_L NOT singlet c_N), same form `s_Z·g_A^SM·(m_Z²/M_KK²)·U†diag(a(c)−a_ref)U`, dimensionless additive zpole convention (no extra g_Z), Hermitian. **a_ref is the SHARED EW-universal `spectrum.a(DEFAULT_A_REF_C=0.65)`** (name/diagnose it EW-universal, not quark-only).
- Contacts: REPLACE 3a's quark-NP×lepton-SM tensor IN-PLACE in the same `neutral_contacts` path with the FULL `(g_q^SM δ_ij+δg_q)(g_l^SM δ_ab+δg_l) − g_q^SM g_l^SM δ_ij δ_ab` for charged leptons AND active neutrinos (GeV^-2). NO double-count (supersede, don't add). Keep `includes_heavy_neutral_exchange=False` (Z'/γ' layer deferred; do NOT claim heavy FULL).
- Wilson bundle: extend with off-diagonal `lfv_llqq` blocks (charged-lepton 3x3 pairs) AND `_nunu` blocks (`b_to_s_nunu`, `s_to_d_nunu`) mapped `X_NP=C/g_SM²` (LH active-ν current), never via `_wilson_prefactor`/second 1/M_KK². (The legit adapter G_F/α/λ factor stays for C9/C10.)

4a GATE TESTS (recompute independently, NO rewiring):
- universal-c (lepton a(c)=a_ref) ⇒ z_delta_g_e=z_delta_g_nu=0 ⇒ LFV contacts/Wilsons=0; family-universal lepton c (≠a_ref) ⇒ diagonal nonzero but off-diagonal (LFV)=0.
- Hermiticity of lepton Z matrices; no g_Z double-mult; no second 1/M_KK²; contacts GeV^-2.
- No double-count with 3a: confirm the same-flavor quark×lepton Wilson from the completed contact equals the 3a value when δg_l=0 (superseded tensor reduces to 3a in that limit).
- νν: Majorana=Dirac for light active νν (universal c_L ⇒ U_nu†(a·I)U_nu=a·I, phases cancel) — assert.
- Independent manual recompute of one lepton z_delta_g and one lfv_llqq / one _nunu Wilson.
- determinism; finite; extras immutable.

CONSTRAINTS: ADD lepton builder + extend rs_ew_couplings/rs_semileptonic_wilsons + new tests; modify ONLY the builder/coupling modules (+ KNOWN_EXTRA_KEYS if needed) — do NOT rewire constraints, do NOT touch the 103 constraint files. `python -m pytest tests/ -q` stays green + new tests pass.

OUTPUT (<=16 lines): short plan; new/extended modules + API; ACTUAL asserted numbers (a sample z_delta_g_e/nu, an lfv_llqq + a _nunu Wilson, the universal-c=0 + family-universal-offdiag=0 + 3a-consistency + Majorana=Dirac checks); pytest counts. End with: P4A-DONE.
