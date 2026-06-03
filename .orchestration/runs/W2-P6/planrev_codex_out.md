1. BLOCKER: none.
2. SHOULD-FIX: none found.
3. NIT: spell Casagrande sums explicitly: `g_L` fermion term uses `|Y_d[2,i]|^2/|Y_d[2,2]|^2` with `c_d[i]`; `g_R` uses `|Y_d[i,2]|^2/|Y_d[2,2]|^2` with `c_Q[i]`.
4. MACHINERY: RESOLVED. Repo lacks normalized fermion-KK profiles and full zero/KK fermion mass diagonalization; ZMA/Casagrande matching is the honest tractable default.
5. Inputs: available, but as arrays/fit fields: `bulk_state.c_Q[2]`, `c_d[2]`, `F_Q/F_d`, `Y_d_bulk_basis`; `m_b=fit_result.masses_down[2]`; `Lambda_IR=spectrum.kk_ew_mass_gev/spectrum.gauge_roots_x[0]` or `spectrum.lambda_ir_gev`.
6. Zbb formula: RESOLVED. Signs/structure match Casagrande ZMA: `delta g_L^b,ferm > 0`, `delta g_R^b,ferm < 0`, with `M_KK` in that paper equal to geometric `Lambda_IR`, not physical first gauge mass.
7. LEP tension: minimal RS sign does not solve `A_FB^b`; it gives the wrong direction for the positive `g_R^b` shift, consistent with Casagrande and Agashe custodial discussion.
8. Top-partner/custodial Zb_L: NEEDS-HUMAN. Deferral is correct; representation/P_LR/exotic spectrum choices are absent and can dominate `Zb_L`.
9. T010/T011 scope: RESOLVED. “FULL” is only defensible as documented `minimal_rs_tree_complete`; custodial/top-partner variants must remain human/model-dependent diagnostics.
10. Higgs-LFV zero-v1: RESOLVED. Current diagonal charged-lepton fit has `Y_E_bar_matrix=diag`, `U_e_L=U_e_R=I`, so tree `Y_h_mass` off-diagonals are exactly zero and non-vetoing.
11. Higgs formula: correct chirally suppressed structure: `Y_ij = -[(m_i/v) delta_E_ij + delta_L_ij (m_j/v)]`, matching Casagrande `Delta g_h`.
12. Symbol check: `higgs_lfv_branching_fraction_with_proxy` exists in the adapter; core name is `h_lfv_branching_fraction_with_proxy`.
13. Missing/new symbols by design: `rs_higgs_yukawas`, `RSHiggsYukawaCouplings`, and `quarkConstraints/rs_higgs_yukawas.py` do not exist yet; add `rs_higgs_yukawas` to `KNOWN_EXTRA_KEYS`.
14. Lepton mass note: charged-lepton masses are in `yukawa.constants.LEPTON_MASSES`, not stored on current `RSLeptonMassBasisCouplings`; add/store or derive explicitly.
15. Sub-step verdict: 6a OK; 6b OK if diagonal-v1 zero and non-diagonal toy nonzero tests are included.
16. Sources checked: local repo files; Casagrande arXiv https://arxiv.org/abs/0807.4937; Agashe et al. https://arxiv.org/abs/hep-ph/0605341.
PLAN-OK