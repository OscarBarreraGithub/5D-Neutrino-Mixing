1. BLOCKER: 3a contacts/Wilsons under-specified: plan builder lists quark inputs only, but design requires `lepton_inputs` and contacts use lepton couplings; otherwise rare rewires undercount heavy neutral/lepton-side terms.
2. BLOCKER: 3a must pin stored `z_delta_g = prefactor * U†diag(a(c)-a_ref)U`; `U†diag(...)U` basis is correct and Hermitian, but raw `a(c)-a_ref` is not the coupling shift.
3. SHOULD-FIX: move full-subtraction SM-limit and Zbb sign/magnitude pins into the “no rewiring” gate, not just end validation.
4. SHOULD-FIX: B017/B018/B019 are LFU ratio proxies, not LFV, but numerator/denominator semantics must be preserved; a lepton-universal Wilson mostly cancels.
5. SHOULD-FIX: B015 also carries an existing C7 proxy; Phase 3 should state only C9/C10/C9p/C10p are replaced and C7 is unchanged.
6. OK: 3b Z-pole structure is sound; T010/T011 must keep PARTIAL/NEEDS-HUMAN until phase 6; T012 covers charm, and T013 has no standalone file.
7. OK: 3c FCNC-Z rewire is sound; SM-zero/direct-width policy stays unchanged.
8. OK with fixes: 3d target scope is right: B005-B008/B015-B019/B021, K006/K008/K009/K010/K012, C004/C005/C007; C006/C008 are LFV and correctly deferred; K013 is radiative and correctly excluded.
9. Missing/new symbols: `build_from_rs_ew_inputs`, `rs_ew_couplings`, `rs_semileptonic_wilsons`, and `rs_ew_spectrum` as a `ParameterPoint` extra do not exist yet.
10. Naming caveat: use adapter names `zpole_sm_couplings`, `zpole_shifted_couplings`, `zpole_evaluate_quark`; core has `sm_couplings`, `shifted_couplings`, `evaluate_quark_pseudo_observables`.
11. Verified symbols: `*_with_proxy`, `zpole_down_fcnc_branching_fraction_from_couplings`, `compute_rare_{b,kaon,charm}_dilepton_wilsons`, rare cores, and quark fit `U_L/U_R/c_*` outputs exist.
12. PLAN-NEEDS-FIXES.