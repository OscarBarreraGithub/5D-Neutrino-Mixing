1. BLOCKER: none.
2. SHOULD-FIX: rename plan `s_W`; it conflicts with standard `sin theta_W`. Use `eta_W`/`s_Wmix` for the W-mixing sign coefficient.
3. SHOULD-FIX: align stored names with consensus: `delta_g_W_ud_L/R`, `delta_g_W_lnu_L/R`, `charged_contact_LL`, `delta_G_F_over_G_F`; plan names `w_delta_g_*` are misnamed.
4. NIT: `x_1~2.405` is only the small-epsilon limit; implementation must use exact `solve_kk("gauge","NN", exact=True)` roots and include the zero mode in W diagonalization.
5. RESOLVED: charged spectrum formula is right for W/W': same gauge NN tower as neutral, brane-Higgs matrix `m_n^2 delta + g2^2 v^2 chi_m(1)chi_n(1)/4`.
6. RESOLVED: CC matching is physically consistent with `z_delta_g` if it uses the same shared `a_ref`; minimal RS has LH W only, RH/scalar/custodial remain out of scope.
7. RESOLVED: `delta_G_F/G_F` scheme is correct if it is the total muon-decay amplitude shift in the same `C_SM` convention and is subtracted exactly once in each consumed epsilon.
8. RESOLVED: K017-vs-K018 arbitration is correct: K017 is `R_K`, K018 is `K_l3 |V_us| f_+(0)`.
9. RESOLVED: EW002 target/status correct: rigorous minimal LH W/W' first-row shift, still SOFT non-veto.
10. RESOLVED: B009 formula is correct for minimal LH W/W' tree amplitude; nonminimal charged-Higgs/RH/scalar remains diagnostic/PARTIAL scope.
11. RESOLVED: B025 must stay PARTIAL; the vector LFU ratio is useful, but scalar/RH WET and form-factor integration are not built.
12. RESOLVED: EW003 should stay data-level unless covariance/scheme inputs exist; add diagnostics only.
13. Missing symbols today: no `quarkConstraints/rs_charged_current.py`, `RSChargedCurrentCouplings`, `rs_charged_current` extra, W-shift/contact/G_F fields.
14. Existing observables confirmed: `EW002`, `EW003`, `K017`, `K018`, `B009`, `B025`; current proxies match the plan’s rewire premise.
15. v1 expectations sound: `(m_W/M_KK)^2~7e-4`, universal lepton pieces mostly absorbed, third-family quark-doublet residuals `1e-4-1e-3`, B025 universal vector cancellation.
16. Sub-step verdict: 5a/5b/5c are independently testable; 5a should state where charged diagonalization lives if not extending `RSEWSpectrum`.
17. PLAN-NEEDS-FIXES.