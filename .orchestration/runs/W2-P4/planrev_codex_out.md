1. BLOCKER: none found for the light-Z Phase 4 path.
2. SHOULD-FIX: RESOLVED-as-EW-universal-a_ref. Leptons should share the same universal subtraction `a_ref = spectrum.a(DEFAULT_A_REF_C=0.65)`, not a separate lepton reference; rename/diagnose it as EW-universal, not quark-only.
3. SHOULD-FIX: RESOLVED-as-light-Z-only-for-Phase-4. No `Z'/gamma'` neutral mass-eigenvector layer exists in current `RSEWSpectrum`; keep `includes_heavy_neutral_exchange=False` unless Phase 4 first builds `neutral_vector_mixing`/heavy couplings. Do not claim heavy-neutral FULL status when false.
4. SHOULD-FIX: pin the `nunu` current convention: existing cores use the LH active-neutrino current, `X_NP=C/g_SM^2`, no `_wilson_prefactor`; Majorana vs Dirac gives no extra rate factor for the light active `nunu` observables, but this needs an explicit validation note/test.
5. NIT: missing/misnamed symbols: none if `c_L,c_E,c_N,M_N` are read from `YukawaResult.params`; they are not top-level dataclass fields. `Y_E_bar/Y_N_bar` are vectors, so Phase 4 must build the diagonal/matrix forms.
6. OK: `compute_all_yukawas`, `YukawaResult`, `Y_N_matrix`, `f_L/f_E/f_N/f_N_UV`, and `neutrinos.neutrinoValues.get_pmns()` exist; charged-lepton fit is diagonal, so `U_e_L=U_e_R=I` and no second PMNS rotation is correct.
7. OK: lepton Z formula is correct: same `s_Z g_A^SM m_Z^2/M_KK^2 U†diag(a-a_ref)U`; charged leptons use `c_L/c_E`, active neutrinos use `c_L` not singlet `c_N`; Hermitian projection is right.
8. OK: contact completion must supersede the 3a `neutral_contacts` tensor. Current 3a is quark-NP times lepton-SM only, so replacing in the same path avoids double-counting.
9. OK: rewire map is sound; note `K019/K020/K021` live under `secondary/kaon`. `L002/L009` and `L003/L004/L005` must remain PARTIAL as stated.
10. Sub-step verdict: 4a/4b/4c/4d are independently testable, with 4a prerequisite for all rewires; add the arbitration/status clarifications above before implementation.

PLAN-NEEDS-FIXES.