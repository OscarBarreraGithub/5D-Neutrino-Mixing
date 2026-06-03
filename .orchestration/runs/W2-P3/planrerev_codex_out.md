1. BLOCKER lepton side: RESOLVED - SM lepton Z couplings are used in light-Z contacts/Wilsons; heavy-vector lepton pieces deferred/flagged.
2. BLOCKER prefactor: RESOLVED - `z_delta_g` is pinned to `s_Z * g_SM * (m_Z^2/M_KK^2) * U†diag(a-a_ref)U`, not raw overlap.
3. SHOULD-FIX gate pins: RESOLVED - full-subtraction SM limit and IR-`b_R` sign/magnitude/scaling are in the 3a no-rewire gate.
4. SHOULD-FIX rare-B: RESOLVED - B017/B018/B019 stay LFU ratios; B015 only replaces C9/C10/C9p/C10p, C7 remains untouched.
5. NITs: RESOLVED - `c_Q/c_u/c_d` come from `bulk_state`; all 5 rare-B adapter files listed; Z-pole export names match.
6. New gaps: none causing wrong/under-constraining rewire; implementation must preserve B015’s separate C7 source and LFU numerator/denominator cancellation.
7. Implementability: ready; scope, flags, gate tests, and validation are concrete.
PLAN-OK