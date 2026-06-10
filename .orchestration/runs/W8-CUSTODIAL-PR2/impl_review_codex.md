No blocking findings.

Verified against the real diff:
- Negative `Delta T` is not fabricated: `top_partner_loop_t_sign=-1` raises unless paired with finite numeric override.
- Carena `delta g_L^b` is added directly to `z_delta_g_L_d[2,2]`; no `s_z/g_Z` factor, no RH loop insertion.
- Core loop rejection is now gated to `minimal_rs`.
- Separate computed/Zbb/T loop flags are present; EW001 uses only the T-loop alias.
- T014 keeps RH FCNC minimal while zeroing LH off-diagonal entries in the all-gen bidoublet mode.
- Minimal/default PR1 behavior, PR1 hash `45e21a07585f7489`, and the 15 TeV custodial survival test are covered by the focused suite.

Independent oracle recomputation:
- singlet `Delta T = 0.15745209098112417`
- singlet `delta g_L^b = 0.00041018530641991833`
- bidoublet vertex `delta g_L^b = -0.0003174965326198926`

Requested test command result: `19 passed in 8.94s`.

VERDICT: APPROVE