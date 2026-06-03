1. BLOCKER: None.
2. SHOULD-FIX: None.
3. NIT: None; builder is opt-in, key registered, lepton inputs fail loudly, RH W arrays are zero/statused, extras are immutable/finite/deterministic, and no constraints were rewired.
4. W mass verdict: `g2=0.651731`, `v=246.2197`, so unmixed `g2*v/2=80.2345 GeV`; diagonalized `m_w=77.5067 GeV`, expected from EWSB-KK mixing, not a `v=174` bug.
5. W' verdict: test uses `Lambda_IR=1224.235123`, `x1=2.450509663813748`, so neutral bare `kk_ew_mass=x1*Lambda_IR=3000.000`; charged eigen-`W'=3071.5336` from the same tower plus mixing, not the `Lambda_IR=3000 => 7351.529` mismatch.
6. `eta_W=-1.0`; eigenvector projection `-6.6792e-4 < 0`; sign test is present and matches neutral `s_z=-1`.
7. Recomputed `delta_g_W_ud_L[0,1]=-8.1699158137e-7+3.4880867618e-6j`; `epsilon[0,1,0]=2.6427715370e-5-1.2332357025e-5j`.
8. `delta_G_F/G_F=-1.8334382815e-5`; epsilon subtracts it exactly once, with explicit tests against double subtraction and second `1/M_KK^2`.
9. Universal-c check: max `delta_g_W_ud_L=0`, max `delta_g_W_lnu_L=0`, max `epsilon=3.3881e-21`.
10. `pytest tests/ -q`: `1671 passed, 1 skipped in 781.33s`.

P5A-OK