1. PASS: offdiag `bs/bd/sd` unchanged: `-2.038738e-4-2.203657e-8j`, `3.850171e-5+7.684502e-6j`, `4.398393e-6-6.072737e-7j`; T014 pred/ratio unchanged `1.3203075908784097e-7` / `4.5527847961324475e-5`.
2. PASS: default custodial diagonals `L_d[0,0..2,2]=0`; residual path with `kappa_b=0.7` gives `b=5.31104695738613e-5`, matching `0.7/L * minimal 0.002620529475486147`; hermiticity max `0.0`.
3. PASS: EW001 plumbing works; custodial S unchanged `0.05052009670510209`, U `0`, T coeff `-0.029189036558757226` vs minimal `71.5131395689552`.
4. PASS: T010/T011 custodial diagnostics active: `custodial_variant_deferred=False`, top-partner deferred `True`, tag rigorous/no `needs_human`.
5. PASS: minimal byte identity confirmed: arrays/metadata equal, T010/T011/T014/EW001 result triples equal, config hash `45e21a07585f7489`, no default `ew_model` payload key.
6. PASS: EW001 manifest has `rs_ew_couplings`; forbidden extras are only lepton/charged-current/higgs/semileptonic slots.
7. PASS: 15 TeV benchmark asserts intended relaxation: minimal T010 fails ratio `1.135208212413654`; custodial T010/T011 pass with `delta_g_L_b=0j`, `delta_g_R_b=0j`.
8. PASS: full requested suite green: `1723 passed, 1 skipped in 700.54s (0:11:40)`.
Defects: none found.
CUSTODIAL-IMPL-REVIEW: APPROVE