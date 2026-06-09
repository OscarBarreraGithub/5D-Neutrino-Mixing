1. Config: use existing `ew_model` kwarg in `point_builder.build_from_rs_ew_inputs` and `rs_ew_builder.build_rs_ew_extras`; add custodial kwargs/default metadata there, not `ScanConfig` first PR.
2. Defaults: `ew_model="minimal_rs"` remains default and must be byte-identical; allow `"custodial_rs_plr"` with `qL_rep`, `tR_rep`, `bR_rep`, `protect_scope`, `bR_strategy`, `kappa_b=0.0`, `custodial_PLR_breaking_residual=False`, `include_top_partner_loops=False`.
3. `rs_ew_spectrum.py`: no `a(c)`/`Omega_n(c)` physics change; keep custodial `model_label` as metadata/validation only.
4. `rs_ew_builder.build_rs_ew_extras`: remove minimal-only rejection, validate supported models, pass custodial options into `build_rs_ew_couplings`, and return existing `rs_ew_couplings`/`kk_ew_mass_gev` slots.
5. `rs_ew_couplings.build_rs_ew_couplings`: build minimal arrays exactly as today first; save `minimal_z_delta_l_d_full` after optional Casagrande add, then in custodial mode set protected mass-basis down-left block `z_delta_l_d[:, :] = 0`.
6. Residual: if `custodial_PLR_breaking_residual` or `kappa_b != 0`, set only `z_delta_l_d[2,2] = kappa_b / L * minimal_z_delta_l_d_full[2,2]`; otherwise exact zero; hermitianize and record minimal/residual values.
7. Fermion mixing: still compute `build_rs_zbb_fermion_kk_mixing` when requested for metadata, but do not add its `delta_g_L_b` in custodial mode except through the saved minimal residual source above.
8. `b_R`: default `bR_strategy="elementary_zero"` sets `z_delta_r_d[2,2]=0`; leave other `R_d` entries minimal; future explicit strategies may add positive `delta_g_R_b`, never automatic.
9. `oblique_stu.py`: keep `minimal_rs_t_coefficient`; add custodial coefficient `-pi/(4*cW^2*L)` and gate `evaluate_rs_oblique_proxy(..., ew_model=...)`; `c_S` unchanged, `U=0`.
10. `physics_adapters/oblique_stu.py` and `EW001.evaluate`: export/pass `ew_model`; read from `rs_ew_couplings.metadata["ew_model"]` when present, else minimal.
11. `T010.py`/`T011.py` diagnostics: copy `ew_model`, reps, `protect_scope`, `bR_strategy`, residual, loop flag, omissions; stop saying custodial deferred when active; ensure harness does not tag custodial tree as minimal.
12. Omissions metadata: `SU2_R_tower=False`, `custodian_spectrum=False`, `exact_NC_mixing=False`, `BKT=False`, `top_partner_loops="deferred"` unless enabled.
13. Tests: add custodial Zbb unit tests in `test_rs_ew_phase6a_zbb_fermion_mixing.py`, oblique coefficient tests in `test_EW001.py`, builder threading/determinism tests, and minimal byte-identical regression.
14. Benchmark test: same ~15 TeV minimal point vetoed by T010 should pass/relax with custodial `z_delta_g_L_d[:, :]==0`, `z_delta_g_R_d[2,2]==0`, while S unchanged and custodial T tiny.
15. Risks: mass-basis protection scope, accidental double-zeroing before residual capture, `b_R` model ambiguity, EW001 fallback for mass-only points, minimal-path numeric regression.
16. PR1 scope: tree-level switch, Zb_L protection, b_R default-zero, oblique T swap, diagnostics/tests; defer numerical one-loop top-partner ΔT and δg_L^b loop formulas to PR2.

CUSTODIAL-PLAN-DONE.