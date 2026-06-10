# W8-CUSTODIAL-PR2 Implementation Summary

Files changed:
- `quarkConstraints/rs_ew_couplings.py`
- `quarkConstraints/oblique_stu.py`
- `flavor_catalog_constraints/primary/top_higgs_ew/EW001.py`
- `flavor_catalog_constraints/primary/top_higgs_ew/T010.py`
- `flavor_catalog_constraints/primary/top_higgs_ew/T011.py`
- `flavor_catalog_constraints/secondary/top_higgs_ew/T014.py`
- `flavor_catalog_constraints/rs_ew_builder.py`
- `flavor_catalog_constraints/point_builder.py`
- `tests/test_rs_ew_custodial_pr2.py`

No change was needed in `flavor_catalog_constraints/physics_adapters/oblique_stu.py`;
PR2 extended existing oblique function signatures without adding a new exported symbol.

Implemented:
- Threaded PR2 custodial knobs through `point_builder`, `rs_ew_builder`, and `build_rs_ew_couplings`, including `top_partner_loop_delta_t_override`.
- Allowed top-partner loops only for `ew_model="custodial_rs_plr"`; `minimal_rs` still raises.
- Added Carena leading singlet and bidoublet-vertex loop proxy with separate Zbb and EW001 T flags.
- Enforced no fabricated negative T: `top_partner_loop_t_sign=-1` without finite numeric override raises.
- Implemented missing-sign deterministic behavior: magnitudes computed, no matrix/EW001 loop application.
- Added Carena `delta g_L^b` directly to `z_delta_g_L_d[2,2]` in the repo additive convention; `z_delta_g_R_d` is untouched.
- Added representation-aware custodial FCNC mode with LH off-diagonal PLR zeroing and RH minimal off-diagonals retained.
- Threaded valid loop `Delta T` into EW001 as a separate `t_loop_prediction`.
- Added T010/T011/T014 diagnostics passthrough for loop and FCNC metadata without changing default PR1 rigorous tagging.

Tests added:
- `test_top_partner_loop_proxy_singlet_numeric_oracle_3tev`
- `test_t010_singlet_zbb_pseudo_observable_oracle`
- `test_bidoublet_vertex_numeric_oracle_no_fake_negative_t`
- `test_negative_t_requires_numeric_override_and_threads_to_ew001`
- `test_top_partner_loop_missing_sign_computes_magnitudes_but_does_not_apply`
- `test_minimal_rs_rejects_top_partner_loops_and_scan_hash_stays_pinned`
- `test_ew001_adds_loop_t_as_separate_term_and_deferred_is_tree_only`
- `test_t010_t011_consume_looped_zbb_matrix_and_report_components`
- `test_custodial_fcnc_all_gen_bidoublet_zero_lh_mode_t014_rh_minimal_oracle`
- `test_custodial_fcnc_residual_kappa_over_l_keeps_rh_minimal`
- `test_honest_omission_flags_present_for_deferred_and_loop_computed_branches`
- `test_top_partner_mass_ratio_metadata_uses_physical_mkk_denominator`

Oracle confirmations:
- Singlet benchmark:
  - `F_Q3 = 0.2656322304046161`
  - `F_u3 = 0.5656854250202850`
  - `Y_t_eff = 3.308347352802294`
  - `Delta T_loop = +0.15745209098112428`
  - `delta g_L_b_loop = +0.00041018530641991833`
  - `delta g_R_b_loop = 0.0`
- Bidoublet vertex:
  - `delta g_L_b_bidoublet_vertex = -0.0003174965326198927`
  - `top_partner_delta_t_loop_applied = 0.0`
  - `top_partner_loop_numerics_included = False`
- T010 pseudo-observable:
  - `g_left = -0.42242314802691344`
  - `g_right = +0.07716666666666666`
  - `R_b^0 = 0.21530246427000885`
  - `A_b = 0.9354140642148572`
  - selected observable `R_b^0`
  - ratio `1.4680590546247065`
- T014 RH-minimal oracle:
  - selected channel `bs`
  - predicted `6.369734635064622e-08`
  - ratio `2.196460218987801e-05`

Preservation checks:
- Minimal/default scan hash remains `45e21a07585f7489`.
- 15 TeV custodial benchmark still survives:
  - T010 passes with ratio `0.9960141559711674`
  - T011 passes with ratio `0.0`
- Focused regressions:
  - `tests/test_rs_ew_custodial_pr2.py`: `12 passed`
  - PR1 + PR2 + EW001 + Zbb focused set: `35 passed`
  - Existing T010/T011/T014 set: `29 passed`
- Full suite command after final validation patch:
  - `source ~/.bashrc && conda activate ising_bootstrap && export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH" && python -m pytest -q 2>&1 | tail -30`
  - Result: `1746 passed, 1 skipped in 856.91s (0:14:16)`

IMPL-READY
