# W8-CUSTODIAL-PR2 implementation plan

This is an implementation plan only. Do not write production code in this run.

## Grounding and current code anchors

- Approved custodial prescription remains unchanged: all-generation `Q_L` bidoublets in `(2,2)_{2/3}` and mass-basis treatment are required by `.orchestration/runs/CUSTODIAL-RESEARCH/CORRECTED_PRESCRIPTION.md:13`; protected `b_L` has `T^3_L=T^3_R=-1/2` at `.orchestration/runs/CUSTODIAL-RESEARCH/CORRECTED_PRESCRIPTION.md:16`; tree `Zb_L` protection zeros both gauge-profile and bottom-mass admixture pieces at `.orchestration/runs/CUSTODIAL-RESEARCH/CORRECTED_PRESCRIPTION.md:18`; oblique custodial tree `T` is `-pi/(4 c_W^2 L)` at `.orchestration/runs/CUSTODIAL-RESEARCH/CORRECTED_PRESCRIPTION.md:27`; one-loop top-partner effects must be flagged/computed honestly at `.orchestration/runs/CUSTODIAL-RESEARCH/CORRECTED_PRESCRIPTION.md:39`; omissions must remain visible at `.orchestration/runs/CUSTODIAL-RESEARCH/CORRECTED_PRESCRIPTION.md:49`.
- Current neutral-current core accepts `include_top_partner_loops` but rejects it at `quarkConstraints/rs_ew_couplings.py:491`. This second rejection site must be replaced, not just the wrapper rejection in `flavor_catalog_constraints/rs_ew_builder.py:78`. The validated loop path is allowed only for `ew_model=="custodial_rs_plr"`; `minimal_rs` keeps raising if top-partner loops are requested.
- Current core builder line anchors: `build_rs_ew_couplings` signature starts at `quarkConstraints/rs_ew_couplings.py:467`; `a(c)` profiles are built at `quarkConstraints/rs_ew_couplings.py:528`; tree `z_delta_l_d` is built at `quarkConstraints/rs_ew_couplings.py:581`; minimal down shifts are snapshotted at `quarkConstraints/rs_ew_couplings.py:595`; `_apply_custodial_rs_plr_proxy` is called at `quarkConstraints/rs_ew_couplings.py:599`.
- PR1 custodial proxy zeroes protected diagonal down-left entries at `quarkConstraints/rs_ew_couplings.py:952`, applies optional `kappa_b/L` residual at `quarkConstraints/rs_ew_couplings.py:958`, zeroes elementary `b_R` at `quarkConstraints/rs_ew_couplings.py:965`, records top-partner loops as deferred at `quarkConstraints/rs_ew_couplings.py:999`, and records FCNC as `deferred_PR2_off_diagonal_kept_minimal` at `quarkConstraints/rs_ew_couplings.py:1004`.
- Current tree `z_delta` metadata says `s_Z * g_A_SM * (m_Z^2/M_KK^2) * U^\dagger diag(a(c)-a_ref) U` at `quarkConstraints/rs_ew_couplings.py:745`; the implementation is `_z_delta = _hermitian(s_z * g_sm * scale * a_mass_basis)` at `quarkConstraints/rs_ew_couplings.py:1137`. Keep this untouched for `minimal_rs`.
- Current `a(c)` convention is the exact tower sum `sum_n (M_KK^2/m_n^2) chi_n(1) Omega_n(c)` in `quarkConstraints/rs_ew_spectrum.py:7`, implemented by `a_terms` at `quarkConstraints/rs_ew_spectrum.py:817`, with subtraction in `a_with_diagnostics` at `quarkConstraints/rs_ew_spectrum.py:890`.
- Existing oblique formulas are documented at `quarkConstraints/oblique_stu.py:18`; `minimal_rs_t_coefficient` and `custodial_rs_plr_t_coefficient` live at `quarkConstraints/oblique_stu.py:151` and `quarkConstraints/oblique_stu.py:165`; `rs_minimal_oblique_proxy` builds `S,T,U` at `quarkConstraints/oblique_stu.py:198`; `evaluate_rs_oblique_proxy` compares the point to the fit at `quarkConstraints/oblique_stu.py:250`.
- EW001 resolves `ew_model` from `rs_ew_couplings.metadata` at `flavor_catalog_constraints/primary/top_higgs_ew/EW001.py:349`, calls `evaluate_rs_oblique_proxy` inside a `try/except ValueError` at `flavor_catalog_constraints/primary/top_higgs_ew/EW001.py:412`, and currently passes no loop information.
- T010 and T011 consume final additive `z_delta_g_L/R_d[2,2]` entries directly: T010 at `flavor_catalog_constraints/primary/top_higgs_ew/T010.py:584`, T011 at `flavor_catalog_constraints/primary/top_higgs_ew/T011.py:761`. `quarkConstraints/zpole.py:202` defines SM `g_L=T3-Q s_eff^2`, and `quarkConstraints/zpole.py:220` adds `delta_g_left` directly.
- T014 consumes both off-diagonal `z_delta_g_L_d` and `z_delta_g_R_d` entries at `flavor_catalog_constraints/secondary/top_higgs_ew/T014.py:450` and `flavor_catalog_constraints/secondary/top_higgs_ew/T014.py:456`; any PR2 FCNC test must account for RH down FCNC unless an explicit RH-zero mode is added.
- Carena, Ponton, Santiago, Wagner, Phys. Rev. D76:035006, arXiv:hep-ph/0701055 is the loop source. The arXiv abstract states that loop corrections to both `T` and `Zb_L b_L` are calculable in custodial models, and the TeX source confirms the leading singlet `Delta T` equation, leading singlet and bidoublet `delta g_{b_L}` equations, and full Appendix formulas `T:eq` / `Zbb:eq`.

## Physics formulas and repo conventions

Keep the repo's two electroweak vev conventions explicit:

- Oblique tree proxy uses `v_246 = quarkConstraints.oblique_stu.DEFAULT_HIGGS_VEV_GEV = 246.21965`.
- Fermion mass/Yukawa overlap convention in the quark fit uses `V_EWSB = 174.0` from `warpConfig/baseParams.py:3`; tests already encode `m_q[i] = 2 * V_EWSB * F_Q[i] * Y[i,i] * F_q[i]` in `tests/test_diagnostics.py:33`.

Keep the mass convention explicit:

- `spectrum.kk_ew_mass_gev` is the physical first electroweak gauge KK mass `M_KK = m_1`.
- `spectrum.lambda_ir_gev` is the geometric IR scale `Lambda_IR`.
- The physical gauge-root convention is `M_KK ~= 2.450509663813736 * Lambda_IR`, already used by `tests/rs_ew_phase3b_helpers.py:14`.
- `L = spectrum.warp_log = ln(k/Lambda_IR)`.

Existing tree formulas to preserve:

```text
F_IR(c)^2 = (0.5 - c) / (1 - epsilon^(1 - 2c))
a_raw(c) = sum_n (M_KK^2 / m_n^2) * chi_n(1) * Omega_n(c)
a_sub(c) = a_raw(c) - a_ref
delta g_A^tree = s_Z * g_A^SM * (m_Z^2 / M_KK^2)
                 * U_A^\dagger diag(a_sub(c_i)) U_A

Delta S_tree = c_S * v_246^2 / M_KK^2
Delta T_tree(minimal) = [pi L / (2 c_W^2)] * v_246^2 / M_KK^2
Delta T_tree(custodial) = [-pi / (4 c_W^2 L)] * v_246^2 / M_KK^2
Delta U_tree = 0
```

Explicit Carena-to-repo `Zbb` normalization/sign map:

- The repo stores dimensionless additive shifts in the same normalized coupling basis consumed by `zpole.shifted_couplings`: `g_L^b = T3_b - Q_b s_eff^2 + z_delta_g_L_d[2,2]`, with no extra factor of `g_Z`. T010/T011 read this final additive entry directly.
- The tree gauge-profile term already includes the repo sign parameter `s_z=-1` and SM chiral factor through `_z_delta`; this tree map is not reused for fermion-loop `Zbb`.
- Carena's leading formulas for `delta g_{b_L}^s` and `delta g_{b_L}^q + delta g_{b_L}^\chi` are dimensionless vertex shifts in the normalized `Z b_L b_L` coupling. Therefore the PR2 loop insertion is exactly:

```text
z_delta_g_L_d[2,2] += delta_g_bL_Carena
```

- Do not multiply the Carena loop by repo `s_z`, by `g_Z`, or by `g_L^SM`; do not divide by `g_Z`. If a future source provides only a plotted relative value `delta g_{b_L}/g_{b_L}`, convert it first as `delta_g_additive = relative_value * g_L^SM_repo`; the analytic equations used in PR2 are already additive.
- With the synthetic singlet oracle below, adding `+0.00041018530641991833` changes T010's shifted bottom coupling to `g_left=-0.42242314802691344`, `g_right=0.07716666666666666`, decreases `R_b^0` to `0.21530246427000885`, and selects the `R_b^0` pull. This pseudo-observable oracle is required to guard the convention.

Loop formulas to implement as a tagged Carena leading proxy, not as a full custodian-spectrum calculation:

```text
s_W^2 = inputs.sin2_theta_w
c_W^2 = 1 - s_W^2
m_W = m_Z * c_W
alpha = inputs.alpha_em_mz
m_t = quark_fit_result.masses_up[2]
F_Q3 = bulk_state.F_Q[2] if present, else f_IR(c_Q[2], spectrum.epsilon)
F_u3 = bulk_state.F_u[2] if present, else f_IR(c_u[2], spectrum.epsilon)
Y_t_eff = m_t / (2 * V_EWSB * F_Q3 * F_u3)
```

The one-mode proxy exposes explicit mass/mixing knobs:

```text
M_t   = rho_t   * M_KK        (default rho_t = 1.0)
M_q   = rho_q   * M_KK        (default rho_q = 1.0)
M_chi = rho_chi * M_KK        (default rho_chi = 1.0)

m_q0t_t   = xi_s   * 2 * V_EWSB * Y_t_eff * F_Q3 = xi_s   * m_t / F_u3
m_qt_t    = xi_q   * 2 * V_EWSB * Y_t_eff * F_u3 = xi_q   * m_t / F_Q3
m_chid_t  = xi_chi * 2 * V_EWSB * Y_t_eff * F_u3 = xi_chi * m_t / F_Q3
```

Use `rho_*` relative to physical `M_KK` by default. If a custodian mass tied to the geometric IR scale is intended, users must pass `rho_* = Lambda_IR/M_KK ~= 1/2.450509663813736`; record this in metadata instead of changing the convention silently.

Carena leading singlet contribution:

```text
T_top = 3 m_t^2 / (16 pi s_W^2 c_W^2 m_Z^2)

Delta T_s =
  T_top * (2 m_q0t_t^2 / M_t^2)
  * [ ln(M_t^2 / m_t^2) - 1 + m_q0t_t^2 / (2 m_t^2) ]

delta g_L_b^s =
  alpha / (16 pi s_W^2 m_W^2)
  * (m_q0t_t^4 / M_t^2)
  * [ 1 + 2 m_t^2 / m_q0t_t^2 * (ln(M_t^2 / m_t^2) - 1) ]
```

Carena leading bidoublet vertex contribution:

```text
delta g_L_b^(q+chi) =
  alpha / (32 pi s_W^2 m_W^2) * m_t^2
  * [ (m_qt_t^2 / M_q^2) * ln(M_q^2 / m_t^2)
      - (m_chid_t^2 / M_chi^2) * ln(M_chi^2 / m_t^2) ]
```

Full one-loop `T` sign and magnitude for bidoublet custodians require the mass-basis `V^L,V^R,U^L,U^R,M` matrices in Carena Appendix `T:eq`; PR2 must not invent those matrices.

## Loop component contract

Separate "computed", "Zbb applied", and "EW001 T applied" semantics:

- `top_partner_loop_magnitudes_computed`: pure helper computed finite Carena leading numbers.
- `top_partner_zbb_loop_numerics_included`: a loop term was added to `z_delta_g_L_d[2,2]`.
- `top_partner_t_loop_numerics_included`: a loop term was added to EW001 `Delta T`.
- `top_partner_loop_numerics_included`: backward-compatible EW001-facing alias for `top_partner_t_loop_numerics_included`, not a generic "some loop happened" flag.

Component rules:

- `top_partner_loop_components="defer"` is the default and preserves PR1 tree behavior byte-for-byte.
- `singlet` applies both `Delta T_s` and `delta g_L_b^s` only when the `T` application is valid. Valid PR2 `T` application is either `top_partner_loop_t_sign=+1` for the positive singlet formula or a finite numeric `top_partner_loop_delta_t_override`.
- `bidoublet_vertex` may apply only `delta g_L_b^(q+chi)`. It must set `top_partner_t_loop_numerics_included=False` and `top_partner_delta_t_loop_applied=0.0` unless a finite numeric `top_partner_loop_delta_t_override` is supplied. This is the supported PR2 route for the model-dependent negative vertex case.
- `singlet_plus_bidoublet_vertex` applies the sum of the two vertex terms and a valid `T` term only if the singlet `T` rule above is satisfied. If the singlet `T` input is missing, compute magnitudes but do not apply the combined loop update.

Deterministic missing-sign behavior:

- Choose "compute magnitudes but do not apply" for missing `T` sign/override. Do not raise solely because `top_partner_loop_t_sign` is missing.
- For a sign-missing singlet request, metadata must say `top_partner_loop_mode="computed_not_applied_missing_t_input"`, `top_partner_loop_magnitudes_computed=True`, `top_partner_zbb_loop_numerics_included=False`, `top_partner_loop_numerics_included=False`, `top_partner_delta_t_loop_applied=0.0`, and the final matrices/EW001 result must equal PR1 tree-only.
- Still raise `ValueError` for nonfinite inputs, nonpositive masses/ratios, unknown component names, `top_partner_loop_t_sign=-1` without an override, or attempting loops in `minimal_rs`.

No fabricated negative `Delta T`:

- `top_partner_loop_t_sign=-1` alone is invalid and must raise. It must never apply `-Delta T_singlet_magnitude`.
- Negative `Delta T` in PR2 requires `top_partner_loop_delta_t_override=<finite negative number>` with metadata `top_partner_loop_t_source="explicit_numeric_override"`. A future separately validated bidoublet/custodian-spectrum `T` calculation may add another source, but PR2 must not label the leading vertex proxy as that calculation.
- Sign metadata remains useful: `+1` means `explicit_singlet_positive`; `override` means `explicit_numeric_override`; `-1` is accepted only as metadata consistent with a numeric override, not as an instruction to flip the singlet magnitude.

## Implementation steps

1. Extend APIs without changing defaults.
   - In `flavor_catalog_constraints/point_builder.py:105`, add and forward optional PR2 kwargs: `top_partner_loop_t_sign`, `top_partner_loop_delta_t_override`, `top_partner_loop_components`, `top_partner_loop_mass_ratios`, `top_partner_loop_mixing_scales`, `custodial_fcnc_mode`, and `kappa_fcnc`.
   - In `flavor_catalog_constraints/rs_ew_builder.py:39`, add the same kwargs. Replace the unconditional PR1 rejection at `flavor_catalog_constraints/rs_ew_builder.py:78` with: allow `include_top_partner_loops=True` only when `ew_model=="custodial_rs_plr"`; keep raising for `minimal_rs`; forward every PR2 kwarg, including `top_partner_loop_delta_t_override`, to `build_rs_ew_couplings`.
   - In `quarkConstraints/rs_ew_couplings.py:467`, add the same kwargs to the core `build_rs_ew_couplings` signature with defaults preserving PR1 behavior: `include_top_partner_loops=False`, `top_partner_loop_t_sign=None`, `top_partner_loop_delta_t_override=None`, `top_partner_loop_components="defer"`, `top_partner_loop_mass_ratios=None`, `top_partner_loop_mixing_scales=None`, `custodial_fcnc_mode="pr1_minimal_offdiag"`, and `kappa_fcnc=0.0`.
   - At `quarkConstraints/rs_ew_couplings.py:491`, replace the hard `include_top_partner_loops=True is deferred to PR2` raise with validated dispatch: raise only if `ew_model=="minimal_rs"`; otherwise continue into the custodial loop path.

2. Add a pure helper for loop numerics.
   - New dataclass in `quarkConstraints/rs_ew_couplings.py`, e.g. `RSCustodialTopPartnerLoopProxy`, holding inputs, unsigned magnitudes, requested/applied components, applied `Delta T`, applied `delta g_L_b`, booleans listed in the loop contract, and metadata.
   - New helper `_build_custodial_top_partner_loop_proxy(quark_fit_result, spectrum, inputs, ...)`.
   - Validate finite positive masses/ratios/mixing scales; require `quark_fit_result.masses_up[2]` for computed loops. Missing `top_partner_loop_t_sign` is not an input error; it produces the deterministic computed-not-applied mode above.
   - Use `bulk_state.F_Q/F_u` when present; otherwise compute from `c_Q/c_u` and `spectrum.epsilon`.
   - Do not import or touch lepton-sector code.

3. Layer loop terms on the PR1 tree result.
   - Build the loop proxy before `_apply_custodial_rs_plr_proxy`, after minimal left/right down snapshots exist.
   - Change `_apply_custodial_rs_plr_proxy` at `quarkConstraints/rs_ew_couplings.py:926` to accept `top_partner_loop_proxy` and FCNC options.
   - Preserve PR1 order: zero protected diagonal, apply optional `kappa_b/L` residual, zero elementary `b_R`.
   - If `top_partner_zbb_loop_numerics_included=True`, add `top_partner_delta_g_L_b_loop_applied` to `left[2,2]` after PR1 tree/residual handling and before `_hermitian`.
   - Store `z_delta_g_L_d_tree_pr1_b_before_top_partner` and `top_partner_delta_g_L_b_loop_applied` so T010/T011 can show what changed.
   - Do not add any top-partner term to `z_delta_g_R_d` in PR2; store `top_partner_delta_g_R_b_loop=0.0` and `top_partner_delta_g_R_b_loop_reason="not present in leading Carena Zb_L proxy"`.

4. Implement representation-aware custodial-FCNC modeling as an explicit mode.
   - Keep `custodial_fcnc_mode="pr1_minimal_offdiag"` as the default so PR1 tree behavior and T014 are byte-identical unless the new mode is requested.
   - Add `custodial_fcnc_mode="all_gen_bidoublet_mass_basis_proxy"`:
     - Interpret all generations as `(2,2)_{2/3}` in the mass basis.
     - The leading PLR-protected down-left gauge-profile FCNC matrix is zero.
     - Residual proxy is `left[i,j] = (kappa_fcnc / L) * minimal_z_delta_l_d_full[i,j]` for `i != j`; ideal PLR default is `kappa_fcnc=0.0`.
     - Right-handed down FCNC stays minimal. Do not zero `z_delta_g_R_d` off-diagonals in PR2; metadata must say RH FCNC is not protected by this proxy and remains minimal/deferred.
   - Add T014 diagnostics passthrough for `custodial_fcnc_modeling`, `custodial_fcnc_mode`, `custodial_fcnc_basis`, `custodial_fcnc_residual_source`, `custodial_fcnc_residual_applied`, `custodial_fcnc_leading_PLR_zeroed`, `custodial_fcnc_rh_status`, and `kappa_fcnc`, because current notes still describe minimal-RS off-diagonal entries at `flavor_catalog_constraints/secondary/top_higgs_ew/T014.py:539`.

5. Thread loop `Delta T` into EW001 without changing minimal defaults.
   - In `quarkConstraints/oblique_stu.py`, extend `rs_minimal_oblique_proxy` and `evaluate_rs_oblique_proxy` with optional `delta_t_loop=0.0` and `loop_metadata=None`.
   - Validate `delta_t_loop` is finite. If loop metadata says `top_partner_loop_numerics_included=True` but `delta_t_loop` is missing/nonfinite, raise `ValueError`; EW001 already catches this at `flavor_catalog_constraints/primary/top_higgs_ew/EW001.py:419`.
   - Keep `s` and `u` unchanged. Compute `t_tree = coeff_t * (v_246/M_KK)^2`, then `t = t_tree + delta_t_loop`.
   - Only emit new diagnostics keys when `delta_t_loop != 0.0` or loop metadata is provided, to avoid default diagnostic churn.
   - Export any new helper/dataclass through `flavor_catalog_constraints/physics_adapters/oblique_stu.py:5`.
   - In `flavor_catalog_constraints/primary/top_higgs_ew/EW001.py:384`, read `rs_ew_couplings.metadata`. Pass `delta_t_loop=top_partner_delta_t_loop_applied` only when `ew_model=="custodial_rs_plr"` and `top_partner_loop_numerics_included=True`; otherwise pass `0.0` and no loop metadata. A loop-deferred custodial point must yield exactly the PR1 tree-only `T`.
   - Diagnostics must include `t_tree_prediction`, `t_loop_prediction`, `top_partner_loop_numerics_included`, `top_partner_t_loop_numerics_included`, and `top_partner_loop_sign_convention` when loop metadata is present.

6. Update T010/T011 diagnostics and notes only.
   - Numeric predictions already use the final matrices, so no extra observable code is needed beyond the builder.
   - Extend metadata passthrough key lists at `flavor_catalog_constraints/primary/top_higgs_ew/T010.py:494` and `flavor_catalog_constraints/primary/top_higgs_ew/T011.py:674` with all new top-partner and custodial-FCNC keys listed below.
   - Change the active note text so it no longer says numerical top-partner loops remain deferred when `top_partner_zbb_loop_numerics_included=True`.
   - Keep `custodial_toppartner_zbL_deferred=True` only when no `Zbb` loop was applied. If a bidoublet vertex proxy was applied but full custodian spectrum is still absent, expose that through explicit proxy/omission metadata rather than the old blanket deferred flag.

## New metadata keys

Top-partner loop keys:

- `top_partner_loop_mode`: `deferred`, `computed_not_applied_missing_t_input`, or `carena_leading_proxy_applied`.
- `include_top_partner_loops`: existing bool, preserved.
- `top_partner_loop_magnitudes_computed`: bool.
- `top_partner_zbb_loop_numerics_included`: bool.
- `top_partner_t_loop_numerics_included`: bool.
- `top_partner_loop_numerics_included`: EW001-facing alias for `top_partner_t_loop_numerics_included`.
- `top_partner_loop_source`: `Carena-Ponton-Santiago-Wagner PhysRevD76 035006, arXiv:hep-ph/0701055`.
- `top_partner_loop_formula_set`: `leading_singlet_main_text_plus_bidoublet_vertex; full_Teq_Zbbeq_not_reconstructed`.
- `top_partner_loop_components_requested`: `defer`, `singlet`, `bidoublet_vertex`, or `singlet_plus_bidoublet_vertex`.
- `top_partner_loop_components_applied`: list.
- `top_partner_loop_components_deferred`: list.
- `top_partner_loop_t_sign`: `+1`, `-1`, `None`, or `override`.
- `top_partner_loop_t_source`: `explicit_singlet_positive`, `explicit_numeric_override`, `missing_t_input_not_applied`, or `not_applicable_vertex_only`.
- `top_partner_loop_sign_convention`: human-readable string; required when any loop is applied.
- `top_partner_loop_exact_Teq_Zbbeq_included`: bool, default `False`.
- `top_partner_loop_proxy_status`: `proxy_not_full_custodian_spectrum`.
- `top_partner_delta_t_singlet_magnitude`.
- `top_partner_delta_t_override`.
- `top_partner_delta_t_loop_applied`.
- `top_partner_delta_g_L_b_singlet`.
- `top_partner_delta_g_L_b_bidoublet_vertex`.
- `top_partner_delta_g_L_b_loop_applied`.
- `top_partner_delta_g_R_b_loop`.
- `top_partner_delta_g_R_b_loop_reason`.
- `top_partner_loop_applied_to_z_delta_g_L_d_22`.
- `z_delta_g_L_d_tree_pr1_b_before_top_partner`.
- `top_partner_loop_inputs`: dict with `m_top_gev`, `c_Q3`, `F_Q3`, `c_u3`, `F_u3`, `Y_t_eff`, `V_EWSB_gev`, `m_Z_gev`, `m_W_gev`, `alpha_em_mz`, `sin2_theta_w`, `cW2`, `M_KK_gev`, `lambda_ir_gev`, `warp_log`, and `M_KK_convention`.
- `top_partner_loop_mass_ratios`: dict with `rho_t`, `rho_q`, `rho_chi`, and `ratio_denominator="physical_M_KK"`.
- `top_partner_loop_mixing_scales`: dict with `xi_s`, `xi_q`, `xi_chi`.

Custodial-FCNC keys:

- `custodial_fcnc_modeling`: replace PR1 string only when the new FCNC mode is used; otherwise keep `deferred_PR2_off_diagonal_kept_minimal` for byte identity.
- `custodial_fcnc_mode`.
- `custodial_fcnc_basis`: `all_gen_bidoublet_mass_basis`.
- `custodial_fcnc_leading_PLR_zeroed`: bool.
- `custodial_fcnc_residual_source`: `kappa_fcnc*(1/L)*minimal_z_delta_l_d_full[i,j]`.
- `custodial_fcnc_residual_applied`: bool.
- `custodial_fcnc_rh_status`: `right_handed_down_offdiagonal_kept_minimal`.
- `kappa_fcnc`.
- `minimal_z_delta_l_d_full_offdiag`.
- `minimal_z_delta_r_d_full_offdiag`.
- `custodial_z_delta_l_d_offdiag_before_after`.

Honest omission keys:

- Keep existing `custodial_omissions`, `SU2_R_tower`, `custodian_spectrum`, `exact_NC_mixing`, and `BKT` fields for compatibility.
- Add `custodial_omission_details`: dict with explicit status strings, because the existing booleans are ambiguous.
- Add `custodian_spectrum_inferred=False`, `exact_top_partner_mass_matrix_included=False`, `brane_kinetic_terms_included=False`, `full_Teq_Zbbeq_loop_matching_included=False`.

## Numeric oracles for tests

Use a deterministic synthetic benchmark matching the existing test convention:

```text
M_KK = 3000.0 GeV, Lambda_IR = M_KK / 2.450509663813736
epsilon = 1e-15
c_Q3 = 0.43
c_u3 = 0.18
m_t = 173.0 GeV
alpha_em_mz = 1/127.952
sin2_theta_w = 0.23122
m_Z = 91.1876 GeV
rho_t = 1.0
xi_s = 1.0
top_partner_loop_components = singlet
top_partner_loop_t_sign = +1
```

Expected intermediate values:

```text
F_Q3 = 0.2656322304046161
F_u3 = 0.5656854250202850
Y_t_eff = 3.308347352802294
m_q0t_t = 305.8236828247721 GeV
T_top = 1.2084941744793902
```

Expected applied loop values for the singlet proxy with `M_t=M_KK`:

```text
Delta T_loop = +0.15745209098112428
delta g_L_b_loop = +0.00041018530641991833
delta g_R_b_loop = 0.0
```

T010 pseudo-observable oracle for this additive `delta g_L_b_loop` with otherwise zero `Zbb` shifts:

```text
shifted g_left = -0.42242314802691344
shifted g_right = +0.07716666666666666
R_b^0 predicted = 0.21530246427000885
A_b predicted = 0.9354140642148572
T010 selected_observable = R_b^0
T010 ratio = 1.4680590546247065
```

Bidoublet-vertex oracle using the same benchmark with `top_partner_loop_components="bidoublet_vertex"`, `rho_q=rho_chi=1.0`, `xi_q=0.5`, and `xi_chi=1.0`:

```text
m_qt_t = 325.63819483893786 GeV
m_chid_t = 651.2763896778757 GeV
delta g_L_b_bidoublet_vertex = -0.0003174965326198927
top_partner_delta_t_loop_applied = 0.0
top_partner_zbb_loop_numerics_included = True
top_partner_loop_numerics_included = False
```

This bidoublet test must assert that no negative `Delta T` is fabricated and metadata says the full bidoublet/custodian-spectrum `T` calculation is not included.

T014 RH-minimal oracle from the existing sample point:

```text
minimal sample bs offdiagonal:
  z_delta_g_L_d[1,2] = -0.00020387380225097264 - 2.2036569174674154e-08j
  z_delta_g_R_d[1,2] = -0.00019683644415246384 + 0j

With LH offdiagonals zeroed and RH offdiagonals kept minimal:
  T014 selected_channel = bs
  predicted = 6.369734635064622e-08
  ratio = 2.196460218987801e-05
```

Add a separate test with `rho_t=1/2.450509663813736` only to confirm the metadata convention, not as the default oracle, because using `Lambda_IR` as the partner mass gives much larger loop values.

## Test plan

- `tests/test_rs_ew_custodial_pr2.py::test_top_partner_loop_proxy_singlet_numeric_oracle_3tev`: build the synthetic benchmark above, assert intermediate and applied values to tight `pytest.approx` tolerances, assert `z_delta_g_L_d[2,2]` equals PR1 tree/residual plus `delta_g_L_b_loop`, and assert `top_partner_loop_numerics_included=True` means the EW001 `T` loop was applied.
- `test_t010_singlet_zbb_pseudo_observable_oracle`: inject only the singlet additive `delta_g_L_b_loop` into final `z_delta_g_L_d[2,2]`; assert shifted `g_left/g_right`, `R_b^0`, `A_b`, selected observable, and ratio exactly match the oracle above. This is the normalization/sign-map guard.
- `test_bidoublet_vertex_numeric_oracle_no_fake_negative_t`: request `bidoublet_vertex` with `xi_q=0.5`, `xi_chi=1.0`; assert `delta_g_L_b_bidoublet_vertex=-0.0003174965326198927`, matrix shift equals that value, `top_partner_delta_t_loop_applied=0.0`, `top_partner_loop_numerics_included=False`, and metadata does not claim a bidoublet `T` calculation.
- `test_negative_t_requires_numeric_override`: `top_partner_loop_t_sign=-1` without `top_partner_loop_delta_t_override` raises `ValueError`; with `top_partner_loop_delta_t_override=-0.05`, EW001 receives `-0.05`, diagnostics say `explicit_numeric_override`, and no source text claims a full bidoublet spectrum was computed.
- `test_top_partner_loop_missing_sign_computes_magnitudes_but_does_not_apply`: `include_top_partner_loops=True`, `top_partner_loop_components="singlet"`, and no sign/override computes magnitudes but leaves matrices and EW001 identical to PR1 tree-only; assert the computed-not-applied metadata keys.
- `test_minimal_rs_byte_identity_and_scan_hash`: keep existing PR1 test at `tests/test_rs_ew_custodial_pr1.py:118`; default minimal outputs and `scripts/run_full_catalog_scan.py` config hash `45e21a07585f7489` must not change.
- `test_custodial_pr1_tree_byte_identity_when_loops_and_fcnc_deferred`: with defaults (`include_top_partner_loops=False`, `custodial_fcnc_mode="pr1_minimal_offdiag"`), arrays, neutral contacts, and existing metadata keys match PR1 outputs.
- `test_15tev_point_survives_custodial_pr1_tree`: preserve `tests/test_rs_ew_custodial_pr1.py:273`; add a loop-enabled variant only if the selected benchmark remains below T010/T011 budgets.
- `test_ew001_adds_loop_t_as_separate_term`: assert `t_prediction = t_tree_prediction + t_loop_prediction`, `s_prediction` unchanged, `u_prediction=0`, and diagnostics include loop metadata only for custodial loop-enabled points.
- `test_ew001_loop_deferred_exactly_pr1_tree_only`: build a custodial point in computed-not-applied mode; assert EW001 diagnostics and numeric `t_prediction` exactly equal the PR1 custodial tree value.
- `test_t010_t011_consume_looped_zbb_matrix_and_report_components`: assert T010/T011 predictions shift exactly by reading final `z_delta_g_L_d[2,2]`, and diagnostics contain `z_delta_g_L_d_tree_pr1_b_before_top_partner`, `top_partner_delta_g_L_b_loop_applied`, `top_partner_zbb_loop_numerics_included`, and `top_partner_t_loop_numerics_included`.
- `test_custodial_fcnc_all_gen_bidoublet_zero_lh_mode`: with `custodial_fcnc_mode="all_gen_bidoublet_mass_basis_proxy", kappa_fcnc=0`, down-left off-diagonal entries are zero, right-handed down off-diagonals equal the minimal point, and T014 equals the RH-minimal oracle above. Do not assert total T014 rates are zero.
- `test_custodial_fcnc_residual_kappa_over_L`: with `kappa_fcnc=1`, off-diagonal left entries equal `minimal_z_delta_l_d_full[i,j] / L`; T014 changes according to the sum of the residual LH contribution plus unchanged RH contribution, not a pure `1/L^2` scaling unless RH is separately zeroed.
- `test_honest_omission_flags_present`: both loop-deferred and loop-computed custodial branches expose `custodial_omission_details`, `full_Teq_Zbbeq_loop_matching_included=False`, and `custodian_spectrum_inferred=False`.
- Run focused tests: `pytest tests/test_rs_ew_custodial_pr1.py tests/test_rs_ew_custodial_pr2.py tests/constraints/primary/top_higgs_ew/test_EW001.py tests/test_rs_ew_phase6a_zbb_fermion_mixing.py`.

## Isolation and overlap

- Files to touch in PR2:
  - `quarkConstraints/rs_ew_couplings.py`
  - `quarkConstraints/oblique_stu.py`
  - `flavor_catalog_constraints/physics_adapters/oblique_stu.py`
  - `flavor_catalog_constraints/primary/top_higgs_ew/EW001.py`
  - `flavor_catalog_constraints/primary/top_higgs_ew/T010.py`
  - `flavor_catalog_constraints/primary/top_higgs_ew/T011.py`
  - `flavor_catalog_constraints/secondary/top_higgs_ew/T014.py` diagnostics only
  - `flavor_catalog_constraints/rs_ew_builder.py`
  - `flavor_catalog_constraints/point_builder.py`
  - New/updated tests under `tests/`
- Do not touch W7 lepton-sector production files. The loop/FCNC work must not alter `z_delta_g_L/R_e`, neutrino metadata, LFV adapters, or charged-current lepton wiring.
- Avoid W9 scan-harness changes in PR2. Do not add default `ScanConfig` fields or CLI defaults that change `_config_payload`/`_config_hash`; W9 can opt into the new kwargs later. If a scan knob is unavoidable, omit false/default PR2 fields from the payload exactly like the existing `ew_model` omission asserted at `tests/test_rs_ew_custodial_pr1.py:137`.
- Determinism: all loop numerics are pure functions of point inputs and explicit kwargs; no random draws, no external data access, no mutable global state.
- Graceful degradation: default behavior remains PR1. Loop-enabled mode must either compute requested tagged values honestly under the component contract or fail loudly on invalid inputs; it must never silently label a zero as a computed loop.

PLAN-READY
