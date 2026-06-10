# W8 plan REVISION — address dual-review fixes (codex + Opus, both NEEDS-FIXES)

Revise your plan `.orchestration/runs/W8-CUSTODIAL-PR2/plan_codex.md` IN PLACE to resolve ALL
items below. Keep the parts both reviewers validated (the singlet Carena oracle is correct; sign-
as-metadata is correct; PR1 byte-identity approach is correct; M_KK convention explicit is correct).
Re-verify each fix against the REAL repo files. End the revised plan with `PLAN-READY`.

## From codex review (NEEDS-FIXES)
1. **No fabricated negative ΔT.** Do NOT allow `top_partner_loop_t_sign=-1` to simply apply
   `-ΔT_singlet_magnitude`. A negative ΔT must require EITHER an explicit numeric
   `top_partner_loop_delta_t_override` OR a separately-validated bidoublet/custodian-spectrum T
   calculation. Sign-flipping the singlet magnitude and calling it a bidoublet proxy = fabrication.
2. **Coherent loop-component semantics.** `singlet` applies BOTH ΔT_s and δg_bL^s; `bidoublet_vertex`
   may apply ONLY the vertex proxy unless a real/override ΔT is supplied. Don't let one boolean
   (`top_partner_loop_numerics_included`) mean both "some Zbb loop applied" AND "EW001 T loop applied".
3. **Fix the T014 FCNC test.** T014 reads BOTH `z_delta_g_L_d` and `z_delta_g_R_d` off-diagonals
   (`flavor_catalog_constraints/secondary/top_higgs_ew/T014.py:450`); RH down off-diagonals are
   nonzero in the sample point. With `custodial_fcnc_mode="all_gen_bidoublet_mass_basis_proxy",
   kappa_fcnc=0`: assert LH off-diagonals are zero and T014 equals the RH-minimal contribution —
   do NOT claim total T014 rates are zero while RH FCNC stays minimal. (Or introduce + justify an
   explicit RH-FCNC-zeroing mode.)
4. **Thread `top_partner_loop_delta_t_override`** through EVERY API layer including the
   `build_rs_ew_couplings` signature (currently omitted there).
5. **Deterministic missing-sign behavior.** Pick ONE: "compute magnitudes but do not apply" OR
   "raise ValueError". Make all tests match the single chosen behavior.

## From Opus review (NEEDS-FIXES)
6. **Normalization/sign map (CORRECTNESS-BLOCKING).** State the EXPLICIT map from Carena's δg_L^b
   (plan lines ~83-95) to the repo's additive `(t3 − Q·s²)` Z-pole convention that T010/T011
   consume (`zpole.py:202,220`, `T011.py:761`; tree entries carry `s_z·g_L^SM·(m_Z/M_KK)²`,
   `rs_ew_couplings.py:1137`), including any `s_z` sign / `g_Z` factor. Add a T010 PSEUDO-OBSERVABLE
   oracle (hand-checked Zbb observable), not just a raw-number assertion.
7. **Second rejection site (CORRECTNESS-BLOCKING).** There is a second hard `raise` at
   `quarkConstraints/rs_ew_couplings.py:491-492` (`include_top_partner_loops=True is deferred to PR2`).
   Explicitly replace it with the validated loop path gated on `ew_model=="custodial_rs_plr"`, and
   keep it raising for `minimal_rs`. (Step 1 only mentioned `rs_ew_builder.py:78`.)
8. **EW001 loop contract.** EW001.evaluate wraps the oblique helper in try/except ValueError
   (`EW001.py:417`). Specify: EW001 passes `delta_t_loop` ONLY when
   `top_partner_loop_numerics_included=True` (else 0.0); helper raises ValueError on inconsistent
   loop inputs; add a test that a loop-DEFERRED custodial point yields EXACTLY the PR1 tree-only T.
9. **Bidoublet-vertex oracle.** Add a deterministic numeric oracle for the `bidoublet_vertex`
   δg_L^b component (the negative-T model-dependent case), not only the singlet.

End the revised plan with `PLAN-READY`.
