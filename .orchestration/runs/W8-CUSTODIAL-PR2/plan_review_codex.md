**Findings**

The leading Carena formulas in the plan are transcribed correctly for the positive singlet `Delta T`, singlet `delta g_bL`, and bidoublet vertex `delta g_bL` proxy. The numeric oracle also checks out. Source: [arXiv:hep-ph/0701055](https://arxiv.org/abs/hep-ph/0701055).

The plan is not approvable as written because it allows a negative `Delta T_loop` by sign-flipping the singlet magnitude and labeling it a bidoublet proxy. Since the plan explicitly says the full bidoublet `T` magnitude needs the unreconstructed Carena mass-basis matrices, this would fabricate a numerical `Delta T` model.

There is also a concrete FCNC test contradiction: T014 reads both `z_delta_g_L_d` and `z_delta_g_R_d` off-diagonals in [T014.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/secondary/top_higgs_ew/T014.py:450). The plan says right-handed down FCNC stays minimal, but then expects zero T014 branching fractions when left off-diagonals are zero. In the real sample point, right-handed off-diagonals are nonzero, so that test cannot pass without also zeroing/changing right-handed FCNC.

**Numbered Fixes**

1. Do not allow `top_partner_loop_t_sign=-1` to apply `-Delta T_singlet_magnitude`. Negative `Delta T` must require either an explicit numeric `top_partner_loop_delta_t_override` or an implemented, separately validated bidoublet/custodian-spectrum `T` calculation.

2. Make loop component semantics coherent: `singlet` should apply both `Delta T_s` and `delta g_bL^s`; `bidoublet_vertex` may apply only the vertex proxy unless a real/override `Delta T` is supplied. Avoid one boolean like `top_partner_loop_numerics_included` meaning both “some Zbb loop applied” and “EW001 T loop applied”.

3. Fix the FCNC test plan: with `custodial_fcnc_mode="all_gen_bidoublet_mass_basis_proxy", kappa_fcnc=0`, assert left off-diagonals are zero and T014 equals the right-handed-minimal contribution, or explicitly introduce and justify a right-handed FCNC-zeroing mode. Do not claim total T014 rates are zero while right-handed FCNC remains minimal.

4. Thread `top_partner_loop_delta_t_override` through every planned API layer. It is listed for `point_builder`, but omitted from the proposed `build_rs_ew_couplings` signature despite being required by the physics section.

5. Resolve the inconsistent missing-sign behavior: the plan alternates between “compute magnitudes but do not apply” and “raise `ValueError`”. Pick one deterministic behavior and make the tests match it.

VERDICT: NEEDS-FIXES