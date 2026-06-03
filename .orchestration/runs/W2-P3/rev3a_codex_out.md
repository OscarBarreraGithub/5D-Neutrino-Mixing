1. BLOCKER: none. `quarkConstraints/rs_ew_couplings.py:224`/`:259`/`:430` use `a(c)-a_ref`, correct rotations, and coupling-shift scaling; Hermitian gate at `:130`.
2. Recomputed `z_delta_g_R^b`: 3TeV `-1.660422757e-3`, 6TeV `-4.151056893e-4`, ratio `0.25`, sign/magnitude OK.
3. SHOULD-FIX: none. Contacts at `quarkConstraints/rs_ew_couplings.py:461` carry `g_Z^2/m_Z^2`; Wilsons at `quarkConstraints/rs_semileptonic_wilsons.py:272` do not use old `_wilson_prefactor`/second `1/M_KK^2`.
4. Recomputed `b->s mu mu`: `C9=0.3046272275+0.0058031605j`, `C10=-4.055208034-0.0772518704j`, matching stored extras.
5. `c->u mu mu C10=-20.44077450-31.38594017j`: convention-physical from tiny `lambda_b` (`|lambda_b|=1.5267e-4`), not a contact/mass-suppression bug.
6. SHOULD-FIX: none. SM-limit recompute: max `z=0`, contacts `0`, Wilsons `0`; family-universal max offdiag `3.83e-20`, diagonal `z_Rd33=-1.660422757e-3`.
7. SHOULD-FIX: none. Schema/fail-loud/immutability OK: `flavor_catalog_constraints/point_builder.py:61`, `:74`; `ParameterPoint.extras` mapping proxy at `flavor_catalog_constraints/base.py:110`.
8. NIT: none. Isolation OK: primary/secondary tree has zero refs to new extras/builder; registered catalog count remains `103`; tests manually recompute at `tests/test_rs_ew_phase3a.py:342`.
9. `pytest tests/ -q`: `1652 passed, 1 skipped in 811.46s`.

P3A-OK