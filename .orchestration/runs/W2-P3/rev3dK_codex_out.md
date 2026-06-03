1. BLOCKER: none; mapping is physically consistent: `C10/C10p -> Y_NP` at `flavor_catalog_constraints/physics_adapters/rare_kaon_dilepton.py:206`, `y7V/y7A` at `:242`.
2. SHOULD-FIX: none; no second `1/M_KK^2` or old proxy prefactor in the RS path (`wilson_prefactor_reused=False`, `second_mkk_suppression_applied=False`).
3. NIT: `git status` also shows untracked orchestration notes plus `tests/rare_kaon_phase3d_helpers.py:1`; include/drop before final commit.
4. K006 recompute: sample `predicted=1.0056093658836135e-09`; independent core same.
5. SM-limit: universal-c/SM point gives K006 `8.209801430715005e-10 == sm_prediction`, matching committed SM.
6. K012-Im: preserved at `flavor_catalog_constraints/physics_adapters/rare_kaon_dilepton_kshort_mumu.py:154`; diag source `Im[-lambda_c Y_c + lambda_t C10]^2`, not real-rescaled.
7. K009/K010-LD: K009 keeps `interference_sign_fixed=False`/vector-only interference at `rare_kaon_dilepton_muon.py:320`; K010 keeps NEEDS-HUMAN/a_S sign envelope at `K010.py:610`.
8. Degradation: all K006/K008/K009/K010/K012 old-style points are non-vetoing: `passes=True`, `predicted=None`, `ratio=None`, `evaluated=False`, `missing_extra=rs_semileptonic_wilsons`.
9. Test-count: unchanged vs HEAD: K006 10->10, K008 10->10, K009 10->10, K010 9->9, K012 11->11.
10. Isolation/fields: tracked diff only scoped adapters/constraints/tests; evaluated `ConstraintResult` numeric fields are real floats, complex values stay in diagnostics.
11. Pytest: `pytest tests/ -q` => `1646 passed, 1 skipped in 788.01s`.
12. P3DK-OK