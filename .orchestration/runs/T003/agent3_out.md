1. BLOCKER `tests/constraints/primary/top_higgs_ew/test_T003.py:73`: numerical cross-check is a hand-coded formula, not a direct recomputation through `quarkConstraints.top_fcnc`. Fix by importing the core evaluator, e.g. `quarkConstraints.top_fcnc.photon_dipole_branching_fraction`, computing the dipoles independently from the couplings, and comparing that core result to `Constraint.evaluate()`.

2. Evidence: T003 isolation passes. `T003.py` imports physics only via `flavor_catalog_constraints.physics_adapters.top_fcnc`; adapter AST check showed no changed/removed pre-existing symbols, only added `TopPhotonDipoleProxyCouplings`, helper validators, `top_photon_dipole_proxy_from_couplings`, and `t_to_q_gamma_from_couplings`.

3. Contract/anchor/determinism checks passed: numeric result fields are floats; complex amplitudes stay in diagnostics; missing extra returns non-crashing HARD/pass with `predicted=None`; missing/wrong anchors raise `AnchorError`; repeated evaluation is identical and does not mutate coupling matrices.

4. Independent core numbers I ran: mixed point `BR=6.0684210348288108e-06`, core `BR=6.0684210348288108e-06`, diff `0`, ratio `0.40188218773700735`; safe left=1 gives `5.2654412449707681e-06`, ratio `0.34870471820998467`, pass `True`; excluded left=2 gives `2.1061764979883072e-05`, ratio `1.3948188728399387`, pass `False`.

5. Tests: `python -m pytest tests/constraints/primary/top_higgs_ew/test_T003.py -q` -> `9 passed`; `python -m pytest tests/constraints/ -q` -> `352 passed`; registry smoke -> `T003` registered, `34` total constraints, `0` import failures.

CODE-NEEDS-FIXES