1. BLOCKER [T011.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/top_higgs_ew/T011.py:170): missing T011 anchors do not fail loudly; the code follows the historical `canonical_home` stub and loads T010 (`source_process_id T010`, `fallback True`). Fix: require T011 `pdg_or_equivalent` anchors for `A_FB^{0,b}` and `A_b`, and update tests to reject fallback.
2. BLOCKER [test_T011.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/top_higgs_ew/test_T011.py:31): tests explicitly accept the same T010 fallback, so they mask the catalog-source violation. Fix: assert `source_process_id == "T011"` and `canonical_home_fallback is False`.

Cross-check: SM manual/result `A_b=0.935535274431`, `A_FB=0.103278860422`; matches T011 validation target.
Cross-check: excluded point `delta_gR=0.0277172613125`, manual/result `A_b=0.884074865239`, `A_FB=0.0975978641376`, `ratio=1.62848107829`, `passes=False`.
Safe/missing/determinism: safe ratio `0`, `passes=True`; missing-extra ratio `0`, `passes=True`; repeated evaluate equal `True`.
Isolation/contract: T011 imports physics only through `physics_adapters.zpole`; adapter/core tracked unchanged; numeric result fields are real floats; registry smoke passed.
Pytest: T011 targeted `11 passed`; full `python -m pytest tests/constraints/ -q` => `862 passed in 18.01s`.

CODE-NEEDS-FIXES