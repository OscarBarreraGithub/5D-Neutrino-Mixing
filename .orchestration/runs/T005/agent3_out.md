1. BLOCKER: none.
2. SHOULD-FIX: none.
3. NIT: none.

Isolation: OK. `T005.py` imports only scaffold/base/registry plus `physics_adapters.top_fcnc`; adapter diff is append-only at `flavor_catalog_constraints/physics_adapters/top_fcnc.py:373`, with no existing function/class body changes by AST check. `quarkConstraints/top_fcnc.py` and `deltaf2.py` unchanged.

Contract/anchor/tests: OK. Missing couplings returns non-crashing `passes=True` with `predicted=None`; numeric result fields are real floats. Anchors route through `load_anchor`; probes raised `AnchorError` for missing value_id and mismatched selected block.

Numerical cross-check: YAML limit `3.7e-4`, SM anchor `4.6e-12`. Probe `(left=1+0.25j,right=0.3j,M_KK=3000)` gives constraint BR `1.118111915877719e-4`, core BR `1.118111915877719e-4`, delta `0`, ratio `0.3021924096966808`. Safe `left=1`: BR `9.701621829741596e-5`, ratio `0.2622059953984215`, passes. Excluded `left=2`: BR `3.8806487318966383e-4`, ratio `1.048823981593686`, fails.

Determinism: same input produced equal `ConstraintResult`; input coupling arrays unchanged in tests.

Pytest: targeted T005 `9 passed`; full `python -m pytest tests/constraints/ -q` `463 passed in 13.23s`. Registry smoke: `44` constraints, T005 registered as `top_higgs_ew/HARD`.

CODE-OK