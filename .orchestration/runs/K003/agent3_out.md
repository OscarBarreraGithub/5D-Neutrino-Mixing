1. No BLOCKER / SHOULD-FIX / NIT findings for K003 code/tests.

Isolation: `K003.py` imports no sibling constraints and no `quarkConstraints`; physics access is only through `flavor_catalog_constraints.physics_adapters.kaon_direct_cp`. `deltaf2.py` is tracked and unmodified; K003 constraint/test/adapter are new files.

Numerical cross-check, recomputed from `K003.yaml` without calling the constraint adapter: experimental `0.00166`, budget `0.00166`, ratio `1.0`, SM context `0.00217`, SM context uncertainty `0.0008378544026261365`; all matched result fields exactly.

Adapter pass/fail probe: safe room ratio `1.0` passes; overfilled room ratio `2.0` fails. Constraint result remains INFO/non-vetoing and missing-couplings empty point returns `K003` without crashing.

Anchor probes: mismatched `load_anchor` raises `AnchorError`; missing candidate raises `AnchorError`.

Determinism: repeated `evaluate(empty_point)` results equal; point extras unchanged.

Registry smoke: `K003 kaon INFO`, registry total `84`, import failures `0`.

Pytest: `python -m pytest tests/constraints/primary/kaon/test_K003.py -q` -> `9 passed`; `python -m pytest tests/constraints/ -q` -> `862 passed in 16.64s`.

CODE-OK