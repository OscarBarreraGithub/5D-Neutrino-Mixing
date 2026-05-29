1. BLOCKER: none.
2. SHOULD-FIX: none.
3. NIT: none.
4. Isolation: B018 imports no other constraint and reaches physics via `flavor_catalog_constraints.physics_adapters.rare_b_meson`; existing rare-b adapter/core files show no diff.
5. Anchor probes: missing observable raises `AnchorError`; mismatched `load_anchor` block raises `AnchorError`.
6. Numerical cross-check: anchor `R_K=0.947`, sigma `0.047`, budget `0.100`; q2 `(1.1, 6.0)`.
7. SM cross-check: independent manual BR `1.851950636963255649e-07`; core BR `1.851950637006692028e-07`; constraint denominator same; `R_K=1`, ratio `0.53`, pass.
8. NP cross-check: direct core BR `1.796818113915925792e-07`; expected/constraint `R_K=0.970230025580014`, ratio `0.23230025580014`, pass.
9. Safe/excluded: safe `R_K=0.977059271664665`, ratio `0.300592716646654`, pass; large NP `R_K=0.792391917396395`, ratio `1.54608082603605`, fail.
10. Determinism/missing-couplings: repeated evaluations equal; missing extra returns non-crashing pass with `predicted=None`, `ratio=None`.
11. Tests: `python -m pytest tests/constraints/primary/beauty/test_B018.py -q` -> `11 passed`; `python -m pytest tests/constraints/ -q` -> `407 passed`.
12. Registry smoke: `B018 beauty PRIMARY HARD R_K central-q2`.
CODE-OK