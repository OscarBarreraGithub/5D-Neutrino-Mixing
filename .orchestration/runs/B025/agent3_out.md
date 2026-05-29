1. BLOCKER: none.
2. SHOULD-FIX: none.
3. NIT: none.
4. Isolation: B025 reaches physics only via adapter at `flavor_catalog_constraints/primary/beauty/B025.py:52`; adapter wraps `quarkConstraints.semileptonic_lfu` at `flavor_catalog_constraints/physics_adapters/semileptonic_lfu.py:12`. B025-added files are new only; no pre-existing function-body edits found.
5. Contract: real float fields returned at `flavor_catalog_constraints/primary/beauty/B025.py:396`; complex amplitudes remain diagnostics at `quarkConstraints/semileptonic_lfu.py:221`. Missing extra returns `passes=True,predicted=None,ratio=None` at `B025.py:347`.
6. Anchor probe: `load_anchor` calls were `[("canonical_average",), ("sm_reference_used_by_hflav",)]`; forced missing and mismatched candidate both raised `AnchorError`.
7. Anchor numbers: exp `0.358`, exp sigma `0.024`, SM `0.296`, SM sigma `0.004`, budget `0.08633105012119288`.
8. Independent cross-check: SM pred `0.296`, ratio `0.7181657110965687`, pass `True`; constraint/direct core/manual formula agree.
9. Independent cross-check: complex NP pred `0.35539830853917825`, ratio `0.030136219322821116`, scalar `(0.09562630065761932+0.01652182157958868j)`; all agree.
10. Safe/excluded: safe `bottom_right=3.8` pred `0.3549928491555573`, ratio `0.034832784267319936`, pass `True`; excluded `bottom_right=20` pred `0.666587328299304`, ratio `3.574465130055806`, pass `False`.
11. Determinism: repeated `evaluate(point)` equality is `True`; B025 test also checks coupling matrices unchanged at `tests/constraints/primary/beauty/test_B025.py:290`.
12. Tests: `python -m pytest tests/constraints/ -q` -> `626 passed in 14.30s`.
13. Registry smoke: `fcc.get("B025")` -> `flavor_catalog_constraints.primary.beauty.B025`, family `beauty`, severity `HARD`; `import_failures_count=0`.

CODE-OK