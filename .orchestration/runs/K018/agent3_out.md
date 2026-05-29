1. SHOULD-FIX `tests/constraints/primary/kaon/test_K018.py:163`: the test recomputes YAML arithmetic by hand but does not call the designated `quarkConstraints.ckm_extraction` core; add a direct core recomputation via `extract_vus_from_kl3` + `vus_consistency_pull` and compare to `Constraint.evaluate()`.
2. Evidence: isolation passes otherwise. K018 adds only its new constraint/test/adapter/core paths; AST imports show no other constraint imports, and physics is reached via `flavor_catalog_constraints.physics_adapters.ckm_extraction` at `K018.py:47`.
3. Contract: `K018.py:336-344` returned real floats; empty/missing-coupling point non-crashed with `passes=True`, `predicted=0.22330377397401527`, `ratio=0.0071207056891883`, severity SOFT.
4. Anchor probes: `K018.py:202-208` routes through `load_anchor`; missing value ID and forced mismatched block both raised `AnchorError`.
5. Numerical verifier: core `Vus=0.22330377397401527`, `sigma=0.00053839741927841009`, `pull=0.0071207056891882999`; constraint matched exactly, deltas 0.
6. Safe/excluded probes: safe `predicted=0.2`, ref `0.2`, ratio `0.0`, pass True; excluded ref `0.16`, ratio `4.000000000000001`, pass False.
7. Determinism: repeated `evaluate()` equal; point extras preserved as `{'kk_ew_mass_gev': 6000.0}`; empty-point result equal.
8. Tests: `python -m pytest tests/constraints/ -q` => `509 passed in 13.98s`; registry smoke `K018_registered=True`, `import_failures_count=0`.

CODE-NEEDS-FIXES