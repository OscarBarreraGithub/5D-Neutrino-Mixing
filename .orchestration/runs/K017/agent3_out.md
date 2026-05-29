1. BLOCKER tests/constraints/primary/kaon/test_K017.py:179: the “independent” numerical cross-check does not call the underlying `quarkConstraints.leptonic_tree` core; it manually recomputes from `constraint.ratio_inputs`, which are produced through the constraint path. Fix: build core inputs from the YAML SM anchor and compare to `quarkConstraints.leptonic_tree.evaluate_leptonic_lfu_ratio(...)`.

Isolation OK: K017 imports only `physics_adapters.leptonic_tree`; adapter imports `quarkConstraints.leptonic_tree`; AST check found no modified/removed existing top-level funcs/classes, only additions.

Contract OK: numeric result fields are floats, complex amplitudes stay in diagnostics, severity is HARD, missing `kk_ew_mass_gev` returns non-crashing SM result.

Anchor OK: exp=2.4880000000000002e-05, SM=2.4770000000000002e-05, budget=2.005538513813748e-07; missing/mismatched anchor probes raise `AnchorError`.

Actual core cross-check: SM pred=2.4770000000000002e-05, core=2.4770000000000002e-05, diff=0, passes=True.

Actual core cross-check: MKK=20 GeV pred=2.480019354356144e-05, core=2.480019354356144e-05, ratio=0.15055080395350695, passes=True.

Actual core cross-check: MKK=3000 GeV pred=2.4770001341526599e-05, core=2.4770001341526599e-05, ratio=6.689109125945149e-06, passes=True.

Actual core cross-check: MKK=3 GeV pred=2.6129690623206301e-05, core=2.6129690623206301e-05, ratio=6.7796784446722036, passes=False.

Determinism OK: repeated evaluate on same point equal; point extra unchanged.

Pytest: `test_K017.py` 10 passed; `python -m pytest tests/constraints/ -q` 724 passed in 15.66s. Registry smoke OK: `K017 kaon HARD`.

CODE-NEEDS-FIXES