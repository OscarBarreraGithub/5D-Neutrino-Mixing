1. NIT `tests/constraints/primary/collider_rs/test_CR001.py:111`: this test uses `build_from_quark_couplings`, which auto-adds `kk_gluon_mass_gev`; add a direct `make_point(quark_mass_basis_couplings=SimpleNamespace(M_KK=...))` case to pin the fallback path.

Isolation: OK. CR001 imports no other constraint, physics access is via `physics_adapters.collider_resonance`; relevant implementation files are new/untracked, with no pre-existing function body modified.

Contract/anchor: OK. Numeric fields are real floats; missing input returns non-crashing result. Anchor value loaded from `CR001.yaml` through `load_anchor`; missing value-id and wrong-unit probes both raise `AnchorError`.

Cross-check: `M_KK=6000 GeV` gives result/core mass `6.0 TeV`, ratio `0.916666666667`, pass `True`. `6500 GeV`: ratio `0.846153846154`, pass `True`. `3000 GeV`: ratio `1.83333333333`, pass `False`.

Determinism: OK. Same input gives identical `ConstraintResult`; point extras were not mutated. Registry smoke: CR001 present, `import_failures=0`.

Pytest: `python -m pytest tests/constraints/ -q` -> `509 passed in 15.33s`.

CODE-OK