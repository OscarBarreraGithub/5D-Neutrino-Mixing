1. NIT tests/constraints/primary/collider_rs/test_CR007.py:53 - missing-anchor is tested, but the CR007.py:277 mismatched-`load_anchor` guard is only covered by my manual probe. Fix: add a monkeypatch test returning a wrong `block_key`.
2. NIT flavor_catalog_constraints/physics_adapters/kk_graviton_resonance.py:76 - uses the KK-gluon-named converter only as a generic GeV->TeV validator. Fix: optional neutral alias in `collider_resonance`.

Isolation: CR007 imports no other constraint; reaches physics via `physics_adapters.kk_graviton_resonance` -> `collider_resonance`; existing tracked `collider_resonance.py`, `deltaf2.py`, and `quarkConstraints/collider_resonance.py` have no diff.
Contract: numeric result fields are real floats; missing inputs return non-crashing `passes=True predicted=None ratio=None`; invalid negative mass returns non-crashing failure.
Anchor: active value loaded from YAML via `load_anchor`; missing ID raises `AnchorError`; mismatch probe raised `AnchorError` for wrong `block_key`.

Cross-check: YAML active limit = 1.4 TeV.
Cross-check: m=2000 GeV -> CR007 ratio 0.7, core ratio 0.7, pass=True.
Cross-check: safe m=2400 GeV -> CR007/core ratio 0.5833333333333334, pass=True.
Cross-check: excluded m=900 GeV -> CR007/core ratio 1.5555555555555554, pass=False.
Determinism: repeated evaluate equal; point extras not mutated.
Registry smoke: `fcc.get("CR007")` OK; `import_failures=0`.
Pytest: `tests/constraints/primary/collider_rs/test_CR007.py -q` = 8 passed.
Pytest: `tests/constraints/ -q` = 810 passed in 17.45s.
CODE-OK