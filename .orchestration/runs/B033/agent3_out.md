1. NIT `tests/constraints/primary/beauty/test_B033.py:171` - hard-coded ratio literal is off in the last digits: direct YAML recompute gives `0.2489562261407833`, literal is `0.24895622877980005`. Fix: update/remove the duplicate literal; current `pytest.approx` still passes.
2. Cross-check: YAML direct `S_phiKs=0.74`, `sigma=0.12`; `sin2beta=0.710`, `sigma=0.011`; `DeltaS=0.030000000000000027`, `sigmaDelta=0.12050311199301036`, `ratio=0.2489562261407833`. `evaluate()` matches these, `passes=True`, `severity=INFO`, `predicted=None`.
3. Adapter safe/excluded: safe `ratio=0.2489562261407833 pass=True`; excluded advisory `Delta=1`, `sigma=0.014142135623730951`, `ratio=70.710678118654755 pass=False`.
4. Isolation/contract/anchor OK: no other constraint imports; no `deltaf2.py` diff; physics path only via `charmless_b_phiks`; empty point non-crashes; numeric fields are float/None; missing/mismatched anchors raise `AnchorError`.
5. Registry smoke OK: `B033 beauty INFO S_phiK_S`; `import_failures=0`.
6. Pytest: `python -m pytest tests/constraints/ -q` -> `810 passed in 16.25s`.

CODE-OK