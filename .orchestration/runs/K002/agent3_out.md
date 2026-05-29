1. NIT `tests/constraints/primary/kaon/test_K002.py:105`: redundant hardcoded budget literal; fix by relying on the YAML-derived assertion at lines 102-104, or label this as a regression pin.
2. No BLOCKER/SHOULD-FIX findings. K002 imports no other constraint; `git diff quarkConstraints/` is empty; adapter AST check shows no removed/changed existing function signatures, only added wrappers.
3. Contract OK: numeric result fields are real floats; complex couplings/Wilsons stay in diagnostics; severity is HARD; missing extra returns non-crashing `passes=True`, `predicted=None`, `ratio=None`.
4. Anchor OK: `load_anchor("K002", family="kaon", candidates=("pdg_fit_assuming_cpt",))`; missing-key probe raised `AnchorError`; unit-mismatch probe raised `AnchorError`.
5. Cross-check safe point: predicted/direct `1.627683060065448e-16`; ratio `9.343985976273966e-02` vs direct-core `9.343760390731619e-02`; budget `1.741957944070575e-15`; pass `True`.
6. Cross-check excluded point: predicted/direct `1.627683060065448e-14`; ratio `9.343985976273967e+00` vs direct-core `9.343760390731619e+00`; budget `1.741957944070575e-15`; pass `False`.
7. Determinism OK: repeated evaluate results equal; coupling arrays unchanged.
8. Registry smoke OK: `fcc.get("K002") -> K002 kaon HARD`, registry contains K002.
9. Pytest: `python -m pytest tests/constraints/ -q` -> `51 passed in 2.31s`.

CODE-OK