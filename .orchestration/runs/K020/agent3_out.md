1. BLOCKER: none.
2. SHOULD-FIX: none.
3. NIT: none.
4. Isolation OK: K020-owned files are new additions; no other constraint imports; targeted diff on reused rare-kaon/K009/K010/K019 files was empty.
5. Contract OK: missing-input smoke returned `passes=True`, `predicted=None`, `ratio=None`, `experimental=1.3e-11`, `budget=1.3e-11`; evaluated numeric fields are floats, complex Wilsons stay in diagnostics.
6. Anchor OK: value comes through YAML/load_anchor; missing value_id probe raised `AnchorError`; focused test also probes mismatched CL.
7. Numerical cross-check: adapter `1.32388300111444405e-12`; independent `quarkConstraints` recomputation `1.32388300111444405e-12`; relative diff `0.0`.
8. Pass/fail points: safe `predicted=1.00875927346318871e-12`, `ratio=0.07759686718947606`, passes; excluded `predicted=2.52189818365797150e-11`, `ratio=1.939921679736901`, fails.
9. Determinism OK: repeated `evaluate()` equality and input non-mutation covered by test.
10. Tests: `python -m pytest tests/constraints/secondary/kaon/test_K020.py -q` -> `11 passed`; `python -m pytest tests/constraints/ -q` -> `1002 passed in 26.69s`.
11. Registry smoke OK: `fcc.get("K020")` resolves to `SECONDARY kaon HARD`.

CODE-OK