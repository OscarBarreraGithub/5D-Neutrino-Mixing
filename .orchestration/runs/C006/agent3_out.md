1. BLOCKER: none.
2. SHOULD-FIX: none.
3. NIT: none.
4. Isolation OK: C006 adds new C006/test/adapter/core only; no tracked diff to C004/C005/C007 or existing rare_charm/deltaf2 code; AST primary-constraint imports empty.
5. Contract OK: severity HARD, numeric result fields are real floats, complex Wilsons stay in diagnostics, missing/invalid couplings return non-crashing unevaluated results.
6. Anchor OK: C006.py uses load_anchor; YAML value = 1.3e-08 from canonical_current_limit. Missing candidate/value and bad units probes raised AnchorError.
7. Cross-check: constraint BR = 3.752990413820220023e-13, ratio = 2.886915702938630754e-05.
8. Independent numbers: quarkConstraints core = 3.752990413820220023e-13, manual formula = 3.752990413820220023e-13, both diff = 0.
9. Safe point: BR = 5.428606395969134911e-13, ratio = 4.175851073822411545e-05, pass=True.
10. Excluded point: BR = 3.392878997480708271e-04, ratio = 2.609906921139006226e+04, pass=False.
11. Registry smoke: C006 charm HARD BR(D0 -> e+- mu-+) LFV.
12. Pytest: `python -m pytest tests/constraints/ -q` -> 407 passed in 7.12s.
CODE-OK