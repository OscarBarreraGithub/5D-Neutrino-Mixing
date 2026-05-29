1. SHOULD-FIX tests/constraints/primary/top_higgs_ew/test_T018.py:183 - numeric test recomputes by local manual formula, not the required `quarkConstraints.higgs_lfv` core evaluator. Fix: add a direct core comparison.

2. NIT none - isolation/contract/anchor checks pass: T018 reaches physics only via adapter, result scalars are real floats, complex Yukawas stay in diagnostics, missing input is non-crashing, missing/mismatched anchors raise `AnchorError`.

3. NUMBERS: anchor=1.5e-3; cross point BR=6.42839817364032e-05, core BR=6.42839817364032e-05, ratio=0.0428559878242688, width=2.61635805667161e-07.

4. PASS/FAIL probes: safe `Y_mu_tau=2e-4` BR=4.89782717991643e-05 ratio=0.0326521811994429 passes=True; excluded `Y_mu_tau=1e-2` BR=0.122445679497911 ratio=81.6304529986073 passes=False.

5. Determinism: repeated `evaluate()` results equal; input extra unchanged. Registry smoke: T018 registered, import_failures=0.

6. Pytest: `python -m pytest tests/constraints/ -q` -> 568 passed in 14.42s.

CODE-NEEDS-FIXES