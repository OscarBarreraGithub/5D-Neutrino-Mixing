1. BLOCKER tests/constraints/primary/charm/test_C005.py:237: stale golden value is masked by `pytest.approx` default `abs=1e-12`. Actual for the tested point is `1.28111287014103257e-17`, but the test expects `1.333296784847654e-17` and still passes. Fix the golden and use tight relative tolerance with `abs=0.0`; also tighten small-value checks at lines 236 and 328.

Cross-check numbers: probe `left=2.0e-2+0.3e-2j,right=0.4e-2j`: constraint/core/manual all `1.28111287014103257e-17`, ratio `1.62166186093801586e-10`.

Safe/excluded: `left=1` gives `3.19479518738412005e-14`, ratio `4.04404454099255696e-07`, pass; `left=5000` gives `7.98698796846029900e-07`, ratio `1.01101113524813897e+01`, fail.

Isolation/contract/anchors/determinism: OK. C005 imports no other constraint, reaches physics via adapter, pre-existing rare-charm/C004 function AST hashes unchanged, missing input returns clean result, mismatch/missing anchor probes raise `AnchorError`, repeated evaluation is equal.

Pytest: `tests/constraints/primary/charm/test_C005.py -q` = 12 passed; `python -m pytest tests/constraints/ -q` = 296 passed. Registry smoke: C005 registered as charm PRIMARY HARD.

CODE-NEEDS-FIXES