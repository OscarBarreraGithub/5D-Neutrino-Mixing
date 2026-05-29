1. BLOCKER: none found.
2. SHOULD-FIX: none found.
3. NIT: none found.

Evidence: C007 reaches physics via adapters in `flavor_catalog_constraints/primary/charm/C007.py:49` and `:52`; no tracked diffs in C004/C005, `rare_charm_dilepton`, or `deltaf2`.

Independent cross-check: manual quarkConstraints-core integration = `9.132644202844895e-12`; C007 smooth diagnostic = `9.132644202844895e-12`; diff = `0.000e+00`. Predicted total = `1.0132644202844895e-11`, ratio = `1.5123349556484916e-04`.

Safe/excluded: safe predicted `2.03409648170265e-08`, ratio `0.3035964898063657`, passes `True`; excluded predicted `5.085001204256627e-07`, ratio `7.589554036203921`, passes `False`.

Contract/probes: numeric result fields are real `float`; complex values stay in diagnostics; missing couplings returns non-crashing `passes=True`, `predicted=None`, `ratio=None`; deterministic equality and no input mutation confirmed. Missing anchor and mismatched-unit anchor both raise `AnchorError`.

Tests: `python -m pytest tests/constraints/primary/charm/test_C007.py -q` -> `11 passed`; `python -m pytest tests/constraints/ -q` -> `352 passed in 5.21s`. Registry smoke: `C007 charm HARD`, anchors `6.7e-08`, `1e-12`.

CODE-OK