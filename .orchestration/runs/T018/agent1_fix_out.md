Change: `test_T018.py` now cross-checks `Constraint.evaluate()` against `quarkConstraints.higgs_lfv` core BR/width evaluator.
Cross-check: eval/core BR `6.4283981736403212e-05`, width `2.6163580566716109e-07`, norm `5.2500000000000007e-08`, ratio `0.042855987824268804`.
Pytest: `test_T018.py -q` -> 12 passed; `tests/constraints/ -q` -> 568 passed.