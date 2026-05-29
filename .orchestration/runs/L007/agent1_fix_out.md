Change: updated `test_L007.py` so the τ→μγ numeric cross-check uses `check_mu_to_e_gamma_raw(..., pmns[[1,2,0], :])` and compares BR/ratio/core diagnostics to `Constraint.evaluate()`.

Pytest: `test_L007.py -q` = 11 passed; `tests/constraints/ -q` = 862 passed.