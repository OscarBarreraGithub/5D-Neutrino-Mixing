Option (a): sign-envelope implemented; K010 now evaluates `s_a=±1` and uses the least-excluded branch for the HARD verdict.

Both-sign diagnostic point `left_sd=1e-5`: `+` branch `BR=5.792892672e-9`, ratio `0.002961`; `-` branch `BR=5.807111685e-9`, ratio `0.002452`, selected `-`.

Large-NP check `left_sd=1e-1`: `+` ratio `50.61`; `-` ratio `99.64`, selected `+`, still fails.

Pytest: `tests/constraints/primary/kaon/test_K010.py -q` = `10 passed`; `tests/constraints/ -q` = `407 passed`.