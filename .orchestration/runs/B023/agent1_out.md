Plan implemented: reused B022 rare-B machinery, added only a K* adapter wrapper, then added B023 + tests.

Physics: `BR(B -> K* nu nubar) = BR_SM * r_kstar`, with `r_kstar=(1+1.31 eta) epsilon^2` from shared `b -> s nu nubar`.

Budget: HARD total-BR limit `2.7e-5` from `B023.yaml`; SM anchor `9.2e-6`, diagnostic limit-minus-SM room `1.78e-5`.

Validation: SM point gives `9.2e-6` (`r_kstar=1`, ratio `0.3407`), matching the expected ~`9.0e-6`.

Gap: flagged `NEEDS-HUMAN-PHYSICS`; full RS EW KK/Z/Z' and neutrino matching are absent, so this uses the documented Z-like proxy.

Files changed: `physics_adapters/rare_b_nunu.py`, `primary/beauty/B023.py`, `test_B023.py`. B023 YAML source is Belle 2017, not Belle II; I followed the sidecar.

Tests: `python -m pytest tests/constraints/ -q` -> `193 passed`.