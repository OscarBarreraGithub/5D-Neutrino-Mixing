Plan executed: studied B005/scaffold/B006 yaml, implemented B006, added focused tests, ran full suite.

Physics: reused Buras `b -> q ll` C10 formula with `lambda_t^d = V_tb V_td*`; `B_d` has `y_d=0`, so SM time factor is `1`.

Machinery: reused `quarkConstraints/rare_b_dilepton.py` and `physics_adapters/rare_b_meson.py`; no shared core changes needed.

Budget: HARD upward-NP room is `1.5e-10 - 1.06e-10 = 4.4e-11` from B006.yaml CMS/PDG vs Bobeth SM; LHCb `2.3e-10` kept as comparison.

SM validation: independent recomputation gives `BR(B0 -> mu mu) = 1.0466846288738303e-10`.

NEEDS-HUMAN-PHYSICS: RS piece remains the documented C9/C10 proxy; full EW/lepton/scalar matching is absent from `ParameterPoint`.

Files changed by me: [B006.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/beauty/B006.py), [test_B006.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/beauty/test_B006.py).

Tests: `test_B006.py` 12 passed; `python -m pytest tests/constraints/ -q` 193 passed.