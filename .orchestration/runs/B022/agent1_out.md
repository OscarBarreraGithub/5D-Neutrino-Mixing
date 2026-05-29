Plan executed: studied K001/K004 + scaffold, built shared b→sνν core/adapter, added B022 + tests.

Physics: `X_t` Inami-Lim short-distance response, `C_L^SM=-X_t/sin^2θ_W`; `B+→K+νν` uses `R_K=|C_L+C_R|^2/|C_L^SM|^2`.

Built: `quarkConstraints/rare_b_nunu.py` and adapter `flavor_catalog_constraints/physics_adapters/rare_b_nunu.py`; B023 can reuse `C_L,C_R,epsilon,eta`.

Budget: Belle II vs HPQCD from `B022.yaml`; central residual `1.742e-5`, direction sigmas `+7.080741486595878e-6/-6.413805422680049e-6`, HARD.

SM validation: `X_t=1.4812287076698927`, `BR_SM=5.58e-6`.

RS NP: documented Z-like proxy, flagged `NEEDS-HUMAN-PHYSICS` in docstrings and diagnostics.

Files added: `B022.py`, `test_B022.py`, `rare_b_nunu.py`, adapter `rare_b_nunu.py`.

Tests: `test_B022.py` 10 passed; full `python -m pytest tests/constraints/ -q` 157 passed.