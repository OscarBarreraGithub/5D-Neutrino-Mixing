Plan completed: studied scaffold/examples/`B025.yaml`, implemented, tested.

Physics: `R_D = R_D^SM |1 + C_tau^proxy|^2`, with `C_tau^proxy ∝ xi_cb m_b m_tau / M_KK^2`.

Sources: `B025.yaml` HFLAV `R_D = 0.358 ± 0.024`, SM `R_D = 0.296 ± 0.004`.

Budget: HARD envelope `0.08633105012119288` from `|exp - SM| + sqrt(sigma_exp^2 + sigma_SM^2)`.

Built: `quarkConstraints/semileptonic_lfu.py` plus adapter `flavor_catalog_constraints/physics_adapters/semileptonic_lfu.py`.

Constraint/test: `flavor_catalog_constraints/primary/beauty/B025.py`, `tests/constraints/primary/beauty/test_B025.py`.

Validation: SM prediction loads/evaluates to `0.296`; SM point ratio to budget is `0.7181657110965687`.

Gap: `NEEDS-HUMAN-PHYSICS` flagged for full charged-current RS/W'/charged-Higgs/leptoquark/lepton-neutrino/form-factor matching.

Tests: `python -m pytest tests/constraints/ -q` -> `617 passed in 14.30s`.

Note: unrelated pre-existing dirty/untracked worktree files were left untouched.