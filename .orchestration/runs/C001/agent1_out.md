Plan executed: mirror `K001`, add an append-only D0 adapter wrapper, implement `C001`, then test direct running-evaluator agreement.

Files changed: `flavor_catalog_constraints/primary/charm/C001.py`, `tests/constraints/primary/charm/test_C001.py`, and additive wrappers in `flavor_catalog_constraints/physics_adapters/deltaf2.py`.

Evaluator used: `quarkConstraints.deltaf2.evaluate_d0_mixing_with_running`, via adapter, with `mu_had=2.0`.

Budget: `Delta_m_D.value_GeV / 2 = 6.562e-15 / 2 = 3.281e-15 GeV`, from `flavor_catalog/processes/charm/C001.yaml:145-157`; convention cited from `quarkConstraints/deltaf2.py:959-961` and `docs/audits/bag_param_inventory.md:53`.

Numerical cross-check: safe point `|M12|=1.1447248065797324e-16`, ratio `3.4889509496486816e-02`, passes. Excluded point `|M12|=4.5788992263189290e-14`, ratio `1.3955803798594724e+01`, fails. Direct running evaluator matched `|M12|`.

Tests: `python -m pytest tests/constraints/primary/charm/test_C001.py -q` -> 7 passed. `python -m pytest tests/constraints/ -q` -> 51 passed.

Diff/stat for C001 files: `C001.py` 226 insertions, `test_C001.py` 223 insertions; adapter is modified additively. Workspace also contains unrelated untracked B001/B003/K002 files from other agents, so global status is not C001-only.

Reviewer scrutiny: charm mixing policy is intentionally conservative and theory-dirty; plugin uses the catalog 2025 budget, while the core evaluator’s embedded default budget is older (`3.125e-15 GeV`).