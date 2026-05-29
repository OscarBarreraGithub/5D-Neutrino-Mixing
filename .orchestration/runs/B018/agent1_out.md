Implemented B018.

- Plan: scoped B018 to the central-q2 `R_K` YAML row, reused B016/B017 `rare_b_dilepton` machinery, added focused tests, ran full constraint suite.
- Physics: `R_K = BR(B+ -> K+ mu mu) / BR(B+ -> K+ e e)` for `1.1 < q^2 < 6.0`, using the existing Buras C9/C10 + BCL-like `B -> K` q2 integral.
- Machinery: reused `flavor_catalog_constraints.physics_adapters.rare_b_meson`; no new adapter or direct physics import in the constraint.
- Budget: HARD, `0.10000000000000005`, from `B018.yaml` central row: `|0.947 - 1.0| + 0.047`.
- Validation: SM LFU proxy `R_K = 1.0`; same-bin SM branching proxy `1.851950637006692e-07`.
- NEEDS-HUMAN-PHYSICS: full electron-vs-muon RS matching, QED/radiative handling, C7/nonlocal charm, scalar/tensor, and covariance are missing; diagnostics document the muon-only-over-SM-electron proxy and lepton-universal cancellation caveat.
- Files changed: [B018.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/beauty/B018.py:1), [test_B018.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/beauty/test_B018.py:1).
- Tests: `test_B018.py` 11 passed; `python -m pytest tests/constraints/ -q` 385 passed.