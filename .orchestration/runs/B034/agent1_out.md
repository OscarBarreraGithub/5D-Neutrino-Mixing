Implemented B034 as an INFO-only non-vetoing penguin-CP stub.

- Physics: records measured `phi_s^{s sbar s}` for `B_s -> phi phi`; no QCDF/SCET SM amplitude and no RS `Delta B=1` penguin matching.
- NEEDS-HUMAN-PHYSICS: flagged separately for SM hadronic-penguin physics and RS penguin matching.
- Machinery: added new adapter [charmless_b_phiphi.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/physics_adapters/charmless_b_phiphi.py); constraint imports physics only through the adapter.
- Anchor/budget: `phi_s = -0.074 ± 0.069 rad`, `|lambda| = 1.009 ± 0.030`; non-vetoing budget is full observed `|phi_s| = 0.074 rad` from `B034.yaml`.
- Validation: independent convention check `-2 * 3.7e-2 = -0.074 rad`; result has `predicted=None`, `sm_prediction=None`, `ratio=1.0`.
- Files changed: [B034.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/beauty/B034.py), adapter above, [test_B034.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/beauty/test_B034.py).
- Tests: `test_B034.py` 9 passed; full `python -m pytest tests/constraints/ -q` 898 passed.
- Workspace also has unrelated dirty/untracked files from other tasks; I left them untouched.