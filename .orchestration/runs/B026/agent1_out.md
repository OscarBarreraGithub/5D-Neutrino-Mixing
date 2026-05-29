Plan executed: added isolated B026 constraint plus tests; B025 left untouched.

Physics: `R_D* = R_D*^SM |1 + C_tau^proxy|^2`, with HFLAV CKM 2025 anchors from `B026.yaml`.

Reused machinery: `quarkConstraints/semileptonic_lfu.py` via `physics_adapters.semileptonic_lfu`; no new physics core added.

Budget: `0.0390830459735946`, from `B026.yaml` HFLAV `0.281±0.011` vs SM `0.254±0.005`.

SM validation: SM-limit prediction is `0.254`, matching the YAML SM anchor.

Gap: RS charged-current NP remains `NEEDS-HUMAN-PHYSICS` and is flagged in diagnostics.

Files changed: [B026.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/beauty/B026.py), [test_B026.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/beauty/test_B026.py).

Tests: `python -m pytest tests/constraints/primary/beauty/test_B026.py -q` → 11 passed; `python -m pytest tests/constraints/ -q` → 664 passed.