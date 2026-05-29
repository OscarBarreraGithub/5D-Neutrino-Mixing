Plan completed: studied scaffold/examples and `B033.yaml`, added the INFO stub + adapter, added focused tests, ran the suite.

Physics: `S_phiK_S = 0.74 +/- 0.12` vs local HFLAV `sin(2 beta) = 0.710 +/- 0.011`; reports `Delta S` only, non-vetoing.

Reused vs built: reused scaffold `load_anchor`; built new adapter `charmless_b_phiks.py`. Did not reuse B002 mixing as a full `S_phiK_S` prediction.

Budget: `0.12050311199301038`, from quadrature of the B033 YAML HFLAV uncertainties.

Validation: independent YAML recompute gives `Delta S = 0.030000000000000027`, ratio `0.24895622877980005`.

NEEDS-HUMAN-PHYSICS: SM QCDF/SCET hadronic-penguin treatment and RS `Delta B=1 b -> s sbar s` penguin matching are both flagged.

Files changed: `flavor_catalog_constraints/primary/beauty/B033.py`, `flavor_catalog_constraints/physics_adapters/charmless_b_phiks.py`, `tests/constraints/primary/beauty/test_B033.py`.

Tests: `test_B033.py` 9 passed; full `python -m pytest tests/constraints/ -q` 810 passed.

Note: unrelated untracked CR/E/orchestration files are present in the worktree; I left them untouched.