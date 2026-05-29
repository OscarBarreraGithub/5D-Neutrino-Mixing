Plan: implemented the C003 INFO stub, adapter, and tests.

Physics: no SM/RS penguin calculation; LHCb `Delta A_CP` anchor is loaded from `C003.yaml`.

Built: `physics_adapters/charm_direct_cp.py` stub adapter only. Reused: scaffold `load_anchor`.

Budget: non-vetoing NP room = full observed `|Delta A_CP| = 1.54e-3`, source `lhcb2019_discovery_average`.

Validation: YAML recompute gives experimental `-1.54e-3`, ratio `1.0`, naive LHCb significance `5.31 sigma`.

Gaps flagged: `NEEDS-HUMAN-PHYSICS` for SM long-distance penguins and RS Delta C=1 penguin matching.

Files changed: `flavor_catalog_constraints/physics_adapters/charm_direct_cp.py`, `flavor_catalog_constraints/primary/charm/C003.py`, `tests/constraints/primary/charm/test_C003.py`.

Tests: `python -m pytest tests/constraints/primary/charm/test_C003.py -q` -> 8 passed; `python -m pytest tests/constraints/ -q` -> 705 passed.