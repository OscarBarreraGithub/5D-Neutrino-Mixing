Plan: inspected scaffold/examples, added adapter + K013 INFO stub, added focused tests, ran full suite.

Physics: ChPT O(p^4/p^6) chiral-loop/vector-exchange dominated; sources are K013.yaml PDG 2025 average with KTeV 2008 + NA48 2002 inputs.

Machinery: reused registry/load_anchor pattern and C003-style non-vetoing stub semantics; built only `physics_adapters/radiative_kaon.py`.

Budget: non-vetoing full observed BR = `1.273e-6` from `flavor_catalog/processes/kaon/K013.yaml:pdg_or_equivalent`.

Validation: independent YAML recompute checks BR `1.273e-6`, uncertainty `3.3e-8`, ratio `1.0`; no SM prediction is returned by design.

NEEDS-HUMAN-PHYSICS: flagged both SM ChPT calculation and RS Delta S=1 radiative/ChPT matching.

Files changed: [radiative_kaon.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/physics_adapters/radiative_kaon.py), [K013.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/kaon/K013.py), [test_K013.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/kaon/test_K013.py).

Pytest: `test_K013.py` 9 passed; latest `python -m pytest tests/constraints/ -q` 753 passed. Note: unrelated concurrent dirty files remain in the worktree.