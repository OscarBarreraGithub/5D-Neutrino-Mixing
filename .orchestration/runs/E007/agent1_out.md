Plan executed: inspected E006/K001/K004 plus `E007.yaml`, implemented the stub, added tests, ran the suite.

Physics: INFO-only direct-limit bookkeeping for `d_Ra, d_Xe`; no SM or RS atomic EDM prediction. Sources: Bishof 2016, Sachdeva 2019, Allmendinger 2019 from `E007.yaml`.

Built/reused: reused scaffold anchors/registry and the E006 stub pattern; added `atomic_edm.py` as a small adapter.

Budget: strongest current direct Ra/Xe limit is `|d_A(129Xe)| < 1.4e-27 e cm` at 95% CL from Sachdeva et al. E007 also records `|d(225Ra)| < 1.4e-23 e cm`.

Numerical validation: independent YAML recompute gives Xe ratio `1.4e-28 / 1.4e-27 = 0.1`; Xe combined uncertainty `6.90e-28 e cm`.

NEEDS-HUMAN-PHYSICS: flagged both nuclear/atomic Schiff + atomic-structure inputs and RS CP-odd matching.

Files changed: `flavor_catalog_constraints/physics_adapters/atomic_edm.py`, `flavor_catalog_constraints/primary/edm_neutrino/E007.py`, `tests/constraints/primary/edm_neutrino/test_E007.py`.

Tests: `test_E007.py` 7 passed; full `python -m pytest tests/constraints/ -q` 898 passed.