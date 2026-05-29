Plan executed: studied K001/K004 plus C003/K013 stub patterns, added K003, then ran focused and full constraint tests.

Physics: K003 is an INFO-only non-vetoing stub for `Re(epsilon'/epsilon)`. It loads the PDG/KTeV/NA48 anchor from `K003.yaml`, records RBC/UKQCD 2020 as SM context, and does no penguin calculation.

Reused: scaffold `load_anchor` via a flat-anchor virtual adapter pattern. Built: new append-only `kaon_direct_cp` room-comparison adapter.

Budget: `0.00166`, sourced from `flavor_catalog/processes/kaon/K003.yaml pdg_or_equivalent`, as full observed non-vetoing NP room.

Validation: experimental `0.00166 +/- 0.00023`; RBC/UKQCD SM context `0.00217 +/- 0.0008378544`; KTeV `0.00192`, NA48 `0.00147`.

NEEDS-HUMAN-PHYSICS flagged on both SM penguin/lattice cancellation and RS `Delta S=1` penguin/chromomagnetic matching.

Files changed: `flavor_catalog_constraints/primary/kaon/K003.py`, `flavor_catalog_constraints/physics_adapters/kaon_direct_cp.py`, `tests/constraints/primary/kaon/test_K003.py`.

Tests: K003 focused `9 passed`; full `python -m pytest tests/constraints/ -q` `830 passed`; `git diff --check` clean.