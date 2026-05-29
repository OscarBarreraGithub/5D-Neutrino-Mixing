Plan executed: studied scaffold/K001/K004/B005.yaml, built shared core + adapter, implemented B005, added tests.

Physics: Buras `b->s l l` Hamiltonian, C10-dominant `B_s->mumu`; anchors from `B005.yaml` PDG2026 exp and Czaja-Misiak 2024 SM.

Reused/built: reused `QuarkMassBasisCouplings`, CKM target, B-meson constants; built `quarkConstraints/rare_b_dilepton.py` and adapter `physics_adapters/rare_b_meson.py`.

Budget: HARD NP-shift budget `5.95465734053883e-10` from YAML exp-vs-SM loose-edge construction.

SM validation: independent formula gives `BR(B_s->mumu)=3.6483515801441554e-9`, consistent with `3.64(12)e-9`.

Gap: RS C9/C10 Z/KK-penguin proxy is flagged `NEEDS-HUMAN-PHYSICS` in docstrings and diagnostics.

Files added: `quarkConstraints/rare_b_dilepton.py`, `flavor_catalog_constraints/physics_adapters/rare_b_meson.py`, `flavor_catalog_constraints/primary/beauty/B005.py`, `tests/constraints/primary/beauty/test_B005.py`.

Tests: `python -m pytest tests/constraints/ -q` -> `157 passed`.