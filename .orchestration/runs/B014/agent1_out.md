Plan completed: studied B012/B013 + B014 YAML, added B014, tested focused and full suite.

Physics: C7/C7' exclusive scaling with `BR_norm^(b->d)` from B014 YAML PDG rows; CKM diagnostic `|Vtd/Vts|^2 = 0.0456810801`.

Machinery: reused `quarkConstraints.bsgamma` through `physics_adapters/bsgamma.py`; appended b→d wrappers selecting down-sector `(0, 2)` and preserving LL C7/C8 running. No core bsgamma functions changed.

Budget: measurement-consistency bands, not theory room: `B0->rho0 gamma` hardest at `1.5e-7`; charged `2.5e-7`; omega `1.8e-7`, all from `B014.yaml`.

Validation: no-NP `B0->rho0 gamma = 8.6e-7`; tests independently recompute C7/RG numerics and verify b-s slot is ignored.

NEEDS-HUMAN-PHYSICS: RS C7/C7' proxy and exclusive B→rho/omega form factors, annihilation/spectator/helicity likelihood.

Files changed: [B014.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/secondary/beauty/B014.py), [test_B014.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/secondary/beauty/test_B014.py), [bsgamma.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/physics_adapters/bsgamma.py).

Tests: `13 passed` focused B014; `1054 passed` for `python -m pytest tests/constraints/ -q`.

Note: unrelated pre-existing CR014/collider worktree changes remain untouched.