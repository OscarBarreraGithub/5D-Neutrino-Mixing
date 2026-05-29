Implemented E001.

- Plan: load E001 YAML anchor, build EDM core + adapter, add E001 constraint, add tests, verify.
- Physics: pure-NP `|d_e| = |c_CPodd| * hbar*c`; SM reported as `0.0` at catalog precision.
- Built: [quarkConstraints/edm.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/edm.py), [physics_adapters/edm.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/physics_adapters/edm.py), [E001.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/edm_neutrino/E001.py), [test_E001.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/edm_neutrino/test_E001.py).
- Budget: `|d_e| < 4.1e-30 e cm` at 90% CL from `E001.yaml` `canonical_limit` (PDG Live 2026 / Roussy 2023 HfF+; ACME 2018 retained as superseded).
- Validation: `hbar*c = 1.973269804e-14 GeV cm`; coefficient budget `2.077769594249e-16 GeV^-1`; `1e-16 -> 1.973269804e-30 e cm`, ratio `0.481285`.
- NEEDS-HUMAN-PHYSICS: full RS loop matching needs complex lepton/KK/Higgs/EW couplings not on `ParameterPoint`; every E001 result flags this in diagnostics.
- Tests: `test_E001.py + test_contract.py`: `24 passed`; full suite excluding unrelated untracked B009: `500 passed`.
- Raw `python -m pytest tests/constraints/ -q` is currently blocked by unrelated untracked `B009`: `490 passed, 1 failed` on a B009 budget constant mismatch.