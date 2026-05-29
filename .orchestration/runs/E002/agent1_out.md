Plan executed: add isolated E002 constraint, load YAML anchor loudly, reuse EDM adapter with `lepton="mu"`, add focused tests, run full constraint suite.

Physics: pure-NP `|d_mu| <= limit`, with `d_mu[e cm] = c_CPodd[GeV^-1] * hbar*c`; SM reported as `0.0`.

Reused machinery: `quarkConstraints/edm.py` through `flavor_catalog_constraints.physics_adapters.edm`; no new physics module.

Budget: `1.8e-19 e cm`, from `E002.yaml` `canonical_direct_limit` / PDG Live 2026 S004EDM, Bennett 2009, 95% CL.

Validation: independent hbar-c recomputation for `c_CPodd=2.0e-6 GeV^-1` gives `|d_mu|=3.946539608e-20 e cm`, ratio `0.21925220044444443`.

NEEDS-HUMAN-PHYSICS: full RS one-loop lepton/KK/Higgs/EW matching is not on `ParameterPoint`; explicit proxy is flagged in diagnostics.

Files changed: [E002.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/edm_neutrino/E002.py), [test_E002.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/edm_neutrino/test_E002.py).

Pytest: `tests/constraints/primary/edm_neutrino/test_E002.py -q` -> 12 passed; `tests/constraints/ -q` -> 526 passed.