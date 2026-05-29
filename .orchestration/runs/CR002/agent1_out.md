Plan: studied CR001, K001/K004, scaffold contracts, CR002 YAML; implemented adapter proxy, CR002 constraint, and tests; ran targeted and full suites.

Physics: CR002 uses `m_VLQ = M_KK` as the documented mass proxy and compares it to the CR002 YAML pair-production mass lower limit.

Reused/built: reused `quarkConstraints.collider_resonance` through the adapter; appended reusable VLQ pair helpers only.

Budget: active HARD limit is `1.46 TeV` (`1460 GeV`) from `PDG2026:CR002:ATLAS2023:X53_pair_Wt`; the `1.60 TeV` pair+single coupling-dependent entry is diagnostic only.

Numerical validation: for `M_KK = 1.6 TeV`, ratio is `1.46 / 1.6 = 0.9125`, checked against the core collider-resonance comparison.

NEEDS-HUMAN-PHYSICS: full `sigma*BR`, branching fractions, VLQ spectrum, widths, acceptance, and mass-dependent limit curves are not available on `ParameterPoint`.

Files changed: [collider_resonance.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/physics_adapters/collider_resonance.py), [CR002.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/collider_rs/CR002.py), [test_CR002.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/collider_rs/test_CR002.py).

Tests: collider targeted `17 passed`; full `python -m pytest tests/constraints/ -q` -> `534 passed`. Unrelated untracked worktree files were left untouched.