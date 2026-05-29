Plan: read CR004/scaffold patterns, append adapter helper, implement CR004/test, run suite.

Physics: HARD mass proxy `m_B = M_KK`, ratio `m_limit / m_B`, sourced from `CR004.yaml` PDG2025/CMS bottom-partner pair limits.

Reused: `quarkConstraints.collider_resonance` via adapter. Built: append-only CR004 B-pair adapter helpers; no CR001/CR002/CR003 changes.

Budget: active `B(B -> H b)=1` limit `m_B > 1570 GeV = 1.57 TeV`; `bZ=1.54 TeV` and `tW=1.56 TeV` kept in diagnostics.

Numerical validation: independent core recompute at `M_KK=1800 GeV` gives ratio `1.57/1.8 = 0.872222`, pass. No SM prediction for this direct collider bound.

NEEDS-HUMAN-PHYSICS: full `sigma*BR`, width, branching-simplex, acceptance, and mass-curve recast is still a documented proxy gap.

Files: [collider_resonance.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/physics_adapters/collider_resonance.py), [CR004.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/collider_rs/CR004.py), [test_CR004.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/collider_rs/test_CR004.py).

Tests: `test_CR004.py` 8 passed; `python -m pytest tests/constraints/ -q` 678 passed; `git diff --check` clean.

Workspace note: unrelated untracked B026/L004/T006/T020/orchestration files were already present and left untouched.