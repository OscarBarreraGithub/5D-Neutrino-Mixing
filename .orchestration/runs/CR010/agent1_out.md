Plan completed: read scaffold/gold standards, reused collider mass-limit machinery, verified CR010 sidecar, ran focused and full tests.

Physics: CR010 compares `m_T = m_B = M_KK` to the simultaneous ATLAS/PDG 2025 `(T,B)` doublet pair-production lower bound from `CR010.yaml`.

Reused/built: reused `collider_resonance` mass-vs-limit recast; added CR010-specific append-only adapter helpers for the `(T,B)` mass proxy. No CR001-CR009 files changed.

Budget: active HARD bound is `1.37 TeV`, from `PDG2025:CR010:TB_doublet_T_mass_lower_limit` and `...B_mass_lower_limit`.

Validation: independent core recomputation at `M_KK = 1600 GeV` gives ratio `1.37 / 1.6 = 0.85625`, passing.

Gap: `sigma x BR`, branching-surface, T/B spectrum, widths, splittings, acceptance, and limit-surface recast are flagged `NEEDS-HUMAN-PHYSICS`.

Files: `flavor_catalog_constraints/primary/collider_rs/CR010.py`, `tests/constraints/primary/collider_rs/test_CR010.py`, `flavor_catalog_constraints/physics_adapters/collider_resonance.py`.

Tests: `test_CR010.py` `8 passed`; full `tests/constraints/` `937 passed`.