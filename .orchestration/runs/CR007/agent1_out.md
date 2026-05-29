Plan: add CR007 only, reuse collider resonance recast, test registration/anchors/math/pass-fail/purity.

Physics: HARD mass lower-bound proxy, `ratio = m_limit / m_GKK`; source values come from `CR007.yaml`.

Reused/built: reused `quarkConstraints.collider_resonance`; added a small spin-2 adapter wrapper plus CR007 constraint.

Budget: `1.4 TeV` from `PDG2025:CR007:bulk_graviton_diboson_mass_kMpl_0p5` (PDG 2025 / CMS 2023 WW/ZZ all-jets bulk-RS benchmark).

Validation: independent core check at `m_GKK = 2.0 TeV` gives `ratio = 0.7`, pass.

NEEDS-HUMAN-PHYSICS: full `sigma*BR`, width, branching fractions, `k/Mbar_Pl`, acceptance, and limit-curve recast remain proxy-only and flagged in diagnostics.

Files changed: `physics_adapters/kk_graviton_resonance.py`, `primary/collider_rs/CR007.py`, `tests/constraints/primary/collider_rs/test_CR007.py`.

Tests: `test_CR007.py` 8 passed; full `python -m pytest tests/constraints/ -q` 803 passed.