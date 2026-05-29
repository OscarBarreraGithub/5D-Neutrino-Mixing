Implemented CR006.

- Physics: `pp -> W_KK^(1) -> ell nu, tb`, using `m_WKK = kk_ew_mass_gev` or fallback `quark_mass_basis_couplings.M_KK`; HARD ratio is `m_limit / m_WKK`.
- Budget/source: active `PDG2025:CR006:Wprime_SSM_enu_mass_lower_bound`, `6.0 TeV`, PDG 2025 quoting ATLAS `arXiv:1906.05609`. CMS `5.6/5.7 TeV` and `tb` `4.3/3.9 TeV` kept in diagnostics.
- Reuse: reused `quarkConstraints.collider_resonance` through the adapter; appended only the charged-current proxy helper. No CR001-CR005 funcs touched.
- Validation: independent core recomputation at `m_WKK = 6.0 TeV` gives ratio `1.0` and pass.
- Gap: `sigma*BR`, width, interference, acceptance, and limit-curve recast flagged `NEEDS-HUMAN-PHYSICS`.

Files changed: `flavor_catalog_constraints/primary/collider_rs/CR006.py`, `flavor_catalog_constraints/physics_adapters/collider_resonance.py`, `tests/constraints/primary/collider_rs/test_CR006.py`.

Tests: `test_CR006.py` 8 passed; `python -m pytest tests/constraints/ -q` 769 passed.