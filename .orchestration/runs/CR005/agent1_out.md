Plan: read CR005/scaffold/examples, append adapter wrapper, add constraint/tests, run suite.

Physics: CR005 uses a HARD mass-lower-bound proxy for `pp -> (gamma^(1), Z^(1))_KK -> l+l-`, sourced from `CR005.yaml` CMS2021/ATLAS2019/PDG2024 entries.

Reused vs built: reused `quarkConstraints.collider_resonance` via `physics_adapters.collider_resonance`; appended only an EW KK mass resolver/prediction wrapper.

Budget: active CMS2021 SSM-Z′ dilepton limit, `m > 5.15 TeV`, from CMS JHEP 07 (2021) 208 / arXiv:2103.02708.

Validation: independent core recomputation at `m_KK_EW = 6.0 TeV` gives ratio `5.15 / 6.0 = 0.858333`, pass.

NEEDS-HUMAN-PHYSICS: sigma*BR recast still needs light-quark/lepton couplings, widths, BRs, interference, acceptance, and mass-dependent limit curves; xsec limits are diagnostics only.

Files changed: `flavor_catalog_constraints/physics_adapters/collider_resonance.py`, `flavor_catalog_constraints/primary/collider_rs/CR005.py`, `tests/constraints/primary/collider_rs/test_CR005.py`.

Tests: `test_CR005.py` 8 passed; collider_rs 41 passed; full `python -m pytest tests/constraints/ -q` 697 passed.