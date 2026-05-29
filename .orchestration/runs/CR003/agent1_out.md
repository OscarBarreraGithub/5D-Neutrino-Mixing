Plan: read scaffold/YAML, implement CR003, add focused tests, run constraints suite.

Physics: CR003 is VLQ charge-2/3 `T Tbar`; proxy is `m_T = M_KK`, compared to mass lower bound from `CR003.yaml`.

Reused vs built: reused `quarkConstraints/collider_resonance.py`; appended adapter-only charge-2/3 VLQ proxy wrappers.

Budget: active HARD limit is `1.48 TeV` from `CMS2023:CR003:T_all_third_generation_decay_mixtures` / CMS arXiv:2209.07327. Pure `Wb/Zt/Ht` limits stay diagnostic.

Numerical validation: independent core check at `M_KK=2000 GeV` gives `m_T=2.0 TeV`, ratio `1.48/2.0 = 0.74`, passes.

NEEDS-HUMAN-PHYSICS: full `sigma*BR`, branching-fraction simplex, widths, acceptance, and mass-dependent limit curve are not available on `ParameterPoint`.

Files changed: `flavor_catalog_constraints/primary/collider_rs/CR003.py`, `tests/constraints/primary/collider_rs/test_CR003.py`, adapter append in `physics_adapters/collider_resonance.py`.

Tests: focused `21 passed`; full `python -m pytest tests/constraints/ -q` -> `590 passed`. `ruff` was unavailable in this env; compile check passed.