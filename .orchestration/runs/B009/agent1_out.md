Plan: studied scaffold/examples/YAML, built reusable core+adapter, added B009 constraint/tests, ran full constraints suite.

Physics: tree-level `B+ -> tau+ nu` formula; anchors from `B009.yaml` HFLAV Dec 2025 average and UTfit Summer 2024 SM prediction.

Built: `quarkConstraints/leptonic_tree.py` and `physics_adapters/leptonic_tree.py`. Reused scaffold `load_anchor`, registry, and `ParameterPoint`.

Budget: `4.493733520830467e-05`, from `|BR_exp - BR_SM| + sqrt(sigma_exp^2 + sigma_SM^2)`.

Validation: SM tree formula gives `8.630796690372665e-05`, consistent with YAML SM anchor `8.65e-05`.

NEEDS-HUMAN-PHYSICS: RS uses documented charged-current proxy `A/A_SM = 1 + m_B^2/M_KK^2`; full W/W'/charged-Higgs/lepton/neutrino matching is absent. `B009.yaml` also lacks standalone `f_B` and `|V_ub|` blocks, so diagnostics flag that source gap.

Files changed: `quarkConstraints/leptonic_tree.py`, `flavor_catalog_constraints/physics_adapters/leptonic_tree.py`, `flavor_catalog_constraints/primary/beauty/B009.py`, `tests/constraints/primary/beauty/test_B009.py`.

Tests: `python -m pytest tests/constraints/ -q` -> `509 passed in 14.82s`.