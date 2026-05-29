Plan: load CR001 YAML anchor, build reusable collider resonance core + adapter, wire CR001, test pass/fail and purity.

Physics: CR001 uses mass-limit proxy `ratio = m_limit / m_KK`; source is `CR001.yaml` CMS-B2G-25-009 / arXiv:2603.23454.

Built: `quarkConstraints/collider_resonance.py`, adapter, `primary/collider_rs/CR001.py`, and `test_CR001.py`.

Reused: scaffold `load_anchor`, registry, `ParameterPoint` extras, and point-builder `M_KK` handling.

Budget: HARD `5.5 TeV` CMS 95% CL KK-gluon `ttbar` benchmark exclusion.

Numerical validation: `M_KK=6000 GeV -> 6.0 TeV`, ratio `5.5/6.0 = 0.9166667`, checked against the new core, not the adapter.

NEEDS-HUMAN-PHYSICS: full `sigma(pp -> g_KK) * BR(g_KK -> ttbar)` recast, width, interference, acceptance, and limit curve are not implemented.

Tests: `test_CR001.py` 8 passed; `python -m pytest tests/constraints/ -q` 509 passed.