Plan executed: read CR008 YAML/examples, implemented one isolated constraint, added tests, ran full suite.

Physics: `m_T = M_KK` mass-limit proxy for `pp -> T Tbar`, singlet `B(Wb:Ht:Zt)=1/2:1/4:1/4`; source `CR008.yaml`.

Reused machinery: existing `collider_resonance` adapter/core mass-vs-limit recast; no new physics machinery and no CR001-CR007 edits.

Budget: HARD lower limit `m_T > 1.36 TeV` from ATLAS 2024 `arXiv:2401.17165`; CMS `1.48 TeV` all-decay envelope kept in diagnostics.

Numerical validation: `M_KK=1600 GeV` gives predicted `1.600 TeV`, ratio `1.36/1.60 = 0.85`, passing and matching the core recomputation.

Gap flagged: σ×BR/BR-mixture/width/acceptance recast is `NEEDS-HUMAN-PHYSICS` in docstring and diagnostics.

Files changed: added `flavor_catalog_constraints/primary/collider_rs/CR008.py` and `tests/constraints/primary/collider_rs/test_CR008.py`.

Tests: `test_CR008.py` 8 passed; contract 13 passed; full `python -m pytest tests/constraints/ -q` 821 passed.