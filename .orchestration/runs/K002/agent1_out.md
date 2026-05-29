Plan done: studied K001/scaffold/core/YAML, added K002, added tests, ran constraints suite.

Files changed by me:
- `flavor_catalog_constraints/primary/kaon/K002.py`
- `tests/constraints/primary/kaon/test_K002.py`
- additive wrappers in `flavor_catalog_constraints/physics_adapters/deltaf2.py`

Running evaluator used: `quarkConstraints.deltaf2.evaluate_delta_mk_with_running`, via adapter, with `mu_had=2.0 GeV`.

Budget: `Delta m_K^exp / 2 = 1.7419579440705747e-15 GeV`, from YAML `pdg_fit_assuming_cpt` `0.5293 x 10^10 hbar s^-1 -> 3.4839158881411494e-15 GeV`. No SM subtraction. Sources: `docs/quark_scan_assumptions_compact.tex:464-466`, `quarkConstraints/deltaf2.py:809-825`, audit cross-check `docs/audits/bag_param_inventory.md:20`.

Numerical cross-check at reference safe point:
`predicted=1.6276830600654479e-16`, plugin `ratio=9.343985976273966e-02`, direct core `ratio=9.343760390731619e-02`, direct inferred core budget `1.742e-15`.

Pytest:
`test_K002.py`: 8 passed; kaon constraints: 16 passed; full `tests/constraints/ -q`: 51 passed.

Diff stat for K002 files:
`K002.py` 236 insertions; `test_K002.py` 217 insertions. Adapter file is dirty and its current `git diff --stat` includes parallel/pre-existing additions; my K002 changes there are append-only wrappers. Current status also shows unrelated untracked B001/B003/C001 and `.orchestration` files that I did not touch.