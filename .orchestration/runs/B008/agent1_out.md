Plan executed: mirror B007, add tau adapter, add B008 constraint, add focused tests, run full suite.

Physics: Buras `b -> q l l` C10-dominant leptonic formula via `rare_b_dilepton`, with `m_l = m_tau = 1.77686 GeV`.

Reused/built: reused `quarkConstraints/rare_b_dilepton.py`; added append-only adapter `flavor_catalog_constraints/physics_adapters/rare_b_tauonic.py`.

Budgets: `BR(B_s -> tau tau) < 6.8e-3`, `BR(B0 -> tau tau) < 2.1e-3`, both 95% CL from `B008.yaml`.

SM validation: `BR(B_s)=7.737985313074278e-07`, `BR(B0)=2.1909728194059746e-08`, independently recomputed in tests and within sidecar SM uncertainties.

Gap flagged: RS C9/C10 is a documented proxy, with `NEEDS-HUMAN-PHYSICS` in docstring and diagnostics.

Files added: `flavor_catalog_constraints/secondary/beauty/B008.py`, `tests/constraints/secondary/beauty/test_B008.py`, `flavor_catalog_constraints/physics_adapters/rare_b_tauonic.py`.

Pytest: `tests/constraints/secondary/beauty/test_B008.py -q` -> 12 passed; `tests/constraints/ -q` -> 980 passed.

Note: the worktree also contains unrelated pre-existing changes/untracked files outside B008, which I left untouched.