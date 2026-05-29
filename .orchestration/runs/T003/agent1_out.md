Plan executed: inspected examples/T001/T002/T003 YAML, appended adapter proxy, added T003 constraint/test, ran suite.

Physics: `BR(t -> c gamma)` photon dipole, `Gamma = 0.5 * alpha_em * m_t * (|lambda_L|^2 + |lambda_R|^2)`, pure-NP vs collider bound.

Reused/built: reused `top_fcnc` photon-dipole rate; added adapter-only RS proxy `lambda_ct=(g_ct/g_s)*(m_t^2/M_KK^2)`. No Z-FCNC functions changed.

Budget: HARD `CMS2024:T003:tcgamma = 1.51e-5` from `T003.yaml` / CMS-TOP-21-013.

Validation: SM anchor `4.6e-14`; independent manual test-point recomputation `BR_NP = 6.068421034828811e-6`.

NEEDS-HUMAN-PHYSICS: full RS EM dipole matching needs loop/Yukawa/EW-KK/collider-recast inputs; flagged in docstring and diagnostics.

Files changed: `flavor_catalog_constraints/physics_adapters/top_fcnc.py`, `flavor_catalog_constraints/primary/top_higgs_ew/T003.py`, `tests/constraints/primary/top_higgs_ew/test_T003.py`.

Tests: T003 targeted `9 passed`; full `python -m pytest tests/constraints/ -q` -> `307 passed`.

Existing unrelated dirty/untracked worktree entries were left untouched.