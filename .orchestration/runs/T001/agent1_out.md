Plan executed: studied scaffold/examples, built shared core+adapter, added T001, added tests, ran suite.

Physics: `t -> c Z` uses vector FCNC convention `g/(2 c_W) Z_mu cbar gamma^mu (X_L P_L + X_R P_R)t`; RS proxy is `X_ct=(g_ct/g_s) m_Z^2/M_KK^2`.

Built machinery: new `quarkConstraints/top_fcnc.py` plus adapter `flavor_catalog_constraints/physics_adapters/top_fcnc.py`; no direct physics import from constraint.

Budget: HARD limit `BR(t -> c Z) < 1.2e-4`, from `T001.yaml` PDG2025 right-handed benchmark; left/right and ATLAS/CMS limits retained in diagnostics.

SM/validation: SM anchor `1.0e-14`; independent manual proxy check gives `4.817352363171973e-7` for the controlled test point.

Gap: `NEEDS-HUMAN-PHYSICS` flagged in docstring and diagnostics for full RS EW KK/Z/Z'/SMEFT/collider matching.

Files changed: `quarkConstraints/top_fcnc.py`, `flavor_catalog_constraints/physics_adapters/top_fcnc.py`, `flavor_catalog_constraints/primary/top_higgs_ew/T001.py`, `tests/constraints/primary/top_higgs_ew/test_T001.py`.

Tests: `python -m pytest tests/constraints/ -q` -> `213 passed`.