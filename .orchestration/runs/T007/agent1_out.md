Plan executed: studied scaffold/examples, reused `top_fcnc`, appended adapter proxy, added T007 constraint/test, ran suite.

Physics: `BR(t -> H c)` from effective Yukawa `L=-H cbar(y_L P_L+y_R P_R)t`; width reused from `quarkConstraints.top_fcnc`.

Built: append-only Higgs proxy `y_ct=(g_ct/g_s)*(m_t^2/M_KK^2)` in [top_fcnc.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/physics_adapters/top_fcnc.py:509), plus [T007.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/top_higgs_ew/T007.py:1).

Budget: `3.4e-4`, loaded from `T007.yaml` value `PDG2026:T007:tHc_combined`; HARD pure-NP collider bound.

Validation: zero-Yukawa rate `0.0`; independent core recomputation test gives `BR=3.4762570502439738e-6` for the reference proxy point.

NEEDS-HUMAN-PHYSICS: full RS top-Higgs/Yukawa-misalignment and collider recast matching is absent; diagnostics flag this explicitly.

Tests: `test_T007.py` `9 passed`; full `python -m pytest tests/constraints/ -q` `724 passed`.

Files changed: [test_T007.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/top_higgs_ew/test_T007.py:1). Note: the workspace also has unrelated modified/untracked files I did not edit.