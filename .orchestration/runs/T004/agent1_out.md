Plan executed: T004 is implemented as the T003 up-quark analogue with YAML-loaded anchors, a pure-NP HARD ratio, and tests.

Physics: reused the shared top-FCNC photon dipole via `physics_adapters.top_fcnc`; proxy is `lambda_ut=(g_ut/g_s)*(m_t^2/M_KK^2)`.

Budget: `8.5e-06` from `ATLAS2023:T004:tugamma_left` in `T004.yaml`; CMS/PDG/right-handed limits are diagnostics context.

SM/validation: no numeric SM central value in T004.yaml, so `sm_prediction=None`; independent core check gives BR `6.068421034828811e-06`.

NEEDS-HUMAN-PHYSICS remains flagged for full RS electromagnetic dipole matching and collider recast.

Files changed: [T004.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/top_higgs_ew/T004.py), [test_T004.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/top_higgs_ew/test_T004.py).

Tests: `python -m pytest tests/constraints/ -q` -> 361 passed.