Plan: studied K001/K004, C004/C005/C007, scaffold, and C006 YAML; implemented isolated C006 path; verified tests.

Physics: pure-NP `BR(D0 -> e+- mu-+)` HARD bound with SM set to `0.0`; unequal-lepton two-body `C9-C9'` / `C10-C10'` rate.

Sources: `C006.yaml` PDG Live/API S032.40, LHCb 2016, 90% CL.

Reused vs built: reused `rare_charm_dilepton` inputs/CKM/quark Wilson proxy; built new LFV wrapper + adapter for explicit `e-mu` lepton proxy.

Budget: `1.3e-8` branching fraction from `canonical_current_limit`.

Validation: pinned proxy point gives `BR = 3.752990413820220e-13`, independently recomputed in test.

NEEDS-HUMAN-PHYSICS: rigorous RS `e-mu` lepton neutral-current coupling is not standard on `ParameterPoint`; explicit proxy is required.

Files changed: `quarkConstraints/rare_charm_lfv_dilepton.py`, `flavor_catalog_constraints/physics_adapters/rare_charm_lfv_dilepton.py`, `flavor_catalog_constraints/primary/charm/C006.py`, `tests/constraints/primary/charm/test_C006.py`.

Pytest: `tests/constraints/primary/charm/test_C006.py -q` -> 12 passed; `tests/constraints/ -q` -> 407 passed.