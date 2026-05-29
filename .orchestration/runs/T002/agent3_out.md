1. SHOULD-FIX tests/constraints/primary/top_higgs_ew/test_T002.py:187: the numerical cross-check uses a local manual formula, not `quarkConstraints.top_fcnc` core as required; fix by comparing against `compute_top_z_fcnc_proxy` + `z_vector_branching_fraction` or another direct core evaluator.

Isolation/contract evidence: `T002.py` imports physics only via `flavor_catalog_constraints.physics_adapters.top_fcnc`; scoped diff/status shows only new `T002.py` and `test_T002.py`, with adapter/core unchanged. Missing-extra path returns non-crashing real fields.

Actual core cross-check: constraint `predicted=4.817352363171973e-07`, core `BR=4.817352363171973e-07`, ratio `0.007769923166406408`, diff `0.0`; budget `6.2e-05`, SM `8e-17`.

Safe/excluded: left=1 gives `BR=4.1799152825787175e-07` pass; left=25 gives `BR=0.00026124470516116985`, ratio `4.213624276793062`, fail.

Tests: `python -m pytest tests/constraints/primary/top_higgs_ew/test_T002.py -q` -> 9 passed; `python -m pytest tests/constraints/ -q` -> 296 passed. Registry smoke: `fcc.get("T002") -> top_higgs_ew/HARD`.

CODE-NEEDS-FIXES