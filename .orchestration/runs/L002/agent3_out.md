1. BLOCKER `flavor_catalog_constraints/physics_adapters/lfv_three_body.py:171` + `quarkConstraints/lfv_three_body.py:445`: caller `initial_flavor`/`final_flavor` can override L002’s mu/e defaults. Probe accepted tau->mu input while reporting `BR(mu -> 3e)`. Fix: pin/reject non-`mu,e` in the adapter or L002 and add a regression.

2. SHOULD-FIX `tests/constraints/primary/charged_lepton/test_L002.py:52`: numerical test uses a local copied formula, not the `quarkConstraints.lfv_three_body` core evaluator required by the review contract. Fix: compute expected via `lfv_three_body_contact_amplitudes` + `lfv_three_body_from_components` and direct L001 dipole core, not the adapter.

Isolation: L002 imports only anchors/base/registry/new adapter; adapter is the only path to `quarkConstraints.lfv_three_body`. Relevant tracked diff empty.

Contract/anchors: real numeric fields OK, complex values stay in diagnostics; missing couplings returns `predicted=None`, `ratio=None`, `passes=True`; missing and mismatched anchors raised `AnchorError`.

Direct-core cross-check: predicted `1.95364112289741326e-14`, core `1.95364112289741326e-14`, diff `0`; ratio `1.95364112289741314e-02`. Components: dipole `6.94135861145652961e-16`, contact `1.88422753678284811e-14`.

Safe point: predicted/core `6.75879658008492711e-15`, ratio `6.75879658008492722e-03`, pass `True`. Excluded point: predicted/core `2.00000000000000033e-10`, ratio `200.00000000000003`, pass `False`.

Determinism: same input produced equal results; point extras unchanged. Registry smoke: `L002 charged_lepton HARD BR(mu -> 3e)`.

Pytest: `test_L002.py` 10 passed; `python -m pytest tests/constraints/ -q` 509 passed.

CODE-NEEDS-FIXES