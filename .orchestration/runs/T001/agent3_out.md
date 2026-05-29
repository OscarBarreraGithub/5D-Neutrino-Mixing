1. NIT [flavor_catalog_constraints/primary/top_higgs_ew/T001.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/top_higgs_ew/T001.py:258): `_load_scaffold_value_anchor` temporarily mutates `anchor_scaffold.load_pdg_block`; acceptable here, but a future scaffold helper for list-entry anchors would avoid global patching.

2. Isolation OK: T001 imports no sibling constraint and reaches physics only via [top_fcnc adapter](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/top_higgs_ew/T001.py:48); T001 files are new additions, no pre-existing function bodies modified.

3. Contract OK: `severity=HARD`; missing couplings returns non-crashing `passes=True, predicted=None, ratio=None`; numeric fields are `float`, complex values stay in diagnostics.

4. Anchor OK: budget loaded from YAML via `load_anchor`; active limit `PDG2025:T001:tZc_right = 1.2e-4`, SM `1e-14`; missing and mismatched anchor probes both raised `AnchorError`.

5. Numerical cross-check: manual formula BR `4.817352363171973e-07`, constraint BR `4.817352363171973e-07`, diff `0`, ratio `0.004014460302643311`.

6. Pass/fail points: safe BR `4.179915282578718e-07`, ratio `0.003483262735482264`, passes `True`; excluded BR `0.0002612447051611698`, ratio `2.177039209676416`, passes `False`.

7. Determinism OK: test confirms repeated evaluate equality and unchanged coupling matrices.

8. Registry smoke OK: imported `flavor_catalog_constraints.primary.top_higgs_ew.T001`; registered protocol `True`; registry count `24`; contains `T001=True`.

9. Pytest: `python -m pytest tests/constraints/ -q` -> `244 passed in 4.01s`.

CODE-OK