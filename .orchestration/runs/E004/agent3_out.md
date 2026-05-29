1. BLOCKER: none.
2. SHOULD-FIX: none.
3. NIT: none.
4. Isolation: `flavor_catalog_constraints/primary/edm_neutrino/E004.py:51` reaches physics only through `physics_adapters.neutron_edm`; no other constraint imports.
5. Contract: `E004.py:398` returns INFO with real floats: `experimental=1.8e-26`, `ratio=0.0`, `budget=1.8e-26`; empty-point path does not crash.
6. Anchor: `E004.py:214` uses `load_anchor`; missing-anchor probe raised `AnchorError`; mismatched primary limit probe raised `AnchorError`.
7. Cross-check: direct YAML recompute gives `limit=1.80000000000000007e-26`, `measured=0`, `stat=1.10000000000000019e-26`, `syst=2.00000000000000008e-27`, `total=1.11803398874989503e-26`, `ratio=0`; `evaluate()` agreed.
8. Pass/fail edge: adapter safe case `ratio=0`, `passes=True`; over-limit case `ratio=2`, `passes=False`; constraint remains non-vetoing via `Severity.INFO`.
9. Determinism: repeated `evaluate()` results equal; point extras unchanged.
10. Tests: `python -m pytest tests/constraints/primary/edm_neutrino/test_E004.py -q` -> 8 passed; `python -m pytest tests/constraints/ -q` -> 769 passed.
11. Registry smoke: `fcc.get("E004")` resolved `PRIMARY edm_neutrino INFO |d_n|`; registry size 74; import failures 0.
12. CODE-OK