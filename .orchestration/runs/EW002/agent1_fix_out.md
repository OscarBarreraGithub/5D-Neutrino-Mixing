Implemented EW002 fixes:
- Routed EW002 list-entry anchors through scaffold `load_anchor` with loud missing/mismatch `AnchorError` probes.
- Reworked EW002 tests to spy on `load_anchor` and validate pass/fail via `Constraint.evaluate()` using hand-computed CKM arithmetic.
- Kept `SOFT` severity and `NEEDS-HUMAN-PHYSICS`; `ckm_unitarity` adapter unchanged.

Pytest:
- `test_EW002.py -q`: 8 passed
- `tests/constraints/ -q`: 127 passed