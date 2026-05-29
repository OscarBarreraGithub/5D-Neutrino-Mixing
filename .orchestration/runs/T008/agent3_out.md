1. CODE-OK: No BLOCKER / SHOULD-FIX / NIT findings for T008. Isolation verified: no T007/top_fcnc/deltaf2 diff; T008 reaches physics only via `physics_adapters.top_fcnc` at `T008.py:48`.
2. Contract OK: real float fields confirmed; complex Yukawas only in diagnostics; missing extra returns non-crashing pass result: `predicted=None`, `ratio=None`, `experimental=0.00019`, `sm=0.0`.
3. Anchor OK: active YAML/load_anchor limit `PDG2026:T008:tHu_headline = 1.9e-4`; missing anchor raises `AnchorError`; mismatched load_anchor block raises `AnchorError`.
4. Independent cross-check OK: suppression `0.0033162241777777773`; core BR `3.4762570502439738e-06`; constraint BR `3.4762570502439738e-06`; delta `0`; ratio `0.018296089738126176`.
5. Safe/excluded OK: left `7.0` gives BR `1.4779747979345306e-04`, ratio `0.7778814725971214`, pass; left `8.0` gives BR `1.930416062608367e-04`, ratio `1.0160084540044036`, fail.
6. Determinism OK: repeated evaluate results equal; input coupling matrices unchanged.
7. Tests: `python -m pytest tests/constraints/primary/top_higgs_ew/test_T008.py -q` -> `9 passed`; `python -m pytest tests/constraints/ -q` -> `769 passed in 16.32s`.
8. Registry smoke: `T008` registered as `top_higgs_ew`, `HARD`; import failures `0`.

CODE-OK