1. BLOCKER: none.
2. SHOULD-FIX: none.
3. NIT: none.
4. Isolation: CR009 adds `CR009.py`, `test_CR009.py`, and an append-only collider adapter block; `collider_resonance.py:381` onward only. No CR001-CR008 function bodies modified; CR009 reaches physics via the adapter.
5. Contract/anchor: registry smoke OK, severity `HARD`; missing input returns non-crashing `passes=True`, `predicted=None`, `ratio=None`, `experimental=budget=36.4` floats. Missing/mismatched anchor probes raised `AnchorError`.
6. Numerical cross-check: active YAML limit is `PDG2025:CR009:CMS_RR_constructive_range_endpoint = 36.4 TeV`. For `kk_ew_mass_gev=40000`, constraint gives `predicted=40.0`, `ratio=0.91`, pass; independent `quarkConstraints.collider_resonance` gives the same.
7. Safe/excluded: `37.4 TeV -> ratio=0.9732620321`, pass; `35.4 TeV -> ratio=1.0282485876`, fail.
8. Determinism: same input produced identical results and left `point.extras` unchanged; numeric result fields were real floats.
9. Tests: `python -m pytest tests/constraints/ -q` -> `910 passed in 17.28s`.
CODE-OK