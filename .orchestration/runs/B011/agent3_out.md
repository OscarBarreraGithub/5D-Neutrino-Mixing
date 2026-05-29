1. SHOULD-FIX tests/constraints/primary/beauty/test_B011.py:165 - numerical cross-check is manual formula-only; add a test that imports `quarkConstraints.bsgamma` and recomputes via `compute_bsgamma_wilsons` + `branching_fraction_from_c7` or `evaluate_inclusive_bsgamma`.

2. Cross-check numbers: B011 constraint pred `3.18080418975069e-4`; direct core pred `3.18080418975069e-4`; abs delta `0`; C7_NP `(0.01-0.003j)`, C7p_NP `(0.003+0.002j)`, ratio_to_SM `0.935530644044321`.

3. Pass/fail probe: SM `3.40000000000000e-4`, ratio `0`, pass; safe point `3.17999480609418e-4`, ratio `0.637786843398116`, pass; excluded point `3.97922437673130e-5`, ratio `8.70291077279034`, fail.

4. Contract/isolation: B011 reaches physics only through `flavor_catalog_constraints.physics_adapters.bsgamma`; B011 files are new only; result numeric fields are real floats/None and complex values stay in diagnostics. Missing extra returns non-crashing HARD result with `predicted=None`, `ratio=None`, `passes=True`.

5. Anchor probes: YAML anchors loaded via `load_anchor`; missing sidecar raises `AnchorError`; mismatched units raises `AnchorError`. Registry smoke: `B011 beauty PRIMARY HARD`.

6. Tests run: `python -m pytest tests/constraints/primary/beauty/test_B011.py -q` -> `10 passed`; `python -m pytest tests/constraints/ -q` -> `157 passed`.

CODE-NEEDS-FIXES