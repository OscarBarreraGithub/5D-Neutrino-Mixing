1. BLOCKER: none.
2. SHOULD-FIX: none.
3. NIT: none.

Evidence: CR006 reaches physics via `flavor_catalog_constraints.physics_adapters.collider_resonance`; adapter AST check showed no modified/removed existing functions, only added `kk_charged_current_prediction_from_m_kk_gev`.

Independent numerical cross-check vs `quarkConstraints.collider_resonance`: 7 TeV ratio `0.857142857143` pass; 6 TeV ratio `1.0` pass; 5 TeV ratio `1.2` fail. Constraint and direct core values matched.

Anchor probes: missing value id raised `AnchorError`; mismatched `load_anchor` block raised `AnchorError`.

Tests: CR006 file `8 passed in 4.43s`; requested full slice `python -m pytest tests/constraints/ -q` was `769 passed in 15.78s`. Registry smoke: `registry_count=74`, `import_failures=0`, CR006 registered as `collider_rs` / `HARD`.

CODE-OK