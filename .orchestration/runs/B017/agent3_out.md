1. NIT flavor_catalog_constraints/primary/beauty/B017.py:313: `_load_scaffold_list_anchor` temporarily patches `anchor_scaffold.load_pdg_block`; fix later with a scaffold-native list-anchor loader. Non-blocking: restored in `finally`, probes pass.
2. ISOLATION: OK. B017 imports no other constraint and reaches physics only via `physics_adapters.rare_b_meson`; adapter diff is append-only at rare_b_meson.py:329-345, no pre-existing function bodies changed.
3. CONTRACT: OK. Numeric `ConstraintResult` fields are real floats/None; complex C9/C10/amplitudes are diagnostics-only; severity HARD; missing couplings returns non-crashing pass with `predicted=None`, `ratio=None`.
4. ANCHOR: OK. B017 values route through scaffold `load_anchor`; missing observable probe raised `AnchorError`; mismatched block probe raised `AnchorError`.
5. NUMERICS: anchor `R_K=0.949`, q2 `[1.1,6.0]`, budget `0.3979713194305231`; SM BR/core denominator `1.851950637006692e-07`, proxy pred `1.0`, ratio `0.12814993822413753`, pass.
6. NUMERICS: NP check constraint pred `0.9702300255800139`, direct `quarkConstraints.rare_b_dilepton.evaluate_b_to_k_mumu` expected `0.9702300255800139`, diff `0.0`, BR `1.7968181139159258e-07`, pass.
7. NUMERICS: safe point pred/core `0.9770592716646653`, ratio `0.07050576334198348`, pass; excluded point pred/core `0.4463862703096126`, ratio `1.2629395766750287`, fail.
8. DETERMINISM: OK. Repeated evaluate on same point returned equal results; coupling arrays unchanged.
9. TESTS: `python -m pytest tests/constraints/primary/beauty/test_B017.py -q` -> 11 passed.
10. TESTS: `python -m pytest tests/constraints/ -q` -> 352 passed in 5.20s.
11. REGISTRY: smoke OK; 34 constraints discovered, B017 present in beauty, import failures 0.
CODE-OK