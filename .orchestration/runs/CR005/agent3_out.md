1. BLOCKER: none found.
2. SHOULD-FIX: none found.
3. NIT: none found.
4. Isolation OK: no CR001-CR004 diffs; adapter AST check changed_existing_funcs=[], removed=[], added `resolve_kk_ew_mass_gev`, `kk_ew_dilepton_prediction_from_m_kk_gev`.
5. Contract OK: HARD severity, real float numeric fields, missing extras returns non-crashing result.
6. Anchor OK: active YAML limit is 5.15 TeV; `load_anchor` path used; missing-value probe and mismatched-block probe both raised `AnchorError`.
7. Independent cross-check: 6000 GeV -> constraint 6.0 TeV, core 6.0 TeV, exp/budget 5.15 TeV, ratio 0.8583333333333334, pass=True; agrees.
8. Safe/excluded: 6150 GeV -> ratio 0.8373983739837398 pass=True; 4150 GeV -> ratio 1.2409638554216866 pass=False.
9. Determinism: 6500 GeV repeated result equal, extras unchanged, ratio 0.7923076923076924.
10. Tests: `test_CR005.py -q` -> 8 passed in 4.40s; `python -m pytest tests/constraints/ -q` -> 724 passed in 15.78s.
11. Registry smoke: `fcc.get("CR005")` -> CR005 collider_rs HARD; registry count 69, contains CR005=True.
CODE-OK