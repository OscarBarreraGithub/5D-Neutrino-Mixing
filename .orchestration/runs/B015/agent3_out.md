1. No BLOCKER / SHOULD-FIX / NIT findings for B015; no fix required.
2. Isolation: B015 reaches physics via `physics_adapters/rare_b_meson.py`; no other constraint imports. AST check found 0 pre-existing changed/removed defs in `rare_b_meson.py` and `quarkConstraints/rare_b_dilepton.py`.
3. Contract: real float result fields verified; complex Wilsons stay in diagnostics; missing couplings return non-crashing `passes=True`, `predicted=None`, `ratio=None`.
4. Anchor: YAML via `load_anchor` verified; missing anchor and mismatched block both raise `AnchorError`.
5. Cross-check: anchor exp=`1.58e-06`, SM=`1.62e-06`, budget=`9.0678865529319545e-07`.
6. Cross-check NP: manual core/Wilson+kernel rebuild=`1.6077698253698167e-06`; constraint=`1.6077698253698176e-06`; core evaluator=`1.6077698253698176e-06`; abs diff manual=`8.47e-22`.
7. SM/safety: SM pred=`1.62e-06`, pass=True; safe pred=`1.5768349708177282e-06`, ratio=`0.0034903713933742337`, pass=True; excluded pred=`3.1949083244342025e-06`, ratio=`1.7809092725272881`, pass=False.
8. Determinism: repeated evaluate equal=True; input coupling matrix mutated=False.
9. Tests: `python -m pytest tests/constraints/primary/beauty/test_B015.py -q` -> 11 passed.
10. Tests: `python -m pytest tests/constraints/ -q` -> 296 passed.
11. Registry smoke: `fcc.get("B015")` -> `B015:HARD`.

CODE-OK