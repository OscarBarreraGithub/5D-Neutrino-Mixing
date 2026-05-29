1. SHOULD-FIX `flavor_catalog_constraints/primary/beauty/B032.py:42`: B032 bypasses `load_anchor`, using `load_pdg_block`/`build_anchor` directly. Fix by virtualizing the selected list entries through `load_anchor` as in other list-anchor constraints, and update `tests/constraints/primary/beauty/test_B032.py:108` to spy/mismatch `load_anchor`.

Numerical cross-check: YAML gives `A_CP(K+ pi0)=2.7%`, `A_CP(K+ pi-)=-8.31%`, so `Delta A_CP=11.01% = 0.1101`; `sigma=sqrt(1.2^2+0.31^2)% = 1.23939501370628% = 0.0123939501370628`. Constraint returned `experimental=0.1101`, `budget=0.1101`, `ratio=1`, `passes=True`, `severity=INFO`; agreement `True`.

Adapter pass/fail probe: safe room `ratio=1`, `passes=True`; overfilled room `ratio=2`, `passes=False`. Missing/mismatched anchor probes raised `AnchorError`.

B032 tests: `8 passed in 4.33s`. Full requested run: `python -m pytest tests/constraints/ -q` -> `769 passed in 15.75s`. Registry smoke: `B032 beauty PRIMARY INFO`.

CODE-NEEDS-FIXES