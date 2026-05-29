1. BLOCKER: none.
2. SHOULD-FIX: none.
3. NIT: none.
4. Isolation: PASS. AST diff shows only new B→K functions/classes added to `quarkConstraints/rare_b_dilepton.py` and append-only wrappers in `rare_b_meson.py`; no pre-existing bodies removed/changed. B016 reaches physics via `physics_adapters.rare_b_meson`.
5. Contract/anchors: PASS. Missing couplings returns non-crashing `passes=True, predicted=None, ratio=None`; numeric fields are real floats; complex Wilsons stay in diagnostics. Missing/mismatched anchor probes raise `AnchorError`.
6. Independent numerical check: SM full `constraint=5.750185255422402e-07`, direct trapezoid core recompute `5.751444525930657e-07`; low-q² `constraint=1.851950637006692e-07`, direct `1.851950637006572e-07`.
7. Safe/fail check: safe `left=1e-2` gives `constraint=5.618271817599917e-07`, direct `5.619502199525533e-07`, ratio `0.3458347564975091`, PASS. Excluded `left=1e-1` gives `constraint=4.5564003199286436e-07`, direct `4.5573981557011934e-07`, ratio `2.936936008273504`, FAIL.
8. Pytest: `tests/constraints/primary/beauty/test_B016.py -q` -> `10 passed`; `python -m pytest tests/constraints/ -q` -> `244 passed`.
9. Registry smoke: `fcc.get("B016")` -> `B016 beauty HARD BR(B+ -> K+ ell+ ell-)`.
CODE-OK