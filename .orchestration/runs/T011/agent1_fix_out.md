T011 is now standalone, with `A_FB^{0,b}=0.0992 +/- 0.0016` and `A_b=0.923 +/- 0.020` anchors in `T011.yaml` plus T011-local PDG/LEP-SLC snapshots.

`T011.py` now requires T011 anchors through `load_anchor`; legacy T010 `canonical_home` fallback raises `AnchorError`.

SM validation: `A_b=0.935535274431`, `A_FB^{0,b}=0.103278860422`.

Tests: `test_T011.py -q` -> 12 passed.  
Tests: `tests/constraints/ -q` -> 863 passed.