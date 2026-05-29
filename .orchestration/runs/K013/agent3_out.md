1. BLOCKER: none.
2. SHOULD-FIX: none.
3. NIT: none.

Isolation OK: K013 imports no other constraint and reaches physics only through `flavor_catalog_constraints.physics_adapters.radiative_kaon`; `deltaf2.py`/`base.py` diffs are empty for this review.

Numerical check: YAML BR `1.273e-06`, uncertainty `3.300e-08`; `evaluate()` returns experimental `1.273e-06`, budget `1.273e-06`, ratio `1.0`. Naive BR/uncertainty `38.5757575758`. Synthetic adapter check: full-room ratio `1.0` passes, half-room ratio `2.0` fails. Values agree.

Anchor probes: missing value and mismatched selected block both raise `AnchorError`. Registry smoke: `fcc.get("K013")` imports/evaluates as `kaon INFO`.

Tests: targeted K013 `9 passed`; full `python -m pytest tests/constraints/ -q` ran `769 passed in 18.12s`.

CODE-OK