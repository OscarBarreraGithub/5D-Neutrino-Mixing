1. BLOCKER: none found.
2. SHOULD-FIX: none found.
3. NIT: none found.

Evidence: B008 adds only `B008.py`, `rare_b_tauonic.py`, and `test_B008.py`; `rare_b_electronic.py` and `quarkConstraints/rare_b_dilepton.py` are unchanged. Constraint reaches physics only via `rare_b_tauonic`.

Anchor probe: `load_anchor` calls were `values[0]`, `values[1]`, `auxiliary_theory_inputs[0]`, `auxiliary_theory_inputs[1]`; forced mismatch and missing row both raised `AnchorError`.

Numerical cross-checks:
SM/core/constraint agree: `Bs=7.737985313074278e-07`, `Bd=2.190972819405975e-08`.
Safe point `bs_left=bd_left=10`: `Bs=3.543240132909352e-04`, `Bd=2.448980107876779e-04`, ratio `0.1166075289465133`, passes.
Excluded point `bs_left=bd_left=40`: `Bs=6.073477595293012e-03`, `Bd=3.867132575881579e-03`, ratio `1.8414811313721804`, fails.

Contract/determinism: public numeric fields are real floats, complex Wilsons stay in diagnostics, missing couplings returns non-crashing HARD/pass result, repeated evaluation is identical.

Tests: `python -m pytest tests/constraints/secondary/beauty/test_B008.py -q` -> `12 passed`; `python -m pytest tests/constraints/ -q` -> `991 passed in 21.96s`. Registry smoke: `B008`, `SECONDARY`, `beauty`.

CODE-OK