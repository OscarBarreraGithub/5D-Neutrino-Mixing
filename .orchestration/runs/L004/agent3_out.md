1. SHOULD-FIX: [test_L004.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/charged_lepton/test_L004.py:189) recomputes with a hand-copied formula, not the `quarkConstraints.mu_e_conversion` core requested by the review contract. Fix by comparing `Constraint.evaluate(...)` to `quarkConstraints.mu_e_conversion.mu_e_conversion_from_components(...)` with Au inputs and an independently computed L001 dipole BR.

2. Cross-check I ran directly against `quarkConstraints`: constraint CR = `2.79769806116240085e-18`, core CR = `2.79769806116240085e-18`, abs diff = `0.0`; ratio = `3.99671151594628631e-06`; upper envelope = `3.09690586791658228e-18`.

3. Safe/excluded probes: `g_lv_p=1.0e-13` passes with CR `3.00044530722683948e-21`, ratio `4.28635043889548443e-09`; `g_lv_p=3.0e-9` fails with CR `2.70040077650415537e-12`, ratio `3.85771539500593619`.

4. Other checks passed: L004 imports physics only via `physics_adapters.mu_e_conversion`; tracked adapter/core files are text/AST-equal to `HEAD`; anchor missing/mismatched probes raise `AnchorError`; numeric result fields are real floats; complex amplitudes stay in diagnostics; deterministic evaluate confirmed.

5. Tests: `python -m pytest tests/constraints/ -q` -> `678 passed in 15.27s`; targeted L004 -> `9 passed`; registry smoke -> `64` constraints, L004 present, `0` import failures.

CODE-NEEDS-FIXES