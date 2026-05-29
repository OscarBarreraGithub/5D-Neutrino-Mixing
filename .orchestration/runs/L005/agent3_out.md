1. SHOULD-FIX `tests/constraints/primary/charged_lepton/test_L005.py:75`: numerical test hand-reimplements the formula instead of calling the underlying `quarkConstraints.mu_e_conversion` core required by the review contract. Fix: compute expected via `mu_e_conversion_from_components(..., nuclear_inputs=titanium_nuclear_inputs())`.

2. Cross-check: constraint predicted `2.893383442977269e-18`; direct core predicted `2.893383442977269e-18`; abs/rel diff `0.0/0.0`. Upper interval `3.205936297453907e-18` agrees.

3. Safe/fail: safe `g_lv_p=1e-13` predicted `2.502845239523481e-21`, ratio `5.820570324473211e-10`, pass `True`; fail `g_lv_p=1e-8` predicted `2.502845239523480e-11`, ratio `5.820570324473211`, pass `False`.

4. Contract/isolation evidence: L005 imports only anchors/base/registry and `physics_adapters.mu_e_conversion`; AST body check vs `HEAD` for adapter/core/deltaf2: added/removed/changed `0/0/0`; missing input returns non-crashing unevaluated result.

5. Anchor evidence: runtime anchor from `L005.yaml` is `4.3e-12`; missing anchor and mismatched target probes both raise `AnchorError`.

6. Determinism: repeated `evaluate()` results equal; input extras unchanged.

7. Registry smoke: `fcc.get("L005")` resolves `flavor_catalog_constraints.primary.charged_lepton.L005`, `PRIMARY/HARD`.

8. Pytest: `tests/constraints/primary/charged_lepton/test_L005.py -q` = `9 passed`; `python -m pytest tests/constraints/ -q` = `724 passed`.

CODE-NEEDS-FIXES