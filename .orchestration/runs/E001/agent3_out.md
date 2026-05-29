1. BLOCKER `flavor_catalog_constraints/primary/edm_neutrino/E001.py:106`: mismatched-anchor probe accepted `superseded_acme_2018` with `1.1e-29 e cm`. Fix: assert `experimental.block_key == "canonical_limit"` and add a monkeypatch test.
2. SHOULD-FIX `tests/constraints/primary/edm_neutrino/test_E001.py:41`: numerical expected value is a hardcoded manual `HBARC` multiply, not direct `quarkConstraints.edm` core. Fix: recompute expected via `quarkConstraints.edm.edm_e_cm_from_cp_odd_dipole` or `evaluate_charged_lepton_edm`.
3. NIT `flavor_catalog_constraints/primary/edm_neutrino/E001.py:42`: unused `HBARC_GEV_CM` import. Fix: remove it.

Isolation: E001 imports no sibling constraint; physics is reached only through `flavor_catalog_constraints.physics_adapters.edm`; AST/import check found new E001/adapter/core functions only, no pre-existing function-body edits.

Cross-check: YAML/load_anchor limit `4.10000000000000034e-30`; missing-input path `passes=True predicted=None ratio=None evaluated=False`.

Cross-check: `c=1.0e-16` => constraint/core/direct `1.97326980399999987e-30`, ratio `0.48128531804878044`, passes `True`, point unchanged `True`.

Cross-check: `c=1.0e-15` => constraint/core/direct `1.97326980399999994e-29`, ratio `4.81285318048780475`, passes `False`, point unchanged `True`.

Cross-check: chiral `(3e-16+1e-16j)` => core Im `9.99999999999999979e-17`, prediction `1.97326980399999987e-30`, complex kept in diagnostics `True`.

Determinism: repeated evaluate on same numeric input returned equal results.

Registry smoke: `E001_registered=True`, family `edm_neutrino`, import failures `0`, total constraints `49`.

Pytest: focused E001 `11 passed`; requested `python -m pytest tests/constraints/ -q` => `509 passed in 13.64s`.

CODE-NEEDS-FIXES