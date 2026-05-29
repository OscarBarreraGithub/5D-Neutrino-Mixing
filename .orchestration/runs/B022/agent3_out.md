1. BLOCKER [flavor_catalog_constraints/primary/beauty/B022.py:35](</n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/beauty/B022.py:35>) / [test_B022.py:135](</n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/beauty/test_B022.py:135>): B022 uses `load_pdg_block`, not `load_anchor`. Probe: patching `anchors.load_anchor` to raise still let `B022.Constraint()` construct with `sm=5.58e-06`, `exp=2.3e-05`. Fix: route value-bearing anchors through `load_anchor` and add a test that spies/forces mismatched `load_anchor`.

2. Numerical cross-check, constraint vs independent core: SM `5.58e-06` vs `5.58e-06`, ratio `2.7160162886`, pass `False`; safe left=-0.7 `2.23318031698e-05`, ratio `0.1041810261`, pass `True`; excluded left=-2.0 `8.30870917308e-05`, ratio `8.4859886277`, pass `False`; `X_t=1.48122870767`.

3. Contract/isolation: result numeric fields are real floats, complex amplitudes are diagnostics, missing couplings returns non-crashing pass result; B022 imports adapter only and no sibling constraint. B022-scoped files are new additions, no pre-existing function body modified.

4. Determinism: repeated `evaluate()` matched exactly and did not mutate coupling matrices.

5. Tests: `python -m pytest tests/constraints/primary/beauty/test_B022.py -q` -> `10 passed`; `python -m pytest tests/constraints/ -q` -> `157 passed`. Registry smoke: `B022 beauty HARD` imported.

CODE-NEEDS-FIXES