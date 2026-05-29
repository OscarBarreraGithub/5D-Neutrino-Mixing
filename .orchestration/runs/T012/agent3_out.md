1. SHOULD-FIX [test_T012.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/top_higgs_ew/test_T012.py:212) - safe/excluded tests only assert pass/fail through the constraint; the independent recomputation at [test_T012.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/top_higgs_ew/test_T012.py:125) is SM-only and does not exercise shifted Zcc via `quarkConstraints.zpole`. Fix: add a test that manually computes `dgL/dgR`, calls core `shifted_couplings` + `evaluate_quark_pseudo_observables`, and asserts `R_c/A_c/ratio`.

Isolation: OK. T012 imports physics only via `physics_adapters.zpole_charm`; no tracked diff to `T010.py`, `physics_adapters/zpole.py`, `quarkConstraints/zpole.py`, or `deltaf2.py`.

Contract/anchor: OK. Numeric fields are real `float`; severity `HARD`; missing couplings return SM pass; YAML anchors route through `load_anchor`; mismatch probe raises `AnchorError`.

Independent numbers: SM `R_c=0.1721`, `A_c=0.667577445962`; YAML `A_c=0.670`; SM selected `A_c`, ratio `0.0897242236336`, pass `True`.

Independent safe point: `R_c=0.1721`, `A_c=0.667577445962`, `dgL=0`, `dgR=0`, ratio `0.0897242236336`; constraint agrees, pass `True`.

Independent excluded point: `R_c=0.191452504897`, `A_c=0.708164612402`, `dgL=0.0277172613125`, `dgR=0`, selected `R_c^0`, ratio `6.45083496578`; constraint agrees, pass `False`.

Determinism: OK, repeated results equal; numeric result field types all `float`.

Tests: `test_T012.py` 10 passed; `python -m pytest tests/constraints/ -q` 296 passed. Registry smoke: `T012 top_higgs_ew PRIMARY HARD`.

CODE-NEEDS-FIXES