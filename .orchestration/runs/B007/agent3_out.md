1. BLOCKER [B007.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/secondary/beauty/B007.py:257): mismatched-anchor probe does not fail loudly. Monkeypatched `load_anchor` to return `block_key='wrong_block'`; `_limit_anchor_from_value_row(...)` accepted it. Fix: assign the returned anchor, check `anchor.block_key == block_key`, and raise `AnchorError` on mismatch; add a regression test near [test_B007.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/secondary/beauty/test_B007.py:208).

2. NIT [B007.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/secondary/beauty/B007.py:241): virtual-anchor loading temporarily monkeypatches `anchor_scaffold.load_pdg_block`; works and restores, but is brittle. Fix later with a scaffold helper for nested/list anchors.

Isolation/contract otherwise OK: no other-constraint imports; `rare_b_dilepton.py` has no diff vs `HEAD`; only B007 constraint/test + electron adapter added; numeric result fields are real floats, complex Wilsons stay in diagnostics; missing couplings returns non-crashing HARD pass.

Cross-check numbers: SM manual/core `Bs=8.540154879838497e-14`, `Bd=2.450169556007537e-15`.

Safe point `bs_left=bd_left=10`: direct `Bs=3.910555304412083e-11` ratio `0.004151080111076685`; direct `Bd=2.7386996545281628e-11` ratio `0.010953806618112651`; result pass `True`.

Excluded point `bs_left=bd_left=100`: direct `Bs=4.2463620109893e-09` ratio `0.4517315543605638`; direct `Bd=2.695760074366185e-09` ratio `1.078303037746474`; result pass `False`, active `b_d`.

Determinism probe: repeated `evaluate()` equal; coupling arrays unchanged. Registry smoke: `B007 beauty SECONDARY HARD registered=True import_failures=0`.

`python -m pytest tests/constraints/ -q`: `959 passed in 19.36s`.

CODE-NEEDS-FIXES