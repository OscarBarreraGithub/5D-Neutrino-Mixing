1. BLOCKER: none.
2. SHOULD-FIX: none.
3. NIT: none.

Isolation OK: T007 imports no other constraint; physics path is via [T007.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/top_higgs_ew/T007.py:48). AST diff on `physics_adapters/top_fcnc.py` shows no removed/modified defs, only added Higgs proxy defs at lines 520, 547, 589.

Numerical cross-check: YAML budget `3.4000000000000002e-04`; sample `g_ct=1+0.25j`, `g'_ct=0.3j` gives constraint BR `3.4762570502439738e-06`, direct `quarkConstraints.top_fcnc` core BR `3.4762570502439738e-06`, delta `0`, ratio `0.01022428544189404`.

Safe/excluded: safe BR `3.0162750978255732e-04`, ratio `0.8871397346545803`, passes `True`; excluded BR `3.6496928683689444e-04`, ratio `1.0734390789320425`, passes `False`.

Contract/determinism OK: missing-couplings path returns non-crashing pass with `predicted=None`, `ratio=None`; numeric fields are real floats; complex amplitudes stay in diagnostics; repeated evaluate results equal and coupling arrays unchanged.

Anchor probes OK: missing value_id raises `AnchorError`; forced mismatched `load_anchor` block raises `AnchorError`.

Pytest: targeted T007 `9 passed`; required `python -m pytest tests/constraints/ -q` -> `724 passed in 15.64s`. Registry smoke OK: `T007 top_higgs_ew HARD BR(t -> H c)`, registry contains T007.

CODE-OK