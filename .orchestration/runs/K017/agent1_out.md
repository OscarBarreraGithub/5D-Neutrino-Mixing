Plan: studied K001/K004/B009 + K017.yaml, extended leptonic_tree, added K017 + tests, ran constraints suite.
Physics: R_K uses helicity-suppressed K_l2 tree ratio with K017.yaml SM radiative normalization.
Proxy: electron-only amplitude shift `1 + m_K^2/M_KK^2`, muon SM; flagged NEEDS-HUMAN-PHYSICS.
Reused/built: reused and extended `quarkConstraints/leptonic_tree.py`; exposed through leptonic_tree adapter.
Budget: `2.005538513813748e-07` from K017.yaml `canonical_ratio + sm_prediction`.
Validation: tree ratio `2.5689629165499624e-05`, radiative multiplier `0.9642023183917868`, SM `2.477e-05`.
Files changed: `quarkConstraints/leptonic_tree.py`, `physics_adapters/leptonic_tree.py`, `primary/kaon/K017.py`, `test_K017.py`.
Tests: `python -m pytest tests/constraints/ -q` -> `724 passed in 15.09s`.

Unrelated worktree changes were present/appeared for other constraint IDs and adapters; I left them untouched.