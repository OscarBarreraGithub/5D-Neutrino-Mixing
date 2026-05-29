Plan: studied K001/K004, scaffold, and B011.yaml; built shared C7 core + adapter; added B011 and tests.

Physics: inclusive C7-normalized `BR(Bbar -> X_s gamma)` at `E_gamma > 1.6 GeV`; anchors from B011.yaml HFLAV 2024 and Misiak-Rehman-Steinhauser 2020.

Built: `quarkConstraints/bsgamma.py` and `flavor_catalog_constraints/physics_adapters/bsgamma.py`; B011 imports only the adapter, so B012 can reuse the C7 machinery.

Budget: HARD NP-shift room `3.449509756796392e-05`, from `|3.49e-4 - 3.40e-4| + sqrt((0.19e-4)^2 + (0.17e-4)^2)` in B011.yaml.

SM validation: SM-limit prediction returns `3.40e-4`.

NEEDS-HUMAN-PHYSICS: RS NP is a documented b-s overlap proxy for `C7/C7p`; full KK fermion/Higgs/charged-current/C8/RG matching is not available on `ParameterPoint`.

Files changed: `quarkConstraints/bsgamma.py`, `flavor_catalog_constraints/physics_adapters/bsgamma.py`, `flavor_catalog_constraints/primary/beauty/B011.py`, `tests/constraints/primary/beauty/test_B011.py`.

Tests: `test_B011.py` 10 passed; `python -m pytest tests/constraints/ -q` 137 passed.

Worktree note: unrelated untracked B005/B022/rare_b* and `.orchestration` files are present; I left them untouched.