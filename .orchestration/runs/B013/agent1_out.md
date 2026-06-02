Plan: mirror B011/B012, add B013, append a B_s→φγ adapter wrapper, and test anchors/numerics/pass-fail/purity.

Physics: C7/C7' exclusive BR scaling uses B013.yaml HFLAV 2024 exp and PDG 2025 3.4e-5 SM/proxy rows; S_phiγ/A_Delta rows are diagnostics only.

Reused/built: reused `quarkConstraints.bsgamma` via `physics_adapters.bsgamma`; no core bsgamma/B011/B012 changes.

Budget: `5.122499e-06` from `|exp-SM| + sqrt(sigma_exp^2 + sigma_SM^2)`; band `2.887750e-05..3.912250e-05`.

SM validation: zero-coupling B013 predicts `3.4e-5`; independent manual LL C7/C8 recomputation matches.

NEEDS-HUMAN-PHYSICS: RS C7/C7' proxy; full exclusive matching/form factors/spectator amplitudes and A_Delta/S_phiγ helicity likelihood are not calculated.

Files changed: [bsgamma.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/physics_adapters/bsgamma.py), [B013.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/secondary/beauty/B013.py), [test_B013.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/secondary/beauty/test_B013.py).

Pytest: B013 `12 passed`; B011/B012/B013 `35 passed`; full `tests/constraints/` `1033 passed in 34.85s`.

Note: `git status` also shows unrelated pre-existing CR013/K021 and adapter changes I did not touch.