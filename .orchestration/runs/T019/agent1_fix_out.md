Decision: use sqrt-norm, not rename-only.

Changed shared `higgs_lfv` so `yukawa_norm = sqrt(|Y_ij|^2+|Y_ji|^2)` and added explicit `yukawa_norm_squared`.

T018/T019 diagnostics now report `effective_yukawa_norm` consistently, plus `effective_yukawa_norm_squared`.

Verified T019 test point: `2.291e-4 < 1.278e-3`.

Tests: `pytest tests/constraints/ -q` -> `627 passed`.