Plan executed: studied scaffold/K001/K004/B015/B016/B017, implemented B017, added focused tests, ran full constraints.

Physics: active HARD scalar is B017 `R_K central-q2` charged-mode proxy, `1.1 < q^2 < 6.0`, using B016 `B+ -> K+ mu mu` C9/C10 form-factor machinery normalized to SM.

Reused: `rare_b_dilepton` exclusive machinery; B017 loads inclusive/K*/RK* YAML rows as diagnostics. Built: only B017 constraint/test plus append-only adapter export for the 30% proxy constants.

Budget: `0.3979713194305231`, from `B017.yaml` `R_K central-q2` uncertainty plus `|0.949 - 1| + 30%*R_K^SM`.

SM validation: independent trapezoid integral gives `BR(B+ -> K+ mu mu)[1.1,6.0] = 1.851950637006692e-07`; SM `R_K proxy = 1.0`.

NEEDS-HUMAN-PHYSICS: full EW KK/Z/Z', lepton non-universality, K* angular/form-factor, nonlocal charm, and covariance matching remain flagged.

Files changed: `flavor_catalog_constraints/primary/beauty/B017.py`, `tests/constraints/primary/beauty/test_B017.py`, `flavor_catalog_constraints/physics_adapters/rare_b_meson.py`.

Tests: `test_B017.py` 11 passed; `python -m pytest tests/constraints/ -q` 330 passed. Unrelated dirty/untracked files already present were left untouched.