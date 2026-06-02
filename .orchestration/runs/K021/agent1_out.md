Plan executed: studied K020/adapter/YAML, appended neutral helper, added K021, added tests, ran full suite.

Physics: pure-NP HARD bound for `K_L -> pi0 e+- mu-+`, SM catalog rate `0.0`; uses K_L mass/lifetime, pi0 mass, neutral K_l3 `f_+(0)`, charge-summed factor `2`.

Reused vs built: reused K020 LFV semileptonic machinery; appended neutral-mode adapter helpers only.

Budget: `BR < 7.6e-11` from `PDG2025:K021:KL_pi0emu_limit` in `K021.yaml`; KTeV row retained as provenance.

Validation point: independent recomputation matched adapter; `BR = 1.219313981243497e-11`, ratio `0.16043605016361803`.

Gap flagged: `NEEDS-HUMAN-PHYSICS` for explicit e-mu proxy and full K_L charge-orientation CP matching.

Files changed: `flavor_catalog_constraints/physics_adapters/rare_kaon_lfv_semileptonic.py`, `flavor_catalog_constraints/secondary/kaon/K021.py`, `tests/constraints/secondary/kaon/test_K021.py`.

Tests: `python -m pytest tests/constraints/secondary/kaon/test_K021.py -q` -> 11 passed; `python -m pytest tests/constraints/ -q` -> 1033 passed.