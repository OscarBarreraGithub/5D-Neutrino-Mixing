Plan: studied K001/K004/scaffold/K005 YAML, appended neutral rare-kaon evaluator, added K005 constraint, added focused tests, ran full constraint suite.

Files changed for K005:
`quarkConstraints/rare_kaon_snd.py`, `flavor_catalog_constraints/physics_adapters/rare_kaon.py`, `flavor_catalog_constraints/primary/kaon/K005.py`, `tests/constraints/primary/kaon/test_K005.py`.

Evaluator used: `evaluate_klong_to_pi0_nunu` via `klong_pi0_nunu_from_couplings`; QCD `_with_running` is not applicable for this Delta-S=1 semileptonic `kappa_L` observable, and diagnostics mark `qcd_running_applied=False`.

Budget: `2.2e-9`, the KOTO 2025 90% CL total-BR upper limit from `K005.yaml experimental_inputs[0]`; ratio is `BR_total / limit`, no SM subtraction. SM anchor: Buras-Venturini `2.94e-11 +/- 0.15e-11`.

Cross-checks against direct core call:
SM `BR=2.9537598934305903e-11`; ref point `BR=2.9707328694509226e-11`, ratio `0.013503331224776921`; excluded point `BR=5.206881891456275e-08`, ratio `23.66764496116489`.

Tests:
`test_K005.py`: `9 passed`.
`python -m pytest tests/constraints/ -q`: `89 passed`.

Scoped diff stat: `rare_kaon.py` +38, `rare_kaon_snd.py` +141; plus new `K005.py` and `test_K005.py`. Worktree also contains unrelated/concurrent `deltaf2`, B002, and C002 changes not made for K005.

Reviewer note: `K005.yaml` has no standalone `kappa_L` field, so the core uses the standard Buras `kappa_L` constant in the same style as K004 constants; YAML supplies the KOTO limit and SM validation anchors.