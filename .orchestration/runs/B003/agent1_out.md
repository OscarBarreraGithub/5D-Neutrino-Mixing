Plan: studied K001/scaffold/B003/core, added B003 + tests, verified direct running-evaluator behavior and full constraints suite.

Files changed: `flavor_catalog_constraints/primary/beauty/B003.py`, `tests/constraints/primary/beauty/test_B003.py`, append-only Bs wrappers in `flavor_catalog_constraints/physics_adapters/deltaf2.py`.

Running evaluator: `quarkConstraints.deltaf2.evaluate_bs_mixing_with_running`, via adapter `bs_mixing_from_wilsons_with_running`.

Budget: `2.635167648629676e-13 GeV` for `|M12^NP|`.
Construction: `(|Delta m_exp - Delta m_SM(core)| + sqrt(sigma_exp^2 + sigma_SM^2))/2`, with `sigma_SM = 2*(sigma_fBs√B/fBs√B)*Delta m_SM`.
Sources: `B003.yaml` `canonical_hflav_recommended` + `flag_2024_bmixing`; `docs/audits/bag_param_inventory.md:38-45`.

Cross-check safe point `1e-3`: predicted `4.760385525086220e-14`, ratio `0.1806482986978717`, pass `True`; direct core ratio/budget `0.008145765785568481 / 5.844e-12`.
Cross-check large point `1e-2`: predicted `4.760385525086221e-12`, ratio `18.06482986978718`, pass `False`.

Tests: `python -m pytest tests/constraints/ -q` -> `51 passed in 2.36s`.

Diff stat scoped to B003 new files: `B003.py | 457 +`, `test_B003.py | 290 +`; adapter currently `227 insertions, 2 deletions` in shared worktree.

Reviewer focus: B003 intentionally overrides the core’s legacy generous Bs budget (`5.844e-12 GeV`) with the uncertainty-aware SM-vs-HFLAV room required by the prompt; core legacy budget/ratio are preserved in diagnostics.