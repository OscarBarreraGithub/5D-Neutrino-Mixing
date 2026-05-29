Plan: mirror K001, add the Bd adapter wrapper, load B001 YAML anchors, apply the uncertainty-aware budget, and test against the running core evaluator.

Files changed: `flavor_catalog_constraints/primary/beauty/B001.py`, `tests/constraints/primary/beauty/test_B001.py`, `flavor_catalog_constraints/physics_adapters/deltaf2.py`.

Evaluator used: `quarkConstraints.deltaf2.evaluate_bd_mixing_with_running`, via an append-only adapter wrapper, with `mu_had=2.0 GeV`.

Budget: `3.36989977945937e-14 GeV` for `|M12^NP|`, constructed as `(|Delta m_exp - Delta m_SM| + sqrt(sigma_exp^2 + sigma_SM^2)) / 2`. Sources: B001 YAML HFLAV/PDG `0.5069 +/- 0.0019 ps^-1`, core Bd SM central, HPQCD 2019 Eq. 17 `arXiv:1907.01025`, and `docs/audits/bag_param_inventory.md:36-37`.

Cross-check: safe point `|M12^NP|=2.1778918692437868e-14`, plugin ratio `0.6462779345898485`, core ratio `0.1306473826780916`; excluded point ratio `16.156948364746213`.

Tests: `python -m pytest tests/constraints/ -q` -> `51 passed`.

Diff stat: B001 constraint `428 insertions`, B001 test `247 insertions`, adapter diff currently `244 insertions, 2 deletions` in the allowed adapter file. Review note: the adapter file also contains concurrent append-only wrappers for other constraints already present in the worktree; I did not revert them.