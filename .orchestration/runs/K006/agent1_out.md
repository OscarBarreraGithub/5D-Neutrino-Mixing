Plan executed: built shared core, adapter, K006 constraint, focused tests, then ran suite checks.

Physics: `BR(K_L -> mu+ mu-)_SD` uses Buras/Isidori `Y_eff` short-distance parametrization; sources: https://arxiv.org/abs/1507.08672 and https://arxiv.org/abs/hep-ph/0311084.

Built: `quarkConstraints/rare_kaon_dilepton.py` and `flavor_catalog_constraints/physics_adapters/rare_kaon_dilepton.py`; did not touch `rare_kaon_snd.py`.

Constraint: `flavor_catalog_constraints/primary/kaon/K006.py`, severity `SOFT` because K006 is long-distance dominated and only the conservative SD allowance is tested.

Budget: `2.5e-9` from `K006.yaml` `supporting_numeric_values[2]`.

SM validation: `BR_SD = 8.209801430715005e-10`, below the K006 SD bound and consistent with the cited `~0.8e-9` SM expectation.

NP gap: diagnostics and docstrings flag `NEEDS-HUMAN-PHYSICS` for full RS EW KK/Z/Z'/muon axial matching; current NP is a documented Z/penguin proxy.

Tests added: `tests/constraints/primary/kaon/test_K006.py`.

Pytest: K006 `11 passed`; kaon constraints `46 passed`; full `tests/constraints/` ran but is blocked by unrelated untracked B016 failure (`243 passed, 1 failed`). With that B016 test ignored: `234 passed`.