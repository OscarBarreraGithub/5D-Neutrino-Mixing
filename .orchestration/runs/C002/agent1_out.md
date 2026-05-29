Plan: load C002 anchors, evaluate QCD-run D0 `M12^NP`, compare CP-odd fraction to HFLAV 95% CPV room, test against the direct running evaluator.

Files changed: `flavor_catalog_constraints/primary/charm/C002.py`, `tests/constraints/primary/charm/test_C002.py`, plus an additive D0 complex-amplitude helper in `flavor_catalog_constraints/physics_adapters/deltaf2.py`. Core `quarkConstraints/deltaf2.py` untouched.

Running evaluator used: `d0_mixing_from_wilsons_with_running(..., mu_had=2.0, m12_np_budget=3.281e-15 GeV)`.

Budget: `0.045`, built as `min(max||q/p|-1|_95=0.045, max|sin(phi_D)|_95=0.063313)`. Sources: `C002.yaml` HFLAV anchors, `C001.yaml` `Delta_m_D`, and `docs/quark_scan_assumptions_compact.tex:464-470`.

Numerical cross-check: reference complex point gives `predicted=0.0234065082`, `ratio=0.5201446261`; excluded complex point gives `predicted=9.36260327`, `ratio=208.0578504`. Direct running `abs_m12` matches the complex helper.

Tests: `tests/constraints/primary/charm/test_C002.py`: `7 passed`; full `python -m pytest tests/constraints/ -q`: `89 passed`.

Scoped status: `M flavor_catalog_constraints/physics_adapters/deltaf2.py`, `?? flavor_catalog_constraints/primary/charm/C002.py`, `?? tests/constraints/primary/charm/test_C002.py`. Full worktree also has unrelated dirty B002/K005/rare-kaon files from other work.

Reviewer focus: the `NEEDS-HUMAN-PHYSICS` caveat is intentional: no full HFLAV `|q/p|-phi_D` contour, `Gamma12`, or grounded SM long-distance charm phase is available.