Plan executed: mirror K001/B001, reuse running Bd Delta F=2 path, compute `phi_d^NP = arg(1 + M12_NP/M12_SM)`, compare total `S_psiK_S` to HFLAV.

Files changed: `flavor_catalog_constraints/primary/beauty/B002.py`, `tests/constraints/primary/beauty/test_B002.py`, and append-only adapter support in `flavor_catalog_constraints/physics_adapters/deltaf2.py`.

Running evaluator: QCD-running complex `M12_NP` wrapper at `mu_had=2 GeV`, cross-checked against existing `bd_mixing_from_wilsons_with_running` magnitude.

Budget: `0.01769290183565696 = sqrt(0.011^2 + 0.011056690249896537^2 + 0.00835394374436627^2)`, from `B002.yaml` HFLAV exp, beta solution, and penguin phase bound.

`M12_SM = DELTA_M_BD_SM / 2 = 1.8e-13 GeV`; `2beta` is flagged `NEEDS-HUMAN-PHYSICS` because no core CKM phase exists, using YAML beta `22.63 deg`.

Numerical cross-check point: `left_db=(1+1j)e-3`, `right_db=0`, `M_KK=3 TeV` gives `M12_NP=2.9852989857064115e-15j`, `phi_d^NP=0.950163068823887 deg`, `S=0.7218829884163276`, `ratio=0.6716246168494265`.

Tests: `python -m pytest tests/constraints/ -q` -> `89 passed`.

Scoped diff/stat for this task: adapter diff `74 insertions`; `B002.py` and `test_B002.py` are new. Global worktree also contains unrelated concurrent C002/K005/rare-kaon changes that I did not touch.