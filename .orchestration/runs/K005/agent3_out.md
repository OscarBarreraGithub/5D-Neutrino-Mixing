1. BLOCKER `quarkConstraints/rare_kaon_snd.py:58`: isolation check fails; `git diff -- quarkConstraints/` is nonempty and adds the K_L machinery. Fix by moving this out of `quarkConstraints/` per the reviewer checklist, or clarify that the checklist allows the requested core extension.

2. SHOULD-FIX `flavor_catalog_constraints/primary/kaon/K005.py:41`: anchors are loaded via `load_full_yaml` plus custom list scanning, not `load_anchor` as requested. Values do come from YAML and my missing/mismatch probes raised `AnchorError`; fix by using/extending the scaffold anchor loader for this list-shaped sidecar.

3. SHOULD-FIX `flavor_catalog_constraints/physics_adapters/deltaf2.py:20`: unrelated deltaf2 adapter additions are in the reviewed worktree, while K005 uses `physics_adapters/rare_kaon.py`. Fix by splitting unrelated adapter changes from K005 or documenting why they are part of this constraint.

Numerical cross-check: SM predicted/direct `2.95375989343059029e-11`, ratio `1.34261813337754107e-02`, budget `2.19999999999999985e-09`, pass `True`.

Safe point: predicted/direct `2.97802188744418121e-11`, ratio `1.35364631247462796e-02`, pass `True`. Excluded point: predicted/direct `5.20688189145627510e-08`, ratio `2.36676449611648891e+01`, pass `False`.

Contract/determinism probes passed: numeric fields are real floats, complex values stay in diagnostics, missing couplings do not crash, repeated evaluate gives identical result and no coupling mutation.

`python -m pytest tests/constraints/ -q`: `89 passed in 2.73s`. Registry smoke: `K005 kaon HARD BR(K_L -> pi0 nu nubar)`.

CODE-NEEDS-FIXES