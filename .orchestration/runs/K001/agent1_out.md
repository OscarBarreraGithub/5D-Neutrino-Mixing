**Plan**
Used the audited Delta F=2 adapter path, loaded K001 exp/SM/FLAG anchors from YAML, applied a HARD central-value NP budget, returned only real numeric result fields, and covered registration, anchors, pass/fail physics, and determinism.

**Files Changed**
- [K001.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/kaon/K001.py:1): production constraint.
- [deltaf2.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/physics_adapters/deltaf2.py:62): additive adapter wrappers.
- [test_K001.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/kaon/test_K001.py:1): production tests.
- Deleted stale `tests/constraints/test_K001_example.py`.

**Formula/Convention**
`predicted = |epsilon_K^NP|` from audited `evaluate_epsilon_k`, with sign-sensitive `im_m12_np_gev` kept in diagnostics. Budget is `|epsilon_exp - epsilon_SM| = |0.002228 - 0.002161| = 6.7e-5`, loaded from `K001.yaml`; `ratio = predicted / budget`; severity is `HARD`.

**Pytest**
- `python -m pytest tests/constraints/ -q`: `20 passed in 3.63s`
- `python -m pytest tests/test_epsilon_k_physics.py -q`: `20 passed in 2.27s`
- `git diff --check`: clean

**Reviewer Notes**
The implementation deliberately uses the audited non-running `evaluate_epsilon_k` path, matching the existing adapter convention and requested numerical cross-check. Missing couplings return a non-vetoing, unevaluated result with diagnostics. Unrelated untracked `.orchestration/runs/` was left untouched.