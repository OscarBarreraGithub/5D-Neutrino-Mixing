You are agent1 again. Two independent reviewers checked your K001 (epsilon_K) implementation and found issues, including 2 BLOCKERS. Fix them. Read their full findings first:
- Physics review: .orchestration/runs/K001/agent2_out.md  (verdict PHYSICS-NEEDS-FIXES)
- Code review:    .orchestration/runs/K001/agent3_out.md  (verdict CODE-OK, minor)

Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing (you are -C'd in)

## Required fixes (BLOCKERS — must fix):
1. **QCD running.** K001 currently wraps the NON-running evaluate_epsilon_k, applying Wilsons matched at M_KK directly onto 2-GeV matrix elements. This is wrong. Switch to the RG-evolved path: quarkConstraints/deltaf2.py:evaluate_epsilon_k_with_running (evolve Wilsons to mu_had=2 GeV before the 2-GeV matrix elements). Add an additive (append-only) wrapper in flavor_catalog_constraints/physics_adapters/deltaf2.py for the running path; do NOT remove the existing wrapper. Verify the change: on the test couplings the ratio should move from ~0.266 (unrun) to ~2.845 (run).
2. **Uncertainty-aware budget.** The HARD budget must NOT be the bare central residual 6.7e-5. READ docs/audits/epsilon_k_sm_decision.md (esp. lines ~34-37, ~96-100) and implement the repo's DOCUMENTED uncertainty-aware budget band for epsilon_K (the loose edge is ~3e-4). Use that policy for the HARD veto budget. Keep severity HARD. Cite the doc in the code/docstring. If the doc defines a specific band construction (e.g. central residual plus N-sigma of SM+exp+kappa uncertainties), implement THAT, do not invent your own.

## SHOULD-FIX:
3. Propagate kappa_epsilon uncertainty (BGS: 0.94(2)) and the BGS grouped SM uncertainties into the budget band / diagnostics, consistent with the audit doc. If the doc already folds these into the band, just make sure you use the doc's prescription.
4. Fix misleading wording in the K001 docstring/notes: state the budget is an uncertainty-aware band (not "central-value room").

## Code-review item:
5. Ensure the new test file tests/constraints/primary/kaon/test_K001.py is part of your change (it will be committed). Keep changes ISOLATED to K001.py, its test, and the additive adapter — do not touch core physics modules or other constraints.

## After fixing:
- Update the tests to reflect the running path + new budget (the safe/excluded thresholds change; cross-check `predicted` and `ratio` against a DIRECT call to evaluate_epsilon_k_with_running with the same couplings; assert numeric fields real+finite; one safe point passes, one clearly-excluded point fails).
- Run `python -m pytest tests/constraints/ -q` and `python -m pytest tests/test_epsilon_k_physics.py -q`; fix until green.

## Output (stdout, concise, <=20 lines):
What you changed for each numbered item, the budget value + band construction you used (with the doc citation), the new ratio numbers on the test couplings (unrun vs run), the pytest counts, and confirm changes are isolated (git diff --stat). Do not paste full files.
