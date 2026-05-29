# Implement C005 — BR(D⁰→e⁺e⁻). family=charm. REUSE rare_charm_dilepton.
REUSE quarkConstraints/rare_charm_dilepton.py (built for C004) + adapter. D0→ee is the electron mode (strong helicity suppression ∝ m_e²); reuse the c→uℓℓ machinery with the electron mass. SM SD≈0, LD-dominated; HARD bound vs PDG limit (from C005.yaml). RS proxy reused → NEEDS-HUMAN-PHYSICS. Append-only if any addition needed (do not modify C004's muon functions).

## You are agent3 — CODE reviewer + numerical verifier (codex). Review ONLY the code/tests of the constraint named above. Do NOT rewrite code; produce a precise numbered findings list and run things yourself. Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing.

Files: flavor_catalog_constraints/primary/<family>/<ID>.py, additive wrapper in physics_adapters/deltaf2.py, tests/constraints/primary/<family>/test_<ID>.py. Scaffold contract in flavor_catalog_constraints/base.py (numeric ConstraintResult fields MUST be real floats — __post_init__ rejects complex; complex amplitudes go in diagnostics).

Checks (do them, with evidence):
1. ISOLATION (judge by what THIS constraint ADDS, not global worktree cleanliness): imports no other constraint; adds ONLY its own constraint file + test + its intended append-only physics additions (a new adapter wrapper, and for machinery-building constraints a NEW/designated physics module). NOTE: this runs in a SHARED worktree during parallel waves, so `git status`/`git diff` WILL show unrelated wave-mate files and a non-empty `quarkConstraints/` diff — that is EXPECTED, NOT a violation. Verify instead that (a) this constraint did not MODIFY/remove any pre-existing function (function-body/AST check on touched files), and (b) it reaches physics only via an adapter. Do NOT flag "worktree dirty" / "git diff quarkConstraints non-empty" as a blocker.
2. CONTRACT: result numeric fields real floats; complex in diagnostics; severity set; missing-couplings path returns a non-crashing result.
3. ANCHOR: values come from the yaml via load_anchor (NOT hardcoded); mismatched/missing anchor fails loudly (probe).
4. TEST ADEQUACY + NUMERICAL VERIFICATION: tests must include a real numerical cross-check that is INDEPENDENT of the adapter the constraint calls — recompute the expected value directly from the underlying core in quarkConstraints/ (e.g. _evolve_wilsons + compute_m12_np, or the relevant core evaluator), NOT by re-calling the same adapter wrapper the constraint uses (that would be circular). Report the ACTUAL numbers and whether they agree. A safe point must pass and a clearly-excluded point must fail.
5. DETERMINISM: evaluate() is pure (same input->same output, no mutation).
6. Run `python -m pytest tests/constraints/ -q`; report actual counts. Registry import smoke.

OUTPUT (stdout, <=22 lines): numbered findings tagged BLOCKER / SHOULD-FIX / NIT with file:line + fix, the ACTUAL cross-check numbers, pytest counts. End with: CODE-OK or CODE-NEEDS-FIXES.
