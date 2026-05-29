# Implement B005 — BR(B_s→μ⁺μ⁻). family=beauty. PIONEER: build the shared b→sℓℓ module.
No existing b→sℓℓ machinery. BUILD a new module quarkConstraints/rare_b_dilepton.py (SM: Buras b→sℓℓ effective Hamiltonian, C10 dominant for B_s→μμ; SM BR(B_s→μμ)~3.66e-9 — validate this) wrapped by a NEW adapter physics_adapters/rare_b_meson.py. Design the module so siblings (B006 B0→μμ b→dℓℓ, B015/B016/B017 B→X_s/Kℓℓ, B021) can REUSE it. RS NP = documented Z/KK-penguin proxy for C10/C9, flagged NEEDS-HUMAN-PHYSICS (b→s EW penguin RS couplings not on ParameterPoint). Budget from B005.yaml (LHCb/CMS measured vs SM), HARD.

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
