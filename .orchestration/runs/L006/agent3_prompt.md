# Implement L006 — muonium-antimuonium (M–M̄) conversion. family=charged_lepton (PRIMARY). small new module / proxy.
Muonium-antimuonium conversion is a ΔL=2 process (lepton-number-violating four-fermion μ⁺e⁻↔μ⁻e⁺). SM≈0 → pure-NP bound vs the PSI MACS limit on the conversion probability G_MM̄ (from L006.yaml), HARD. The RS ΔL=2 effective operator is NOT on ParameterPoint → documented proxy NEEDS-HUMAN-PHYSICS. Build a small adapter; compare the predicted G_MM̄ (proxy) to the MACS limit. Document the proxy.

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
