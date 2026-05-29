# Implement B004 — φ_s (CP phase in B_s→J/ψφ). family=beauty.
REUSE the existing running B_s ΔF=2 complex M12 (physics_adapters/deltaf2.py bs-mixing complex helper, analogous to B002's bd helper; add an append-only bs phase helper if needed). NP phase shift φ_s^NP = arg(1 + M12^NP/M12^SM); compare to measured φ_s with an uncertainty-aware budget. exp + SM from B004.yaml. CP observable → use Im/phase of M12^NP. If the SM φ_s reference phase isn't in the core, take it from the yaml and flag NEEDS-HUMAN-PHYSICS (like B002's 2β). Do NOT modify deltaf2 core.

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
