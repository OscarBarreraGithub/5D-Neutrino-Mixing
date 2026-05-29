# Implement B002 — S_ψKs / sin(2β), CP phase in B_d→J/ψK_S. family=beauty.
RIGOROUS (NOT a proxy): NP enters as the PHASE of the B_d ΔF=2 mixing amplitude. REUSE the EXISTING adapter `bd_mixing_from_wilsons_with_running` (already in physics_adapters/deltaf2.py; returns complex M12^NP) — do NOT build new ΔF=2 machinery and do NOT modify deltaf2.py. Compute NP phase shift φ_d^NP = arg(1 + M12^NP / M12^SM) and S_ψKs = sin(2β + φ_d^NP); compare to the measured value with an uncertainty-aware budget. exp (sin2β≈0.69) + SM from B002.yaml. This is a CP observable: use the IMAGINARY part / phase of M12^NP, not |M12|. If M12^SM (or 2β) is not obtainable from the core/CKM, flag NEEDS-HUMAN-PHYSICS for that piece and use the documented best value from the yaml.

## You are agent3 — CODE reviewer + numerical verifier (codex). Review ONLY the code/tests of the constraint named above. Do NOT rewrite code; produce a precise numbered findings list and run things yourself. Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing.

Files: flavor_catalog_constraints/primary/<family>/<ID>.py, additive wrapper in physics_adapters/deltaf2.py, tests/constraints/primary/<family>/test_<ID>.py. Scaffold contract in flavor_catalog_constraints/base.py (numeric ConstraintResult fields MUST be real floats — __post_init__ rejects complex; complex amplitudes go in diagnostics).

Checks (do them, with evidence):
1. ISOLATION: constraint imports no other constraint; changes limited to its file + test + ADDITIVE adapter wrapper; `git diff quarkConstraints/` is empty; the adapter wrapper is append-only and doesn't alter existing functions.
2. CONTRACT: result numeric fields real floats; complex in diagnostics; severity set; missing-couplings path returns a non-crashing result.
3. ANCHOR: values come from the yaml via load_anchor (NOT hardcoded); mismatched/missing anchor fails loudly (probe).
4. TEST ADEQUACY + NUMERICAL VERIFICATION: tests must include a real numerical cross-check — independently construct a QuarkMassBasisCouplings, run the constraint, and compare `predicted`/`ratio`/`budget` to a DIRECT call to the *_with_running evaluator in quarkConstraints/deltaf2.py. Report the ACTUAL numbers and whether they agree. A safe point must pass and a clearly-excluded point must fail.
5. DETERMINISM: evaluate() is pure (same input->same output, no mutation).
6. Run `python -m pytest tests/constraints/ -q`; report actual counts. Registry import smoke.

OUTPUT (stdout, <=22 lines): numbered findings tagged BLOCKER / SHOULD-FIX / NIT with file:line + fix, the ACTUAL cross-check numbers, pytest counts. End with: CODE-OK or CODE-NEEDS-FIXES.
