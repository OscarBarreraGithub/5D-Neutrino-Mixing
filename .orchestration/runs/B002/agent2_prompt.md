# Implement B002 — S_ψKs / sin(2β), CP phase in B_d→J/ψK_S. family=beauty.
RIGOROUS (NOT a proxy): NP enters as the PHASE of the B_d ΔF=2 mixing amplitude. REUSE the EXISTING adapter `bd_mixing_from_wilsons_with_running` (already in physics_adapters/deltaf2.py; returns complex M12^NP) — do NOT build new ΔF=2 machinery and do NOT modify deltaf2.py. Compute NP phase shift φ_d^NP = arg(1 + M12^NP / M12^SM) and S_ψKs = sin(2β + φ_d^NP); compare to the measured value with an uncertainty-aware budget. exp (sin2β≈0.69) + SM from B002.yaml. This is a CP observable: use the IMAGINARY part / phase of M12^NP, not |M12|. If M12^SM (or 2β) is not obtainable from the core/CKM, flag NEEDS-HUMAN-PHYSICS for that piece and use the documented best value from the yaml.

## You are agent2 — PHYSICS fact-checker (codex). Review ONLY the physics of the constraint named above. Do NOT rewrite code; produce a precise numbered findings list. Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing.

Files: flavor_catalog_constraints/primary/<family>/<ID>.py, its additive wrapper in flavor_catalog_constraints/physics_adapters/deltaf2.py, tests/constraints/primary/<family>/test_<ID>.py. Source of truth = the catalog yaml flavor_catalog/processes/<family>/<ID>.yaml + its reference snapshots; physics core = quarkConstraints/deltaf2.py; audit docs in docs/audits/.

MANDATORY physics checks (verdict each):
1. Correct amplitude: mass-difference observables (Δm) use |M_12^NP| (magnitude); CP observables use the IMAGINARY part. Confirm the code uses the right one (not the wrong part).
2. **QCD running ACTUALLY used**: the constraint must call the *_with_running evaluator (Wilsons evolved to mu_had=2 GeV). Confirm it does NOT use the non-running path for its verdict. Quantify the running effect if you can.
3. **Budget defensible & uncertainty-aware**: NOT a bare central residual. For long-distance-dominated Δm (kaon, charm) the accepted convention is |M_12^NP| <= Δm^exp/2; for Bd/Bs an SM-vs-exp room with uncertainties. Check the value against the yaml + core + any audit doc. State whether the budget is right, too tight, or too loose, with numbers.
4. Anchor numbers loaded match PDG/HFLAV/the yaml snapshots (spot-check actual values).
5. Severity appropriate; units consistent (GeV); any wrong/misleading physics in docstring/notes/diagnostics.

OUTPUT (stdout, <=22 lines): numbered findings, each tagged BLOCKER / SHOULD-FIX / NIT with the precise issue + correct physics + file:line + numbers. End with: PHYSICS-OK or PHYSICS-NEEDS-FIXES.
