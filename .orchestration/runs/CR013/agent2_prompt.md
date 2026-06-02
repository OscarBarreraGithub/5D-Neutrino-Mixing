# Implement CR013 — diphoton high-mass resonance (spin-0/spin-2): pp → X → γγ. family=collider_rs (PRIMARY). REUSE flavor_catalog_constraints/physics_adapters/kk_graviton_resonance.py (the RS spin-2 KK-graviton resonance machinery; first spin-2 Bessel root ≈ 3.8317 × Λ_IR, as built for CR007); STUDY committed primary/collider_rs/CR007.py and CR010.py for the house style. Compare the lightest RS KK graviton mass (spin-2, from the graviton resonance core) to the CMS 2024 RS-graviton diphoton observed lower limit from CR013.yaml (value 4.8 TeV at k̃=0.1, 95% CL; use the strongest applicable γγ resonance-mass limit row in the yaml). σ×BR / branching-surface / coupling (k̃=k/M_Pl) / width / acceptance recast is a documented proxy → NEEDS-HUMAN-PHYSICS (flag in docstring + diagnostics). HARD. Append-only: add CR013-specific helpers to the adapter; do NOT modify pre-existing CR001–CR012 functions or the kk_graviton_resonance core (used by CR007).

## You are agent2 — PHYSICS fact-checker (codex). Review ONLY the physics of the constraint named above. Do NOT rewrite code; produce a precise numbered findings list. Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing.

Files: flavor_catalog_constraints/primary/<family>/<ID>.py, its additive wrapper in flavor_catalog_constraints/physics_adapters/deltaf2.py, tests/constraints/primary/<family>/test_<ID>.py. Source of truth = the catalog yaml flavor_catalog/processes/<family>/<ID>.yaml + its reference snapshots; physics core = quarkConstraints/deltaf2.py; audit docs in docs/audits/.

MANDATORY physics checks (verdict each):
1. Correct amplitude: mass-difference observables (Δm) use |M_12^NP| (magnitude); CP observables use the IMAGINARY part. Confirm the code uses the right one (not the wrong part).
2. **QCD running ACTUALLY used**: the constraint must call the *_with_running evaluator (Wilsons evolved to mu_had=2 GeV). Confirm it does NOT use the non-running path for its verdict. Quantify the running effect if you can.
3. **Budget defensible & uncertainty-aware**: NOT a bare central residual. For long-distance-dominated Δm (kaon, charm) the accepted convention is |M_12^NP| <= Δm^exp/2; for Bd/Bs an SM-vs-exp room with uncertainties. Check the value against the yaml + core + any audit doc. State whether the budget is right, too tight, or too loose, with numbers.
4. Anchor numbers loaded match PDG/HFLAV/the yaml snapshots (spot-check actual values).
5. Severity appropriate; units consistent (GeV); any wrong/misleading physics in docstring/notes/diagnostics.

OUTPUT (stdout, <=22 lines): numbered findings, each tagged BLOCKER / SHOULD-FIX / NIT with the precise issue + correct physics + file:line + numbers. End with: PHYSICS-OK or PHYSICS-NEEDS-FIXES.
