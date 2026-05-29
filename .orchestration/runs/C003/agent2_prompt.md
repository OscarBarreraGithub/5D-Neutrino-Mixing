# Implement C003 — ΔA_CP(D⁰→K⁺K⁻, π⁺π⁻) direct CP violation. family=charm. HARD STUB (dual NEEDS-HUMAN).
Charm DIRECT CP violation is long-distance penguin / non-perturbative dominated — the SM prediction is NOT reliably computable from first principles AND the RS NP needs penguin matching not on ParameterPoint. Implement an honest STUB: load the measured ΔA_CP anchor from C003.yaml (LHCb), record the observable, and return a NON-VETOING result that compares the documented NP-room to the measurement. Flag NEEDS-HUMAN-PHYSICS on BOTH the SM side (long-distance penguin) AND the NP side (RS penguin matching). Severity INFO or SOFT (do NOT hard-veto on an observable whose SM is incalculable). Document thoroughly. Wrap via a small adapter; do not fake a penguin calculation.

## You are agent2 — PHYSICS fact-checker (codex). Review ONLY the physics of the constraint named above. Do NOT rewrite code; produce a precise numbered findings list. Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing.

Files: flavor_catalog_constraints/primary/<family>/<ID>.py, its additive wrapper in flavor_catalog_constraints/physics_adapters/deltaf2.py, tests/constraints/primary/<family>/test_<ID>.py. Source of truth = the catalog yaml flavor_catalog/processes/<family>/<ID>.yaml + its reference snapshots; physics core = quarkConstraints/deltaf2.py; audit docs in docs/audits/.

MANDATORY physics checks (verdict each):
1. Correct amplitude: mass-difference observables (Δm) use |M_12^NP| (magnitude); CP observables use the IMAGINARY part. Confirm the code uses the right one (not the wrong part).
2. **QCD running ACTUALLY used**: the constraint must call the *_with_running evaluator (Wilsons evolved to mu_had=2 GeV). Confirm it does NOT use the non-running path for its verdict. Quantify the running effect if you can.
3. **Budget defensible & uncertainty-aware**: NOT a bare central residual. For long-distance-dominated Δm (kaon, charm) the accepted convention is |M_12^NP| <= Δm^exp/2; for Bd/Bs an SM-vs-exp room with uncertainties. Check the value against the yaml + core + any audit doc. State whether the budget is right, too tight, or too loose, with numbers.
4. Anchor numbers loaded match PDG/HFLAV/the yaml snapshots (spot-check actual values).
5. Severity appropriate; units consistent (GeV); any wrong/misleading physics in docstring/notes/diagnostics.

OUTPUT (stdout, <=22 lines): numbered findings, each tagged BLOCKER / SHOULD-FIX / NIT with the precise issue + correct physics + file:line + numbers. End with: PHYSICS-OK or PHYSICS-NEEDS-FIXES.
