# Implement B016 — BR(B→Kℓℓ) (exclusive b→sℓℓ). family=beauty. EXTEND the shared rare_b_dilepton module.
EXTEND quarkConstraints/rare_b_dilepton.py (built for B005) with exclusive B→K form-factor machinery (q² integration / C9,C10 with form factors), APPEND-ONLY (do not modify B005's leptonic functions). Wrap via the existing rare_b_meson adapter (append-only). Design so B015 (inclusive B→X_sℓℓ), B017 (B→K(*)ℓℓ), B018/B019 (R_K/R_K*), B021 (Λ_b→Λℓℓ) can reuse the exclusive/inclusive machinery. Validate SM BR(B→Kμμ) in a reference q² bin from B016.yaml. RS C9/C10 proxy → NEEDS-HUMAN-PHYSICS. Budget from B016.yaml, HARD.

## You are agent2 — PHYSICS fact-checker (codex). Review ONLY the physics of the constraint named above. Do NOT rewrite code; produce a precise numbered findings list. Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing.

Files: flavor_catalog_constraints/primary/<family>/<ID>.py, its additive wrapper in flavor_catalog_constraints/physics_adapters/deltaf2.py, tests/constraints/primary/<family>/test_<ID>.py. Source of truth = the catalog yaml flavor_catalog/processes/<family>/<ID>.yaml + its reference snapshots; physics core = quarkConstraints/deltaf2.py; audit docs in docs/audits/.

MANDATORY physics checks (verdict each):
1. Correct amplitude: mass-difference observables (Δm) use |M_12^NP| (magnitude); CP observables use the IMAGINARY part. Confirm the code uses the right one (not the wrong part).
2. **QCD running ACTUALLY used**: the constraint must call the *_with_running evaluator (Wilsons evolved to mu_had=2 GeV). Confirm it does NOT use the non-running path for its verdict. Quantify the running effect if you can.
3. **Budget defensible & uncertainty-aware**: NOT a bare central residual. For long-distance-dominated Δm (kaon, charm) the accepted convention is |M_12^NP| <= Δm^exp/2; for Bd/Bs an SM-vs-exp room with uncertainties. Check the value against the yaml + core + any audit doc. State whether the budget is right, too tight, or too loose, with numbers.
4. Anchor numbers loaded match PDG/HFLAV/the yaml snapshots (spot-check actual values).
5. Severity appropriate; units consistent (GeV); any wrong/misleading physics in docstring/notes/diagnostics.

OUTPUT (stdout, <=22 lines): numbered findings, each tagged BLOCKER / SHOULD-FIX / NIT with the precise issue + correct physics + file:line + numbers. End with: PHYSICS-OK or PHYSICS-NEEDS-FIXES.
