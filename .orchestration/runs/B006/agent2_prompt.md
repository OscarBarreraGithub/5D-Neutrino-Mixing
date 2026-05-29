# Implement B006 — BR(B⁰→μ⁺μ⁻). family=beauty. REUSE the shared b→sℓℓ module.
REUSE quarkConstraints/rare_b_dilepton.py + physics_adapters/rare_b_meson.py (built for B005). B0→μμ is the b→DOWN-quark analog: use the b→d CKM factor (λ_t^d = V_tb V_td*) and the b→dℓℓ entry point the module already exposes (the module was designed with bq/bs/bd entrypoints). SM BR(B0→μμ)~1.0e-10 — validate. Same A_ΔΓ time-integration treatment as B005 (B_d has y_d≈0). RS C9/C10 proxy reused → NEEDS-HUMAN-PHYSICS. Budget from B006.yaml (LHCb/CMS vs SM), HARD. If a small b→d helper must be added to the module, keep it APPEND-ONLY and do not modify B005's b→s functions.

## You are agent2 — PHYSICS fact-checker (codex). Review ONLY the physics of the constraint named above. Do NOT rewrite code; produce a precise numbered findings list. Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing.

Files: flavor_catalog_constraints/primary/<family>/<ID>.py, its additive wrapper in flavor_catalog_constraints/physics_adapters/deltaf2.py, tests/constraints/primary/<family>/test_<ID>.py. Source of truth = the catalog yaml flavor_catalog/processes/<family>/<ID>.yaml + its reference snapshots; physics core = quarkConstraints/deltaf2.py; audit docs in docs/audits/.

MANDATORY physics checks (verdict each):
1. Correct amplitude: mass-difference observables (Δm) use |M_12^NP| (magnitude); CP observables use the IMAGINARY part. Confirm the code uses the right one (not the wrong part).
2. **QCD running ACTUALLY used**: the constraint must call the *_with_running evaluator (Wilsons evolved to mu_had=2 GeV). Confirm it does NOT use the non-running path for its verdict. Quantify the running effect if you can.
3. **Budget defensible & uncertainty-aware**: NOT a bare central residual. For long-distance-dominated Δm (kaon, charm) the accepted convention is |M_12^NP| <= Δm^exp/2; for Bd/Bs an SM-vs-exp room with uncertainties. Check the value against the yaml + core + any audit doc. State whether the budget is right, too tight, or too loose, with numbers.
4. Anchor numbers loaded match PDG/HFLAV/the yaml snapshots (spot-check actual values).
5. Severity appropriate; units consistent (GeV); any wrong/misleading physics in docstring/notes/diagnostics.

OUTPUT (stdout, <=22 lines): numbered findings, each tagged BLOCKER / SHOULD-FIX / NIT with the precise issue + correct physics + file:line + numbers. End with: PHYSICS-OK or PHYSICS-NEEDS-FIXES.
