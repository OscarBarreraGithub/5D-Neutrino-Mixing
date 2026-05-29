# Implement K005 — BR(K_L→π0νν̄), the purely-CP-violating sibling of K004. family=kaon.
REUSE/EXTEND the K004 machinery: quarkConstraints/rare_kaon_snd.py and physics_adapters/rare_kaon.py (APPEND-ONLY — add a K_L→π0νν function; do NOT alter K004's K+ functions). K_L→π0νν is purely CP-violating (top/imaginary part only): BR_SM = kappa_L * (Im(λ_t)/λ^5 · X_t)^2. Reuse X_t and take kappa_L + exp limit (KOTO) + SM from K005.yaml. RS NP: SAME documented Z-like proxy approach as K004 (the imaginary part); flag NEEDS-HUMAN-PHYSICS for the full EW KK/Z/Z′+ν matching exactly as K004 did. Budget from K005.yaml. HARD. SM-limit validation: reproduce BR_SM(K_L)~3e-11.

## You are agent2 — PHYSICS fact-checker (codex). Review ONLY the physics of the constraint named above. Do NOT rewrite code; produce a precise numbered findings list. Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing.

Files: flavor_catalog_constraints/primary/<family>/<ID>.py, its additive wrapper in flavor_catalog_constraints/physics_adapters/deltaf2.py, tests/constraints/primary/<family>/test_<ID>.py. Source of truth = the catalog yaml flavor_catalog/processes/<family>/<ID>.yaml + its reference snapshots; physics core = quarkConstraints/deltaf2.py; audit docs in docs/audits/.

MANDATORY physics checks (verdict each):
1. Correct amplitude: mass-difference observables (Δm) use |M_12^NP| (magnitude); CP observables use the IMAGINARY part. Confirm the code uses the right one (not the wrong part).
2. **QCD running ACTUALLY used**: the constraint must call the *_with_running evaluator (Wilsons evolved to mu_had=2 GeV). Confirm it does NOT use the non-running path for its verdict. Quantify the running effect if you can.
3. **Budget defensible & uncertainty-aware**: NOT a bare central residual. For long-distance-dominated Δm (kaon, charm) the accepted convention is |M_12^NP| <= Δm^exp/2; for Bd/Bs an SM-vs-exp room with uncertainties. Check the value against the yaml + core + any audit doc. State whether the budget is right, too tight, or too loose, with numbers.
4. Anchor numbers loaded match PDG/HFLAV/the yaml snapshots (spot-check actual values).
5. Severity appropriate; units consistent (GeV); any wrong/misleading physics in docstring/notes/diagnostics.

OUTPUT (stdout, <=22 lines): numbered findings, each tagged BLOCKER / SHOULD-FIX / NIT with the precise issue + correct physics + file:line + numbers. End with: PHYSICS-OK or PHYSICS-NEEDS-FIXES.
