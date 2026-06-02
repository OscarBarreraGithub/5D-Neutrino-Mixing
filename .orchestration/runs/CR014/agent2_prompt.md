# Implement CR014 — four-top-quark production (pp → t t̄ t t̄); BSM top-philic vector mediator search. family=collider_rs (PRIMARY). REUSE flavor_catalog_constraints/physics_adapters/collider_resonance.py (mass-vs-limit recast); STUDY committed primary/collider_rs/CR010.py and CR012.py for the house style. Compare the relevant RS top-philic vector resonance mass (M_KK proxy) to the CMS four-top two-lepton exclusion from CR014.yaml: a top-philic vector mediator (Z') is excluded up to m(Z')=850 GeV (observed) at Γ/m=50%, 95% CL (use the strongest applicable mass-exclusion row in the yaml). NOTE this is a WIDTH-DEPENDENT exclusion (Γ/m=50% benchmark) and a relatively low mass — the σ×BR / width-dependence / four-top acceptance / SM-four-top-background recast is a documented proxy → NEEDS-HUMAN-PHYSICS (flag in docstring + diagnostics). HARD. Append-only: add CR014-specific helpers to the adapter; do NOT modify pre-existing CR001–CR013 functions or the collider_resonance core.

## You are agent2 — PHYSICS fact-checker (codex). Review ONLY the physics of the constraint named above. Do NOT rewrite code; produce a precise numbered findings list. Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing.

Files: flavor_catalog_constraints/primary/<family>/<ID>.py, its additive wrapper in flavor_catalog_constraints/physics_adapters/deltaf2.py, tests/constraints/primary/<family>/test_<ID>.py. Source of truth = the catalog yaml flavor_catalog/processes/<family>/<ID>.yaml + its reference snapshots; physics core = quarkConstraints/deltaf2.py; audit docs in docs/audits/.

MANDATORY physics checks (verdict each):
1. Correct amplitude: mass-difference observables (Δm) use |M_12^NP| (magnitude); CP observables use the IMAGINARY part. Confirm the code uses the right one (not the wrong part).
2. **QCD running ACTUALLY used**: the constraint must call the *_with_running evaluator (Wilsons evolved to mu_had=2 GeV). Confirm it does NOT use the non-running path for its verdict. Quantify the running effect if you can.
3. **Budget defensible & uncertainty-aware**: NOT a bare central residual. For long-distance-dominated Δm (kaon, charm) the accepted convention is |M_12^NP| <= Δm^exp/2; for Bd/Bs an SM-vs-exp room with uncertainties. Check the value against the yaml + core + any audit doc. State whether the budget is right, too tight, or too loose, with numbers.
4. Anchor numbers loaded match PDG/HFLAV/the yaml snapshots (spot-check actual values).
5. Severity appropriate; units consistent (GeV); any wrong/misleading physics in docstring/notes/diagnostics.

OUTPUT (stdout, <=22 lines): numbered findings, each tagged BLOCKER / SHOULD-FIX / NIT with the precise issue + correct physics + file:line + numbers. End with: PHYSICS-OK or PHYSICS-NEEDS-FIXES.
