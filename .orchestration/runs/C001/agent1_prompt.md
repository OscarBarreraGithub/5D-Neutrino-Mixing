# Implement constraint C001 — Neutral charm mixing (x_D, y_D, Delta m_D), family=charm.
Observable: D0-D0bar mixing. RUNNING evaluator: quarkConstraints/deltaf2.py:evaluate_d0_mixing_with_running. Catalog: flavor_catalog/processes/charm/C001.yaml. Charm mixing is long-distance dominated and theoretically dirty; SM short-distance is small. Use the core's documented D0 convention (likely full measured mixing as NP room, analogous to Delta m_K). Document the budget choice clearly for reviewers; flag if the physics is uncertain.

## Shared pipeline rules (this is agent1; codex implementer)
You implement ONE constraint into the plugin scaffold at flavor_catalog_constraints/. Two codex reviewers (physics + code) and an Opus reviewer will check you afterward, so be correct and defensible. First write a SHORT plan, then implement. Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing.

STUDY FIRST (read, don't modify): the committed worked example flavor_catalog_constraints/primary/kaon/K001.py and its test tests/constraints/primary/kaon/test_K001.py (this is the GOLD STANDARD — mirror its structure), the scaffold (base.py, anchors.py, registry.py, point_builder.py, TEMPLATE.py, physics_adapters/deltaf2.py), and your constraint's catalog sidecar yaml (the source of truth for experimental/SM values + provenance).

MANDATORY CONVENTIONS (learned from K001 review — non-negotiable):
1. **Use the QCD-RUNNING evaluator** for your observable (the *_with_running variant) that evolves Wilsons to mu_had=2 GeV before 2-GeV matrix elements. The non-running path is WRONG. Add an additive (append-only) wrapper to physics_adapters/deltaf2.py if needed; never remove existing wrappers.
2. **Budget must be uncertainty-aware**, NOT a bare central residual. Check docs/audits/ for a documented band/decision for your observable and the convention already used in quarkConstraints/deltaf2.py for your evaluator. For long-distance-dominated observables (Δm_K, D0 mixing) the established core convention is |M_12^NP| <= Δm^exp/2 (full measured splitting as NP room) — use it if that's what the core does. For Bd/Bs use the SM-vs-exp room with uncertainties. Cite your source. If unsure which band, pick the physically defensible one and document it for the reviewers.
3. Numeric ConstraintResult fields MUST be real floats; complex amplitudes go in diagnostics. Severity: HARD for an observed bound where NP must fit the measurement (typical for these mixing observables). Missing-couplings -> non-crashing result.
4. Anchor loaded from the yaml via the scaffold loader (NO hardcoded numbers, fail loud on mismatch).
5. Keep ALL changes ISOLATED to your constraint file, its test, and (if needed) an additive adapter. Do NOT touch core physics modules or other constraints.

TESTS: write tests/constraints/primary/<family>/test_<ID>.py: registration, anchor matches yaml, evaluate runs end-to-end on a real QuarkMassBasisCouplings with a NUMERICAL cross-check against a DIRECT call to the running evaluator (predicted/ratio/budget real+finite), a safe point passes and a clearly-excluded large-NP point fails, and evaluate is pure. Run `python -m pytest tests/constraints/ -q` until green.

OUTPUT (stdout, concise, <=18 lines): your plan, files changed, the running evaluator used, the budget value + construction + source citation, the numerical cross-check numbers, pytest counts, git diff --stat confirming isolation, and anything reviewers should scrutinize. Do not paste full files.
