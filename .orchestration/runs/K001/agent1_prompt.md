You are agent1 (implementer) in a strict per-constraint pipeline. Implement the flavor-physics constraint K001 (epsilon_K, indirect CP violation in neutral-kaon mixing) into the freshly-built plugin scaffold. After you, two independent codex reviewers (physics + code) and a final Opus reviewer will check your work, so be correct and defensible. First write a SHORT plan, then implement.

Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing  (you are -C'd into it)

## Study first (read, do not modify):
- The scaffold: flavor_catalog_constraints/{base.py,anchors.py,registry.py,point_builder.py,TEMPLATE.py,README.md} and flavor_catalog_constraints/physics_adapters/{__init__.py,deltaf2.py}. Understand the ConstraintProtocol contract, ConstraintResult fields (numeric fields MUST be real floats; complex go in diagnostics), how @register + auto-discovery works, how load_anchor (schema-flex, typed, fails loudly) works, and how point_builder.build_from_quark_couplings produces a ParameterPoint.
- The EXISTING deletable example at flavor_catalog_constraints/primary/kaon/K001.py — you will replace it with the production implementation. Note how it wraps physics and loads its anchor.
- The catalog metadata sidecar: flavor_catalog/processes/kaon/K001.yaml — this is the SOURCE OF TRUTH for the experimental value (|epsilon|=2.228e-3), SM reference (BGS2020 2.161e-3), FLAG B_K inputs, and provenance. Your anchor MUST load from this yaml via the scaffold anchor loader, not hardcode numbers.
- The audited physics: quarkConstraints/deltaf2.py functions evaluate_epsilon_k / evaluate_epsilon_k_with_running and compute_delta_f2_wilsons, and the existing tests tests/test_epsilon_k_physics.py (see how a QuarkMassBasisCouplings input is constructed, e.g. the _sd_couplings helper).

## Physics (epsilon_K):
epsilon_K measures indirect CP violation in K0-K0bar mixing. The NP contribution is epsilon_K^NP = kappa_epsilon * Im(M12^NP) / (sqrt(2) * Delta m_K), with M12^NP from the kaon DeltaF=2 Wilson coefficients (KK-gluon + EW exchange in RS). The constraint: the NP contribution must fit within the room between the experimental |epsilon_K| and the SM prediction (HARD severity — observed bound). Use the repo's audited evaluate_epsilon_k machinery via the physics-adapter boundary; do NOT reinvent the DeltaF=2 matching.

## Implement:
1. Write a SHORT plan (to stdout, ~10-15 lines): the formula, what budget you use (exp vs SM-subtracted), severity, which adapter function, what the anchor candidates are, and the test plan.
2. Implement flavor_catalog_constraints/primary/kaon/K001.py as the PRODUCTION constraint (replacing the example): proper module docstring (physics + severity rationale + catalog sidecar ref + physics-core ref), a typed anchor loaded from K001.yaml via the scaffold loader, evaluate() that reaches physics ONLY through physics_adapters (add an epsilon_k wrapper to physics_adapters/deltaf2.py if one isn't there — append-only, re-export the upstream result dataclass unchanged), returns a ConstraintResult with real numeric fields (predicted=NP epsilon_K, sm_prediction, experimental, ratio=|pred|/budget, budget, severity=HARD, informative notes, diagnostics with any complex amplitudes). Handle the missing-couplings case gracefully (non-crashing result).
3. Write/upgrade tests/constraints/primary/kaon/test_K001.py: registration, anchor matches yaml, evaluate runs end-to-end on a real QuarkMassBasisCouplings (numeric fields real + finite), a known-safe point passes and a clearly-excluded large-NP point fails (numerical verification against the audited evaluate_epsilon_k), and evaluate is pure/deterministic.
4. Keep ALL changes isolated to: K001.py, its test, and (if needed) an append to physics_adapters/deltaf2.py. Do NOT touch other constraints or core physics modules.
5. Run `python -m pytest tests/constraints/ -q` and `python -m pytest tests/test_epsilon_k_physics.py -q` and confirm green. Fix until green.

## Output (to stdout, concise):
Your plan, the list of files changed, the key formula/convention you used, the pytest results (actual counts), and anything the reviewers should scrutinize. Do not paste full files.
