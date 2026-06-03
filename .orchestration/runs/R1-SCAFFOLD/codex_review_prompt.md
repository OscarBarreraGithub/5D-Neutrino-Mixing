# RETROACTIVE CODE REVIEW (codex, gpt-5.x xhigh). Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. Do NOT modify code — produce a precise findings list + verdict. This is the constraint-framework SCAFFOLD (commit 02e2424). It was built by one Claude/Opus agent and reviewed by another Claude/Opus agent, but was NEVER reviewed by codex. All 103 catalog constraints depend on it, so scrutinize the contract hard.

REVIEW these files:
- flavor_catalog_constraints/base.py — `ConstraintResult` (numeric fields MUST be real floats; `__post_init__` must reject complex; complex → diagnostics), `ConstraintProtocol`, `Severity`, `ConstraintLevel`.
- flavor_catalog_constraints/anchors.py — `load_anchor` (schema-flex, typed, MUST fail loud on missing/mismatched anchor: value_id, block_key, units, CL).
- flavor_catalog_constraints/registry.py — `@register`/auto-`discover`, `all_constraints`, `evaluate_all` (per-constraint try/except so one failure can't abort a batch), `by_family`/level helpers, `import_failures`.
- flavor_catalog_constraints/point_builder.py — `ParameterPoint` (frozen container), `build_from_quark_couplings`, KNOWN_EXTRA_KEYS.
- flavor_catalog_constraints/TEMPLATE.py and a couple of representative constraints to confirm the contract is actually enforced in practice.
- tests/constraints/ scaffold-level tests.

CHECK (with evidence, run things):
1. CONTRACT INTEGRITY: does `ConstraintResult.__post_init__` truly reject complex/NaN/Inf in numeric fields? Probe it directly (construct a result with a complex/NaN field → must raise). Are severity/level enums sound?
2. ANCHOR FAIL-LOUD: does `load_anchor` raise (not silently default) on: missing value_id, mismatched block_key, wrong units, wrong CL? Probe each. Any silent-accept path is a BLOCKER (a sibling constraint had exactly this bug).
3. REGISTRY ROBUSTNESS: does `evaluate_all` isolate per-constraint exceptions (one crash → that constraint's result carries the error, batch continues)? Probe by registering a deliberately-throwing dummy. Is auto-discovery deterministic? Does `import_failures` surface broken modules?
4. POINT_BUILDER / FROZEN CONTRACT: is `ParameterPoint` truly immutable; does it fail loud on unknown keys; is `build_from_quark_couplings` deterministic and pure?
5. DETERMINISM + PURITY of the framework; any global mutable state that could corrupt a 100M-point loop (caches keyed wrong, shared mutables, RNG).
6. Run `python -m pytest tests/constraints/ -q`; report counts. Registry import smoke.

OUTPUT (stdout, <=22 lines): numbered findings tagged BLOCKER/SHOULD-FIX/NIT with file:line + fix, the ACTUAL probe results (what raised / what didn't), pytest counts. End with: SCAFFOLD-OK or SCAFFOLD-NEEDS-FIXES.
