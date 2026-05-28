# Phase 3 Scaffolding Plan — Review R2

**Reviewer:** Opus scaffolding reviewer (round 2)
**Date:** 2026-05-28
**Subject:** `.orchestration/PHASE3_SCAFFOLDING_PLAN.md` (revised) + `PHASE3_SCAFFOLDING_CHANGES_R2.md`
**Prior:** `.orchestration/PHASE3_SCAFFOLDING_REVIEW_R1.md`

---

## Verdict

**APPROVE.** Ready for execution (scaffold code commit). All three R1
critical issues are resolved and all four R1 recommended changes are
folded in. The minor R1 nits and "things missed" items (N-1..N-7,
M-1..M-6) are also addressed. No new blockers; one orchestration note
below (non-blocking).

---

## Critical Issues: 3/3 Resolved

### CR-1 — Schema-flex anchor loader: RESOLVED

§E rewritten. `anchors.load_raw(process_id, family, tier)` returns the
parsed `pdg_or_equivalent` mapping verbatim (no schema assumptions
beyond presence of the top-level key). Each constraint defines its
own typed view inline via `_Anchor` dataclass + `_build_anchor(raw)`
function. Worked patterns A-D spot-checked against the catalog:

- Pattern A → `kaon/K001.yaml` uses `canonical_experimental_value:` (confirmed)
- Pattern B → `beauty/B005.yaml` uses flat `source:` + `canonical_experimental_average:` (confirmed)
- Pattern C → `edm_neutrino/E001.yaml` uses `canonical_limit:` (confirmed)
- Pattern D → `collider_rs/CR001.yaml` uses `values:` list (confirmed)

All four files exist on disk; the documented top-level keys match. The
102-yaml frequency tally and the secondary/beauty/B007.yaml example
add useful color. `test_anchor_matches_yaml` is now per-constraint and
pins the specific fields each constraint extracts (not a global
schema). The R1 blocker is dead.

### CR-2 — ParameterPoint.extras owner: RESOLVED

- `flavor_catalog_constraints/point_builder.py` is present in §A
  directory layout and fully sketched in §F as the single producer of
  `ParameterPoint.extras`.
- `ParameterPointExtras` TypedDict is in §B's `base.py` with
  `total=False`, comment marking it as the frozen registry, and four
  initial keys (`quark_mass_basis_couplings`, `lepton_mass_basis_couplings`,
  `kk_gluon_mass_gev`, `kk_w_mass_gev`, `kk_z_mass_gev`, `deltaf2_wilsons`).
- §J Open Question 2 is closed authoritatively.
- §D global contract test adds `test_extras_keys_are_declared`
  (static-scan check) so a constraint cannot invent a stray key.

### CR-3 — K001 walkthrough: RESOLVED

§E's `deltaf2_adapter.py` sketch now imports `EpsilonKResult` and
returns it unchanged from `evaluate_epsilon_k(couplings)`. §H K001
walkthrough reads `result.passes`, `result.epsilon_k_np`,
`result.ratio_to_budget`, `result.epsilon_k_np_budget`,
`result.im_m12_np` by field name — no tuple unpacking. The R1 4-vs-2
unpack mismatch is gone. §E adds the append-only adapter convention
(M-6) and pins the dependency direction.

---

## Recommended Changes: 4/4 Resolved

- **REC-1 lazy `discover()`:** PRESENT. §C introduces `_DISCOVERED`
  flag and removes auto-discovery from
  `flavor_catalog_constraints/__init__.py`. §D adds session-scoped
  autouse `_discover_constraints` fixture. Lookup methods trigger
  discovery on first use.
- **REC-2 path-derived `level` and `family`:** PRESENT. §C
  `_derive_tier_and_family_from_module()` derives both from
  `cls.__module__` and pins them as instance attributes via
  `object.__setattr__`. §B removes them from the Protocol. §G
  "Promote SECONDARY → PRIMARY" is now `git mv` with zero code edits.
- **REC-3 yaml-sidecar existence check at registration:** PRESENT.
  `register_constraint` calls `yaml_path_for(...).is_file()` and
  raises if missing. §D adds `test_snapshot_files_exist` to the
  global contract test.
- **REC-4 ID/family regex enforced:** PRESENT. `_ID_RE =
  re.compile(r"^[A-Z]+[0-9]+$")` and `_VALID_FAMILIES = {...}` both
  enforced in the decorator with loud `RuntimeError`s.

---

## Minor R1 Items (N-1..N-7, M-1..M-6): All Folded

- N-1..N-7: confirmed in §J and §K. Inventory count corrected to
  **102** (was "103"). `complex` dropped from
  `ConstraintResult.predicted`. `Severity.INFO` policy documented in
  §I.6 with a sunset clause.
- M-1..M-6: `test_severity_documented`, schema-bijection test,
  extras-key static scan, exception isolation in `evaluate_all`,
  `_IMPORT_FAILURES` capture, append-only adapter convention — all
  present.

---

## New Issues / Risks

None blocking. Two soft notes for the orchestrator:

1. `point_builder.build_parameter_point()` signature uses
   `build_quark_mass_basis_couplings(scan_row)` (line 933) — this
   function does not yet exist in `quarkConstraints.couplings`.
   Acceptable as a scaffold placeholder; the scan-driver integration
   step pins it later (acknowledged in §K open items).
2. The `from quarkConstraints.couplings import ...` at the top of
   `point_builder.py` is a Day-1 hard dependency. If
   `flavor_catalog_constraints` is ever imported in a context where
   `quarkConstraints/` is unavailable, the package will fail. Wrap
   that import inside the function body or guard it. Non-blocking.

---

## Friction Check

Re-confirmed: a codex agent handed one constraint ID writes exactly
two files (`<tier>/<family>/<ID>.py` + `tests/.../test_<ID>.py`) plus
at most one *append* to a `physics_adapters/*.py` file. No edits to
`__init__.py`, `base.py`, `registry.py`, or `anchors.py` unless the
constraint demands a new extras key (one-line edit to
`ParameterPointExtras`, orchestrator-routed). Meets the directive.

---

## Recommendation

Execute. Scaffold-code commit can begin. The orchestrator should:

1. Land `flavor_catalog_constraints/{base,registry,anchors,point_builder}.py`
   + `physics_adapters/deltaf2_adapter.py` + empty package
   `__init__.py` files as the first commit.
2. Land `tests/constraints/conftest.py` + `test_registry_contract.py`
   as the second commit (registry will report all 102 IDs missing —
   expected).
3. Begin per-constraint codex loop with `K001` as the dogfooding case.

