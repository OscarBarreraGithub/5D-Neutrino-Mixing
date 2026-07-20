# Phase 3 Scaffolding Plan — Round 2 Changelog

**Author:** Opus scaffolding architect (planner R2)
**Date:** 2026-05-28
**Source:** `.orchestration/PHASE3_SCAFFOLDING_PLAN.md` (rewritten in place)
**Driver:** `.orchestration/PHASE3_SCAFFOLDING_REVIEW_R1.md` critique

This file documents what changed between R1 and R2 of the scaffolding
plan and why. Each entry maps a reviewer finding to the plan section
that responds to it.

---

## Critical Fixes Applied

### CR-1 (blocker) — Schema-flex anchor loader

**Reviewer finding:** R1's `anchors.load()` hard-codes
`data["pdg_or_equivalent"]["canonical_experimental_value"]`, a path
present in exactly **1 of 102 yamls** (K001 only). Spot-check of 6
yamls across `kaon/`, `beauty/`, `edm_neutrino/`, `collider_rs/`,
`secondary/beauty/` confirmed ~30 distinct top-level keys under
`pdg_or_equivalent`.

**Plan change (§E rewritten):**
- Replaced the typed `ExperimentalAnchor` dataclass at the framework
  level with `anchors.load_raw(process_id, family, tier)` that returns
  the parsed `pdg_or_equivalent` mapping **verbatim**.
- Removed the typed `anchor: ExperimentalAnchor` field from
  `ConstraintBase`.
- Added a per-constraint typed-view pattern: each `<ID>.py` declares
  its own `@dataclass _Anchor` plus a `_build_anchor(raw)` function
  that picks the keys *that constraint* needs.
- Added four worked patterns (A: `canonical_experimental_value`,
  B: `canonical_experimental_average` flat, C: `canonical_limit` with
  limit_operator, D: `values:` list) covering the major schema
  families observed.
- Added an aggregated frequency table summarising the schema variety
  across all 102 yamls.
- Updated `test_anchor_matches_yaml` (§D) to be per-constraint: each
  test re-runs `_build_anchor(yaml.safe_load(...))` and pins the
  fields the constraint actually extracts.

**Spot-checked yamls** (confirms the variety):
- `flavor_catalog/processes/kaon/K001.yaml` → `canonical_experimental_value`
- `flavor_catalog/processes/beauty/B005.yaml` → flat `source:` + `canonical_experimental_average`
- `flavor_catalog/processes/edm_neutrino/E001.yaml` → `canonical_limit`
- `flavor_catalog/processes/edm_neutrino/E007.yaml` → `ra225_current_direct_limit`
- `flavor_catalog/processes/collider_rs/CR001.yaml` → `values:` (list)
- `flavor_catalog/processes/secondary/beauty/B007.yaml` → flat `source:` + `values:` (different shape from CR001)

### CR-2 (high) — `ParameterPoint.extras` owner

**Reviewer finding:** R1 left the owner ambiguous between
`scanParams/` and a new `point_builder.py`. The reviewer's
disposition: new file `flavor_catalog_constraints/point_builder.py`
plus a `ParameterPointExtras` TypedDict in `base.py`.

**Plan change:**
- Added `flavor_catalog_constraints/point_builder.py` to §A directory
  layout.
- Added a fully-spec'd `ParameterPointExtras` TypedDict to §B's
  `base.py` sketch, with `total=False` and a comment stating that
  adding a key is a one-line edit to this file.
- Added `point_builder.py` sketch to §F showing
  `build_parameter_point(scan_row) -> ParameterPoint` as the single
  producer of `extras`.
- Added `test_extras_keys_are_declared` to the global contract test
  in §D: a static-scan check that any `point.extras["<key>"]` literal
  in any constraint file matches a key in `ParameterPointExtras`.
- Updated the K001 worked example (§H) to consume
  `point.extras["quark_mass_basis_couplings"]` as a declared
  TypedDict key.
- Resolved Open Question 2 (§J) confirming `scanParams/` is the wrong
  home.

### CR-3 (medium) — Adapter signature mismatch in K001 walkthrough

**Reviewer finding:** R1 had
`epsilon_k_np_ratio() -> tuple[float, bool]` but the K001 walkthrough
unpacked it into four names. The reviewer recommended changing the
adapter to return the existing `EpsilonKResult` dataclass unchanged.

**Plan change (§E + §H):**
- Adapter convention: `physics_adapters/*` returns the **existing
  physics-module dataclasses unchanged**. The adapter's only job is
  the import boundary; no signature repackaging.
- §E `deltaf2_adapter.py` sketch updated:
  `evaluate_epsilon_k(couplings) -> EpsilonKResult` (matches
  `quarkConstraints/deltaf2.py:764` signature; verified by reading
  the source).
- §H K001 walkthrough rewritten to call the adapter and read fields
  by name from the `EpsilonKResult` dataclass: `result.passes`,
  `result.epsilon_k_np`, `result.ratio_to_budget`,
  `result.epsilon_k_np_budget`, `result.im_m12_np`. No tuple
  unpacking.

---

## Recommended Changes Applied

### REC-1 — Lazy discovery

**Plan change (§C):**
- `flavor_catalog_constraints/__init__.py` no longer calls
  `ConstraintRegistry.discover()` at package import.
- Added an idempotent `_DISCOVERED` flag inside `ConstraintRegistry`.
- All lookup methods (`get`, `all`, `filter`, `evaluate_all`) call
  `discover()` themselves so external callers can stay lazy.
- Added a session-scoped `autouse` `_discover_constraints` fixture to
  `tests/constraints/conftest.py` so per-constraint pytest runs only
  pay the ~100 ms discovery cost once.

### REC-2 — Snapshot-file existence check / provenance pin

**Plan change (§C + §D):**
- `register_constraint` decorator now calls `yaml_path_for(...)` and
  refuses to register if the yaml sidecar is missing (registration-time
  failure, not import-time silent acceptance).
- Added `test_snapshot_files_exist` to the global contract test:
  every path in `c.references` must resolve to an existing file under
  `flavor_catalog/references/`.

### REC-3 — Derive `level` from filesystem path

**Plan change (§B + §C):**
- Removed the `level` class attribute from the `ConstraintBase`
  Protocol. The decorator derives `level` from the module name
  (`flavor_catalog_constraints.primary.<family>.<ID>` →
  `ConstraintLevel.PRIMARY`).
- The decorator pins `level` as an instance attribute via
  `object.__setattr__`.
- `test_path_metadata_consistency` (global) asserts the derived
  `level` matches the filesystem tier.
- **Promotion SECONDARY → PRIMARY is now literally `git mv`** with
  zero code edits (was: "one-line edit" in R1).

### REC-4 — Derive `family` from filesystem path

**Plan change (§B + §C):**
- Same treatment as `level`: removed from the Protocol, derived by
  the decorator from `module.__name__`.
- Decorator validates `family` is in `_VALID_FAMILIES`
  (`beauty, charged_lepton, charm, collider_rs, edm_neutrino, kaon,
  top_higgs_ew`) and refuses unknown families.
- Three-way pin in `test_path_metadata_consistency`:
  `instance.family == yaml["family"] == module_dir_name`.
- Decorator also pins `process_id` validation against
  `^[A-Z]+[0-9]+$` (covers the REC-4 ID schema enforcement).

---

## Minor R1 Findings Folded In

- **N-1** Package name `flavor_catalog_constraints/`: kept.
- **N-2** No K016 anomaly: confirmed in §J; preamble corrected to
  "**102** catalogued constraints" (was "103 catalogued").
- **N-3** All seven family slugs are snake_case and match dir names:
  confirmed.
- **N-4** DEFERRED excluded from `evaluate_all` by default; flag
  `include_deferred=True` opts in. Pinned in §F driver API.
- **N-5** Dropped `complex` from `ConstraintResult.predicted` union
  (now `float | None`). Justified in §B.
- **N-6** Added §I.6 documenting the policy that at least one
  SECONDARY constraint must use `Severity.INFO` by end of the codex
  loop, or the enum value is dropped.
- **N-7** README will state the adapter-only import rule; flagged in
  §K as an orchestrator-level lint enforcement TODO.

## R1 "Things Missed" Folded In

- **M-1 severity provenance** — Added `test_severity_documented` to
  the global contract test; severity rationale must appear in the
  constraint module docstring. Flagged in §K that yaml-level
  provenance for severity is a future cleanup-wave decision.
- **M-2 schema-validation test** — Already covered by §C's
  `test_bijection_with_yaml`; the global test list in §D now
  references it explicitly alongside the 5 per-constraint mandatory
  tests.
- **M-3 invented extras key** — Fixed by CR-2 (TypedDict + global
  `test_extras_keys_are_declared`).
- **M-4 exception isolation** — `evaluate_all` now wraps each call
  in try/except, returning a `passes=False` ConstraintResult with a
  diagnostic note. Matches the existing per-draw isolation pattern
  in `scripts/run_rs_anarchy.py`.
- **M-5 import-failure tolerance** — `ConstraintRegistry.discover()`
  catches per-module exceptions into `_IMPORT_FAILURES`;
  `test_no_import_failures` flags any breakage by ID without
  poisoning the rest of the registry.
- **M-6 append-only adapters** — Documented in §E as a hard
  convention; README must say so. Modifying an existing adapter
  signature is a separate change unit outside the per-constraint
  codex loop.

---

## R1 Open Questions Resolved (§J)

- **Q1 (package name):** Confirmed `flavor_catalog_constraints/`.
- **Q2 (extras owner):** `point_builder.py` + `ParameterPointExtras`
  TypedDict. Both new artifacts spec'd in §A, §B, §F.
- **Q3 (ID / family slugs):** All 102 IDs match
  `^[A-Z]+[0-9]+$`. No K016 anomaly. Preamble inventory count
  corrected from "103" to "102".

---

## Items Declined or Deferred

- **Catalog yaml schema normalisation** (R1 CR-1 option C) — declined.
  The directive forbids touching catalog yaml/tex. The schema-flex
  loader (option A) achieves the same isolation guarantee without
  catalog edits.
- **Per-family resolver adapter** (R1 CR-1 option B) — declined in
  favour of per-constraint typed views. The catalog variety is not
  cleanly partitioned by family (B005 and B007 are both `beauty` but
  use different schemas; E001 and E007 are both `edm_neutrino` but
  use different keys). Per-constraint typed views match the actual
  variance.
- **Adding `severity:` field to the yaml schema** (R1 M-1) — deferred
  to a future cleanup wave. For now, severity is documented in the
  constraint module docstring and pinned by `test_severity_documented`.
- **Lint-rule enforcement for adapter-only imports** (R1 N-7) —
  flagged in §K as an orchestrator-level item; not load-bearing on
  the scaffold itself.

---

## Verdict

**Executable: yes.** All three critical issues are resolved, all
four recommended changes are folded in, and the R1 dispositions of
the three open questions are made authoritative in §J. The scaffold
remains a sibling-package + Protocol + decorator + mirrored-tests
design; the major R2 deltas are (a) the schema-flex anchor pattern,
(b) the `ParameterPointExtras` TypedDict + `point_builder.py`, and
(c) path-derived `level`/`family`. None of these require code
written before the orchestrator's per-constraint loop begins; they
are scaffolding directives the codex implementations consume on
day 1.
