# Phase 3 Scaffolding Plan — Review R1

**Reviewer:** Opus scaffolding reviewer
**Date:** 2026-05-28
**Subject:** `.orchestration/PHASE3_SCAFFOLDING_PLAN.md` (773 lines)

---

## Verdict

**APPROVE-WITH-CHANGES.** The scaffold design (sibling package, Protocol +
decorator, mirrored test tree, adapter moat) is structurally sound and
matches the directive's isolation requirements. The friction check, the
plug-pull invariant, and the worked K001 example all hold up. However,
one **critical** issue and several **high-priority** issues must be
addressed before code is written. They are concentrated in the
yaml-anchor contract (§E of the plan), which is built on a schema
assumption that does not hold across the 102 catalog entries.

---

## Critical Issues (Must Fix Before Implementation)

### CR-1. Anchor-loader schema mismatch (catalog-wide)

**Severity:** Blocker.

The plan's `anchors.py` (lines 444–478) hard-codes the path
`data["pdg_or_equivalent"]["canonical_experimental_value"]["value"]`
plus a fixed set of fields (`value`, `uncertainty`, `units`, `source`,
`source_url`, `snapshot_path`, `sha256_of_local_snapshot`). I audited
all 102 yamls. **Exactly one yaml** (`K001.yaml`) uses that key. The
remaining 101 use a heterogeneous menagerie:

```
21 yamls start `pdg_or_equivalent.year:` (flat, no canonical subkey)
19 yamls start `pdg_or_equivalent.source:` (flat)
10 yamls have an empty/null `pdg_or_equivalent`
 8 yamls use `primary_current_limit`
 5 yamls use `canonical_source`
 4 yamls use `canonical_current_limit`
 3 yamls use `canonical_experimental_average`
 3 yamls use `canonical_average`
 2 yamls use `canonical_limit`
 2 yamls use `canonical_direct_limit`
 1 yaml each: `canonical_experimental_value` (K001 only), `canonical_hflav_recommended`,
  `canonical_hflav_average`, `canonical_ratio`, `measured_experimental_anchor`,
  `experimental_anchors`, `branching_fraction_pdg_2025`, `values`, `summary`, ...
```

CR001 (a collider entry) holds a **list** under `values:` with
per-experiment exclusion limits — the loader would not even index into
it. B005, B006, B015–B019, etc. use a flat schema where
`pdg_or_equivalent` is essentially a single source-record at the top
level with no nested observable record.

If the anchor loader is built as currently sketched, **101 of 102
constraints will crash at import time** the moment auto-discovery walks
them. The orchestration loop will block at constraint #2.

**Fix.** One of:
- **(A) Schema-flex loader.** `anchors.load()` returns a thin
  ``RawAnchor`` that is the parsed `pdg_or_equivalent` dict verbatim,
  plus a per-constraint helper inside each `<ID>.py` that picks the
  right path. The Protocol carries `raw_anchor: Mapping[str, Any]`
  instead of a typed `ExperimentalAnchor`, and individual constraints
  define their own typed view. This is honest about the schema reality
  and lets implementation proceed.
- **(B) Per-family adapter.** Add `anchors/{family}.py` with one
  resolver per family that knows where the canonical observable lives
  for that yaml subset (e.g. `kaon` uses `canonical_experimental_value`,
  `collider_rs` uses `values[0]`, `b->X` decays use
  `canonical_experimental_average` or `canonical_average`, EDM uses
  `canonical_limit`). Five-to-eight resolver functions instead of one.
- **(C) Schema-normalize the yamls first** in a separate cleanup wave
  before scaffolding begins. The directive explicitly says yaml/tex are
  documentation-of-record and stay untouched; **this option violates
  the directive** and I do not recommend it.

I recommend **(A) + per-ID typed view inside each constraint module**.
It preserves the directive ("yaml is the source of truth"), keeps the
isolation invariant ("each constraint owns its anchor extraction"),
and avoids touching the catalog. The contract test in
`test_registry_contract.py` then asserts only the **path equality**
(yaml exists, sha matches), not field shape.

This change also makes `test_anchor_matches_yaml` (Test #2 in §D)
per-constraint: each constraint tests the fields *it* actually
extracts, not a one-size-fits-all schema.

---

### CR-2. `ParameterPoint.extras` owner is genuinely ambiguous

**Severity:** High.

Open question #2 names two candidates (`scanParams/` vs new
`flavor_catalog_constraints/point_builder.py`) but doesn't pick one.
The scan driver `scripts/run_rs_anarchy.py` lines 484–544 builds
`QuarkMassBasisCouplings` *inline* via a local
`_build_kk_gluon_couplings()` helper — there is no canonical
"`ParameterPoint` builder" anywhere in the repo today.

If the owner is left ambiguous:
- Codex agent #1 (implementing K001) will reach for
  `point.extras["quark_mass_basis_couplings"]` (as the worked example
  shows on line 626).
- Codex agent #2 (implementing some lepton constraint, e.g. L001)
  will reach for `point.extras["lepton_mass_basis_couplings"]` or
  `point.extras["c_lepton"]` or whatever the LFV machinery needs.
- These keys are invented per-constraint, never registered or
  documented. By constraint #20 the keys will collide / shadow / be
  misspelt and the orchestrator will be debugging key-typos instead
  of physics.

**Fix.** The plan must pick an owner *and* a key registry. Two-line
addition:

> `flavor_catalog_constraints/point_builder.py` owns
> `build_parameter_point(scan_row) -> ParameterPoint`.
> The keys of `extras` are declared in
> `flavor_catalog_constraints/base.py` as a frozen `ExtrasKey` Enum or
> a TypedDict (`ParameterPointExtras`). Adding a new key requires a
> one-line edit to that file, which is the only legitimate
> cross-constraint coordination point.

This actually *strengthens* isolation: a constraint cannot silently
require a new pre-compute; it must declare the key in the shared
TypedDict, which the registry contract test pins. It also keeps
`scripts/run_rs_anarchy.py` simple — the scan driver calls
`build_parameter_point()` and forgets about the constraint zoo.

`scanParams/` is the wrong home: it is currently a c-value /
delta-degeneracy scan helper, not a parameter-point factory. The
sibling `flavor_catalog_constraints/point_builder.py` is the right
home.

---

### CR-3. Adapter signature mismatch in the K001 worked example

**Severity:** Medium-but-trivial.

`physics_adapters/deltaf2_adapter.py` (lines 412, 417–419):
```python
def epsilon_k_np_ratio(couplings) -> tuple[float, bool]: ...
```
K001 example (line 628):
```python
ratio, passes, im_m12_np, eps_k_np = epsilon_k_np_ratio(couplings)
```
A 2-tuple cannot unpack into 4 names. Either:
- broaden the adapter return to a 4-tuple `(ratio, passes, im_m12_np,
  epsilon_k_np)`, or
- return the existing `EpsilonKResult` dataclass directly (the latter
  is what I recommend — the dataclass already carries every field K001
  needs, including `epsilon_k_np_budget`, and the adapter then becomes
  a 5-line function that just calls `evaluate_epsilon_k`).

The bug is small but it's exactly the kind of thing the codex
implementation agent will silently auto-fix in 30 different ways. Pin
the adapter return type now so the convention is uniform across all
δF=2 consumers (K001, K002 ΔmK, B-mixing, D-mixing).

**Fix.** Specify in §E that adapters return the existing physics-module
dataclasses unchanged. The adapter's only job is the import boundary,
not signature repackaging. K001 then reads:
```python
result = epsilon_k_adapter.evaluate(couplings, ...)  # returns EpsilonKResult
return ConstraintResult(
    process_id=..., passes=result.passes, ratio=result.ratio_to_budget,
    predicted=result.epsilon_k_np, budget=result.epsilon_k_np_budget, ...
)
```

---

## High-Priority Recommended Changes

### REC-1. The Protocol/decorator timing has a sharp edge

The decorator (line 221) instantiates the class at import:
```python
def register_constraint(cls):
    instance = cls()
    ...
```
But the constraint class's `__init__` (line 617) does
`self.anchor = load(self.process_id, family=self.family)`. So
**every import of every constraint module reads a yaml from disk**.

For 102 modules at import time this is:
- ~102 file reads + yaml parses (~100 ms with `lru_cache`)
- 102 instances pinned in `_REGISTRY` whether or not the caller will
  ever evaluate them

This is fine in production but it makes the *test loop* slow.
`pytest tests/constraints/primary/kaon/test_K001.py` will, via
`from flavor_catalog_constraints import ConstraintRegistry`, trigger
`ConstraintRegistry.discover()` (line 240) which walks **all 102
modules**. The codex loop's per-constraint pytest run pays a
~constant 100 ms tax.

Not a blocker. But the plan should add a fixture-level guard:
```python
# tests/constraints/conftest.py
@pytest.fixture(scope="session", autouse=True)
def discover_constraints():
    ConstraintRegistry.discover()
    yield
```
and `flavor_catalog_constraints/__init__.py` should **not**
auto-discover (line 275). Lazy import keeps `import
flavor_catalog_constraints` cheap; the test conftest opts in.
Production scan drivers call `ConstraintRegistry.discover()` once at
startup.

### REC-2. The 5 mandatory tests need a 6th: `test_anchor_provenance_resolves`

`test_anchor_matches_yaml` (Test #2) pins values and sha. But several
of the 102 yamls reference `snapshot_path` files that *might not be
checked in*. If the snapshot is missing, the constraint silently
imports (the loader only reads the yaml) and only crashes when a user
introspects `anchor.snapshot_path`. The contract test should additionally
assert `Path(repo_root/anchor.snapshot_path).is_file()` and (optionally)
that its on-disk sha256 matches `sha256_of_local_snapshot`.

This is the test pattern that catches the R03-I1 class of bug —
a number drifting away from its provenance.

### REC-3. `level` attribute vs directory placement — pick one source of truth

The plan puts the constraint in `primary/.../K001.py` *and* sets
`level = ConstraintLevel.PRIMARY` (line 613). Promoting SECONDARY →
PRIMARY (§G) says "git mv and edit one line."

But the auto-discovery walker doesn't know which tier it's in — it
just imports. If a constraint's directory is `secondary/.../K019.py`
but the class declares `level = PRIMARY`, the registry happily accepts
it. The two facts can disagree.

**Fix.** Either:
- (preferred) Derive `level` from `module.__name__`
  (`"...secondary..."` → SECONDARY) inside the decorator. The class no
  longer needs a `level` attribute. Promotion is **literally one
  `git mv`** with zero code edits.
- Or: contract test asserts `instance.level.value.lower() ==
  module_path_tier`.

This also resolves the planner's own concern that "SECONDARY must be
visible from a bare `ls`" — making the tree position load-bearing
makes that invariant cheap to enforce.

### REC-4. The `family` attribute likewise should be derived

Same argument: deriving `family` from the module's parent package
(`flavor_catalog_constraints.primary.kaon.K001` → `kaon`) means
moving a constraint between families is a `git mv`, and a constraint
cannot be mis-tagged. The codex loop is one less typo to worry about.

The contract test asserts `instance.family == yaml["family"] ==
module_dir_name`. Three-way pin.

---

## Minor Nits

- **N-1.** `flavor_catalog_constraints/` is a fine name. It is verbose
  but unambiguous, and `flavor_catalog/` is already the canonical
  short name for the catalog. Confirmed in §"Open question 1": **keep
  it.** Alternatives `processes_impl/` and `constraints_impl/` both
  break ID-grep affinity with the doc tree. `flavor_catalog/code/`
  must be rejected — see §A of the plan, which already argues this
  correctly.
- **N-2.** No `K016` anomaly exists in the catalog. Inventory pass
  confirms 102 IDs, all matching `^[A-Z]+[0-9]+$`. The planner's
  uncertainty in §"Open question 3" is unfounded.
- **N-3.** `family` values across the 102 yamls match the seven
  family directories exactly (`beauty`, `charged_lepton`, `charm`,
  `collider_rs`, `edm_neutrino`, `kaon`, `top_higgs_ew`). Snake_case
  is uniform. No naming surprises.
- **N-4.** `DEFERRED` constraints (§I, suspicious-flag mechanism)
  should be excluded from `evaluate_all` *by default* but reachable
  via `filter(level=ConstraintLevel.DEFERRED)`. The plan says this in
  passing; lock it down explicitly.
- **N-5.** `ConstraintResult.predicted` is typed `float | complex |
  None`. The downstream `as_ratio_dict()` consumer in `deltaf2.py`
  uses real floats — but for collider mass exclusions a "predicted
  KK-gluon mass" is real, for EDMs it's real, for branching ratios
  it's real. Drop `complex` from the union unless a concrete
  constraint requires it; it complicates `min/max/abs` consumers.
- **N-6.** `Severity.INFO` is declared but never used in the plan.
  Either show one constraint that uses it (e.g. a SECONDARY collider
  resonance) or drop it.
- **N-7.** Mention in the contributor `README.md` that
  `physics_adapters/` is the **only** place constraints may import
  `quarkConstraints.*` / `flavorConstraints.*` from. A ruff
  per-file-ignores or import-rule check enforces this; codex agents
  will need the hint.

---

## Dispositions on Planner's Open Questions

### Open Q1 — Package name
**Disposition:** `flavor_catalog_constraints/`. Verbose but accurate.
`constraints_impl/` loses family-tree affinity; `flavor_catalog/code/`
violates the catalog-untouched invariant. Confirmed: no collision with
existing top-level packages (`warpConfig`, `neutrinos`, `solvers`,
`diagonalization`, `yukawa`, `scanParams`, `quarkConstraints`,
`flavorConstraints`, `qcd`, `flavor_catalog`, `flavorConstraints`).

### Open Q2 — `ParameterPoint.extras` owner
**Disposition:** New file
`flavor_catalog_constraints/point_builder.py`. Keys frozen in
`base.py` as `ParameterPointExtras` TypedDict. `scanParams/` is the
wrong home (it's a scan-config helper, not a point factory).
See **CR-2** for the full prescription.

### Open Q3 — ID / family slug audit
**Disposition:**
- All 102 IDs are well-formed `^[A-Z]+[0-9]+$`. No K016 anomaly.
- All `family` values are snake_case and match directory names.
- Letter prefixes used: `B`, `C`, `CR`, `E`, `EW`, `K`, `L`, `T`. None
  collide. (CR is two letters; the regex above admits multi-letter
  prefixes, which is correct.)
- Total inventory: 94 PRIMARY + 8 SECONDARY = **102** (the planner's
  "103 catalogued" double-counts: the secondary tree contains 8 entries
  total, of which 4 are beauty / 3 kaon / 1 top_higgs_ew, all of which
  the planner's table accounts for. The "+1 extra" in the planner's
  preamble is spurious).

---

## Things the Planner Missed

### M-1. Severity is decoupled from yaml provenance
The yaml does not currently carry a `severity` field
(HARD/SOFT/INFO). The plan introduces severity as a Python attribute
on the constraint class. That's fine, but it means severity becomes a
**new** source of truth that does not appear in the catalog. The
scaffold should either:
- pick severity from a constraint-class attribute (planner's choice)
  AND add a contract test that asserts the choice is documented in
  the constraint module's docstring, **or**
- add a `severity:` field to the yaml schema (catalog-touching, not
  recommended).

Document the choice. Right now severity is silently invented.

### M-2. The 5-mandatory-test contract is missing a schema-validation test

There should be a test (probably in `test_registry_contract.py`, not
per-constraint) that asserts:
- every registered constraint's `process_id` matches a yaml file
- every yaml file in `flavor_catalog/processes/{primary,secondary}`
  has a registered constraint
- the union is exactly the inventory (no orphans either direction).

The plan mentions this in §C ("`test_registry_contract.py` asserts
that the set of registered IDs equals the set of yaml IDs") but the
mandatory per-constraint test list in §D should reference this
explicitly: the 5-test contract is per-constraint, the registry test
is once-globally. Both must exist and the README must say so.

### M-3. Worked-example K001 reads `point.extras["quark_mass_basis_couplings"]` without that key being declared anywhere

Line 626. The string literal is invented in §H. Per **CR-2** above,
the key must be in the frozen `ParameterPointExtras` TypedDict.

### M-4. The plan does not say what happens when `evaluate()` raises

If K001's evaluator divides by zero (M_KK=0), or KeyError on
`extras`, the current `evaluate_all` (line 272) will crash the entire
scan. The directive requires isolation: one constraint failure cannot
poison the others. `evaluate_all` should wrap each call:

```python
def evaluate_all(point):
    out = {}
    for pid, c in _REGISTRY.items():
        try:
            out[pid] = c.evaluate(point)
        except Exception as exc:
            out[pid] = ConstraintResult(
                process_id=pid, passes=False, severity=c.severity,
                predicted=None, sm_prediction=None, experimental=None,
                ratio=None, budget=None,
                notes=f"evaluate() raised {type(exc).__name__}: {exc}",
                diagnostics={"exception_type": type(exc).__name__},
            )
    return out
```
Add this to the plan in §F. The existing scan driver already does
per-draw isolation (`run_rs_anarchy.py` line 541 `except Exception as
err`). The pattern is already in the repo.

### M-5. Module-import side effects mean lint failures cannot be tolerated

If `K005.py` has a syntax error, `pkgutil.walk_packages` will raise on
import, breaking the entire registry. Codex agents *will* commit
broken syntax sometimes. The auto-discovery walker should catch
`ImportError`/`SyntaxError` per module, log it, and continue. The
contract test then **fails** for the offending ID but the rest of the
suite still runs.

```python
@staticmethod
def discover() -> None:
    for pkg in (primary, secondary):
        for _, modname, _ in pkgutil.walk_packages(pkg.__path__, pkg.__name__ + "."):
            try:
                importlib.import_module(modname)
            except Exception as exc:
                _IMPORT_FAILURES[modname] = exc
```
Then `test_registry_contract.py` asserts `_IMPORT_FAILURES == {}`. One
broken constraint flagged; the other 101 work.

### M-6. The codex orchestration loop assumption "agent never touches shared code" is partially false

The plan's claim (§I.4) is that codex working on K001 cannot touch
`quarkConstraints/`. True. But codex will need to add a new function
in `physics_adapters/deltaf2_adapter.py` when implementing K002 (ΔmK)
because the current adapter only exposes `epsilon_k_np_ratio`. That
*is* a cross-constraint touch — multiple kaon constraints share one
adapter file.

Mitigation: each adapter file is **append-only** (codex adds a new
function; never modifies an existing one). The plan should make this
explicit. Existing adapter functions are part of the public ABI and
require a separate change unit to modify.

If two codex agents try to extend the same adapter concurrently (the
directive forbids this, but it's a footgun), git will conflict-merge.
Acceptable.

---

## Friction Check (your question 8)

**Verdict:** With the CR-2 and REC-3/REC-4 fixes applied, **yes** — a
codex agent handed one constraint ID can work entirely within:
```
flavor_catalog_constraints/<tier>/<family>/<ID>.py
tests/constraints/<tier>/<family>/test_<ID>.py
```
plus *append* to one adapter file in
`flavor_catalog_constraints/physics_adapters/` when the constraint
needs a physics import. No edits to:
- `__init__.py` files (auto-discovery)
- `base.py` (Protocol is closed; only `ParameterPointExtras` TypedDict
  ever grows, and only by orchestrator decree)
- `registry.py` (decorator handles everything)
- `anchors.py` (constraint defines its own typed view inline per CR-1)
- `quarkConstraints/*`, `flavorConstraints/*` (firewall)
- `scripts/run_rs_anarchy.py` (it imports the registry and is
  agnostic of which constraints exist)

That meets the directive.

---

## Yaml-Source-of-Truth Check (your question 9)

The plan does the right thing in spirit (constraint reads yaml at
import via `anchors.load`), but the actual loader (line 467,
`exp = data["pdg_or_equivalent"]["canonical_experimental_value"]`)
**duplicates** values for K001 only and **fails to read** for the
other 101.

After CR-1 fix, the contract test (Test #2,
`test_anchor_matches_yaml`) becomes the architectural firewall that
prevents the R03-I1 class of bug ("two surfaces with the same
constant diverged"). Each constraint pins exactly the field(s) it
reads. This is the right pattern; the broken part is the loader, not
the contract.

---

## Whether One More Planner Round Is Needed

**Yes.** CR-1 alone requires a non-trivial rewrite of §E ("How a
Constraint Pulls Physics Inputs"). CR-2 requires a new file
(`point_builder.py`) and a frozen TypedDict in `base.py`. CR-3 is a
30-second fix. REC-1 through REC-4 are clean-up but worth doing
together.

One more planner round, then this is ready to execute.

---
