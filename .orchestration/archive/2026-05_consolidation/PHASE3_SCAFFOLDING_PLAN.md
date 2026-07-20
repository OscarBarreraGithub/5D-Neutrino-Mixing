# Phase 3 Scaffolding Plan — Flavor-Catalog Constraint Implementation

**Author:** Opus scaffolding architect
**Date:** 2026-05-28 (round 2 revision after reviewer critique)
**Repo:** `/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing`
**Inventory:** **102** catalogued constraints under `flavor_catalog/processes/`
(94 PRIMARY across `beauty/`, `charged_lepton/`, `charm/`, `collider_rs/`,
`edm_neutrino/`, `kaon/`, `top_higgs_ew/`; 8 SECONDARY under
`secondary/{beauty,kaon,top_higgs_ew}/`). All IDs match
`^[A-Z]+[0-9]+$` (prefixes used: B, C, CR, E, EW, K, L, T). The
round-1 planner's "+1 extra" and the suspected `K016` snake-case
anomaly are both **non-issues** — confirmed in R1 review.

This plan is **read-only design**; nothing is created yet. The directive
in `.orchestration/PHASE3_CONSTRAINTS_DIRECTIVE.md` requires that each
constraint be self-contained and atomic. Below is the scaffolding that
preserves that invariant while letting the orchestrator drive the
per-constraint codex loop with minimum friction.

**Revision note (R2).** Round 1 made a schema assumption about the yaml
anchor that does not hold across the catalog — exactly one yaml
(`K001.yaml`) uses the path
`pdg_or_equivalent.canonical_experimental_value`; the other 101 use
~30 different top-level keys (e.g. `canonical_average`,
`canonical_limit`, `primary_current_limit`, `values:` list,
`canonical_source`, flat `source:`/`year:`, etc.). Section E has been
**rewritten** to use a schema-flex loader; section B carries a frozen
`ParameterPointExtras` TypedDict; sections C/F/G have lazy-discovery,
path-derived metadata, and exception-isolated bulk evaluation. See
`.orchestration/PHASE3_SCAFFOLDING_CHANGES_R2.md` for the full
diff-against-R1.

---

## A. Directory Structure

### Chosen layout

```
flavor_catalog_constraints/                       # NEW top-level Python package
├── __init__.py                                   # exposes the registry + base types (NO discovery here)
├── base.py                                       # ConstraintBase, ConstraintResult, ParameterPoint, ParameterPointExtras, ConstraintLevel, Severity
├── registry.py                                   # ConstraintRegistry + @register_constraint decorator (lazy discover())
├── anchors.py                                    # Schema-flex YAML sidecar loader (returns raw mapping)
├── point_builder.py                              # NEW: build_parameter_point(scan_row) → ParameterPoint
├── physics_adapters/                             # thin wrappers that proxy quarkConstraints/* etc.
│   ├── __init__.py
│   ├── deltaf2_adapter.py                        # wraps quarkConstraints.deltaf2 calls; returns physics dataclasses unchanged
│   ├── meg_adapter.py                            # wraps flavorConstraints.muToEGamma
│   └── ...                                       # one adapter per shared physics module
├── primary/
│   ├── __init__.py                               # empty (or docstring only)
│   ├── beauty/
│   │   ├── __init__.py                           # empty
│   │   ├── B001.py
│   │   ├── B002.py
│   │   └── ...
│   ├── charged_lepton/
│   ├── charm/
│   ├── collider_rs/
│   ├── edm_neutrino/
│   ├── kaon/
│   │   ├── K001.py
│   │   └── ...
│   └── top_higgs_ew/
├── secondary/
│   ├── __init__.py                               # empty
│   ├── beauty/...
│   ├── kaon/...
│   └── top_higgs_ew/...
└── README.md                                     # one-page contributor guide

tests/constraints/                                # NEW, mirrors implementation tree
├── __init__.py
├── conftest.py                                   # session-scoped discover() fixture + sm_point/excluded_point fixtures
├── test_registry_contract.py                     # global: yaml↔registry bijection, no import failures, snapshot files exist
├── primary/
│   ├── kaon/
│   │   ├── test_K001.py
│   │   └── ...
│   └── ...
└── secondary/
    └── ...
```

The catalog yaml/tex remain in place at
`flavor_catalog/processes/<family>/<id>.{yaml,tex}` — they stay **pure
documentation** and continue to be the canonical anchor source.
The Python implementations live in a **sibling package**
`flavor_catalog_constraints/` that *reads* the yamls at runtime but
does not own them.

### Why this layout (vs alternatives)

| Option | Pros | Cons | Verdict |
|--------|------|------|---------|
| **Chosen:** sibling `flavor_catalog_constraints/{tier}/{family}/{id}.py` | Mirrors the catalog tree exactly, family namespace prevents ID collisions, tier visible from `ls`, one file per constraint = trivial add/delete, leaves the existing PDF/website build untouched, tier/family can be derived from path | Two trees to keep in sync (catalog vs implementation); reconciled by `test_registry_contract.py` that asserts 1:1 mapping | **Selected** |
| Promote each `processes/<family>/<id>/` into a Python package containing yaml+tex+constraint.py | Maximum locality | Breaks current `\input{processes/<family>/<id>}` LaTeX path; forces a flavor_catalog refactor | Rejected — touches finalised catalog |
| Flat `constraints/<id>.py` | Trivial imports | 102 files in one directory; loses family grouping | Rejected |
| Subpackage inside `flavorConstraints/` | Reuses an existing namespace | `flavorConstraints/` currently holds only `muToEGamma`; mixing 102 new modules into a 1-module package muddies semantics | Rejected |

### Updating `pyproject.toml`

Add `flavor_catalog_constraints*` to `tool.setuptools.packages.find.include`.
One-line edit; no other build changes.

---

## B. Constraint Interface

### `base.py` (sketch)

```python
# flavor_catalog_constraints/base.py
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Mapping, Protocol, TypedDict


class ConstraintLevel(str, Enum):
    PRIMARY = "PRIMARY"
    SECONDARY = "SECONDARY"
    DEFERRED = "DEFERRED"


class Severity(str, Enum):
    HARD = "HARD"   # observed bound, NP must fit inside experimental error
    SOFT = "SOFT"   # SM-vs-exp tension / projection / theory-dominated bound
    INFO = "INFO"   # informational only, never fails a point


# ------------------------------------------------------------------ #
# Frozen extras key registry.
# Adding a new key REQUIRES editing this TypedDict. This is the only
# legitimate cross-constraint coordination point in the scaffold.
# `total=False` because not every point pre-computes every helper.
# ------------------------------------------------------------------ #
class ParameterPointExtras(TypedDict, total=False):
    quark_mass_basis_couplings: Any   # QuarkMassBasisCouplings
    lepton_mass_basis_couplings: Any  # LeptonMassBasisCouplings (future)
    kk_gluon_mass_gev: float
    kk_w_mass_gev: float
    kk_z_mass_gev: float
    deltaf2_wilsons: Any              # DeltaF2WilsonCoefficients (cached)
    # ... every new key requires a one-line edit here.


@dataclass(frozen=True)
class ParameterPoint:
    """A single RS parameter point handed to every constraint.

    `raw` is whatever the scan driver produces (warpConfig.Setup +
    scanParams row, opaque to constraints). `extras` is a TypedDict
    mapping (declared in ParameterPointExtras above) of pre-computed
    intermediates the constraint zoo shares.
    """
    raw: Any
    extras: ParameterPointExtras = field(default_factory=dict)


@dataclass(frozen=True)
class ConstraintResult:
    process_id: str
    passes: bool
    predicted: float | None              # NP prediction, observable units (real)
    sm_prediction: float | None          # SM central value where known
    experimental: float | None           # central experimental value (or bound)
    ratio: float | None                  # |predicted| / budget (or similar)
    budget: float | None                 # bound used in the ratio
    severity: Severity
    notes: str = ""
    diagnostics: Mapping[str, Any] = field(default_factory=dict)


class ConstraintBase(Protocol):
    """Every constraint module defines exactly one class implementing this."""

    process_id: str            # "K001"  (REQUIRED class attr)
    severity: Severity         # HARD/SOFT/INFO (REQUIRED class attr)
    observable: str            # "epsilon_K", "BR(mu->e gamma)", ...
    # NOT required as class attrs (derived from module path in the registry):
    #   family, level
    # NOT required as a typed object (varies per yaml schema):
    #   anchor  →  raw mapping accessed via flavor_catalog_constraints.anchors.load_raw()

    def evaluate(self, point: ParameterPoint) -> ConstraintResult: ...
```

### Why a Protocol-plus-decorator instead of an ABC

Each `K001.py` declares a *plain class* — no inheritance. A typo in one
constraint cannot poison a shared parent. The Protocol still gives full
static type checking via mypy.

### Why `ParameterPointExtras` is a TypedDict (R1 CR-2 fix)

Round 1 left the owner of `extras` ambiguous. R1 review (CR-2) flagged
this as high-priority. **R2 resolution:**
- The TypedDict above is the single registry of valid extras keys.
- `point_builder.py` (new file, §E) is the only producer.
- Codex agents implementing a constraint must look up the key in
  `base.py` and use it verbatim — typos are caught by mypy and by
  `test_registry_contract.py` (which asserts every key referenced in
  any constraint's `evaluate()` source matches a TypedDict key).
- A constraint that needs a *new* extras key edits **this one file**
  to add it; the orchestrator catches that diff and routes it through
  the same sign-off pipeline as a constraint.

### Why a ratio-based result

`quarkConstraints/deltaf2.py` already uses ratio-to-bound semantics
(`EpsilonKResult.ratio_to_budget`, `passes`). Reusing that shape lets
adapters pass the existing dataclasses straight through and the
constraint just repackages a few fields.

### Why `predicted` drops `complex` (R1 nit N-5)

Every observable in the 102-constraint inventory is real-valued
(branching ratio, mass exclusion, EDM bound, ratio, asymmetry, decay
width). `complex` would force `min/max/abs` consumers to special-case.
If a future constraint needs a complex amplitude, it goes in
`diagnostics`.

---

## C. Discovery / Plugin Registration

### Mechanism: decorator + **lazy** auto-import + path-derived metadata

R1 raised three points addressed below:
- REC-1: don't eager-walk at package import (slow tests).
- REC-3: derive `level` from filesystem path, not class attribute.
- REC-4: derive `family` from filesystem path likewise.
- M-5: catch import errors per-module so one broken constraint doesn't
  break the whole registry.

```python
# flavor_catalog_constraints/registry.py
from __future__ import annotations
import importlib
import pkgutil
import re
from typing import Any, Dict, Iterable

from .base import ConstraintBase, ConstraintLevel, ParameterPoint, ConstraintResult, Severity


_REGISTRY: Dict[str, ConstraintBase] = {}
_IMPORT_FAILURES: Dict[str, BaseException] = {}
_DISCOVERED: bool = False

_ID_RE = re.compile(r"^[A-Z]+[0-9]+$")
_VALID_FAMILIES = {
    "beauty", "charged_lepton", "charm", "collider_rs",
    "edm_neutrino", "kaon", "top_higgs_ew",
}


def _derive_tier_and_family_from_module(module_name: str) -> tuple[ConstraintLevel, str]:
    """
    flavor_catalog_constraints.primary.kaon.K001     -> (PRIMARY, "kaon")
    flavor_catalog_constraints.secondary.beauty.B007 -> (SECONDARY, "beauty")
    """
    parts = module_name.split(".")
    # ["flavor_catalog_constraints", "<tier>", "<family>", "<ID>"]
    if len(parts) < 4:
        raise RuntimeError(f"Cannot derive tier/family from module name {module_name!r}")
    tier_str, family = parts[1], parts[2]
    tier_map = {"primary": ConstraintLevel.PRIMARY, "secondary": ConstraintLevel.SECONDARY}
    if tier_str not in tier_map:
        raise RuntimeError(f"Unknown tier '{tier_str}' in module {module_name!r}")
    if family not in _VALID_FAMILIES:
        raise RuntimeError(f"Unknown family '{family}' in module {module_name!r}")
    return tier_map[tier_str], family


def register_constraint(cls):
    """Decorator: instantiate the constraint, derive tier/family from path, register.

    The decorator runs at module import. It:
    - validates process_id format (REC-4)
    - derives level + family from the module path (REC-3/REC-4)
    - checks the yaml sidecar exists (REC-3 snapshot check)
    - checks for duplicate IDs
    - attaches `level` and `family` as instance attributes
    """
    instance = cls()

    # Validate process_id
    pid = getattr(instance, "process_id", None)
    if not isinstance(pid, str) or not _ID_RE.match(pid):
        raise RuntimeError(
            f"{cls.__module__}: process_id {pid!r} does not match {_ID_RE.pattern}"
        )

    # Derive tier + family from module path; pin as instance attributes
    level, family = _derive_tier_and_family_from_module(cls.__module__)
    object.__setattr__(instance, "level", level)
    object.__setattr__(instance, "family", family)

    # Check yaml sidecar exists (REC-3)
    from .anchors import yaml_path_for
    yaml_path = yaml_path_for(pid, family=family, tier=level)
    if not yaml_path.is_file():
        raise RuntimeError(
            f"{cls.__module__}: yaml sidecar missing at {yaml_path}"
        )

    # Duplicate guard
    if pid in _REGISTRY:
        prev = _REGISTRY[pid].__class__.__module__
        raise RuntimeError(
            f"Duplicate constraint registration for {pid}: "
            f"already provided by {prev}, now also by {cls.__module__}"
        )

    _REGISTRY[pid] = instance
    return cls


class ConstraintRegistry:
    @staticmethod
    def discover() -> None:
        """Walk flavor_catalog_constraints.{primary,secondary} and import all leaf modules.

        Idempotent: subsequent calls are no-ops. Per-module ImportError is
        captured into _IMPORT_FAILURES (M-5) so one broken file doesn't
        poison the whole registry; the contract test then fails for the
        broken ID while the other 101 still work.
        """
        global _DISCOVERED
        if _DISCOVERED:
            return
        import flavor_catalog_constraints.primary as primary
        import flavor_catalog_constraints.secondary as secondary
        for pkg in (primary, secondary):
            for _, modname, ispkg in pkgutil.walk_packages(pkg.__path__, pkg.__name__ + "."):
                if ispkg:
                    continue  # only leaf modules carry constraints
                try:
                    importlib.import_module(modname)
                except Exception as exc:
                    _IMPORT_FAILURES[modname] = exc
        _DISCOVERED = True

    @staticmethod
    def import_failures() -> Dict[str, BaseException]:
        return dict(_IMPORT_FAILURES)

    @staticmethod
    def all() -> Dict[str, ConstraintBase]:
        ConstraintRegistry.discover()
        return dict(_REGISTRY)

    @staticmethod
    def get(process_id: str) -> ConstraintBase:
        ConstraintRegistry.discover()
        return _REGISTRY[process_id]

    @staticmethod
    def filter(
        *,
        family: str | None = None,
        level: ConstraintLevel | None = None,
        severity: Severity | None = None,
    ) -> Dict[str, ConstraintBase]:
        ConstraintRegistry.discover()
        out: Dict[str, ConstraintBase] = {}
        for pid, c in _REGISTRY.items():
            if family is not None and c.family != family: continue
            if level is not None and c.level != level: continue
            if severity is not None and c.severity != severity: continue
            out[pid] = c
        return out

    @staticmethod
    def evaluate_all(
        point: ParameterPoint,
        *,
        include_deferred: bool = False,
    ) -> Dict[str, ConstraintResult]:
        """Evaluate every registered constraint at `point`, isolating exceptions.

        Per directive isolation (R1 M-4): a per-constraint exception is
        captured into a `passes=False` ConstraintResult with a diagnostic
        note; the other 101 still run. DEFERRED constraints are skipped
        by default (R1 N-4).
        """
        ConstraintRegistry.discover()
        out: Dict[str, ConstraintResult] = {}
        for pid, c in _REGISTRY.items():
            if c.level == ConstraintLevel.DEFERRED and not include_deferred:
                continue
            try:
                out[pid] = c.evaluate(point)
            except Exception as exc:
                out[pid] = ConstraintResult(
                    process_id=pid, passes=False,
                    predicted=None, sm_prediction=None, experimental=None,
                    ratio=None, budget=None,
                    severity=getattr(c, "severity", Severity.HARD),
                    notes=f"evaluate() raised {type(exc).__name__}: {exc}",
                    diagnostics={"exception_type": type(exc).__name__},
                )
        return out
```

`flavor_catalog_constraints/__init__.py` exposes the registry types
but **does NOT** call `discover()` (R1 REC-1). Tests opt in via a
session-scoped autouse fixture; production drivers call
`ConstraintRegistry.discover()` once at startup; lazy callers
(`.get()`, `.all()`, `.filter()`, `.evaluate_all()`) trigger discovery
on first use anyway.

```python
# flavor_catalog_constraints/__init__.py
from .base import (
    ConstraintBase, ConstraintLevel, ConstraintResult,
    ParameterPoint, ParameterPointExtras, Severity,
)
from .registry import ConstraintRegistry, register_constraint

__all__ = [
    "ConstraintBase", "ConstraintLevel", "ConstraintResult",
    "ParameterPoint", "ParameterPointExtras", "Severity",
    "ConstraintRegistry", "register_constraint",
]
```

### Why decorator-based wins for the directive

The directive requires "add a new constraint" to be a one-file change.
With this mechanism:
- create `flavor_catalog_constraints/primary/<family>/<NEW_ID>.py`
- write a class declaring **only** `process_id`, `severity`,
  `observable` (and an `evaluate` method); decorate with
  `@register_constraint`
- done — no manifest edit, no `__init__.py` edit, no registry list,
  no `level`/`family` to keep in sync with the directory

Removal is one-file: delete the `.py` file (and its test file).

### Failure modes guarded

- duplicate IDs raise at import (loud failure)
- malformed `process_id` (not matching `^[A-Z]+[0-9]+$`) raises at
  registration
- missing yaml sidecar raises at registration (REC-3)
- module-level syntax error → captured in `_IMPORT_FAILURES`, surfaced
  by the contract test (M-5)
- the per-family `__init__.py` is **empty** so removing a module
  never leaves a dangling import
- `test_registry_contract.py` asserts that:
  - `_IMPORT_FAILURES == {}` (no broken modules)
  - the set of registered IDs equals the set of yaml IDs in
    `flavor_catalog/processes/`
  - for each registered constraint, `c.level` matches its filesystem
    tier and `c.family` matches its filesystem family

---

## D. Test Convention

### Location

Tests mirror the implementation tree under `tests/constraints/`:

```
tests/constraints/primary/kaon/test_K001.py
tests/constraints/primary/beauty/test_B001.py
tests/constraints/secondary/beauty/test_B007.py
```

Tests are **not** colocated with the constraint modules (repo
convention; keeps test code out of distributions; lets the
orchestrator scope per-constraint pytest to
`tests/constraints/<tier>/<family>/test_<id>.py`).

### Session-level conftest

```python
# tests/constraints/conftest.py
import pytest
from flavor_catalog_constraints import ConstraintRegistry, ParameterPoint

@pytest.fixture(scope="session", autouse=True)
def _discover_constraints():
    """Trigger discovery exactly once per pytest session (R1 REC-1)."""
    ConstraintRegistry.discover()

@pytest.fixture
def sm_point() -> ParameterPoint: ...
@pytest.fixture
def excluded_point() -> ParameterPoint: ...
```

`sm_point` and `excluded_point` are backed by `tests/golden/`
(directory already exists). Both fixtures populate
`ParameterPointExtras` with the canonical pre-computed helpers so
constraint tests do not reinvent benchmark points.

### Minimum test contract (every constraint)

Each `test_<ID>.py` must include at least these tests:

1. `test_registered` — `ConstraintRegistry.get("<ID>")` returns an
   instance, and its `process_id`, derived `family`, derived `level`,
   and `severity` are correct.
2. `test_anchor_matches_yaml` — the constraint's **typed view** of its
   yaml (see §E for the per-constraint anchor pattern) equals the
   values pulled directly via `yaml.safe_load`. This pins each field
   the constraint actually reads, against the yaml. Schema-flex
   loaders mean each constraint tests *its own* fields — not a
   one-size-fits-all schema (R1 CR-1).
3. `test_sm_point_passes` — at the SM-like reference point the
   constraint must `passes=True` and `ratio < 1`.
4. `test_extreme_np_point_fails` — at an obviously excluded point
   the constraint must `passes=False`.
5. `test_evaluate_is_pure` — calling `evaluate` twice on the same
   point yields identical `ConstraintResult` (deterministic, no
   hidden state).

### Global contract test (`test_registry_contract.py`)

This single file enforces the cross-cutting invariants once, not
per-constraint:

1. `test_no_import_failures` — `ConstraintRegistry.import_failures()
   == {}` (R1 M-5).
2. `test_bijection_with_yaml` — set of registered IDs equals set of
   yaml IDs under `flavor_catalog/processes/{,secondary/}*/`.
3. `test_path_metadata_consistency` — for each registered constraint:
   - `c.family` matches the yaml's `family:` value
   - `c.family` matches the constraint module's parent dir name
   - `c.level` matches the constraint module's tier dir
4. `test_snapshot_files_exist` — for each constraint, the snapshot
   paths it advertises in `c.references` resolve to existing files
   under `flavor_catalog/references/` (R1 REC-3 + R1 REC-2 — the
   provenance pin).
5. `test_severity_documented` — for each constraint, the choice of
   `severity` is documented in the constraint module's docstring
   (R1 M-1: severity is not in the yaml; this test prevents silent
   invention).
6. `test_extras_keys_are_declared` — static-scan check that any
   `point.extras["<key>"]` literal referenced in
   `flavor_catalog_constraints/{primary,secondary}/**/*.py` is a key
   in `ParameterPointExtras` (R1 M-3).

---

## E. How a Constraint Pulls Physics Inputs (REWRITTEN for R1 CR-1)

### The three input channels

Each constraint module needs three things:

1. **Experimental anchor (raw yaml mapping)** — loaded by
   `anchors.load_raw(process_id, family, tier)`, which returns the
   parsed `pdg_or_equivalent` dict **verbatim**. Each constraint
   defines its own typed view inline (a small dataclass or NamedTuple
   that pulls the fields *that constraint* needs). The yaml stays the
   single source of truth; the loader makes no schema assumptions.
2. **SM prediction** — if present in the yaml (under whatever key the
   constraint expects, e.g. `standard_model_reference`,
   `standard_model_prediction`, etc.), pulled by the same typed view.
   Constraints that need SM helpers go via a physics adapter.
3. **NP contribution from the parameter point** — computed via a
   physics adapter in `physics_adapters/`.

### The schema-flex anchor loader

Round 1 hard-coded `data["pdg_or_equivalent"]["canonical_experimental_value"]`,
a path that exists in **exactly 1 of 102 yamls**. Round 2 replaces it
with a verbatim raw loader and pushes the path knowledge into each
constraint:

```python
# flavor_catalog_constraints/anchors.py
from __future__ import annotations
import functools
from pathlib import Path
from typing import Any, Mapping

import yaml

from .base import ConstraintLevel

REPO_ROOT = Path(__file__).resolve().parents[1]
CATALOG_ROOT = REPO_ROOT / "flavor_catalog" / "processes"


def yaml_path_for(
    process_id: str,
    *,
    family: str,
    tier: ConstraintLevel | str,
) -> Path:
    """Return the absolute path to a constraint's yaml sidecar."""
    tier_str = tier.value.lower() if isinstance(tier, ConstraintLevel) else str(tier).lower()
    if tier_str == "primary":
        return CATALOG_ROOT / family / f"{process_id}.yaml"
    elif tier_str == "secondary":
        return CATALOG_ROOT / "secondary" / family / f"{process_id}.yaml"
    raise ValueError(f"Unknown tier {tier!r}")


@functools.lru_cache(maxsize=None)
def _load_full_yaml(path_str: str) -> Mapping[str, Any]:
    with open(path_str) as f:
        return yaml.safe_load(f)


def load_raw(
    process_id: str,
    *,
    family: str,
    tier: ConstraintLevel | str = ConstraintLevel.PRIMARY,
) -> Mapping[str, Any]:
    """Return the parsed `pdg_or_equivalent` dict verbatim.

    The schema of this dict varies across the catalog. Each constraint
    is responsible for picking the field(s) it needs. This loader makes
    no assumptions beyond the presence of the top-level
    `pdg_or_equivalent` key (verified to be present in all 102 yamls).
    """
    path = yaml_path_for(process_id, family=family, tier=tier)
    data = _load_full_yaml(str(path))
    raw = data.get("pdg_or_equivalent")
    if raw is None:
        raise KeyError(
            f"{path}: missing top-level 'pdg_or_equivalent' key"
        )
    return raw


def load_full(
    process_id: str,
    *,
    family: str,
    tier: ConstraintLevel | str = ConstraintLevel.PRIMARY,
) -> Mapping[str, Any]:
    """Return the full parsed yaml document. Used by the contract test."""
    path = yaml_path_for(process_id, family=family, tier=tier)
    return _load_full_yaml(str(path))
```

**That is the entire loader.** No `ExperimentalAnchor` dataclass at
the framework level; no field schema. Each constraint inlines a
typed view appropriate to *its* yaml subset.

### Per-constraint typed view pattern

Each `<ID>.py` defines a small dataclass or NamedTuple that captures
the fields *this constraint* needs from the raw anchor. Example
patterns the codex loop will produce:

```python
# Pattern A — K001-style (canonical_experimental_value subkey)
from dataclasses import dataclass

@dataclass(frozen=True)
class K001Anchor:
    epsilon_K_exp_value: float
    epsilon_K_exp_uncertainty: float
    epsilon_K_sm_value: float
    snapshot_path_exp: str
    sha256_exp: str

def _build_anchor(raw):
    exp = raw["canonical_experimental_value"]
    sm = raw["standard_model_reference"]
    return K001Anchor(
        epsilon_K_exp_value=float(exp["value"]),
        epsilon_K_exp_uncertainty=float(exp["uncertainty"]),
        epsilon_K_sm_value=float(sm["value"]),
        snapshot_path_exp=exp["snapshot_path"],
        sha256_exp=exp["sha256_of_local_snapshot"],
    )

# Pattern B — B005-style (flat schema with canonical_experimental_average)
@dataclass(frozen=True)
class B005Anchor:
    br_exp_value: float
    br_exp_uncertainty: float
    br_sm_value: float
    snapshot_path_exp: str
    sha256_exp: str

def _build_anchor(raw):
    exp = raw["canonical_experimental_average"]
    sm = raw["standard_model_prediction"]
    return B005Anchor(
        br_exp_value=float(exp["value"]),
        br_exp_uncertainty=float(exp["uncertainty"]),
        br_sm_value=float(sm["value"]),
        snapshot_path_exp=exp["snapshot_path"],
        sha256_exp=exp["sha256_of_text_snapshot"],
    )

# Pattern C — E001-style (canonical_limit, no SM, with limit_operator)
@dataclass(frozen=True)
class E001Anchor:
    de_upper_limit: float    # in e cm
    confidence_level: str    # "90%"
    snapshot_path: str
    sha256: str

def _build_anchor(raw):
    lim = raw["canonical_limit"]
    return E001Anchor(
        de_upper_limit=float(lim["value"]),
        confidence_level=lim["confidence_level"],
        snapshot_path=lim["snapshot_path"],
        sha256=lim["sha256_of_local_snapshot"],
    )

# Pattern D — CR001-style (values: [list of dicts])
@dataclass(frozen=True)
class CR001Anchor:
    excluded_mass_upper_tev: float
    confidence_level: str
    snapshot_path: str
    sha256: str

def _build_anchor(raw):
    # CR001 uses values: [...] with one or more dicts. Take the first
    # (or filter by value_id) per the constraint's design.
    first = raw["values"][0]
    return CR001Anchor(
        excluded_mass_upper_tev=float(first["value"]),
        confidence_level=first["cl"],
        snapshot_path=first["snapshot_path"],
        sha256=first["sha256"],
    )
```

The yaml schema variety is **real and large** — spot-check confirmed:

| Yaml | Top key inside `pdg_or_equivalent` |
|------|-------------------------------------|
| `kaon/K001.yaml` | `canonical_experimental_value:` (nested) |
| `beauty/B005.yaml` | flat `source:` + `canonical_experimental_average:` (and `hflav_historical_average`, `standard_model_prediction`, `input_measurements: [...]`) |
| `edm_neutrino/E001.yaml` | `canonical_limit:` (with `limit_operator`, no `uncertainty`) |
| `edm_neutrino/E007.yaml` | `ra225_current_direct_limit:` (family-prefixed key) |
| `collider_rs/CR001.yaml` | `values:` (list of per-experiment dicts) |
| `secondary/beauty/B007.yaml` | flat `source:` + `values:` (different shape from CR001) |

Aggregated tally across all 102 yamls (first non-comment key under
`pdg_or_equivalent`, sorted by frequency): `year:` (21), flat `source:`
(19), `primary_current_limit:` (8), `canonical_source:` (5),
`canonical_current_limit:` (4), `canonical_experimental_average:` (3),
`canonical_average:` (3), `values:` (3), `canonical_limit:` (2),
`canonical_direct_limit:` (2), and ~15 more singletons including
`canonical_experimental_value:` (K001 only),
`canonical_hflav_recommended:`, `canonical_ratio:`,
`measured_experimental_anchor:`, `experimental_anchors:`,
`branching_fraction_pdg_2025:`, `summary:`, `values:` (list), etc.

The implication: a one-size-fits-all loader is impossible without
catalog rewrites (which the directive forbids). The schema-flex
pattern above is the only correct design.

### Contract test pins the per-constraint view

`test_anchor_matches_yaml` in each `test_<ID>.py` re-runs
`_build_anchor(yaml.safe_load(open(path)).pdg_or_equivalent)` and
compares to the constraint's stored typed view. The contract is "the
constraint extracts what it claims to extract"; the schema details
live where they belong, with the constraint that owns them.

### Physics adapters (the isolation moat)

The directive: modifying one constraint cannot break another. Shared
physics modules are the risk vector. Solution: a thin **adapter
layer** in `flavor_catalog_constraints/physics_adapters/`. Each
adapter returns the **existing physics-module dataclasses unchanged**
(R1 CR-3 fix). The adapter's *only* job is the import boundary; it
does not repackage signatures.

```python
# flavor_catalog_constraints/physics_adapters/deltaf2_adapter.py
"""Import boundary for quarkConstraints.deltaf2.

Constraints get the existing dataclasses (EpsilonKResult, DeltaMKResult,
BMixingResult) unchanged. The adapter exists so a signature change in
quarkConstraints/deltaf2.py touches *one* file, not 20.
"""
from quarkConstraints.deltaf2 import (
    DeltaF2WilsonCoefficients,
    EpsilonKResult,
    DeltaMKResult,
    BMixingResult,
    DEFAULT_DELTA_F2_INPUTS_V1,
    compute_delta_f2_wilsons,
    evaluate_epsilon_k as _evaluate_epsilon_k,
    evaluate_delta_mk as _evaluate_delta_mk,
    evaluate_bd_mixing as _evaluate_bd_mixing,
    evaluate_bs_mixing as _evaluate_bs_mixing,
)
from quarkConstraints.couplings import QuarkMassBasisCouplings


def evaluate_epsilon_k(
    couplings: QuarkMassBasisCouplings,
    *,
    inputs=DEFAULT_DELTA_F2_INPUTS_V1,
    epsilon_k_np_budget_override: float | None = None,
) -> EpsilonKResult:
    """Pass-through to quarkConstraints.deltaf2.evaluate_epsilon_k."""
    wilsons = compute_delta_f2_wilsons(couplings, inputs=inputs, system="epsilon_k")
    return _evaluate_epsilon_k(
        wilsons,
        epsilon_k_np_budget_override=epsilon_k_np_budget_override,
    )

# Analogous pass-throughs for delta_mk, bd_mixing, bs_mixing.
```

Three things this achieves:
- **isolation**: constraints depend on the adapter API, not on
  `quarkConstraints/` internals.
- **single point of update**: if `quarkConstraints/deltaf2.py` changes
  a signature, only the adapter changes — 20 constraint files do not.
- **dependency direction**: `quarkConstraints/` does not import
  `flavor_catalog_constraints/` (would be circular). Confirmed.

**Append-only convention (R1 M-6).** Adapter files are **append-only**
in the codex loop: a constraint may add a new wrapper function; it
may not modify the signature of an existing one. Modifying an
existing adapter wrapper is a cross-constraint change and requires a
separate sign-off pass (the directive forbids two codex agents from
touching the same adapter signature concurrently). The contributor
README must say this explicitly.

One adapter per shared physics module:
- `deltaf2_adapter.py` → `quarkConstraints/deltaf2.py`
- `modern_phenomenology_adapter.py` → `quarkConstraints/modern/phenomenology.py`
- `meg_adapter.py` → `flavorConstraints/muToEGamma.py`
- `(future) edm_adapter.py`, `collider_adapter.py`, ...

Constraints with no shared physics put their physics in the constraint
module and only refactor into an adapter on the second consumer
(YAGNI).

---

## F. Bulk Evaluation

### Driver API

```python
from flavor_catalog_constraints import (
    ConstraintRegistry, ConstraintLevel, Severity,
)

# evaluate everything (DEFERRED excluded by default — R1 N-4)
results = ConstraintRegistry.evaluate_all(point)
n_fail = sum(1 for r in results.values() if not r.passes)

# evaluate kaon family
kaon = ConstraintRegistry.filter(family="kaon")
kaon_results = {pid: c.evaluate(point) for pid, c in kaon.items()}

# evaluate only PRIMARY HARD constraints
hard = ConstraintRegistry.filter(level=ConstraintLevel.PRIMARY, severity=Severity.HARD)

# also evaluate DEFERRED (diagnostics mode)
all_results = ConstraintRegistry.evaluate_all(point, include_deferred=True)
```

### Exception isolation (R1 M-4)

`evaluate_all` wraps each `evaluate(point)` call in a try/except (see
§C). A constraint that raises is reported as `passes=False` with a
diagnostic note; the other constraints proceed normally. This matches
the per-draw isolation already used by `scripts/run_rs_anarchy.py`
line 541.

### Performance notes (informational, not load-bearing)

- 102 constraints, ~1 ms each → ~100 ms per point. Acceptable for
  scans at <1e5 points/run.
- `point_builder.py` (§E) pre-computes shared intermediates
  (`QuarkMassBasisCouplings`, KK masses, cached Wilson coefficients)
  into `ParameterPoint.extras` so adapters do not recompute them per
  constraint.
- Profile hotspots via an `evaluate_all` wrapper that records
  per-id timing — no scaffold changes needed.

### Aggregation / reporting

```python
# flavor_catalog_constraints/report.py  (convenience only)
def summarise(results: dict[str, ConstraintResult]) -> dict[str, Any]:
    fails = [r for r in results.values() if not r.passes]
    return {
        "n_total": len(results),
        "n_fail": len(fails),
        "fail_ids": sorted(r.process_id for r in fails),
        "max_ratio": max(
            (r.ratio for r in results.values() if r.ratio is not None),
            default=0.0,
        ),
    }
```

### `point_builder.py` — the canonical ParameterPoint factory (R1 CR-2)

```python
# flavor_catalog_constraints/point_builder.py
"""Owner of ParameterPoint construction. The scan driver calls this
once per draw; constraints read `point.extras` only via the
ParameterPointExtras TypedDict keys declared in base.py.
"""
from __future__ import annotations
from typing import Any

from .base import ParameterPoint, ParameterPointExtras
from quarkConstraints.couplings import build_quark_mass_basis_couplings  # if available
from quarkConstraints.deltaf2 import compute_delta_f2_wilsons, DEFAULT_DELTA_F2_INPUTS_V1


def build_parameter_point(scan_row: Any) -> ParameterPoint:
    """Construct a ParameterPoint from a scan-driver row.

    `scan_row` is whatever `scripts/run_rs_anarchy.py` produces today
    (warpConfig.Setup + per-draw anarchic-Yukawa instance). The exact
    fields are owned upstream; this builder is the choke point that
    translates them into the typed extras the constraint zoo consumes.
    """
    couplings = build_quark_mass_basis_couplings(scan_row)  # signature TBD by scan driver
    deltaf2_wilsons = compute_delta_f2_wilsons(
        couplings, inputs=DEFAULT_DELTA_F2_INPUTS_V1, system="epsilon_k",
    )
    extras: ParameterPointExtras = {
        "quark_mass_basis_couplings": couplings,
        "kk_gluon_mass_gev": float(scan_row.m_kk_gluon_gev),
        "deltaf2_wilsons": deltaf2_wilsons,
        # ... future keys go here AND in ParameterPointExtras (one-line edit).
    }
    return ParameterPoint(raw=scan_row, extras=extras)
```

This is the *only* place `ParameterPoint.extras` is populated.
`scripts/run_rs_anarchy.py` swaps its inline
`_build_kk_gluon_couplings()` block for one call to this function;
that is its one new dependency on `flavor_catalog_constraints/`.

---

## G. The Plug-Pull Invariant

### Add a new constraint (one-file change, <5 minutes)

1. `flavor_catalog_constraints/primary/<family>/<NEW_ID>.py`:
   ```python
   """K042 — short physics description.

   Severity: HARD (rationale: PDG-bounded, no SM-vs-exp tension).
   """
   from dataclasses import dataclass
   from flavor_catalog_constraints.base import (
       Severity, ConstraintResult, ParameterPoint,
   )
   from flavor_catalog_constraints.registry import register_constraint
   from flavor_catalog_constraints.anchors import load_raw


   @dataclass(frozen=True)
   class _Anchor:
       br_limit: float
       cl: str
       snapshot_path: str

   def _build_anchor(raw):
       lim = raw["primary_current_limit"]  # the key THIS yaml uses
       return _Anchor(
           br_limit=float(lim["value"]),
           cl=lim["confidence_level"],
           snapshot_path=lim["snapshot_path"],
       )


   @register_constraint
   class Constraint:
       process_id = "K042"
       severity = Severity.HARD
       observable = "BR(K_L -> pi0 nu nu)"

       def __init__(self):
           # family and level are derived from module path by the decorator
           raw = load_raw(self.process_id, family="kaon", tier="primary")
           self.anchor = _build_anchor(raw)
           self.references = (self.anchor.snapshot_path,)

       def evaluate(self, point: ParameterPoint) -> ConstraintResult:
           ...  # physics
           return ConstraintResult(...)
   ```
2. `tests/constraints/primary/kaon/test_K042.py` (mirror).
3. Done. No edits to any other file. Auto-discovery picks the module
   up. The decorator validates the ID and pins `family`/`level` from
   the module path; `test_severity_documented` checks the docstring.

### Remove a constraint (one-file delete)

1. `rm flavor_catalog_constraints/primary/<family>/<ID>.py`
2. `rm tests/constraints/primary/<family>/test_<ID>.py`
3. Done. The `test_bijection_with_yaml` test will fail if the yaml is
   still present — either delete the yaml too or move it (catalog
   deletion is out of scope for the constraint codex loop).

### Promote SECONDARY → PRIMARY

`git mv` of the two files. Because `level` is derived from the module
path, **zero code edits** are required (R1 REC-3). The
`test_path_metadata_consistency` test re-pins the new tier
automatically.

---

## H. Worked Example — K001 (`epsilon_K`)

### `flavor_catalog_constraints/primary/kaon/K001.py` (structure sketch)

```python
"""K001 — Indirect CP violation in K0-K0bar mixing (epsilon_K).

Severity: HARD (PDG-bounded, hadronic-uncertainty-dominated budget but
        the experimental uncertainty itself is ~0.5%).

Catalog sidecar: flavor_catalog/processes/kaon/K001.yaml
Physics core:    quarkConstraints/deltaf2.py
"""
from __future__ import annotations

from dataclasses import dataclass

from flavor_catalog_constraints.base import (
    Severity, ConstraintResult, ParameterPoint,
)
from flavor_catalog_constraints.registry import register_constraint
from flavor_catalog_constraints.anchors import load_raw
from flavor_catalog_constraints.physics_adapters.deltaf2_adapter import (
    evaluate_epsilon_k,
)


@dataclass(frozen=True)
class _K001Anchor:
    eps_K_exp_value: float
    eps_K_exp_uncertainty: float
    eps_K_sm_value: float
    snapshot_path_exp: str
    sha256_exp: str
    snapshot_path_sm: str
    sha256_sm: str


def _build_anchor(raw) -> _K001Anchor:
    """K001 uses canonical_experimental_value (the only yaml that does)."""
    exp = raw["canonical_experimental_value"]
    sm = raw["standard_model_reference"]
    return _K001Anchor(
        eps_K_exp_value=float(exp["value"]),
        eps_K_exp_uncertainty=float(exp["uncertainty"]),
        eps_K_sm_value=float(sm["value"]),
        snapshot_path_exp=exp["snapshot_path"],
        sha256_exp=exp["sha256_of_local_snapshot"],
        snapshot_path_sm=sm["snapshot_path"],
        sha256_sm=sm["sha256_of_local_snapshot"],
    )


@register_constraint
class Constraint:
    process_id = "K001"
    severity = Severity.HARD
    observable = "epsilon_K"

    def __init__(self):
        raw = load_raw(self.process_id, family="kaon", tier="primary")
        self.anchor = _build_anchor(raw)
        self.references = (
            self.anchor.snapshot_path_exp,
            self.anchor.snapshot_path_sm,
            "flavor_catalog/references/K001/flag2024_bk_arxiv2411_04268.txt",
            "flavor_catalog/references/K001/cfw2008_rs_flavor_arxiv0804_1954.txt",
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        couplings = point.extras["quark_mass_basis_couplings"]
        # Adapter returns the existing physics dataclass unchanged (R1 CR-3):
        result = evaluate_epsilon_k(couplings)
        return ConstraintResult(
            process_id=self.process_id,
            passes=result.passes,
            predicted=result.epsilon_k_np,
            sm_prediction=self.anchor.eps_K_sm_value,
            experimental=self.anchor.eps_K_exp_value,
            ratio=result.ratio_to_budget,
            budget=result.epsilon_k_np_budget,
            severity=self.severity,
            notes="epsilon_K^NP from Im(M12^NP) via kappa_epsilon convention.",
            diagnostics={"im_m12_np": result.im_m12_np},
        )
```

Note the unpacking is **clean**: `evaluate_epsilon_k(couplings)`
returns an `EpsilonKResult` dataclass; we read the fields by name.
The R1 CR-3 bug ("`(float, bool)` tuple cannot unpack into 4 names")
is gone because the adapter no longer repackages — it passes the
physics dataclass through.

### `tests/constraints/primary/kaon/test_K001.py` (structure sketch)

```python
import yaml
from pathlib import Path

from flavor_catalog_constraints import ConstraintRegistry, ConstraintLevel, Severity

REPO_ROOT = Path(__file__).resolve().parents[3]
YAML_PATH = REPO_ROOT / "flavor_catalog/processes/kaon/K001.yaml"


def test_registered():
    c = ConstraintRegistry.get("K001")
    assert c.process_id == "K001"
    assert c.family == "kaon"            # derived from module path
    assert c.level == ConstraintLevel.PRIMARY  # derived from module path
    assert c.severity == Severity.HARD


def test_anchor_matches_yaml():
    c = ConstraintRegistry.get("K001")
    raw = yaml.safe_load(open(YAML_PATH))["pdg_or_equivalent"]
    exp = raw["canonical_experimental_value"]
    sm = raw["standard_model_reference"]
    assert c.anchor.eps_K_exp_value == float(exp["value"])
    assert c.anchor.eps_K_exp_uncertainty == float(exp["uncertainty"])
    assert c.anchor.eps_K_sm_value == float(sm["value"])
    assert c.anchor.sha256_exp == exp["sha256_of_local_snapshot"]


def test_sm_point_passes(sm_point):
    r = ConstraintRegistry.get("K001").evaluate(sm_point)
    assert r.passes
    assert r.ratio is not None and r.ratio < 1.0


def test_excluded_point_fails(excluded_point):
    r = ConstraintRegistry.get("K001").evaluate(excluded_point)
    assert not r.passes


def test_evaluate_is_pure(sm_point):
    c = ConstraintRegistry.get("K001")
    r1 = c.evaluate(sm_point)
    r2 = c.evaluate(sm_point)
    assert r1 == r2
```

Total new code for K001: ~50 lines of constraint + ~30 lines of test
(slight bump from R1 due to the inline typed anchor view). The
load-bearing demonstration still holds: the scaffold makes
implementation cheap once the adapter exists.

---

## I. Per-Constraint Workflow (codex orchestration loop)

The directive (Step 2) drives this loop, one constraint at a time:

```
for process_id in priority_ordered_list:
    agent1 (codex):  plan + implement
        -> writes flavor_catalog_constraints/<tier>/<family>/<id>.py
           (defines _Anchor dataclass + _build_anchor function + Constraint class)
        -> writes tests/constraints/<tier>/<family>/test_<id>.py
        -> runs `pytest tests/constraints/<tier>/<family>/test_<id>.py`
        -> commits

    agent2 (codex, physics fact-check, parallel with agent3):
        -> reads the yaml sidecar + the constraint module
        -> verifies formula matches the source literature snapshot
        -> writes findings to .orchestration/reviews/<id>_physics.md

    agent3 (codex, code review, parallel with agent2):
        -> runs ruff, mypy, pytest
        -> checks adapter usage, ConstraintResult fields, isolation
        -> writes findings to .orchestration/reviews/<id>_code.md

    agent1: fix based on findings, repeat until both reviews PASS

    opus reviewer: final sign-off
        -> checks per-constraint test passes
        -> checks ConstraintRegistry contract test still passes
        -> updates .orchestration/progress.json
```

### Why the scaffold makes this loop frictionless

1. **One file per constraint**: codex never has to read 100-file
   diffs.
2. **Decorator-only registration**: codex cannot forget to "wire it
   up"; if the file exists, the registry finds it.
3. **Schema-flex anchor loader**: codex defines the typed view that
   matches *this constraint's* yaml shape — no global schema to
   satisfy.
4. **Path-derived tier/family**: codex cannot mis-tag a constraint;
   the decorator pins the tier and family from the module path. Open
   questions about "which `level=` to declare" go away.
5. **Adapter layer**: codex working on `K001` cannot accidentally
   touch `quarkConstraints/deltaf2.py` — they only touch the
   adapter, and only to *append* a new wrapper function, never to
   modify an existing one (M-6 append-only convention).
6. **Mirrored test tree**: the orchestrator runs exactly one test
   file per loop iteration.
7. **Idempotent re-runs**: `pytest tests/constraints/` always
   re-verifies the full inventory; the global
   `test_registry_contract.py` catches any drift since the previous
   commit.
8. **Lazy discovery in tests**: importing
   `flavor_catalog_constraints` is fast (no eager walk);
   per-constraint pytest only pays the discovery cost via the
   session-scoped autouse fixture, amortised across all tests in the
   run.

### Suspicious-constraint flag mechanism

Per directive: "Suspicious/problematic constraints flagged in a
dedicated document."

Add `.orchestration/PHASE3_SUSPICIOUS.md` with a one-row-per-constraint
table:
```
| process_id | flag | reason | recommendation |
```
A constraint can additionally be marked
`severity = Severity.INFO` (downgrades it to advisory) or — for a
genuinely deferred constraint — promoted to `DEFERRED` status by
moving its module to a `deferred/<family>/` subtree (the decorator
will refuse to register it since `deferred` is not in the tier map,
which is the desired outcome). `evaluate_all(..., include_deferred=False)`
is the default; setting it true reports DEFERRED constraints
diagnostically.

### `Severity.INFO` usage (R1 N-6)

Reserved for advisory secondary observables where a failure should
**not** veto a point (e.g. a SECONDARY LHCb tension that's still
within 3σ). At least one SECONDARY constraint must use INFO at the
end of the codex loop, or N-6 says drop the enum value. Tracked as
an open commitment in the readme.

---

## J. Resolved Open Questions (R1 dispositions)

1. **Package name.** **`flavor_catalog_constraints/`** — confirmed.
   Verbose but accurate; no collision with existing top-level
   packages (`warpConfig`, `neutrinos`, `solvers`, `diagonalization`,
   `yukawa`, `scanParams`, `quarkConstraints`, `flavorConstraints`,
   `qcd`, `flavor_catalog`).

2. **`ParameterPoint.extras` owner.** **New file
   `flavor_catalog_constraints/point_builder.py` owns construction;
   the `ParameterPointExtras` TypedDict in `base.py` is the frozen
   key registry.** `scanParams/` is the wrong home (it's a
   scan-config helper, not a point factory). Adding a new extras key
   is a one-line edit to `base.py`, the only legitimate
   cross-constraint coordination point.

3. **Family slugs vs ID prefixes.** **All 102 IDs are well-formed
   `^[A-Z]+[0-9]+$` (B, C, CR, E, EW, K, L, T).** All `family` values
   are snake_case and match directory names exactly. **No K016
   anomaly exists.** The "+1 extra" in the R1 planner's preamble was
   a double-count.

---

## K. Open Items / Notes for the Orchestrator

These do not block code but warrant flagging:

- **N-7 lint enforcement.** Ruff per-file-ignores or an `import-rule`
  check should forbid constraint files from importing
  `quarkConstraints.*` or `flavorConstraints.*` directly. Only
  `flavor_catalog_constraints/physics_adapters/*` may import shared
  physics modules. Document in `flavor_catalog_constraints/README.md`
  and consider an `isort` first-party rule or a pre-commit hook.
- **M-1 severity provenance.** Severity is not in the yaml. The
  scaffold's `test_severity_documented` enforces a docstring rationale,
  which prevents silent invention but does not pin severity to a
  reviewable source. If the orchestrator wants stronger provenance,
  add a `severity:` field to the yaml schema in a future cleanup wave
  (catalog-touching; defer).
- **CR001-style list anchors.** Some yamls (`values:` list with
  multiple per-experiment entries) require the constraint to pick
  *which* entry it uses. The constraint's `_build_anchor` makes this
  choice explicit; the contract test pins it. If the catalog later
  changes the order or contents of the list, the constraint's test
  will fail loudly — desired behaviour.
- **point_builder.py's signature.** The `scan_row` argument's exact
  type is TBD by `scripts/run_rs_anarchy.py`. The builder lives in
  the constraint package so the scan driver only depends on the
  builder; the constraint zoo only depends on `ParameterPointExtras`.
  This is a soft coupling that the orchestrator should pin once
  during the scan-driver integration step.
- **Future families.** If the catalog gains a new family slug, both
  `_VALID_FAMILIES` in `registry.py` and the
  `tests/constraints/<tier>/<new_family>/` directory must be created.
  The registry refuses to register modules in unknown families
  (loud failure) — desired.
