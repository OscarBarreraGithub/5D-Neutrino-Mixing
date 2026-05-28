"""Global contract test for ``flavor_catalog_constraints``.

Cross-cutting invariants enforced once here, not per-constraint:

- No constraint module fails to import (R1 M-5).
- Yaml inventory <-> registered IDs reconciliation. During Phase 3
  scaffolding only K001 is implemented; the other ~101 yamls are
  marked *pending implementation* so the test still passes today but
  will start to detect missing-implementation drift once the codex
  loop begins landing constraints.
- ``family``/``level`` instance attributes derived from the module
  path match the constraint's filesystem location.
- Process IDs match the canonical regex ``^[A-Z]+[0-9]+$``.
- All extras-keys referenced by constraint modules are declared in
  :class:`ParameterPointExtras` (R1 M-3).
- The yaml sidecar referenced by each registered constraint exists.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path
from typing import Iterable, Set, Tuple

import yaml

from flavor_catalog_constraints import (
    ConstraintLevel,
    ConstraintRegistry,
    ParameterPointExtras,
)
from flavor_catalog_constraints.anchors import CATALOG_ROOT, yaml_path_for

REPO_ROOT = Path(__file__).resolve().parents[2]
PACKAGE_ROOT = REPO_ROOT / "flavor_catalog_constraints"

_ID_RE = re.compile(r"^[A-Z]+[0-9]+$")

# IDs implemented in the Phase 3 scaffold commit. Every other catalogued
# yaml is expected to be *pending* — the contract test reports them but
# does not fail the suite. As the codex loop lands additional
# constraints, they will register themselves and the bijection
# assertions below will tighten automatically.
IMPLEMENTED_IDS_AT_SCAFFOLD_TIME = {
    "K001",
    "K002",
    "K003",
    "K004",
    "K005",
    "K006",
    "K008",
    "K009",
    "K010",
    "K012",
    "K013",
    "K017",
    "K018",
    "E001",
    "E002",
    "E004",
    "E006",
    "E007",
    "E008",
    "E009",
    "L001", "L002", "L003", "L004", "L005", "L006",
    "L007", "L008", "L009", "L010", "L023",
    "C001", "C002", "C003", "C004", "C005", "C006", "C007", "C008",
    "B001", "B002", "B003", "B004", "B005", "B006", "B009", "B011",
    "B012", "B015", "B016", "B017", "B018", "B019", "B021", "B022",
    "B023", "B025", "B026", "B032", "B033", "B034",
    "EW001", "EW002", "EW003",
    "T001", "T002", "T003", "T004", "T005", "T006", "T007", "T008",
    "T010", "T011", "T012", "T015", "T016", "T017", "T018", "T019", "T020",
}


# --------------------------------------------------------------------- #
# Yaml inventory helpers
# --------------------------------------------------------------------- #


def _yaml_inventory() -> Tuple[Set[str], dict]:
    """Walk the catalog tree and return ``(set_of_ids, id -> (family, tier))``."""
    ids: Set[str] = set()
    meta: dict[str, Tuple[str, ConstraintLevel]] = {}
    # PRIMARY tier
    for family_dir in CATALOG_ROOT.iterdir():
        if not family_dir.is_dir():
            continue
        if family_dir.name == "secondary":
            continue
        for ypath in family_dir.glob("*.yaml"):
            pid = ypath.stem
            ids.add(pid)
            meta[pid] = (family_dir.name, ConstraintLevel.PRIMARY)
    # SECONDARY tier
    sec_root = CATALOG_ROOT / "secondary"
    if sec_root.is_dir():
        for family_dir in sec_root.iterdir():
            if not family_dir.is_dir():
                continue
            for ypath in family_dir.glob("*.yaml"):
                pid = ypath.stem
                ids.add(pid)
                meta[pid] = (family_dir.name, ConstraintLevel.SECONDARY)
    return ids, meta


# --------------------------------------------------------------------- #
# Tests
# --------------------------------------------------------------------- #


def test_no_import_failures():
    """Discovery walked every module successfully (R1 M-5)."""
    ConstraintRegistry.discover()
    failures = ConstraintRegistry.import_failures()
    assert failures == {}, (
        "Constraint modules failed to import: "
        + "; ".join(f"{m}: {type(e).__name__}: {e}" for m, e in failures.items())
    )


def test_all_yaml_ids_match_regex():
    """All catalogued yaml stem IDs match ``^[A-Z]+[0-9]+$``."""
    ids, _ = _yaml_inventory()
    bad = sorted(i for i in ids if not _ID_RE.match(i))
    assert bad == [], f"Yaml IDs violating the regex: {bad}"


def test_registered_ids_are_subset_of_yaml_inventory():
    """Every registered constraint has a corresponding yaml sidecar.

    The reverse inclusion (every yaml has an implementation) is
    enforced by :func:`test_implemented_ids_pin` below as a soft
    pending-set check; once all 102 are implemented this assertion
    will become exact.
    """
    registered = set(ConstraintRegistry.all())
    yaml_ids, _ = _yaml_inventory()
    extra = registered - yaml_ids
    assert extra == set(), (
        f"Registered constraint(s) without a yaml sidecar: {sorted(extra)}"
    )


def test_implemented_ids_pin():
    """Pin the implemented-IDs set to the scaffold expectation.

    During scaffolding (this commit) the set is exactly ``{K001}``.
    As the codex loop lands further constraints, contributors update
    :data:`IMPLEMENTED_IDS_AT_SCAFFOLD_TIME` in this file (a one-line
    diff). The test then enforces "no silent registration drift" both
    ways:

    - A constraint module that quietly stops registering itself
      surfaces here as a missing entry.
    - A new module merged without bumping the pinned set surfaces
      here as a stray entry.
    """
    registered = set(ConstraintRegistry.all())
    assert registered == IMPLEMENTED_IDS_AT_SCAFFOLD_TIME, (
        f"Registered={sorted(registered)} "
        f"vs pinned={sorted(IMPLEMENTED_IDS_AT_SCAFFOLD_TIME)}. "
        "If you added a new constraint, also bump "
        "IMPLEMENTED_IDS_AT_SCAFFOLD_TIME in this test file."
    )


def test_pending_ids_report_only():
    """Report (but do not fail on) the ~101 yaml IDs still awaiting code.

    This produces a visible diagnostic in pytest output so reviewers can
    see how many of the catalog's constraints still need an
    implementation, without breaking CI before the codex loop lands
    them.
    """
    yaml_ids, _ = _yaml_inventory()
    pending = sorted(yaml_ids - IMPLEMENTED_IDS_AT_SCAFFOLD_TIME)
    # Soft assertion via print: pytest -v will surface this in the
    # captured-stdout block on success.
    print(
        f"\nPending constraint implementations: {len(pending)} of "
        f"{len(yaml_ids)} catalogued IDs. First 10: {pending[:10]}"
    )
    # Sanity floor: at scaffold time there must be at least one pending
    # implementation (otherwise this test file is stale).
    assert pending, (
        "No pending constraints — has every catalog yaml been "
        "implemented? Update IMPLEMENTED_IDS_AT_SCAFFOLD_TIME to "
        "match the full inventory and delete this assertion."
    )


def test_registered_constraints_have_correct_family_and_level():
    """For each registered constraint, derived ``family``/``level``
    match its filesystem location and its yaml's ``family:`` field.
    """
    _, meta = _yaml_inventory()
    for pid, c in ConstraintRegistry.all().items():
        exp_family, exp_level = meta[pid]
        assert getattr(c, "family", None) == exp_family, (
            f"{pid}: family={c.family!r} (expected {exp_family!r})"
        )
        assert getattr(c, "level", None) == exp_level, (
            f"{pid}: level={c.level!r} (expected {exp_level!r})"
        )

        # Cross-check yaml's family: field.
        ypath = yaml_path_for(pid, family=exp_family, tier=exp_level)
        data = yaml.safe_load(open(ypath))
        yaml_family = data.get("family")
        assert yaml_family == exp_family, (
            f"{pid}: yaml family field={yaml_family!r} vs path family={exp_family!r}"
        )


def test_registered_constraints_have_yaml_sidecar():
    """Every registered constraint's yaml path resolves to a real file (REC-3)."""
    _, meta = _yaml_inventory()
    for pid in ConstraintRegistry.all():
        exp_family, exp_level = meta[pid]
        path = yaml_path_for(pid, family=exp_family, tier=exp_level)
        assert path.is_file(), f"{pid}: yaml sidecar missing at {path}"


# --------------------------------------------------------------------- #
# Static scan: extras keys referenced by constraints must be declared
# --------------------------------------------------------------------- #


def _iter_constraint_module_files() -> Iterable[Path]:
    for tier in ("primary", "secondary"):
        tier_root = PACKAGE_ROOT / tier
        if not tier_root.is_dir():
            continue
        for family_dir in tier_root.iterdir():
            if not family_dir.is_dir():
                continue
            for path in family_dir.glob("*.py"):
                if path.name == "__init__.py":
                    continue
                yield path


def _extract_extras_string_keys(source: str) -> Set[str]:
    """Scan AST for ``point.extras["<key>"]`` and ``extras.get("<key>")``."""
    tree = ast.parse(source)
    found: Set[str] = set()
    for node in ast.walk(tree):
        # point.extras["key"]
        if (
            isinstance(node, ast.Subscript)
            and isinstance(node.value, ast.Attribute)
            and node.value.attr == "extras"
            and isinstance(node.slice, ast.Constant)
            and isinstance(node.slice.value, str)
        ):
            found.add(node.slice.value)
        # .extras.get("key", ...) / .extras.get("key")
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "get"
            and isinstance(node.func.value, ast.Attribute)
            and node.func.value.attr == "extras"
            and node.args
            and isinstance(node.args[0], ast.Constant)
            and isinstance(node.args[0].value, str)
        ):
            found.add(node.args[0].value)
    return found


def test_extras_keys_are_declared():
    """Static-scan check: every ``point.extras["k"]`` literal in a constraint
    module is a key declared in :class:`ParameterPointExtras` (R1 M-3).
    """
    declared = set(ParameterPointExtras.__annotations__.keys())
    offenders: dict[str, Set[str]] = {}
    for path in _iter_constraint_module_files():
        source = path.read_text()
        try:
            referenced = _extract_extras_string_keys(source)
        except SyntaxError as exc:  # pragma: no cover — bad source caught by import test
            raise AssertionError(f"{path}: syntax error {exc}") from exc
        stray = referenced - declared
        if stray:
            offenders[str(path.relative_to(REPO_ROOT))] = stray
    assert not offenders, (
        "Constraint modules reference extras keys not declared in "
        f"ParameterPointExtras: {offenders}. Add the key to "
        "flavor_catalog_constraints/base.py::ParameterPointExtras."
    )
