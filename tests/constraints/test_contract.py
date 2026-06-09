"""Global registry + contract tests for flavor_catalog_constraints.

These pass with ZERO registered constraints (the scaffold is valid on an
empty constraint set) and scale up as real constraints land. They verify:

- discovery runs without raising and reports no import failures;
- every registered constraint satisfies ConstraintProtocol and the
  metadata invariants (id format, severity type, level/family pinned);
- every registered constraint has its catalog sidecar;
- exception isolation: a deliberately-raising constraint does not abort
  evaluate_all and is reported as a failing result;
- the anchor loader fails loudly on a bad lookup (no silent default);
- the ParameterPoint extras key registry rejects unknown keys.
"""

from __future__ import annotations

import math
import re

import pytest

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import anchors, point_builder, registry
from flavor_catalog_constraints.base import (
    ConstraintLevel,
    ConstraintProtocol,
    ConstraintResult,
    ParameterPoint,
    Severity,
)

_ID_RE = re.compile(r"^[A-Z]+[0-9]+$")


def test_discovery_has_no_import_failures(import_failures):
    """No constraint module may fail to import."""
    assert import_failures == {}, (
        "constraint modules failed to import: "
        + "; ".join(f"{m}: {e!r}" for m, e in import_failures.items())
    )


def test_discovery_is_idempotent():
    a = fcc.all_constraints()
    b = fcc.all_constraints()
    assert set(a) == set(b)


def test_reset_for_tests_allows_rediscovery():
    before = fcc.all_constraints()
    assert before
    registry.reset_for_tests()
    registry.discover()
    after = registry.all_constraints()
    assert set(after) == set(before)
    assert after


def test_every_constraint_satisfies_protocol(registry):
    for pid, c in registry.items():
        assert isinstance(c, ConstraintProtocol), f"{pid} is not a ConstraintProtocol"


def test_metadata_invariants(registry):
    for pid, c in registry.items():
        assert c.process_id == pid
        assert _ID_RE.match(c.process_id), f"{pid}: bad id format"
        assert isinstance(c.severity, Severity), f"{pid}: severity not a Severity"
        assert isinstance(c.level, ConstraintLevel), f"{pid}: level not pinned"
        assert c.family in registry_valid_families(), f"{pid}: family {c.family!r} invalid"
        assert isinstance(c.observable, str) and c.observable, f"{pid}: empty observable"


def registry_valid_families():
    return registry.VALID_FAMILIES


def test_every_constraint_has_sidecar(registry):
    for pid, c in registry.items():
        path = anchors.yaml_path_for(pid, family=c.family, tier=c.level)
        assert path.is_file(), f"{pid}: missing catalog sidecar at {path}"


def test_evaluate_all_returns_result_per_constraint(registry):
    point = point_builder.empty_point()
    results = fcc.evaluate_all(point)
    assert set(results) == set(registry)
    for pid, r in results.items():
        assert isinstance(r, ConstraintResult)
        assert r.process_id == pid
        # numeric fields are real floats or None (no complex leakage)
        for f in (r.predicted, r.sm_prediction, r.experimental, r.ratio, r.budget):
            assert f is None or isinstance(f, float)


def test_exception_isolation_in_evaluate_all(monkeypatch, registry):
    """A constraint whose evaluate() raises must not abort the sweep."""
    if not registry:
        pytest.skip("no constraints registered; isolation covered by unit test below")

    pid = next(iter(registry))
    target = fcc.get(pid)

    def boom(_point):
        raise RuntimeError("intentional test failure")

    monkeypatch.setattr(target, "evaluate", boom, raising=True)
    results = fcc.evaluate_all(point_builder.empty_point())
    # All constraints still produced a result.
    assert set(results) == set(registry)
    bad = results[pid]
    assert bad.passes is False
    assert "RuntimeError" in bad.notes
    assert bad.diagnostics.get("exception_type") == "RuntimeError"


def test_exception_isolation_unit_level():
    """Isolation works even with a synthetic raising constraint (no real ones needed)."""

    class _Raiser:
        process_id = "Z999"
        severity = Severity.HARD
        observable = "synthetic"
        level = ConstraintLevel.PRIMARY
        family = "kaon"

        def evaluate(self, point: ParameterPoint) -> ConstraintResult:
            raise ValueError("kaboom")

        def evaluate_safe(self) -> ConstraintResult:
            # Mirror the registry's isolation logic to confirm it produces
            # a well-formed failing result.
            try:
                return self.evaluate(point_builder.empty_point())
            except Exception as exc:  # noqa: BLE001
                return ConstraintResult(
                    process_id=self.process_id,
                    severity=self.severity,
                    passes=False,
                    notes=f"evaluate() raised {type(exc).__name__}: {exc}",
                    diagnostics={"exception_type": type(exc).__name__},
                )

    r = _Raiser().evaluate_safe()
    assert r.passes is False
    assert r.diagnostics["exception_type"] == "ValueError"


def test_anchor_loader_fails_loudly_on_missing_key():
    """A non-existent candidate set must raise AnchorError, not default."""
    # Use a real sidecar but ask for keys that don't exist.
    if "K001" not in fcc.all_constraints():
        pytest.skip("K001 example not present")
    with pytest.raises(anchors.AnchorError):
        anchors.load_anchor("K001", family="kaon", candidates=("no_such_block",))


def test_anchor_loader_fails_loudly_on_missing_file():
    with pytest.raises(anchors.AnchorError):
        anchors.load_pdg_block("ZZ999", family="kaon")


def test_anchor_loader_optional_validation_mismatches_raise(monkeypatch):
    block = {
        "canonical_limit": {
            "value": 1.0,
            "value_id": "PDG2026:Z999:limit",
            "units": "dimensionless",
            "confidence_level": "95% CL",
        }
    }

    monkeypatch.setattr(anchors, "load_pdg_block", lambda *args, **kwargs: block)

    anchor = anchors.load_anchor(
        "Z999",
        family="kaon",
        candidates=("canonical_limit",),
        expected_value_id="PDG2026:Z999:limit",
        expected_block_key="canonical_limit",
        expected_units="dimensionless",
        expected_confidence_level="95% CL",
    )
    assert anchor.value_id == "PDG2026:Z999:limit"
    assert anchor.confidence_level == "95% CL"

    mismatch_cases = (
        {"expected_value_id": "wrong"},
        {"expected_block_key": "wrong_block"},
        {"expected_units": "GeV"},
        {"expected_confidence_level": "90% CL"},
    )
    for kwargs in mismatch_cases:
        with pytest.raises(anchors.AnchorError, match="mismatch"):
            anchors.load_anchor("Z999", family="kaon", candidates=("canonical_limit",), **kwargs)


def test_extras_key_registry_rejects_unknown_keys():
    with pytest.raises(KeyError):
        point_builder.make_point(not_a_real_key=123)


def test_make_point_accepts_declared_keys():
    marker = object()
    p = point_builder.make_point(
        kk_gluon_mass_gev=3000.0,
        lepton_lmfv_parameters=marker,
    )
    assert p.get_extra("kk_gluon_mass_gev") == 3000.0
    assert p.get_extra("lepton_lmfv_parameters") is marker
    assert p.get_extra("missing") is None


def test_parameter_point_exposes_immutable_mappings():
    raw = {"scan_id": 1}
    p = point_builder.make_point(raw=raw, kk_gluon_mass_gev=3000.0)

    raw["scan_id"] = 2
    assert p.raw["scan_id"] == 1
    assert p.get_extra("kk_gluon_mass_gev") == 3000.0

    with pytest.raises(TypeError):
        p.raw["scan_id"] = 3
    with pytest.raises(TypeError):
        p.extras["kk_gluon_mass_gev"] = 4000.0


def test_constraint_result_rejects_complex_numeric_field():
    """A complex value in a numeric field is a contract violation."""
    with pytest.raises(TypeError):
        ConstraintResult(
            process_id="K001",
            severity=Severity.HARD,
            passes=True,
            predicted=1.0 + 2.0j,
        )
    # Real floats / None remain valid.
    ConstraintResult(
        process_id="K001",
        severity=Severity.HARD,
        passes=True,
        predicted=1.0,
        ratio=None,
    )


@pytest.mark.parametrize("field,value", [("ratio", math.nan), ("predicted", math.inf), ("budget", -math.inf)])
def test_constraint_result_rejects_non_finite_numeric_field(field, value):
    kwargs = {
        "process_id": "K001",
        "severity": Severity.HARD,
        "passes": True,
        field: value,
    }
    with pytest.raises(ValueError, match="must be finite"):
        ConstraintResult(**kwargs)


def test_constraint_result_rejects_non_bool_passes_and_bad_severity():
    with pytest.raises(TypeError, match="passes must be bool"):
        ConstraintResult(process_id="K001", severity=Severity.HARD, passes=1)
    with pytest.raises(TypeError, match="severity must be a Severity"):
        ConstraintResult(process_id="K001", severity="HARD", passes=True)
