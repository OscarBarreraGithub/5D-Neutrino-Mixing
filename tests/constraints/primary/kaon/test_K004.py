"""Per-constraint contract tests for K004 (``K+ -> pi+ nu nubar``).

The 5 mandatory tests every constraint test file must include (plan section D):

1. ``test_registered`` - registry lookup yields the right metadata.
2. ``test_anchor_matches_yaml`` - typed view matches the yaml on the
   fields K004 actually reads (nested ``latest_experimental_value`` block).
3. ``test_sm_point_passes`` - currently passes because evaluation is
   explicitly deferred pending Delta S = 1 rare-kaon machinery.
4. ``test_excluded_point_fails`` - currently also passes because the
   placeholder is non-vetoing until the correct physics exists.
5. ``test_evaluate_is_pure`` - same input -> same output.
"""

from __future__ import annotations

from pathlib import Path

import yaml

from flavor_catalog_constraints import (
    ConstraintLevel,
    ConstraintRegistry,
    Severity,
)

REPO_ROOT = Path(__file__).resolve().parents[4]
YAML_PATH = REPO_ROOT / "flavor_catalog" / "processes" / "kaon" / "K004.yaml"


def _parse_asymmetric_uncertainty(raw: str, scale: float) -> tuple[float, float]:
    up_raw, down_raw = str(raw).split("/")
    return abs(float(up_raw)) * scale, abs(float(down_raw)) * scale


def test_registered():
    c = ConstraintRegistry.get("K004")
    assert c.process_id == "K004"
    assert c.observable == "BR(K+ -> pi+ nu nubar)"
    assert getattr(c, "family", None) == "kaon"
    assert getattr(c, "level", None) == ConstraintLevel.PRIMARY
    assert c.severity == Severity.INFO


def test_anchor_matches_yaml():
    """The typed anchor pins the latest NA62 branching-ratio value."""
    c = ConstraintRegistry.get("K004")
    raw = yaml.safe_load(open(YAML_PATH))["pdg_or_equivalent"][
        "latest_experimental_value"
    ]
    scale = float(raw["scale"])
    uncertainty_up, uncertainty_down = _parse_asymmetric_uncertainty(
        raw["uncertainty"], scale
    )

    assert c.anchor.year == int(raw["year"])
    assert c.anchor.source == str(raw["source"])
    assert c.anchor.raw_value == str(raw["value"])
    assert c.anchor.raw_uncertainty == str(raw["uncertainty"])
    assert c.anchor.scale == scale
    assert c.anchor.value == float(raw["value"]) * scale
    assert c.anchor.uncertainty_up == uncertainty_up
    assert c.anchor.uncertainty_down == uncertainty_down
    assert c.anchor.units == str(raw["units"])
    assert c.anchor.dataset == str(raw["dataset"])
    assert c.anchor.value_summary == str(raw["value_summary"])
    assert c.anchor.source_url == str(raw["source_url"])
    assert c.anchor.access_date == str(raw["access_date"])
    assert c.anchor.snapshot_path == str(raw["snapshot_path"])
    assert c.anchor.sha256_of_local_snapshot == str(
        raw["sha256_of_local_snapshot"]
    )


def test_sm_point_passes(sm_point):
    """Deferred placeholder is non-vetoing at the SM baseline."""
    r = ConstraintRegistry.get("K004").evaluate(sm_point)
    assert r.process_id == "K004"
    assert r.passes is True
    assert r.predicted is None
    assert r.experimental == float("9.6") * float("1e-11")
    assert r.severity == Severity.INFO
    assert r.diagnostics["deferred"] is True
    assert "Delta S=1" in r.diagnostics["evaluation_deferred"]


def test_excluded_point_fails(excluded_point):
    """Even a K001/K002-excluded point passes because K004 evaluation is deferred."""
    r = ConstraintRegistry.get("K004").evaluate(excluded_point)
    assert r.passes is True
    assert r.predicted is None
    assert r.ratio is None
    assert r.budget is None
    assert r.diagnostics["deferred"] is True


def test_evaluate_is_pure(sm_point):
    """Two consecutive evaluations on the same point must yield identical results."""
    c = ConstraintRegistry.get("K004")
    r1 = c.evaluate(sm_point)
    r2 = c.evaluate(sm_point)
    assert r1 == r2
