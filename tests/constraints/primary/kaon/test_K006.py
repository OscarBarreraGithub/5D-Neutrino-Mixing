"""Per-constraint contract tests for K006 (``K_L -> mu+ mu-``).

The 5 mandatory tests every constraint test file must include (plan section D):

1. ``test_registered`` - registry lookup yields the right metadata.
2. ``test_anchor_matches_yaml`` - typed view matches the yaml on the
   fields K006 actually reads (flat ``pdg_or_equivalent`` block).
3. ``test_sm_point_passes`` - currently passes because evaluation is
   explicitly deferred pending Delta S = 1 rare-kaon dimuon machinery.
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
YAML_PATH = REPO_ROOT / "flavor_catalog" / "processes" / "kaon" / "K006.yaml"


def test_registered():
    c = ConstraintRegistry.get("K006")
    assert c.process_id == "K006"
    assert c.observable == "BR(K_L -> mu+ mu-)"
    assert getattr(c, "family", None) == "kaon"
    assert getattr(c, "level", None) == ConstraintLevel.PRIMARY
    assert c.severity == Severity.INFO


def test_anchor_matches_yaml():
    """The typed anchor pins the PDG K_L -> mu+ mu- branching fraction."""
    c = ConstraintRegistry.get("K006")
    raw = yaml.safe_load(open(YAML_PATH))["pdg_or_equivalent"]

    assert c.anchor.year == int(raw["year"])
    assert c.anchor.source == str(raw["source"])
    assert c.anchor.observable == str(raw["observable"])
    assert c.anchor.value == float(raw["value"])
    assert c.anchor.uncertainty == float(raw["uncertainty"])
    assert c.anchor.units == str(raw["units"])
    assert c.anchor.confidence_level == raw["confidence_level"]
    assert c.anchor.source_url == str(raw["source_url"])
    assert c.anchor.access_date == str(raw["access_date"])
    assert c.anchor.snapshot_path == str(raw["snapshot_path"])
    assert c.anchor.sha256_of_text_snapshot == str(
        raw["sha256_of_text_snapshot"]
    )


def test_sm_point_passes(sm_point):
    """Deferred placeholder is non-vetoing at an SM-like baseline."""
    r = ConstraintRegistry.get("K006").evaluate(sm_point)
    assert r.process_id == "K006"
    assert r.passes is True
    assert r.predicted is None
    assert r.experimental == 6.84e-9
    assert r.severity == Severity.INFO
    assert r.diagnostics["deferred"] is True
    assert r.diagnostics["short_distance_sm_reference"] == 9.0e-10
    assert r.diagnostics["long_distance_dominated"] is True
    assert "Delta S=1" in r.diagnostics["evaluation_deferred"]
    assert "two-photon" in r.diagnostics["evaluation_deferred"]


def test_excluded_point_fails(excluded_point):
    """Even a K001/K002-excluded point passes because K006 evaluation is deferred."""
    r = ConstraintRegistry.get("K006").evaluate(excluded_point)
    assert r.passes is True
    assert r.predicted is None
    assert r.ratio is None
    assert r.budget is None
    assert r.diagnostics["deferred"] is True


def test_evaluate_is_pure(sm_point):
    """Two consecutive evaluations on the same point must yield identical results."""
    c = ConstraintRegistry.get("K006")
    r1 = c.evaluate(sm_point)
    r2 = c.evaluate(sm_point)
    assert r1 == r2
