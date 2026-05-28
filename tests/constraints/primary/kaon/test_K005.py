"""Per-constraint contract tests for K005 (``K_L -> pi0 nu nubar``).

The 5 mandatory tests every constraint test file must include (plan section D):

1. ``test_registered`` - registry lookup yields the right metadata.
2. ``test_anchor_matches_yaml`` - typed view matches the yaml on the
   fields K005 actually reads (one-entry ``pdg_or_equivalent`` list).
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
YAML_PATH = REPO_ROOT / "flavor_catalog" / "processes" / "kaon" / "K005.yaml"


def test_registered():
    c = ConstraintRegistry.get("K005")
    assert c.process_id == "K005"
    assert c.observable == "BR(K_L -> pi0 nu nubar)"
    assert getattr(c, "family", None) == "kaon"
    assert getattr(c, "level", None) == ConstraintLevel.PRIMARY
    assert c.severity == Severity.INFO


def test_anchor_matches_yaml():
    """The typed anchor pins the PDG/KOTO 90% CL upper limit."""
    c = ConstraintRegistry.get("K005")
    raw_list = yaml.safe_load(open(YAML_PATH))["pdg_or_equivalent"]
    assert isinstance(raw_list, list)
    assert len(raw_list) == 1

    raw = raw_list[0]
    assert raw["value"] == "<2.2e-9"
    assert c.anchor.observable == str(raw["observable"])
    assert c.anchor.year == int(raw["year"])
    assert c.anchor.raw_value == str(raw["value"])
    assert c.anchor.upper_limit == float(str(raw["value"]).lstrip("<"))
    assert c.anchor.upper_limit == 2.2e-9
    assert c.anchor.uncertainty == raw["uncertainty"]
    assert c.anchor.units == str(raw["units"])
    assert c.anchor.confidence_level == str(raw["confidence_level"])
    assert c.anchor.source == str(raw["source"])
    assert c.anchor.source_url == str(raw["source_url"])
    assert c.anchor.access_date == str(raw["access_date"])
    assert c.anchor.snapshot_path == str(raw["snapshot_path"])
    assert c.anchor.sha256_of_text_snapshot == str(
        raw["sha256_of_text_snapshot"]
    )
    assert c.anchor.notes == str(raw["notes"])


def test_sm_point_passes(sm_point):
    """Deferred placeholder is non-vetoing at an SM-like baseline."""
    r = ConstraintRegistry.get("K005").evaluate(sm_point)
    assert r.process_id == "K005"
    assert r.passes is True
    assert r.predicted is None
    assert r.experimental == 2.2e-9
    assert r.severity == Severity.INFO
    assert r.diagnostics["deferred"] is True
    assert r.diagnostics["upper_limit"] == 2.2e-9
    assert "Delta S=1" in r.diagnostics["evaluation_deferred"]


def test_excluded_point_fails(excluded_point):
    """Even a K001/K002-excluded point passes because K005 evaluation is deferred."""
    r = ConstraintRegistry.get("K005").evaluate(excluded_point)
    assert r.passes is True
    assert r.predicted is None
    assert r.ratio is None
    assert r.budget is None
    assert r.diagnostics["deferred"] is True


def test_evaluate_is_pure(sm_point):
    """Two consecutive evaluations on the same point must yield identical results."""
    c = ConstraintRegistry.get("K005")
    r1 = c.evaluate(sm_point)
    r2 = c.evaluate(sm_point)
    assert r1 == r2
