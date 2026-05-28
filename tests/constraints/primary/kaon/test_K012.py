"""Per-constraint contract tests for K012 (``K_S -> mu+ mu-``).

The 5 mandatory tests every constraint test file must include (plan section D):

1. ``test_registered`` - registry lookup yields the right metadata.
2. ``test_anchor_matches_yaml`` - typed view matches the yaml on the
   fields K012 actually reads (``pdg_or_equivalent`` list records).
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
from flavor_catalog_constraints.physics_adapters.rare_kaon_decays import (
    evaluate_kshort_mumu_placeholder,
)

REPO_ROOT = Path(__file__).resolve().parents[4]
YAML_PATH = REPO_ROOT / "flavor_catalog" / "processes" / "kaon" / "K012.yaml"


def _by_value_id(raw, value_id):
    return next(record for record in raw if record["value_id"] == value_id)


def test_registered():
    c = ConstraintRegistry.get("K012")
    assert c.process_id == "K012"
    assert c.observable == "BR(K_S -> mu+ mu-)"
    assert getattr(c, "family", None) == "kaon"
    assert getattr(c, "level", None) == ConstraintLevel.PRIMARY
    assert c.severity == Severity.INFO


def test_anchor_matches_yaml():
    """The typed anchor pins the PDG/LHCb K_S -> mu+ mu- upper limit."""
    c = ConstraintRegistry.get("K012")
    raw = yaml.safe_load(open(YAML_PATH))["pdg_or_equivalent"]
    headline = _by_value_id(raw, "PDG2025:K012:headline_limit")
    sm_estimate = _by_value_id(raw, "Chobanova2018:K012:sm_estimate")

    assert headline["value"] == "<2.1e-10"
    assert sm_estimate["value"] == "approximately 5e-12"
    assert c.anchor.year == int(headline["year"])
    assert c.anchor.source == str(headline["source"])
    assert c.anchor.observable == str(headline["observable"])
    assert c.anchor.raw_value == str(headline["value"])
    assert c.anchor.upper_limit == float(str(headline["value"]).lstrip("<"))
    assert c.anchor.upper_limit == 2.1e-10
    assert c.anchor.uncertainty == headline["uncertainty"]
    assert c.anchor.units == str(headline["units"])
    assert c.anchor.confidence_level == headline["confidence_level"]
    assert c.anchor.source_url == str(headline["source_url"])
    assert c.anchor.access_date == str(headline["access_date"])
    assert c.anchor.snapshot_path == str(headline["snapshot_path"])
    assert c.anchor.sha256_of_text_snapshot == str(headline["sha256"])
    assert c.anchor.sm_observable == str(sm_estimate["observable"])
    assert c.anchor.sm_raw_value == str(sm_estimate["value"])
    assert c.anchor.catalog_sm_estimate == 5.0e-12
    assert c.anchor.sm_source == str(sm_estimate["source"])
    assert c.anchor.sm_source_url == str(sm_estimate["source_url"])


def test_sm_point_passes(sm_point):
    """Deferred placeholder is non-vetoing at an SM-like baseline."""
    r = ConstraintRegistry.get("K012").evaluate(sm_point)
    assert r.process_id == "K012"
    assert r.passes is True
    assert r.predicted is None
    assert r.experimental == 2.1e-10
    assert r.severity == Severity.INFO
    assert r.diagnostics["deferred"] is True
    assert r.diagnostics["upper_limit"] == 2.1e-10
    assert r.diagnostics["sm_branching_ratio_reference"] == 5.1e-12
    assert r.diagnostics["catalog_sm_estimate"] == 5.0e-12
    assert r.diagnostics["long_distance_dominated"] is True
    assert "Delta S=1" in r.diagnostics["evaluation_deferred"]
    assert "s -> d mu+ mu-" in r.diagnostics["evaluation_deferred"]
    assert "K_S" in r.diagnostics["evaluation_deferred"]
    assert "two-photon" in r.diagnostics["evaluation_deferred"]


def test_excluded_point_fails(excluded_point):
    """Even a K001/K002-excluded point passes because K012 evaluation is deferred."""
    r = ConstraintRegistry.get("K012").evaluate(excluded_point)
    assert r.passes is True
    assert r.predicted is None
    assert r.ratio is None
    assert r.budget is None
    assert r.diagnostics["deferred"] is True


def test_evaluate_is_pure(sm_point):
    """Two consecutive evaluations on the same point must yield identical results."""
    c = ConstraintRegistry.get("K012")
    r1 = c.evaluate(sm_point)
    r2 = c.evaluate(sm_point)
    assert r1 == r2


def test_adapter_result_is_marked_deferred():
    """Adapter exposes the explicit deferred bit requested for K012."""
    result = evaluate_kshort_mumu_placeholder()
    assert result.deferred is True
    assert result.passes is True
    assert result.predicted is None
    assert "s -> d mu+ mu-" in result.reason
    assert "K_S" in result.reason
