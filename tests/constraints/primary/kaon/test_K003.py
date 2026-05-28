"""Per-constraint contract tests for K003 (``Re(epsilon'/epsilon)``).

The 5 mandatory tests every constraint test file must include (plan section D):

1. ``test_registered`` - registry lookup yields the right metadata.
2. ``test_anchor_matches_yaml`` - typed view matches the yaml on the
   fields K003 actually reads (flat ``pdg_or_equivalent`` block).
3. ``test_sm_point_passes`` - currently passes because evaluation is
   explicitly deferred pending Delta S = 1 machinery.
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
from flavor_catalog_constraints.physics_adapters.eps_prime import (
    evaluate_eps_prime_placeholder,
)

REPO_ROOT = Path(__file__).resolve().parents[4]
YAML_PATH = REPO_ROOT / "flavor_catalog" / "processes" / "kaon" / "K003.yaml"


def test_registered():
    c = ConstraintRegistry.get("K003")
    assert c.process_id == "K003"
    assert c.observable == "Re(epsilon'/epsilon)"
    assert getattr(c, "family", None) == "kaon"
    assert getattr(c, "level", None) == ConstraintLevel.PRIMARY
    assert c.severity == Severity.INFO


def test_anchor_matches_yaml():
    """The typed anchor pins the same PDG eps'/eps value as the yaml."""
    c = ConstraintRegistry.get("K003")
    raw = yaml.safe_load(open(YAML_PATH))["pdg_or_equivalent"]

    assert c.anchor.year == int(raw["year"])
    assert c.anchor.source == str(raw["source"])
    assert c.anchor.observable == str(raw["observable"])
    assert c.anchor.value == float(raw["value"])
    assert c.anchor.uncertainty == float(raw["uncertainty"])
    assert c.anchor.units == str(raw["units"])
    assert c.anchor.display_value == str(raw["display_value"])
    assert c.anchor.source_url == str(raw["source_url"])
    assert c.anchor.access_date == str(raw["access_date"])
    assert c.anchor.snapshot_path == str(raw["snapshot_path"])
    assert c.anchor.sha256_of_text_snapshot == str(
        raw["sha256_of_text_snapshot"]
    )


def test_sm_point_passes(sm_point):
    """Deferred placeholder is non-vetoing at the SM baseline."""
    r = ConstraintRegistry.get("K003").evaluate(sm_point)
    assert r.process_id == "K003"
    assert r.passes is True
    assert r.predicted is None
    assert r.experimental == 0.00166
    assert r.severity == Severity.INFO
    assert r.diagnostics["deferred"] is True
    assert "Delta S=1" in r.diagnostics["evaluation_deferred"]


def test_excluded_point_fails(excluded_point):
    """Even a K001/K002-excluded point passes because K003 evaluation is deferred."""
    r = ConstraintRegistry.get("K003").evaluate(excluded_point)
    assert r.passes is True
    assert r.predicted is None
    assert r.ratio is None
    assert r.budget is None
    assert r.diagnostics["deferred"] is True


def test_evaluate_is_pure(sm_point):
    """Two consecutive evaluations on the same point must yield identical results."""
    c = ConstraintRegistry.get("K003")
    r1 = c.evaluate(sm_point)
    r2 = c.evaluate(sm_point)
    assert r1 == r2


def test_adapter_result_is_marked_deferred():
    """Adapter exposes the explicit deferred bit requested for K003."""
    result = evaluate_eps_prime_placeholder()
    assert result.deferred is True
    assert result.passes is True
    assert result.predicted is None
