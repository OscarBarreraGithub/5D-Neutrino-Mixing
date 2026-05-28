"""Per-constraint contract tests for K010 (``K_S -> pi0 e+ e-``).

The 5 mandatory tests every constraint test file must include (plan section D):

1. ``test_registered`` - registry lookup yields the right metadata.
2. ``test_anchor_matches_yaml`` - typed view matches the yaml on the
   fields K010 actually reads (nested ``pdg_or_equivalent`` measurement block).
3. ``test_sm_point_passes`` - currently passes because evaluation is
   explicitly deferred pending Delta S = 1 rare-kaon electron-mode machinery.
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
    evaluate_kshort_pi0ee_placeholder,
)

REPO_ROOT = Path(__file__).resolve().parents[4]
YAML_PATH = REPO_ROOT / "flavor_catalog" / "processes" / "kaon" / "K010.yaml"


def test_registered():
    c = ConstraintRegistry.get("K010")
    assert c.process_id == "K010"
    assert c.observable == "BR(K_S -> pi0 e+ e-)"
    assert getattr(c, "family", None) == "kaon"
    assert getattr(c, "level", None) == ConstraintLevel.PRIMARY
    assert c.severity == Severity.INFO


def test_anchor_matches_yaml():
    """The typed anchor pins the PDG K_S -> pi0 e+ e- branching ratio."""
    c = ConstraintRegistry.get("K010")
    raw = yaml.safe_load(open(YAML_PATH))["pdg_or_equivalent"]
    canonical = raw["canonical_measurement"]
    extrapolated = raw["extrapolated_total_rate"]

    assert extrapolated["value"] == 5.8e-9
    assert canonical["value"] == 3.0e-9
    assert c.anchor.year == int(extrapolated["year"])
    assert c.anchor.source == str(raw["source"])
    assert c.anchor.observable == str(extrapolated["observable"])
    assert c.anchor.branching_ratio == float(extrapolated["value"])
    assert c.anchor.branching_ratio == 5.8e-9
    assert c.anchor.uncertainty == extrapolated["uncertainty"]
    assert c.anchor.units == str(extrapolated["units"])
    assert c.anchor.assumptions == str(extrapolated["assumptions"])
    assert c.anchor.confidence_level == canonical["confidence_level"]
    assert c.anchor.source_url == str(extrapolated["source_url"])
    assert c.anchor.access_date == str(extrapolated["access_date"])
    assert c.anchor.snapshot_path == str(extrapolated["snapshot_path"])
    assert c.anchor.sha256_of_text_snapshot == str(
        extrapolated["sha256_of_text_snapshot"]
    )
    assert c.anchor.canonical_observable == str(canonical["observable"])
    assert c.anchor.canonical_branching_ratio == float(canonical["value"])
    assert c.anchor.canonical_branching_ratio == 3.0e-9
    assert c.anchor.canonical_uncertainty == canonical["uncertainty"]


def test_sm_point_passes(sm_point):
    """Deferred placeholder is non-vetoing at an SM-like baseline."""
    r = ConstraintRegistry.get("K010").evaluate(sm_point)
    assert r.process_id == "K010"
    assert r.passes is True
    assert r.predicted is None
    assert r.experimental == 5.8e-9
    assert r.severity == Severity.INFO
    assert r.diagnostics["deferred"] is True
    assert r.diagnostics["branching_ratio"] == 5.8e-9
    assert r.diagnostics["canonical_partial_region_branching_ratio"] == 3.0e-9
    assert r.diagnostics["sm_branching_ratio_order"] == 5.8e-9
    assert (
        r.diagnostics["component_summary"]
        == "short-distance + long-distance form-factor + "
        "phase-space extrapolation"
    )
    assert "Delta S=1" in r.diagnostics["evaluation_deferred"]
    assert "s -> d e+ e-" in r.diagnostics["evaluation_deferred"]
    assert "form-factor" in r.diagnostics["evaluation_deferred"]


def test_excluded_point_fails(excluded_point):
    """Even a K001/K002-excluded point passes because K010 evaluation is deferred."""
    r = ConstraintRegistry.get("K010").evaluate(excluded_point)
    assert r.passes is True
    assert r.predicted is None
    assert r.ratio is None
    assert r.budget is None
    assert r.diagnostics["deferred"] is True


def test_evaluate_is_pure(sm_point):
    """Two consecutive evaluations on the same point must yield identical results."""
    c = ConstraintRegistry.get("K010")
    r1 = c.evaluate(sm_point)
    r2 = c.evaluate(sm_point)
    assert r1 == r2


def test_adapter_result_is_marked_deferred():
    """Adapter exposes the explicit deferred bit requested for K010."""
    result = evaluate_kshort_pi0ee_placeholder()
    assert result.deferred is True
    assert result.passes is True
    assert result.predicted is None
    assert "s -> d e+ e-" in result.reason
