"""Per-constraint contract tests for K018 (``K_l3 |V_us| f_+(0)``).

The 5 mandatory tests every constraint test file must include (plan section D):

1. ``test_registered`` - registry lookup yields the right metadata.
2. ``test_anchor_matches_yaml`` - typed view matches the yaml on the
   fields K018 actually reads (the PDG average record in the list-valued
   ``pdg_or_equivalent`` block).
3. ``test_sm_point_passes`` - currently passes because evaluation is
   explicitly deferred pending source-level K_l3 and charged-current
   BSM machinery.
4. ``test_excluded_point_fails`` - currently also passes because the
   placeholder is non-vetoing until the correct physics exists.
5. ``test_evaluate_is_pure`` - same input -> same output.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from flavor_catalog_constraints import (
    ConstraintLevel,
    ConstraintRegistry,
    Severity,
)
from flavor_catalog_constraints.physics_adapters.rare_kaon_decays import (
    evaluate_kl3_vus_placeholder,
)

REPO_ROOT = Path(__file__).resolve().parents[4]
YAML_PATH = REPO_ROOT / "flavor_catalog" / "processes" / "kaon" / "K018.yaml"
ANCHOR_VALUE_ID = "PDG2025:K018:average_fplus_vus"


def _anchor_record() -> dict[str, object]:
    raw = yaml.safe_load(open(YAML_PATH))["pdg_or_equivalent"]
    for record in raw:
        if record["value_id"] == ANCHOR_VALUE_ID:
            return record
    raise AssertionError(f"missing K018 anchor record {ANCHOR_VALUE_ID!r}")


def test_registered():
    c = ConstraintRegistry.get("K018")
    assert c.process_id == "K018"
    assert c.observable == "K_l3 |V_us| f_+(0)"
    assert getattr(c, "family", None) == "kaon"
    assert getattr(c, "level", None) == ConstraintLevel.PRIMARY
    assert c.severity == Severity.INFO


def test_anchor_matches_yaml():
    """The typed anchor pins the PDG K_l3 ``|V_us| f_+(0)`` average."""
    c = ConstraintRegistry.get("K018")
    raw = _anchor_record()

    assert c.anchor.value_id == str(raw["value_id"])
    assert c.anchor.observable == str(raw["observable"])
    assert c.anchor.year == int(raw["year"])
    assert c.anchor.value == float(raw["value"])
    assert c.anchor.value == pytest.approx(0.21656)
    assert c.anchor.uncertainty == float(raw["uncertainty"])
    assert c.anchor.uncertainty == pytest.approx(0.00035)
    assert c.anchor.units == str(raw["units"])
    assert c.anchor.source == str(raw["source"])
    assert c.anchor.source_url == str(raw["source_url"])
    assert c.anchor.access_date == str(raw["access_date"])
    assert c.anchor.snapshot_path == str(raw["snapshot_path"])
    assert c.anchor.sha256_of_local_snapshot == str(raw["sha256"])
    assert c.anchor.value_summary == str(raw["value_summary"])


def test_sm_point_passes(sm_point):
    """Deferred placeholder is non-vetoing at an SM-like baseline."""
    r = ConstraintRegistry.get("K018").evaluate(sm_point)
    assert r.process_id == "K018"
    assert r.passes is True
    assert r.predicted is None
    assert r.experimental == pytest.approx(0.21656)
    assert r.sm_prediction is None
    assert r.severity == Severity.INFO
    assert r.diagnostics["deferred"] is True
    assert r.diagnostics["tree_level_sm_observable"] is True
    assert r.diagnostics["charged_current_semileptonic"] is True
    assert r.diagnostics["not_delta_s1_rare_decay"] is True
    assert r.diagnostics["radiative_corrections_required"] is True
    assert r.diagnostics["value_id"] == ANCHOR_VALUE_ID
    assert "requires |V_us| extraction" in r.diagnostics["evaluation_deferred"]
    assert "5D RS effects on W-quark vertex" in (
        r.diagnostics["evaluation_deferred"]
    )
    assert "radiative corrections" in r.diagnostics["evaluation_deferred"]
    assert "four-fermion contact operators" in r.diagnostics["np_entry"]
    assert "not a Delta S=1 rare kaon decay" in r.notes


def test_excluded_point_fails(excluded_point):
    """Even a K001/K002-excluded point passes because K018 evaluation is deferred."""
    r = ConstraintRegistry.get("K018").evaluate(excluded_point)
    assert r.passes is True
    assert r.predicted is None
    assert r.ratio is None
    assert r.budget is None
    assert r.diagnostics["deferred"] is True


def test_evaluate_is_pure(sm_point):
    """Two consecutive evaluations on the same point must yield identical results."""
    c = ConstraintRegistry.get("K018")
    r1 = c.evaluate(sm_point)
    r2 = c.evaluate(sm_point)
    assert r1 == r2


def test_adapter_result_is_marked_deferred():
    """Adapter exposes the explicit deferred bit requested for K018."""
    result = evaluate_kl3_vus_placeholder()
    assert result.deferred is True
    assert result.passes is True
    assert result.predicted is None
    assert result.reason == (
        "requires |V_us| extraction + 5D RS effects on W-quark vertex + "
        "radiative corrections"
    )
