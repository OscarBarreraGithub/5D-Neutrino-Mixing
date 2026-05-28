"""Per-constraint contract tests for K017 (charged-kaon ``R_K`` LFU).

The 5 mandatory tests every constraint test file must include (plan section D):

1. ``test_registered`` - registry lookup yields the right metadata.
2. ``test_anchor_matches_yaml`` - typed view matches the yaml on the
   fields K017 actually reads (nested ``canonical_ratio`` and
   ``sm_prediction`` blocks).
3. ``test_sm_point_passes`` - currently passes because evaluation is
   explicitly deferred pending charged-current LFU machinery.
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
    evaluate_rk_lfu_placeholder,
)

REPO_ROOT = Path(__file__).resolve().parents[4]
YAML_PATH = REPO_ROOT / "flavor_catalog" / "processes" / "kaon" / "K017.yaml"


def _scaled(record: dict[str, object], key: str) -> float:
    return float(record[key]) * float(record["scale"])


def test_registered():
    c = ConstraintRegistry.get("K017")
    assert c.process_id == "K017"
    assert c.observable == "R_K = Gamma(K+ -> e+ nu) / Gamma(K+ -> mu+ nu)"
    assert getattr(c, "family", None) == "kaon"
    assert getattr(c, "level", None) == ConstraintLevel.PRIMARY
    assert c.severity == Severity.INFO


def test_anchor_matches_yaml():
    """The typed anchor pins the PDG R_K ratio and precise SM prediction."""
    c = ConstraintRegistry.get("K017")
    raw = yaml.safe_load(open(YAML_PATH))["pdg_or_equivalent"]
    canonical = raw["canonical_ratio"]
    sm_prediction = raw["sm_prediction"]

    assert c.anchor.year == int(canonical["year"])
    assert c.anchor.source == str(canonical["source"])
    assert c.anchor.edition == str(canonical["edition"])
    assert c.anchor.citation == str(canonical["citation"])
    assert c.anchor.raw_value == str(canonical["value"])
    assert c.anchor.raw_uncertainty == str(canonical["uncertainty"])
    assert c.anchor.scale == float(canonical["scale"])
    assert c.anchor.value == _scaled(canonical, "value")
    assert c.anchor.value == pytest.approx(2.488e-5)
    assert c.anchor.uncertainty == _scaled(canonical, "uncertainty")
    assert c.anchor.units == str(canonical["units"])
    assert c.anchor.value_summary == str(canonical["value_summary"])
    assert c.anchor.pdg_identifier == str(canonical["pdg_identifier"])
    assert c.anchor.source_url == str(canonical["source_url"])
    assert c.anchor.access_date == str(canonical["access_date"])
    assert c.anchor.snapshot_path == str(canonical["snapshot_path"])
    assert c.anchor.sha256_of_local_snapshot == str(
        canonical["sha256_of_local_snapshot"]
    )

    assert c.anchor.sm_source == str(sm_prediction["source"])
    assert c.anchor.sm_arxiv_id == str(sm_prediction["arxiv_id"])
    assert c.anchor.sm_year == int(sm_prediction["year"])
    assert c.anchor.sm_raw_value == str(sm_prediction["value"])
    assert c.anchor.sm_raw_uncertainty == str(sm_prediction["uncertainty"])
    assert c.anchor.sm_scale == float(sm_prediction["scale"])
    assert c.anchor.sm_prediction == _scaled(sm_prediction, "value")
    assert c.anchor.sm_prediction == pytest.approx(2.477e-5)
    assert c.anchor.sm_uncertainty == _scaled(sm_prediction, "uncertainty")
    assert c.anchor.sm_units == str(sm_prediction["units"])
    assert c.anchor.sm_value_summary == str(sm_prediction["value_summary"])
    assert c.anchor.sm_source_url == str(sm_prediction["source_url"])
    assert c.anchor.sm_access_date == str(sm_prediction["access_date"])
    assert c.anchor.sm_snapshot_path == str(sm_prediction["snapshot_path"])
    assert c.anchor.sm_sha256_of_local_snapshot == str(
        sm_prediction["sha256_of_local_snapshot"]
    )


def test_sm_point_passes(sm_point):
    """Deferred placeholder is non-vetoing at an SM-like baseline."""
    r = ConstraintRegistry.get("K017").evaluate(sm_point)
    assert r.process_id == "K017"
    assert r.passes is True
    assert r.predicted is None
    assert r.experimental == pytest.approx(2.488e-5)
    assert r.sm_prediction == pytest.approx(2.477e-5)
    assert r.severity == Severity.INFO
    assert r.diagnostics["deferred"] is True
    assert r.diagnostics["tree_level_sm_ratio"] is True
    assert r.diagnostics["not_delta_s1_rare_decay"] is True
    assert r.diagnostics["radiative_qed_corrections_required"] is True
    assert "modified W-lepton coupling treatment in 5D RS" in (
        r.diagnostics["evaluation_deferred"]
    )
    assert "radiative QED corrections" in r.diagnostics["evaluation_deferred"]
    assert "four-fermion contact operators" in r.diagnostics["np_entry"]
    assert "not a Delta S=1 rare kaon decay" in r.notes


def test_excluded_point_fails(excluded_point):
    """Even a K001/K002-excluded point passes because K017 evaluation is deferred."""
    r = ConstraintRegistry.get("K017").evaluate(excluded_point)
    assert r.passes is True
    assert r.predicted is None
    assert r.ratio is None
    assert r.budget is None
    assert r.diagnostics["deferred"] is True


def test_evaluate_is_pure(sm_point):
    """Two consecutive evaluations on the same point must yield identical results."""
    c = ConstraintRegistry.get("K017")
    r1 = c.evaluate(sm_point)
    r2 = c.evaluate(sm_point)
    assert r1 == r2


def test_adapter_result_is_marked_deferred():
    """Adapter exposes the explicit deferred bit requested for K017."""
    result = evaluate_rk_lfu_placeholder()
    assert result.deferred is True
    assert result.passes is True
    assert result.predicted is None
    assert result.reason == (
        "requires modified W-lepton coupling treatment in 5D RS + "
        "radiative QED corrections"
    )
