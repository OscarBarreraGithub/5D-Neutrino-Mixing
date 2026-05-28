"""Per-constraint contract tests for K001 (``epsilon_K``).

The 5 mandatory tests every constraint test file must include (plan §D):

1. ``test_registered`` — registry lookup yields the right metadata.
2. ``test_anchor_matches_yaml`` — typed view matches the yaml on the
   fields K001 actually reads (Pattern A,
   ``canonical_experimental_value``).
3. ``test_sm_point_passes`` — at the SM baseline, ``passes=True`` and
   ``ratio < 1``.
4. ``test_excluded_point_fails`` — at an obviously NP-excluded point,
   ``passes=False``.
5. ``test_evaluate_is_pure`` — same input → same output.
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
YAML_PATH = REPO_ROOT / "flavor_catalog" / "processes" / "kaon" / "K001.yaml"


def test_registered():
    c = ConstraintRegistry.get("K001")
    assert c.process_id == "K001"
    assert c.observable == "epsilon_K"
    # Derived from module path:
    assert getattr(c, "family", None) == "kaon"
    assert getattr(c, "level", None) == ConstraintLevel.PRIMARY
    assert c.severity == Severity.HARD


def test_anchor_matches_yaml():
    """The typed anchor pins the same exp/SM values as the canonical yaml."""
    c = ConstraintRegistry.get("K001")
    raw = yaml.safe_load(open(YAML_PATH))["pdg_or_equivalent"]
    exp = raw["canonical_experimental_value"]
    sm = raw["standard_model_reference"]

    assert c.anchor.eps_K_exp_value == float(exp["value"])
    assert c.anchor.eps_K_exp_uncertainty == float(exp["uncertainty"])
    assert c.anchor.eps_K_sm_value == float(sm["value"])
    assert c.anchor.snapshot_path_exp == str(exp["snapshot_path"])
    assert c.anchor.sha256_exp == str(exp["sha256_of_local_snapshot"])
    assert c.anchor.snapshot_path_sm == str(sm["snapshot_path"])
    assert c.anchor.sha256_sm == str(sm["sha256_of_local_snapshot"])


def test_sm_point_passes(sm_point):
    """At the SM baseline (zero NP couplings) epsilon_K^NP = 0."""
    r = ConstraintRegistry.get("K001").evaluate(sm_point)
    assert r.process_id == "K001"
    assert r.passes is True
    assert r.predicted is not None and r.predicted == 0.0
    assert r.ratio is not None and r.ratio < 1.0
    assert r.budget is not None and r.budget > 0.0


def test_excluded_point_fails(excluded_point):
    """A large imaginary (s,d) coupling at 3 TeV blows epsilon_K^NP past the budget."""
    r = ConstraintRegistry.get("K001").evaluate(excluded_point)
    assert r.passes is False
    assert r.ratio is not None and r.ratio > 1.0


def test_evaluate_is_pure(sm_point):
    """Two consecutive evaluations on the same point must yield identical results."""
    c = ConstraintRegistry.get("K001")
    r1 = c.evaluate(sm_point)
    r2 = c.evaluate(sm_point)
    assert r1 == r2
