"""Per-constraint contract tests for K002 (``Delta m_K``).

The 5 mandatory tests every constraint test file must include (plan §D):

1. ``test_registered`` — registry lookup yields the right metadata.
2. ``test_anchor_matches_yaml`` — typed view matches the yaml on the
   fields K002 actually reads (``pdg_fit_assuming_cpt``).
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
YAML_PATH = REPO_ROOT / "flavor_catalog" / "processes" / "kaon" / "K002.yaml"


def test_registered():
    c = ConstraintRegistry.get("K002")
    assert c.process_id == "K002"
    assert c.observable == "Delta_m_K"
    # Derived from module path:
    assert getattr(c, "family", None) == "kaon"
    assert getattr(c, "level", None) == ConstraintLevel.PRIMARY
    assert c.severity == Severity.HARD


def test_anchor_matches_yaml():
    """The typed anchor pins the same PDG CPT-assuming fit as the yaml."""
    c = ConstraintRegistry.get("K002")
    raw = yaml.safe_load(open(YAML_PATH))["pdg_or_equivalent"]
    fit = raw["pdg_fit_assuming_cpt"]

    assert c.anchor.year == int(fit["year"])
    assert c.anchor.source == str(fit["source"])
    assert c.anchor.observable == str(fit["observable"])
    assert c.anchor.value == float(fit["value"])
    assert c.anchor.uncertainty == float(fit["uncertainty"])
    assert c.anchor.units == str(fit["units"])
    assert c.anchor.display_value == str(fit["display_value"])
    assert c.anchor.assumptions == tuple(str(item) for item in fit["assumptions"])
    assert c.anchor.source_url == str(fit["source_url"])
    assert c.anchor.access_date == str(fit["access_date"])
    assert c.anchor.snapshot_path == str(fit["snapshot_path"])
    assert c.anchor.sha256 == str(fit["sha256"])
    assert c.anchor.sha256_of_text_snapshot == str(
        fit["sha256_of_text_snapshot"]
    )


def test_sm_point_passes(sm_point):
    """At the SM baseline (zero NP couplings) |M_12^NP| = 0."""
    r = ConstraintRegistry.get("K002").evaluate(sm_point)
    assert r.process_id == "K002"
    assert r.passes is True
    assert r.predicted is not None and r.predicted == 0.0
    assert r.ratio is not None and r.ratio < 1.0
    assert r.budget is not None and r.budget > 0.0


def test_excluded_point_fails(excluded_point):
    """A large (s,d) coupling at 3 TeV blows |M_12^NP| past Delta m_K / 2."""
    r = ConstraintRegistry.get("K002").evaluate(excluded_point)
    assert r.passes is False
    assert r.ratio is not None and r.ratio > 1.0


def test_evaluate_is_pure(sm_point):
    """Two consecutive evaluations on the same point must yield identical results."""
    c = ConstraintRegistry.get("K002")
    r1 = c.evaluate(sm_point)
    r2 = c.evaluate(sm_point)
    assert r1 == r2
