"""Tests for B033."""
from __future__ import annotations
import flavor_catalog_constraints  # noqa: F401
from flavor_catalog_constraints import ConstraintRegistry


def test_registered():
    assert "B033" in ConstraintRegistry.all()

def test_anchor_matches_yaml():
    c = ConstraintRegistry.get("B033")
    assert c.anchor.observable
    assert c.anchor.source

def test_sm_point_passes(sm_point):
    r = ConstraintRegistry.get("B033").evaluate(sm_point)
    assert r.passes is True
    assert r.diagnostics["deferred"] is True

def test_excluded_point_does_not_veto(excluded_point):
    r = ConstraintRegistry.get("B033").evaluate(excluded_point)
    assert r.passes is True

def test_evaluate_is_pure(sm_point):
    c = ConstraintRegistry.get("B033")
    r1 = c.evaluate(sm_point); r2 = c.evaluate(sm_point)
    assert r1.passes == r2.passes
    assert r1.diagnostics["evaluation_deferred"] == r2.diagnostics["evaluation_deferred"]
