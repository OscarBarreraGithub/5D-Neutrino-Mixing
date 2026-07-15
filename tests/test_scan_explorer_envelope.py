"""Focused tests for Scan Explorer envelope aggregation."""

from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path

from flavor_catalog.website.scripts import build_scan_explorer as explorer


def test_scan_explorer_joint_veto_floor_keeps_disjoint_subthreshold_constraints(monkeypatch):
    monkeypatch.setattr(explorer, "MODEL_ROOTS", {"minimal": Path("unused")})
    monkeypatch.setattr(explorer, "R_GRID", [0.1])
    monkeypatch.setattr(explorer, "MKK_GRID_TEV", [1, 2, 3])

    evaluated_counts = {
        ("minimal", "0.1", 1): 10,
        ("minimal", "0.1", 2): 10,
        ("minimal", "0.1", 3): 10,
    }
    veto_counts = defaultdict(Counter)
    # Each constraint is individually below the 50% threshold at 1 TeV, but
    # the disjoint union vetoes 80% of draws. The old max-of-marginals envelope
    # would have accepted the 1 TeV cell.
    veto_counts[("minimal", "0.1", 1)]["C1"] = 4
    veto_counts[("minimal", "0.1", 1)]["C2"] = 4
    joint_veto_counts = {
        ("minimal", "0.1", 1): 8,
        ("minimal", "0.1", 2): 0,
        ("minimal", "0.1", 3): 0,
    }
    tag_counts = {
        "C1": Counter({"rigorous": 1}),
        "C2": Counter({"proxy": 1}),
    }

    constraints = explorer._build_constraints(
        evaluated_counts=evaluated_counts,
        veto_counts=veto_counts,
        tag_counts=tag_counts,
    )
    joint = explorer._build_joint_veto_tree(
        evaluated_counts=evaluated_counts,
        joint_veto_counts=joint_veto_counts,
    )

    assert [item["id"] for item in constraints] == ["C1", "C2"]
    assert joint["minimal"]["0.1"] == [0.8, 0.0, 0.0]
