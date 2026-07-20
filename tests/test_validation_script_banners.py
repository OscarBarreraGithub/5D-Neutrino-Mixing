"""Honesty-banner tests for self-consistency validation scripts."""

from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]


def _script_text(relative: str) -> str:
    return (REPO_ROOT / relative).read_text(encoding="utf-8")


def _normalized_text(relative: str) -> str:
    return " ".join(_script_text(relative).split())


def test_benchmark_quark_0710_1869_banner_says_self_consistency_not_paper_validation():
    text = _normalized_text("scripts/benchmark_quark_0710_1869.py")

    assert "SELF-CONSISTENCY" in text
    assert "not physics validation" in text
    assert "no independent paper target numbers" in text
    assert "Lane C" in text and "quarantined" in text
    assert "shrinks to about 1.8x" in text


def test_audit_wilson_rg_banner_says_linear_algebra_not_external_validation():
    text = _normalized_text("scripts/audit_wilson_rg.py")

    assert "SELF-CONSISTENCY / linear-algebra check" in text
    assert "not physics validation against external paper targets" in text
    assert "shares the codebase alpha_s routine" in text
    assert "Lane C" in text and "quarantined" in text
    assert "shrinks to about 1.8x" in text


def test_perez_randall_audit_printed_conclusion_is_convention_dependent():
    text = _normalized_text("scripts/audit_perez_randall_consistency.py")

    assert "convention-dependent, not a robust validation failure" in text
    assert "shrinks to about 1.8x" in text
    assert "not an O(1) convention drift" not in text
