"""Contract tests for the deterministic modern per-point evaluation surface."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from types import MappingProxyType

import pytest

from quarkConstraints.benchmarks import evaluate_default_benchmark
from quarkConstraints.deltaf2 import (
    delta_f2_epsilon_k_budget_policy,
    evaluate_delta_f2_constraints,
)
from quarkConstraints.modern import (
    MODERN_PHENOMENOLOGY_SYSTEM_IDS,
    MODERN_POINT_COUPLINGS_SCHEMA_ID,
    MODERN_POINT_EVALUATION_OBSERVABLE_IDS,
    MODERN_POINT_EVALUATION_SYSTEM_IDS,
    MODERN_POINT_EVALUATION_VERDICT_SCHEMA_ID,
    MODERN_POINT_MATCHING_SCHEMA_ID,
    ModernPointEvaluation,
    ModernPointVerdict,
    default_modern_phenomenology_policy,
    evaluate_modern_point,
)
from quarkConstraints.modern.inputs import (
    MODERN_DEFAULT_INPUT_BUNDLE_ID,
    MODERN_DEFAULT_INPUT_PROVENANCE_ID,
    MODERN_DEFAULT_INPUTS_SCHEMA_ID,
    MODERN_DEFAULT_RESOLUTION_POLICY_ID,
    ModernDefaultInputs,
)
from quarkConstraints.paper_0710_1869.validation import module_has_forbidden_import

REPO_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = REPO_ROOT / "quarkConstraints" / "modern"

FORBIDDEN_SCAN_AND_ARTIFACT_MODULES = {
    "quarkConstraints.modern.artifacts",
    "quarkConstraints.modern.scan",
    "quarkConstraints.paper_0710_1869.artifacts",
    "quarkConstraints.paper_0710_1869.scan",
    "quarkConstraints.scan",
}


def _benchmark_point_evaluation() -> ModernPointEvaluation:
    result = evaluate_default_benchmark().result
    return evaluate_modern_point(result, policy=default_modern_phenomenology_policy())


def test_modern_point_evaluation_emits_explicit_verdicts_for_k_b_d_b_s_and_d0():
    evaluation = _benchmark_point_evaluation()

    assert evaluation.policy.system_ids == MODERN_PHENOMENOLOGY_SYSTEM_IDS
    assert evaluation.system_ids == MODERN_POINT_EVALUATION_SYSTEM_IDS
    assert evaluation.input_bundle_schema_id == MODERN_DEFAULT_INPUTS_SCHEMA_ID
    assert evaluation.input_bundle_id == MODERN_DEFAULT_INPUT_BUNDLE_ID
    assert evaluation.input_provenance_id == MODERN_DEFAULT_INPUT_PROVENANCE_ID
    assert evaluation.input_resolution_policy_id == MODERN_DEFAULT_RESOLUTION_POLICY_ID
    assert evaluation.point_id == evaluation.point_label
    assert evaluation.couplings.schema_id == MODERN_POINT_COUPLINGS_SCHEMA_ID
    assert evaluation.matching.schema_id == MODERN_POINT_MATCHING_SCHEMA_ID
    assert evaluation.couplings.point_id == evaluation.point_id
    assert evaluation.matching.point_id == evaluation.point_id
    assert tuple(verdict.observable_id for verdict in evaluation.verdicts) == (
        MODERN_POINT_EVALUATION_OBSERVABLE_IDS
    )
    assert evaluation.deltaf2.input_bundle_label == MODERN_DEFAULT_INPUT_BUNDLE_ID

    kaon = evaluation.verdict_for("K")
    assert kaon.system_id == "K"
    assert kaon.observable_id == "epsilon_K"
    assert kaon.policy_system_id == "epsilon_K"
    assert kaon.backend_system_id == "K"
    assert kaon.backend_key == "epsilon_k"
    assert isinstance(kaon.weighted_operator_sizes, MappingProxyType)
    assert list(kaon.weighted_operator_sizes) == sorted(kaon.weighted_operator_sizes)
    with pytest.raises(TypeError):
        kaon.weighted_operator_sizes["C1_VLL"] = 0.0

    assert evaluation.verdict_for("B_d").observable_id == "B_d"
    assert evaluation.verdict_for("B_s").observable_id == "B_s"
    assert evaluation.verdict_for("D0").observable_id == "D0"
    assert evaluation.verdict_for("D0").backend_key == "d"
    assert evaluation.verdict_for("D0").backend_system_id == "D"
    assert evaluation.verdict_for("K").schema_id == MODERN_POINT_EVALUATION_VERDICT_SCHEMA_ID


def test_modern_point_evaluation_rejects_unknown_systems_or_missing_verdicts():
    evaluation = _benchmark_point_evaluation()
    policy = evaluation.policy
    deltaf2 = evaluation.deltaf2

    with pytest.raises(ValueError, match="verdicts must contain exactly"):
        ModernPointEvaluation(
            policy=policy,
            point_id=evaluation.point_id,
            point_label=evaluation.point_label,
            M_KK=evaluation.M_KK,
            xi_KK=evaluation.xi_KK,
            couplings=evaluation.couplings,
            matching=evaluation.matching,
            deltaf2=deltaf2,
            verdicts=evaluation.verdicts[:-1],
        )

    unknown = ModernPointVerdict(**evaluation.verdict_for("D0").as_dict())
    object.__setattr__(unknown, "system_id", "X")
    object.__setattr__(unknown, "observable_id", "X")
    object.__setattr__(unknown, "policy_system_id", "X")
    object.__setattr__(unknown, "backend_system_id", "X")
    object.__setattr__(unknown, "backend_key", "x")
    object.__setattr__(unknown, "policy_id", "x")
    object.__setattr__(unknown, "policy_display_name", "X policy")
    object.__setattr__(unknown, "policy_notes", "x")

    with pytest.raises(ValueError, match="unknown systems"):
        ModernPointEvaluation(
            policy=policy,
            point_id=evaluation.point_id,
            point_label=evaluation.point_label,
            M_KK=evaluation.M_KK,
            xi_KK=evaluation.xi_KK,
            couplings=evaluation.couplings,
            matching=evaluation.matching,
            deltaf2=deltaf2,
            verdicts=evaluation.verdicts[:-1] + (unknown,),
        )


@pytest.mark.parametrize(
    ("field_name", "replacement", "expected_message"),
    [
        ("policy_id", "widened-policy", "policy_id"),
        ("policy_display_name", "widened display", "policy_display_name"),
        ("policy_notes", "widened notes", "policy_notes"),
    ],
)
def test_modern_point_evaluation_rejects_verdict_policy_metadata_mutation(
    field_name, replacement, expected_message
):
    evaluation = _benchmark_point_evaluation()
    policy = evaluation.policy
    deltaf2 = evaluation.deltaf2

    verdict = ModernPointVerdict(**evaluation.verdict_for("K").as_dict())
    object.__setattr__(verdict, field_name, replacement)

    with pytest.raises(ValueError, match=expected_message):
        ModernPointEvaluation(
            policy=policy,
            point_id=evaluation.point_id,
            point_label=evaluation.point_label,
            M_KK=evaluation.M_KK,
            xi_KK=evaluation.xi_KK,
            couplings=evaluation.couplings,
            matching=evaluation.matching,
            deltaf2=deltaf2,
            verdicts=(verdict,) + evaluation.verdicts[1:],
        )


def test_modern_point_evaluation_is_deterministic_for_same_point_and_policy():
    evaluation_1 = _benchmark_point_evaluation()
    evaluation_2 = _benchmark_point_evaluation()

    assert evaluation_1.as_dict() == evaluation_2.as_dict()


def test_modern_point_evaluation_bridge_matches_repo_summary():
    result = evaluate_default_benchmark().result
    evaluation = evaluate_modern_point(
        result,
        policy=default_modern_phenomenology_policy(),
        inputs=ModernDefaultInputs(),
    )
    # Match the modern fixed-g_sstar_3 benchmark normalization.
    from quarkConstraints.couplings import compute_quark_kk_gluon_couplings

    repo_summary = evaluate_delta_f2_constraints(
        compute_quark_kk_gluon_couplings(
            result, g_s_star=ModernDefaultInputs().qcd_metadata.g_s_star,
        ),
    )

    assert evaluation.deltaf2.input_bundle_label == MODERN_DEFAULT_INPUT_BUNDLE_ID
    assert evaluation.couplings.input_bundle_id == MODERN_DEFAULT_INPUT_BUNDLE_ID
    assert evaluation.matching.input_bundle_id == MODERN_DEFAULT_INPUT_BUNDLE_ID

    # The modern evaluation uses a different epsilon_K bound (hadronic-based),
    # so only compare the non-kaon systems against the repo backend.
    modern_ratios = evaluation.deltaf2.as_ratio_dict()
    repo_ratios = repo_summary.as_ratio_dict()
    for key in ("b_d_ratio", "b_s_ratio", "d_ratio"):
        assert modern_ratios[key] == pytest.approx(repo_ratios[key])

    epsilon_summary = evaluation.deltaf2.get("epsilon_k")
    expected_budget, expected_direction = (
        delta_f2_epsilon_k_budget_policy().selected_signed_budget(
            epsilon_summary.diagnostics["epsilon_k_np_signed"]
        )
    )
    assert epsilon_summary.diagnostics["epsilon_k_selected_budget_direction"] == expected_direction
    assert epsilon_summary.bound == pytest.approx(expected_budget)
    assert evaluation.verdict_for("K").bound == pytest.approx(expected_budget)


def test_modern_point_evaluation_rejects_xi_kk_mismatch_against_modern_bundle():
    result = evaluate_default_benchmark().result

    with pytest.raises(ValueError, match="xi_KK must match"):
        evaluate_modern_point(
            result,
            policy=default_modern_phenomenology_policy(),
            inputs=ModernDefaultInputs(),
            xi_KK=1.1,
        )


def test_modern_point_evaluation_module_does_not_import_scan_or_artifact_modules():
    assert not module_has_forbidden_import(
        PACKAGE_ROOT / "evaluation.py",
        FORBIDDEN_SCAN_AND_ARTIFACT_MODULES,
    )

    script = """
import importlib
import json
import sys

importlib.import_module("quarkConstraints.modern.evaluation")

forbidden = sorted(
    name
    for name in sys.modules
    if (
        name in {
            "quarkConstraints.modern.artifacts",
            "quarkConstraints.modern.scan",
            "quarkConstraints.paper_0710_1869.artifacts",
            "quarkConstraints.paper_0710_1869.scan",
            "quarkConstraints.scan",
        }
        or name.startswith("quarkConstraints.modern.artifacts.")
        or name.startswith("quarkConstraints.modern.scan.")
        or name.startswith("quarkConstraints.paper_0710_1869.artifacts.")
        or name.startswith("quarkConstraints.paper_0710_1869.scan.")
        or name.startswith("quarkConstraints.scan.")
    )
)
print(json.dumps(forbidden))
"""
    completed = subprocess.run(
        [sys.executable, "-c", script],
        check=True,
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )
    assert json.loads(completed.stdout) == []
