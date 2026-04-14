"""Contract tests for the modern phenomenology policy fence."""

from __future__ import annotations

import importlib
import json
import subprocess
import sys
from pathlib import Path

import pytest

from quarkConstraints.benchmarks import evaluate_default_benchmark
from quarkConstraints.modern.conventions import MODERN_INPUT_REGISTRY_SCHEMA_ID, MODERN_LANE_ID
from quarkConstraints.modern.phenomenology import (
    MODERN_PHENOMENOLOGY_POLICY_ID,
    MODERN_PHENOMENOLOGY_SCHEMA_ID,
    MODERN_PHENOMENOLOGY_SYSTEM_IDS,
    MODERN_PHENOMENOLOGY_SYSTEM_POLICY_IDS,
    MODERN_PHENOMENOLOGY_SYSTEM_SCHEMA_ID,
    DEFAULT_MODERN_PHENOMENOLOGY_SYSTEM_POLICIES,
    ModernPhenomenologyPolicy,
    ModernPointPhenomenologyArtifactV1,
    ModernPhenomenologySystemPolicy,
    build_modern_point_phenomenology_artifact,
    default_modern_phenomenology_policy,
    default_modern_phenomenology_system_policies,
    read_modern_point_phenomenology_artifact,
    write_modern_point_phenomenology_artifact,
    MODERN_POINT_PHENOMENOLOGY_NON_CP_ACCEPTANCE_SYSTEM_IDS,
    MODERN_POINT_PHENOMENOLOGY_SYSTEM_TREATMENT_IDS,
)
from quarkConstraints.modern.bridge_artifacts import build_modern_point_bridge_artifact
from quarkConstraints.modern.evaluation import evaluate_modern_point
from quarkConstraints.modern.verifier import verify_phenomenology_artifact_path
from quarkConstraints.paper_0710_1869.validation import module_has_forbidden_import

REPO_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = REPO_ROOT / "quarkConstraints" / "modern"

FORBIDDEN_REPO_V1_MODULES = {
    "deltaf2",
    "quarkConstraints.benchmarks",
    "quarkConstraints.couplings",
    "quarkConstraints.deltaf2",
    "quarkConstraints.fit",
    "quarkConstraints.model",
    "quarkConstraints.proxies",
    "quarkConstraints.scan",
    "quarkConstraints.scales",
    "quarkConstraints.validation",
}


def _modern_phenomenology_artifact() -> ModernPointPhenomenologyArtifactV1:
    result = evaluate_default_benchmark().result
    evaluation = evaluate_modern_point(result)
    bridge_artifact = build_modern_point_bridge_artifact(evaluation)
    return build_modern_point_phenomenology_artifact(
        bridge_artifact,
        policy=evaluation.policy,
    )


def test_modern_phenomenology_policy_is_closed_and_versioned():
    policy = default_modern_phenomenology_policy()

    assert policy.schema_id == MODERN_PHENOMENOLOGY_SCHEMA_ID
    assert policy.policy_id == MODERN_PHENOMENOLOGY_POLICY_ID
    assert policy.lane_id == MODERN_LANE_ID
    assert policy.registry_schema_id == MODERN_INPUT_REGISTRY_SCHEMA_ID
    assert policy.system_ids == MODERN_PHENOMENOLOGY_SYSTEM_IDS
    assert tuple(system.schema_id for system in policy.systems) == (
        MODERN_PHENOMENOLOGY_SYSTEM_SCHEMA_ID,
    ) * len(MODERN_PHENOMENOLOGY_SYSTEM_IDS)
    assert tuple(system.policy_id for system in policy.systems) == MODERN_PHENOMENOLOGY_SYSTEM_POLICY_IDS
    assert tuple(system.as_dict()["system_id"] for system in policy.systems) == MODERN_PHENOMENOLOGY_SYSTEM_IDS
    assert tuple(system.as_dict()["policy_id"] for system in policy.systems) == MODERN_PHENOMENOLOGY_SYSTEM_POLICY_IDS
    assert policy.system_policy("B_s").system_id == "B_s"
    assert policy.require_system_ids(MODERN_PHENOMENOLOGY_SYSTEM_IDS) == MODERN_PHENOMENOLOGY_SYSTEM_IDS

    round_trip = ModernPhenomenologyPolicy(
        schema_id=policy.schema_id,
        policy_id=policy.policy_id,
        lane_id=policy.lane_id,
        registry_schema_id=policy.registry_schema_id,
        systems=default_modern_phenomenology_system_policies(),
        notes=policy.notes,
    )
    assert round_trip.as_dict() == policy.as_dict()


def test_modern_phenomenology_lists_epsilon_k_k_bd_k_bs_and_d0_explicitly():
    policy = ModernPhenomenologyPolicy()

    expected = (
        ("epsilon_K", "epsilon_K CP-violating kaon-mixing policy"),
        ("K", "K neutral-kaon mixing policy"),
        ("B_d", "B_d neutral-B mixing policy"),
        ("B_s", "B_s neutral-B mixing policy"),
        ("D0", "D0 neutral-D mixing policy"),
    )

    assert tuple((system.system_id, system.display_name) for system in policy.systems) == expected
    assert [system.system_id for system in DEFAULT_MODERN_PHENOMENOLOGY_SYSTEM_POLICIES] == list(
        MODERN_PHENOMENOLOGY_SYSTEM_IDS
    )
    assert [system.policy_id for system in DEFAULT_MODERN_PHENOMENOLOGY_SYSTEM_POLICIES] == list(
        MODERN_PHENOMENOLOGY_SYSTEM_POLICY_IDS
    )


@pytest.mark.parametrize(
    ("system_ids", "expected_message"),
    [
        (("epsilon_K", "K", "B_d", "B_s"), "missing systems"),
        (("epsilon_K", "K", "B_d", "B_s", "X"), "unknown systems"),
    ],
)
def test_modern_phenomenology_rejects_unknown_or_missing_systems(system_ids, expected_message):
    policy = ModernPhenomenologyPolicy()

    with pytest.raises(ValueError, match=expected_message):
        policy.require_system_ids(system_ids)


def test_modern_phenomenology_rejects_widened_policy_system_tuple():
    with pytest.raises(ValueError, match="systems must contain exactly"):
        ModernPhenomenologyPolicy(systems=DEFAULT_MODERN_PHENOMENOLOGY_SYSTEM_POLICIES[:-1])

    unknown = ModernPhenomenologySystemPolicy(
        system_id="D0",
        policy_id=MODERN_PHENOMENOLOGY_SYSTEM_POLICY_IDS[-1],
        display_name="D0 neutral-D mixing policy",
        notes="Explicit D0 neutral-meson mixing policy entry.",
    )
    object.__setattr__(unknown, "system_id", "X")
    object.__setattr__(unknown, "policy_id", "quarkConstraints.modern.phenomenology.system.X.v1")
    object.__setattr__(unknown, "display_name", "X neutral-meson mixing policy")
    object.__setattr__(unknown, "notes", "Explicit unknown-system policy entry.")

    with pytest.raises(ValueError, match="unknown systems"):
        ModernPhenomenologyPolicy(
            systems=DEFAULT_MODERN_PHENOMENOLOGY_SYSTEM_POLICIES[:-1] + (unknown,)
        )


def test_modern_point_phenomenology_sidecar_separates_epsilon_k_from_blocked_k(tmp_path: Path):
    artifact = _modern_phenomenology_artifact()
    path = tmp_path / "phenomenology.json"
    write_modern_point_phenomenology_artifact(artifact, path)
    loaded = read_modern_point_phenomenology_artifact(path)

    assert loaded == artifact
    assert loaded.system_ids == MODERN_PHENOMENOLOGY_SYSTEM_IDS
    assert (
        loaded.non_cp_acceptance_system_ids
        == MODERN_POINT_PHENOMENOLOGY_NON_CP_ACCEPTANCE_SYSTEM_IDS
    )
    assert loaded.kaon_viability_claimed is False

    epsilon_k = loaded.system_result("epsilon_K")
    assert epsilon_k.evaluated_from_bridge is True
    assert epsilon_k.included_in_non_cp_acceptance is False
    assert epsilon_k.bridge_system_id == "K"
    assert epsilon_k.bridge_backend_key == "epsilon_k"
    assert epsilon_k.treatment_id == MODERN_POINT_PHENOMENOLOGY_SYSTEM_TREATMENT_IDS["epsilon_K"]

    kaon = loaded.system_result("K")
    assert kaon.evaluated_from_bridge is False
    assert kaon.included_in_non_cp_acceptance is False
    assert kaon.passes is None
    assert kaon.ratio_to_bound is None
    assert kaon.treatment_id == MODERN_POINT_PHENOMENOLOGY_SYSTEM_TREATMENT_IDS["K"]

    d0 = loaded.system_result("D0")
    assert d0.evaluated_from_bridge is True
    assert d0.included_in_non_cp_acceptance is True
    assert d0.treatment_id == MODERN_POINT_PHENOMENOLOGY_SYSTEM_TREATMENT_IDS["D0"]
    assert "conservative" in d0.note


def test_modern_phenomenology_verifier_rejects_widened_non_cp_acceptance(tmp_path: Path):
    artifact = _modern_phenomenology_artifact()
    payload = artifact.as_dict()
    payload["system_results"][0]["included_in_non_cp_acceptance"] = True
    path = tmp_path / "tampered-phenomenology.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    report = verify_phenomenology_artifact_path(path)
    assert report.ok is False
    assert "phenomenology_epsilon_k_in_non_cp_acceptance_unexpected" in report.issue_codes


def test_modern_phenomenology_module_does_not_import_repo_v1_scan_or_deltaf2_stack():
    assert not module_has_forbidden_import(
        PACKAGE_ROOT / "phenomenology.py",
        FORBIDDEN_REPO_V1_MODULES,
    )


def test_importing_modern_phenomenology_does_not_load_repo_v1_scan_or_deltaf2_stack():
    script = """
import importlib
import json
import sys

importlib.import_module("quarkConstraints.modern.phenomenology")

forbidden = sorted(
    name
    for name in sys.modules
    if (
        name in {
            "deltaf2",
            "quarkConstraints.benchmarks",
            "quarkConstraints.couplings",
            "quarkConstraints.deltaf2",
            "quarkConstraints.fit",
            "quarkConstraints.model",
            "quarkConstraints.proxies",
            "quarkConstraints.scan",
            "quarkConstraints.scales",
            "quarkConstraints.validation",
        }
        or name.startswith("quarkConstraints.benchmarks.")
        or name.startswith("quarkConstraints.couplings.")
        or name.startswith("quarkConstraints.deltaf2.")
        or name.startswith("quarkConstraints.fit.")
        or name.startswith("quarkConstraints.model.")
        or name.startswith("quarkConstraints.proxies.")
        or name.startswith("quarkConstraints.scan.")
        or name.startswith("quarkConstraints.scales.")
        or name.startswith("quarkConstraints.validation.")
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
