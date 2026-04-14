"""Contract tests for the modern QS2 bridge sidecar artifact."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from quarkConstraints.benchmarks import evaluate_default_benchmark
from quarkConstraints.modern.bridge_artifacts import (
    MODERN_POINT_BRIDGE_ARTIFACT_SCHEMA_ID,
    MODERN_POINT_COUPLINGS_SCHEMA_ID,
    MODERN_POINT_MATCHING_SCHEMA_ID,
    build_modern_point_bridge_artifact,
    bridge_artifact_from_dict,
    read_modern_point_bridge_artifact,
    write_modern_point_bridge_artifact,
)
from quarkConstraints.modern.evaluation import evaluate_modern_point
from quarkConstraints.modern.verifier import (
    ArtifactSchemaError,
    verify_bridge_artifact_path,
)
from quarkConstraints.paper_0710_1869.validation import module_has_forbidden_import

REPO_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = REPO_ROOT / "quarkConstraints" / "modern"

FORBIDDEN_BACKEND_AND_SCAN_MODULES = {
    "quarkConstraints.benchmarks",
    "quarkConstraints.couplings",
    "quarkConstraints.deltaf2",
    "quarkConstraints.fit",
    "quarkConstraints.model",
    "quarkConstraints.modern.couplings",
    "quarkConstraints.modern.evaluation",
    "quarkConstraints.modern.matching",
    "quarkConstraints.modern.scan",
    "quarkConstraints.paper_0710_1869",
    "quarkConstraints.scan",
}


def _modern_bridge_artifact():
    result = evaluate_default_benchmark().result
    evaluation = evaluate_modern_point(result)
    return build_modern_point_bridge_artifact(evaluation)


def _run_subprocess_json(script: str) -> dict[str, object]:
    completed = subprocess.run(
        [sys.executable, "-c", script],
        check=True,
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )
    return json.loads(completed.stdout)


def test_modern_bridge_artifact_export_is_deterministic(tmp_path: Path) -> None:
    artifact_a = _modern_bridge_artifact()
    artifact_b = _modern_bridge_artifact()

    path_a = tmp_path / "bridge-a.json"
    path_b = tmp_path / "bridge-b.json"
    write_modern_point_bridge_artifact(artifact_a, path_a)
    write_modern_point_bridge_artifact(artifact_b, path_b)

    assert artifact_a.as_dict() == artifact_b.as_dict()
    assert artifact_a.to_json() == artifact_b.to_json()
    assert path_a.read_text(encoding="utf-8") == path_b.read_text(encoding="utf-8")


def test_modern_bridge_artifact_round_trip_preserves_nested_bridge(tmp_path: Path) -> None:
    artifact = _modern_bridge_artifact()
    path = tmp_path / "bridge.json"
    write_modern_point_bridge_artifact(artifact, path)

    loaded = read_modern_point_bridge_artifact(path)

    assert loaded == artifact
    assert loaded.schema_id == MODERN_POINT_BRIDGE_ARTIFACT_SCHEMA_ID
    assert loaded.couplings["schema_id"] == MODERN_POINT_COUPLINGS_SCHEMA_ID
    assert loaded.matching["schema_id"] == MODERN_POINT_MATCHING_SCHEMA_ID
    assert bridge_artifact_from_dict(artifact.as_dict()) == artifact


def test_modern_bridge_verifier_accepts_clean_export(tmp_path: Path) -> None:
    artifact = _modern_bridge_artifact()
    path = tmp_path / "clean-bridge.json"
    write_modern_point_bridge_artifact(artifact, path)

    script = f"""
import json
from pathlib import Path

from quarkConstraints.modern.verifier import verify_bridge_artifact_path

report = verify_bridge_artifact_path(Path({str(path)!r}))
print(json.dumps({{
    "ok": report.ok,
    "issue_codes": list(report.issue_codes),
    "system_ids": list(report.system_ids),
}}))
"""
    payload = _run_subprocess_json(script)

    assert payload["ok"] is True
    assert payload["issue_codes"] == []
    assert tuple(payload["system_ids"]) == ("K", "B_d", "B_s", "D0")


def test_modern_bridge_verifier_rejects_nested_point_id_tampering(tmp_path: Path) -> None:
    artifact = _modern_bridge_artifact()
    payload = artifact.as_dict()
    payload["matching"]["point_id"] = "tampered"
    path = tmp_path / "tampered-bridge.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    with pytest.raises(ArtifactSchemaError):
        verify_bridge_artifact_path(path)


def test_modern_bridge_artifact_module_has_artifact_only_imports() -> None:
    assert not module_has_forbidden_import(
        PACKAGE_ROOT / "bridge_artifacts.py",
        FORBIDDEN_BACKEND_AND_SCAN_MODULES,
    )

    script = """
import importlib
import json
import sys

importlib.import_module("quarkConstraints.modern.bridge_artifacts")

forbidden = sorted(
    name
    for name in sys.modules
    if (
        name in {
            "quarkConstraints.benchmarks",
            "quarkConstraints.couplings",
            "quarkConstraints.deltaf2",
            "quarkConstraints.fit",
            "quarkConstraints.model",
            "quarkConstraints.modern.couplings",
            "quarkConstraints.modern.evaluation",
            "quarkConstraints.modern.matching",
            "quarkConstraints.modern.scan",
            "quarkConstraints.paper_0710_1869",
            "quarkConstraints.scan",
        }
        or name.startswith("quarkConstraints.benchmarks.")
        or name.startswith("quarkConstraints.couplings.")
        or name.startswith("quarkConstraints.deltaf2.")
        or name.startswith("quarkConstraints.fit.")
        or name.startswith("quarkConstraints.model.")
        or name.startswith("quarkConstraints.modern.couplings.")
        or name.startswith("quarkConstraints.modern.evaluation.")
        or name.startswith("quarkConstraints.modern.matching.")
        or name.startswith("quarkConstraints.modern.scan.")
        or name.startswith("quarkConstraints.paper_0710_1869.")
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
