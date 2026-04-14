"""Contract tests for the modern per-point artifact/verifier slice."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from quarkConstraints.benchmarks import evaluate_default_benchmark
from quarkConstraints.modern import (
    MODERN_POINT_ARTIFACT_REQUIRED_OBSERVABLE_IDS,
    MODERN_POINT_ARTIFACT_REQUIRED_POLICY_SYSTEM_IDS,
    MODERN_POINT_ARTIFACT_REQUIRED_SYSTEM_IDS,
    build_modern_point_artifact,
    default_modern_phenomenology_policy,
)
from quarkConstraints.modern.artifacts import (
    ArtifactSchemaError as ArtifactFormatError,
    artifact_from_dict,
    read_modern_point_artifact,
    write_modern_point_artifact,
)
from quarkConstraints.modern.evaluation import evaluate_modern_point
from quarkConstraints.modern.verifier import (
    ArtifactSchemaError,
    verify_artifact,
    verify_artifact_path,
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
    "quarkConstraints.modern.evaluation",
    "quarkConstraints.modern.scan",
    "quarkConstraints.paper_0710_1869",
    "quarkConstraints.scan",
}


def _modern_point_artifact():
    result = evaluate_default_benchmark().result
    evaluation = evaluate_modern_point(result, policy=default_modern_phenomenology_policy())
    return build_modern_point_artifact(evaluation)


def _run_subprocess_json(script: str) -> dict[str, object]:
    completed = subprocess.run(
        [sys.executable, "-c", script],
        check=True,
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )
    return json.loads(completed.stdout)


def test_modern_point_artifact_export_is_deterministic(tmp_path: Path) -> None:
    artifact_a = _modern_point_artifact()
    artifact_b = _modern_point_artifact()

    path_a = tmp_path / "artifact-a.json"
    path_b = tmp_path / "artifact-b.json"
    write_modern_point_artifact(artifact_a, path_a)
    write_modern_point_artifact(artifact_b, path_b)

    assert artifact_a.as_dict() == artifact_b.as_dict()
    assert artifact_a.to_json() == artifact_b.to_json()
    assert path_a.read_text(encoding="utf-8") == path_b.read_text(encoding="utf-8")


def test_modern_point_artifact_export_round_trips_verdicts_for_k_bd_bs_and_d0(tmp_path: Path) -> None:
    artifact = _modern_point_artifact()
    path = tmp_path / "artifact.json"
    write_modern_point_artifact(artifact, path)

    loaded = read_modern_point_artifact(path)

    assert loaded == artifact
    assert loaded.system_ids == MODERN_POINT_ARTIFACT_REQUIRED_SYSTEM_IDS
    assert loaded.policy.system_ids == MODERN_POINT_ARTIFACT_REQUIRED_POLICY_SYSTEM_IDS
    assert tuple(verdict.observable_id for verdict in loaded.verdicts) == (
        MODERN_POINT_ARTIFACT_REQUIRED_OBSERVABLE_IDS
    )
    assert loaded.verdict_for("K").observable_id == "epsilon_K"
    assert loaded.verdict_for("K").policy_system_id == "epsilon_K"
    assert loaded.verdict_for("B_d").observable_id == "B_d"
    assert loaded.verdict_for("B_s").observable_id == "B_s"
    assert loaded.verdict_for("D0").observable_id == "D0"
    assert artifact_from_dict(artifact.as_dict()) == artifact


def test_modern_verifier_accepts_clean_export_and_rejects_tampering(tmp_path: Path) -> None:
    artifact = _modern_point_artifact()
    path = tmp_path / "clean.json"
    write_modern_point_artifact(artifact, path)

    clean_script = f"""
import json
from pathlib import Path

from quarkConstraints.modern.verifier import verify_artifact_path

report = verify_artifact_path(Path({str(path)!r}))
print(json.dumps({{
    "ok": report.ok,
    "issue_codes": list(report.issue_codes),
    "verdict_system_ids": list(report.verdict_system_ids),
}}))
"""
    payload = _run_subprocess_json(clean_script)
    assert payload["ok"] is True
    assert payload["issue_codes"] == []
    assert tuple(payload["verdict_system_ids"]) == MODERN_POINT_ARTIFACT_REQUIRED_SYSTEM_IDS

    tampered_payload = json.loads(path.read_text(encoding="utf-8"))
    tampered_payload["verdicts"][0]["observable_id"] = "B_d"
    tampered_path = tmp_path / "tampered.json"
    tampered_path.write_text(json.dumps(tampered_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    tampered_script = f"""
import json
from pathlib import Path

from quarkConstraints.modern.verifier import verify_artifact_path

try:
    report = verify_artifact_path(Path({str(tampered_path)!r}))
except Exception as exc:
    print(json.dumps({{
        "raised": type(exc).__name__,
        "message": str(exc),
    }}))
else:
    print(json.dumps({{
        "ok": report.ok,
        "issue_codes": list(report.issue_codes),
    }}))
"""
    tampered_payload = _run_subprocess_json(tampered_script)
    assert tampered_payload.get("ok") is False or "raised" in tampered_payload


def test_modern_verifier_requires_explicit_epsilon_k_plus_k_bd_bs_and_d0_coverage(tmp_path: Path) -> None:
    artifact = _modern_point_artifact()
    payload = artifact.as_dict()
    payload["header"]["policy_system_ids"] = [
        "K",
        "B_d",
        "B_s",
        "D0",
    ]
    path = tmp_path / "missing-epsilon-k.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    with pytest.raises((ArtifactSchemaError, ArtifactFormatError)):
        verify_artifact_path(path)


def test_modern_verifier_has_artifact_only_import_isolation() -> None:
    assert not module_has_forbidden_import(
        PACKAGE_ROOT / "artifacts.py",
        FORBIDDEN_BACKEND_AND_SCAN_MODULES,
    )
    assert not module_has_forbidden_import(
        PACKAGE_ROOT / "verifier.py",
        FORBIDDEN_BACKEND_AND_SCAN_MODULES,
    )

    script = """
import importlib
import json
import sys

importlib.import_module("quarkConstraints.modern.verifier")

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
            "quarkConstraints.modern.evaluation",
            "quarkConstraints.modern.scan",
            "quarkConstraints.paper_0710_1869",
            "quarkConstraints.paper_0710_1869.artifacts",
            "quarkConstraints.paper_0710_1869.benchmarks",
            "quarkConstraints.paper_0710_1869.couplings",
            "quarkConstraints.paper_0710_1869.deltaf2",
            "quarkConstraints.paper_0710_1869.fit",
            "quarkConstraints.paper_0710_1869.model",
            "quarkConstraints.paper_0710_1869.scan",
            "quarkConstraints.scan",
        }
        or name.startswith("quarkConstraints.benchmarks.")
        or name.startswith("quarkConstraints.couplings.")
        or name.startswith("quarkConstraints.deltaf2.")
        or name.startswith("quarkConstraints.fit.")
        or name.startswith("quarkConstraints.model.")
        or name.startswith("quarkConstraints.modern.evaluation.")
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


def test_modern_verifier_rejects_preloaded_forbidden_module_before_import(tmp_path: Path) -> None:
    artifact = _modern_point_artifact()
    path = tmp_path / "clean.json"
    write_modern_point_artifact(artifact, path)

    script = f"""
import json
import sys
import types
from pathlib import Path

sys.modules["quarkConstraints.fit"] = types.ModuleType("quarkConstraints.fit")

from quarkConstraints.modern.verifier import verify_artifact_path

report = verify_artifact_path(Path({str(path)!r}))
print(json.dumps({{
    "ok": report.ok,
    "issue_codes": list(report.issue_codes),
    "forbidden_loaded_modules": list(report.import_isolation.forbidden_loaded_modules),
}}))
"""
    payload = _run_subprocess_json(script)
    assert payload["ok"] is False
    assert "import_isolation_runtime_violation" in payload["issue_codes"]
    assert "import_isolation_failed" in payload["issue_codes"]
    assert "quarkConstraints.fit" in payload["forbidden_loaded_modules"]


def test_modern_verifier_rejects_preloaded_paper_module_before_import(tmp_path: Path) -> None:
    artifact = _modern_point_artifact()
    path = tmp_path / "clean.json"
    write_modern_point_artifact(artifact, path)

    script = f"""
import json
import sys
from pathlib import Path

import quarkConstraints.paper_0710_1869.fit  # noqa: F401

from quarkConstraints.modern.verifier import verify_artifact_path

report = verify_artifact_path(Path({str(path)!r}))
print(json.dumps({{
    "ok": report.ok,
    "issue_codes": list(report.issue_codes),
    "forbidden_loaded_modules": list(report.import_isolation.forbidden_loaded_modules),
}}))
"""
    payload = _run_subprocess_json(script)
    assert payload["ok"] is False
    assert "import_isolation_runtime_violation" in payload["issue_codes"]
    assert "import_isolation_failed" in payload["issue_codes"]
    assert any(
        name == "quarkConstraints.paper_0710_1869.fit"
        or name.startswith("quarkConstraints.paper_0710_1869.")
        for name in payload["forbidden_loaded_modules"]
    )


def test_modern_verifier_rejects_runtime_only_isolation_injection(tmp_path: Path) -> None:
    artifact = _modern_point_artifact()
    path = tmp_path / "clean.json"
    write_modern_point_artifact(artifact, path)

    script = f"""
import importlib
import json
import sys
import types

from pathlib import Path

from quarkConstraints.modern.verifier import verify_artifact_path

artifact_path = Path({str(path)!r})
sys.modules["quarkConstraints.fit"] = types.ModuleType("quarkConstraints.fit")

report = verify_artifact_path(artifact_path)
print(json.dumps({{
    "ok": report.ok,
    "issue_codes": list(report.issue_codes),
    "forbidden_loaded_modules": list(report.import_isolation.forbidden_loaded_modules),
}}))
"""
    completed = subprocess.run(
        [sys.executable, "-c", script],
        check=True,
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )
    payload = json.loads(completed.stdout)
    assert payload["ok"] is False
    assert "import_isolation_runtime_violation" in payload["issue_codes"]
    assert "import_isolation_failed" in payload["issue_codes"]
    assert any(name == "quarkConstraints.fit" for name in payload["forbidden_loaded_modules"])


@pytest.mark.parametrize("finite_value", [float("nan"), float("inf")])
def test_modern_point_artifact_rejects_non_finite_operator_sizes(
    tmp_path: Path, finite_value: float
) -> None:
    artifact = _modern_point_artifact()
    payload = artifact.as_dict()
    payload["verdicts"][0]["weighted_operator_sizes"]["Q1_VLL"] = finite_value
    path = tmp_path / "non_finite.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    script = f"""
import json
from pathlib import Path

from quarkConstraints.modern.verifier import verify_artifact_path

try:
    report = verify_artifact_path(Path({str(path)!r}))
except Exception as exc:
    print(json.dumps({{
        "raised": type(exc).__name__,
        "message": str(exc),
    }}))
else:
    print(json.dumps({{
        "ok": report.ok,
        "issue_codes": list(report.issue_codes),
    }}))
"""
    payload = _run_subprocess_json(script)
    assert payload.get("ok") is False or "raised" in payload


def test_modern_point_artifact_rejects_pass_ratio_mismatch(tmp_path: Path) -> None:
    artifact = _modern_point_artifact()
    payload = artifact.as_dict()
    payload["verdicts"][0]["passes"] = True
    payload["verdicts"][0]["ratio_to_bound"] = 99.0
    path = tmp_path / "ratio_mismatch.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    script = f"""
import json
from pathlib import Path

from quarkConstraints.modern.verifier import verify_artifact_path

report = verify_artifact_path(Path({str(path)!r}))
print(json.dumps({{
    "ok": report.ok,
    "issue_codes": list(report.issue_codes),
}}))
"""
    payload = _run_subprocess_json(script)
    assert payload["ok"] is False
    assert "verdict_pass_ratio_mismatch" in payload["issue_codes"]


def test_modern_point_artifact_rejects_dominant_operator_inconsistency(tmp_path: Path) -> None:
    artifact = _modern_point_artifact()
    payload = artifact.as_dict()
    payload["verdicts"][0]["dominant_operator"] = "Q1_VRR"
    payload["verdicts"][0]["dominant_operator_size"] = 0.0
    path = tmp_path / "dominant_mismatch.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    script = f"""
import json
from pathlib import Path

from quarkConstraints.modern.verifier import verify_artifact_path

try:
    report = verify_artifact_path(Path({str(path)!r}))
except Exception as exc:
    print(json.dumps({{
        "raised": type(exc).__name__,
        "message": str(exc),
    }}))
else:
    print(json.dumps({{
        "ok": report.ok,
        "issue_codes": list(report.issue_codes),
    }}))
"""
    payload = _run_subprocess_json(script)
    assert payload.get("ok") is False or "raised" in payload
