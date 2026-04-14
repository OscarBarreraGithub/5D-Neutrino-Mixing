"""Focused tests for the operational modern scan engine."""

from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path

import pytest

from quarkConstraints.benchmarks import default_spurion_seed, evaluate_default_benchmark
from quarkConstraints.fit import QuarkFitSolution
from quarkConstraints.modern.artifacts import read_modern_point_artifact
from quarkConstraints.modern.bridge_artifacts import read_modern_point_bridge_artifact
from quarkConstraints.modern.phenomenology import (
    ModernPointPhenomenologyArtifactV1,
    read_modern_point_phenomenology_artifact,
)


@pytest.fixture
def modern_scan_module():
    module = importlib.import_module("quarkConstraints.modern.scan")
    try:
        yield module
    finally:
        sys.modules.pop("quarkConstraints.modern.scan", None)
        package = sys.modules.get("quarkConstraints.modern")
        if package is not None and hasattr(package, "scan"):
            delattr(package, "scan")


def _stub_fit_solution() -> QuarkFitSolution:
    benchmark_result = evaluate_default_benchmark().result
    return QuarkFitSolution(
        seed=default_spurion_seed(),
        result=benchmark_result,
        success=True,
        message="stubbed benchmark result",
        nfev=0,
        initial_score=float(benchmark_result.score),
    )


def _load_jsonl(path: Path) -> list[dict[str, object]]:
    return [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def _mutated_phenomenology_builder(
    modern_scan_module,
    *,
    system_id: str,
    ratio_to_bound: float,
    passes: bool,
):
    original_builder = modern_scan_module.build_modern_point_phenomenology_artifact

    def _builder(*args, **kwargs):
        artifact = original_builder(*args, **kwargs)
        payload = artifact.as_dict()
        for result in payload["system_results"]:
            if result["system_id"] != system_id:
                continue
            result["ratio_to_bound"] = ratio_to_bound
            result["passes"] = passes
            break
        if system_id in payload["non_cp_acceptance_system_ids"]:
            payload["failing_non_cp_system_ids"] = [system_id] if not passes else []
            payload["non_cp_passes"] = bool(passes)
        else:
            payload["failing_non_cp_system_ids"] = []
            payload["non_cp_passes"] = True
        return ModernPointPhenomenologyArtifactV1.from_dict(payload)

    return _builder


def test_modern_scan_point_enumeration_is_deterministic(modern_scan_module) -> None:
    config = modern_scan_module.ModernScanConfig(
        r_values=[0.25, 0.10, 0.25],
        overall_scale_values=[2.8],
        Lambda_IR_values=[3500.0, 3000.0],
        record_git_metadata=False,
        max_nfev=1,
    )

    points_a = modern_scan_module.enumerate_modern_scan_points(config)
    points_b = modern_scan_module.enumerate_modern_scan_points(config)

    assert config.r_values.tolist() == [0.10, 0.25]
    assert config.Lambda_IR_values.tolist() == [3000.0, 3500.0]
    assert len(points_a) == 4
    assert [point.as_dict() for point in points_a] == [point.as_dict() for point in points_b]
    assert len({point.point_id for point in points_a}) == len(points_a)
    assert all(point.config_hash == config.config_hash for point in points_a)

    shard_sets = []
    for shard_index in range(3):
        shard_points = modern_scan_module.enumerate_modern_scan_shard_points(
            config,
            shard_index=shard_index,
            shard_count=3,
        )
        shard_sets.append({point.point_id for point in shard_points})

    assert set.union(*shard_sets) == {point.point_id for point in points_a}
    assert shard_sets[0].isdisjoint(shard_sets[1])
    assert shard_sets[0].isdisjoint(shard_sets[2])
    assert shard_sets[1].isdisjoint(shard_sets[2])


def test_modern_scan_resume_is_idempotent(tmp_path: Path, monkeypatch, modern_scan_module) -> None:
    monkeypatch.setattr(
        modern_scan_module,
        "fit_quark_sector",
        lambda *args, **kwargs: _stub_fit_solution(),
    )
    config = modern_scan_module.ModernScanConfig(
        r_values=[0.25],
        overall_scale_values=[2.8],
        Lambda_IR_values=[3000.0],
        record_git_metadata=False,
        max_nfev=1,
    )
    output_root = tmp_path / "resume-run"

    first_manifest = modern_scan_module.run_modern_scan_shard(
        config,
        output_dir=output_root,
        shard_index=0,
        shard_count=1,
    )
    results_path = output_root / first_manifest.results_path
    first_results_text = results_path.read_text(encoding="utf-8")

    rows = _load_jsonl(results_path)
    assert len(rows) == 1
    assert rows[0]["schema_id"] == modern_scan_module.MODERN_SCAN_RESULT_SCHEMA_ID
    assert rows[0]["verifier_ok"] is True
    assert rows[0]["bridge_verifier_ok"] is True
    assert rows[0]["phenomenology_verifier_ok"] is True
    assert rows[0]["phenomenology_release_scope_id"] == (
        modern_scan_module.MODERN_POINT_PHENOMENOLOGY_RELEASE_SCOPE_ID
    )
    assert rows[0]["non_cp_acceptance_system_ids"] == ["epsilon_K", "K", "B_d", "B_s", "D0"]
    assert rows[0]["diagnostic_only_system_ids"] == []
    assert rows[0]["blocked_system_ids"] == []
    assert rows[0]["phenomenology_verifier_ok"] is True

    second_manifest = modern_scan_module.run_modern_scan_shard(
        config,
        output_dir=output_root,
        shard_index=0,
        shard_count=1,
    )
    assert second_manifest.completed_point_count == 1
    assert second_manifest.skipped_existing_point_count == 1
    assert results_path.read_text(encoding="utf-8") == first_results_text


def test_modern_scan_resume_requires_bridge_artifact(
    tmp_path: Path,
    monkeypatch,
    modern_scan_module,
) -> None:
    monkeypatch.setattr(
        modern_scan_module,
        "fit_quark_sector",
        lambda *args, **kwargs: _stub_fit_solution(),
    )
    config = modern_scan_module.ModernScanConfig(
        r_values=[0.25],
        overall_scale_values=[2.8],
        Lambda_IR_values=[3000.0],
        record_git_metadata=False,
        max_nfev=1,
    )
    output_root = tmp_path / "resume-bridge-run"

    modern_scan_module.run_modern_scan_shard(
        config,
        output_dir=output_root,
        shard_index=0,
        shard_count=1,
    )
    rows = _load_jsonl(
        output_root / "shards" / "shard-00000-of-00001" / "results.jsonl"
    )
    bridge_artifact_path = output_root / str(rows[0]["bridge_artifact_path"])
    bridge_artifact_path.unlink()

    with pytest.raises(FileNotFoundError, match="missing bridge artifact"):
        modern_scan_module.run_modern_scan_shard(
            config,
            output_dir=output_root,
            shard_index=0,
            shard_count=1,
        )


def test_modern_scan_resume_requires_phenomenology_artifact(
    tmp_path: Path,
    monkeypatch,
    modern_scan_module,
) -> None:
    monkeypatch.setattr(
        modern_scan_module,
        "fit_quark_sector",
        lambda *args, **kwargs: _stub_fit_solution(),
    )
    config = modern_scan_module.ModernScanConfig(
        r_values=[0.25],
        overall_scale_values=[2.8],
        Lambda_IR_values=[3000.0],
        record_git_metadata=False,
        max_nfev=1,
    )
    output_root = tmp_path / "resume-phenomenology-run"

    modern_scan_module.run_modern_scan_shard(
        config,
        output_dir=output_root,
        shard_index=0,
        shard_count=1,
    )
    rows = _load_jsonl(
        output_root / "shards" / "shard-00000-of-00001" / "results.jsonl"
    )
    phenomenology_artifact_path = output_root / str(
        rows[0]["phenomenology_artifact_path"]
    )
    phenomenology_artifact_path.unlink()

    with pytest.raises(FileNotFoundError, match="missing phenomenology artifact"):
        modern_scan_module.run_modern_scan_shard(
            config,
            output_dir=output_root,
            shard_index=0,
            shard_count=1,
        )


def test_modern_scan_shard_and_merge_smoke_path(
    tmp_path: Path,
    monkeypatch,
    modern_scan_module,
) -> None:
    monkeypatch.setattr(
        modern_scan_module,
        "fit_quark_sector",
        lambda *args, **kwargs: _stub_fit_solution(),
    )
    config = modern_scan_module.ModernScanConfig(
        r_values=[0.10, 0.25, 0.40],
        overall_scale_values=[2.8],
        Lambda_IR_values=[3000.0],
        record_git_metadata=False,
        max_nfev=1,
    )
    output_root = tmp_path / "merge-run"

    manifest_0 = modern_scan_module.run_modern_scan_shard(
        config,
        output_dir=output_root,
        shard_index=0,
        shard_count=2,
    )
    manifest_1 = modern_scan_module.run_modern_scan_shard(
        config,
        output_dir=output_root,
        shard_index=1,
        shard_count=2,
    )

    assert manifest_0.complete is True
    assert manifest_1.complete is True
    assert manifest_0.assigned_point_count + manifest_1.assigned_point_count == 3

    merged_manifest = modern_scan_module.merge_modern_scan_shards(output_root)
    merged_results_path = output_root / merged_manifest.merged_results_path
    merged_rows = _load_jsonl(merged_results_path)

    assert merged_manifest.schema_id == modern_scan_module.MODERN_SCAN_MERGE_MANIFEST_SCHEMA_ID
    assert merged_manifest.complete is True
    assert len(merged_rows) == 3
    assert [row["point_index"] for row in merged_rows] == sorted(
        row["point_index"] for row in merged_rows
    )
    assert all(row["verifier_ok"] is True for row in merged_rows)
    assert all(row["bridge_verifier_ok"] is True for row in merged_rows)
    assert all(row["phenomenology_verifier_ok"] is True for row in merged_rows)
    assert all(
        row["phenomenology_release_scope_id"]
        == modern_scan_module.MODERN_POINT_PHENOMENOLOGY_RELEASE_SCOPE_ID
        for row in merged_rows
    )
    assert all(row["non_cp_acceptance_system_ids"] == ["epsilon_K", "K", "B_d", "B_s", "D0"] for row in merged_rows)
    assert all(row["diagnostic_only_system_ids"] == [] for row in merged_rows)
    assert all(row["blocked_system_ids"] == [] for row in merged_rows)

    artifact_path = output_root / str(merged_rows[0]["artifact_path"])
    bridge_artifact_path = output_root / str(merged_rows[0]["bridge_artifact_path"])
    phenomenology_artifact_path = output_root / str(
        merged_rows[0]["phenomenology_artifact_path"]
    )
    artifact = read_modern_point_artifact(artifact_path)
    bridge_artifact = read_modern_point_bridge_artifact(bridge_artifact_path)
    phenomenology_artifact = read_modern_point_phenomenology_artifact(
        phenomenology_artifact_path
    )
    assert artifact.header.point_id == merged_rows[0]["point_id"]
    assert artifact.header.point_label == merged_rows[0]["point_label"]
    assert bridge_artifact.point_id == merged_rows[0]["point_id"]
    assert bridge_artifact.point_label == merged_rows[0]["point_label"]
    assert phenomenology_artifact.point_id == merged_rows[0]["point_id"]
    assert phenomenology_artifact.point_label == merged_rows[0]["point_label"]

    verification = modern_scan_module.verify_merged_scan(
        config,
        output_root=output_root,
        total_shards=2,
    )
    assert verification["schema_id"] == modern_scan_module.MODERN_SCAN_VERIFICATION_SCHEMA_ID
    assert verification["ok"] is True
    assert verification["bridge_verifier_failed_point_count"] == 0
    assert verification["phenomenology_verifier_failed_point_count"] == 0


def test_modern_scan_config_round_trip_and_presets(tmp_path: Path, modern_scan_module) -> None:
    smoke = modern_scan_module.smoke_scan_config()
    pilot = modern_scan_module.pilot_scan_config()
    config_path = tmp_path / "smoke.json"

    modern_scan_module.write_scan_config(smoke, config_path)
    loaded = modern_scan_module.read_scan_config(config_path)

    assert smoke.schema_id == modern_scan_module.MODERN_SCAN_CONFIG_SCHEMA_ID
    assert loaded.as_dict() == smoke.as_dict()
    assert pilot.total_points > smoke.total_points
    assert modern_scan_module.build_point_id(
        smoke,
        r=float(smoke.r_values[0]),
        overall_scale=float(smoke.overall_scale_values[0]),
        Lambda_IR=float(smoke.Lambda_IR_values[0]),
    ) == modern_scan_module.enumerate_scan_points(smoke)[0].point_id


def test_modern_scan_accepts_preset_config_at_run_root(
    tmp_path: Path,
    monkeypatch,
    modern_scan_module,
) -> None:
    monkeypatch.setattr(
        modern_scan_module,
        "fit_quark_sector",
        lambda *args, **kwargs: _stub_fit_solution(),
    )
    output_root = tmp_path / "preset-run"
    config_path = output_root / "config.json"
    config = modern_scan_module.smoke_scan_config()

    modern_scan_module.write_scan_config(config, config_path)
    manifest = modern_scan_module.run_modern_scan_shard(
        config,
        output_dir=output_root,
        shard_index=0,
        shard_count=1,
    )

    run_config_payload = json.loads(config_path.read_text(encoding="utf-8"))
    assert manifest.complete is True
    assert run_config_payload["schema_id"] == modern_scan_module.MODERN_SCAN_RUN_CONFIG_SCHEMA_ID
    assert run_config_payload["config_hash"] == config.config_hash
    assert run_config_payload["config"]["schema_id"] == modern_scan_module.MODERN_SCAN_CONFIG_SCHEMA_ID


def test_modern_scan_rejects_epsilon_k_acceptance_failure(
    tmp_path: Path,
    monkeypatch,
    modern_scan_module,
) -> None:
    monkeypatch.setattr(
        modern_scan_module,
        "fit_quark_sector",
        lambda *args, **kwargs: _stub_fit_solution(),
    )
    monkeypatch.setattr(
        modern_scan_module,
        "build_modern_point_phenomenology_artifact",
        _mutated_phenomenology_builder(
            modern_scan_module,
            system_id="epsilon_K",
            ratio_to_bound=2.0,
            passes=False,
        ),
    )
    config = modern_scan_module.ModernScanConfig(
        r_values=[0.25],
        overall_scale_values=[2.8],
        Lambda_IR_values=[3000.0],
        record_git_metadata=False,
        max_nfev=1,
    )
    output_root = tmp_path / "epsilon-k-failure-run"

    modern_scan_module.run_modern_scan_shard(
        config,
        output_dir=output_root,
        shard_index=0,
        shard_count=1,
    )
    row = _load_jsonl(
        output_root / "shards" / "shard-00000-of-00001" / "results.jsonl"
    )[0]

    # epsilon_K is now acceptance-bearing, so its failure blocks acceptance
    assert row["accepted"] is False
    assert row["phenomenology_passes"] is False
    assert row["failing_non_cp_system_ids"] == ["epsilon_K"]
    assert row["diagnostic_failing_system_ids"] == []
    assert "epsilon_K" in row["non_cp_ratio_to_bound_by_system"]


def test_modern_scan_rejects_non_cp_acceptance_failure(
    tmp_path: Path,
    monkeypatch,
    modern_scan_module,
) -> None:
    monkeypatch.setattr(
        modern_scan_module,
        "fit_quark_sector",
        lambda *args, **kwargs: _stub_fit_solution(),
    )
    monkeypatch.setattr(
        modern_scan_module,
        "build_modern_point_phenomenology_artifact",
        _mutated_phenomenology_builder(
            modern_scan_module,
            system_id="B_d",
            ratio_to_bound=1.5,
            passes=False,
        ),
    )
    config = modern_scan_module.ModernScanConfig(
        r_values=[0.25],
        overall_scale_values=[2.8],
        Lambda_IR_values=[3000.0],
        record_git_metadata=False,
        max_nfev=1,
    )
    output_root = tmp_path / "non-cp-failure-run"

    modern_scan_module.run_modern_scan_shard(
        config,
        output_dir=output_root,
        shard_index=0,
        shard_count=1,
    )
    row = _load_jsonl(
        output_root / "shards" / "shard-00000-of-00001" / "results.jsonl"
    )[0]

    assert row["accepted"] is False
    assert row["phenomenology_passes"] is False
    assert row["failing_non_cp_system_ids"] == ["B_d"]
    assert row["diagnostic_failing_system_ids"] == []
    assert row["blocked_system_ids"] == []
