"""Registry and phenomenology surface for the modern quark lane."""

from __future__ import annotations

from importlib import import_module

_EXPORTS = {
    "MODERN_BUNDLE_ID_TEMPLATE": (".conventions", "MODERN_BUNDLE_ID_TEMPLATE"),
    "MODERN_CONVENTIONS_SCHEMA_ID": (".conventions", "MODERN_CONVENTIONS_SCHEMA_ID"),
    "MODERN_DEFAULT_BUNDLE_FAMILY_ID": (".conventions", "MODERN_DEFAULT_BUNDLE_FAMILY_ID"),
    "MODERN_DEFAULT_INPUT_BUNDLE_ID": (".inputs", "MODERN_DEFAULT_INPUT_BUNDLE_ID"),
    "MODERN_DEFAULT_INPUT_PROVENANCE_ID": (
        ".inputs",
        "MODERN_DEFAULT_INPUT_PROVENANCE_ID",
    ),
    "MODERN_DEFAULT_INPUTS_SCHEMA_ID": (".inputs", "MODERN_DEFAULT_INPUTS_SCHEMA_ID"),
    "MODERN_PHENOMENOLOGY_POLICY_ID": (".phenomenology", "MODERN_PHENOMENOLOGY_POLICY_ID"),
    "MODERN_PHENOMENOLOGY_SCHEMA_ID": (".phenomenology", "MODERN_PHENOMENOLOGY_SCHEMA_ID"),
    "MODERN_POINT_PHENOMENOLOGY_ARTIFACT_SCHEMA_ID": (
        ".phenomenology",
        "MODERN_POINT_PHENOMENOLOGY_ARTIFACT_SCHEMA_ID",
    ),
    "MODERN_POINT_PHENOMENOLOGY_ARTIFACT_SCHEMA_VERSION": (
        ".phenomenology",
        "MODERN_POINT_PHENOMENOLOGY_ARTIFACT_SCHEMA_VERSION",
    ),
    "MODERN_POINT_PHENOMENOLOGY_NON_CP_ACCEPTANCE_SYSTEM_IDS": (
        ".phenomenology",
        "MODERN_POINT_PHENOMENOLOGY_NON_CP_ACCEPTANCE_SYSTEM_IDS",
    ),
    "MODERN_POINT_PHENOMENOLOGY_RELEASE_SCOPE_ID": (
        ".phenomenology",
        "MODERN_POINT_PHENOMENOLOGY_RELEASE_SCOPE_ID",
    ),
    "MODERN_POINT_PHENOMENOLOGY_SYSTEM_RESULT_SCHEMA_ID": (
        ".phenomenology",
        "MODERN_POINT_PHENOMENOLOGY_SYSTEM_RESULT_SCHEMA_ID",
    ),
    "MODERN_POINT_PHENOMENOLOGY_SYSTEM_TREATMENT_IDS": (
        ".phenomenology",
        "MODERN_POINT_PHENOMENOLOGY_SYSTEM_TREATMENT_IDS",
    ),
    "MODERN_PHENOMENOLOGY_SYSTEM_IDS": (
        ".phenomenology",
        "MODERN_PHENOMENOLOGY_SYSTEM_IDS",
    ),
    "MODERN_PHENOMENOLOGY_SYSTEM_POLICY_ID_TEMPLATE": (
        ".phenomenology",
        "MODERN_PHENOMENOLOGY_SYSTEM_POLICY_ID_TEMPLATE",
    ),
    "MODERN_PHENOMENOLOGY_SYSTEM_POLICY_IDS": (
        ".phenomenology",
        "MODERN_PHENOMENOLOGY_SYSTEM_POLICY_IDS",
    ),
    "MODERN_PHENOMENOLOGY_SYSTEM_SCHEMA_ID": (
        ".phenomenology",
        "MODERN_PHENOMENOLOGY_SYSTEM_SCHEMA_ID",
    ),
    "MODERN_POINT_ARTIFACT_BACKEND_KEYS": (".artifacts", "MODERN_POINT_ARTIFACT_BACKEND_KEYS"),
    "MODERN_POINT_ARTIFACT_BACKEND_SYSTEM_IDS": (
        ".artifacts",
        "MODERN_POINT_ARTIFACT_BACKEND_SYSTEM_IDS",
    ),
    "MODERN_POINT_ARTIFACT_EVALUATION_SCHEMA_ID": (
        ".artifacts",
        "MODERN_POINT_ARTIFACT_EVALUATION_SCHEMA_ID",
    ),
    "MODERN_POINT_ARTIFACT_HEADER_SCHEMA_ID": (
        ".artifacts",
        "MODERN_POINT_ARTIFACT_HEADER_SCHEMA_ID",
    ),
    "MODERN_POINT_ARTIFACT_REQUIRED_OBSERVABLE_IDS": (
        ".artifacts",
        "MODERN_POINT_ARTIFACT_REQUIRED_OBSERVABLE_IDS",
    ),
    "MODERN_POINT_ARTIFACT_REQUIRED_POLICY_SYSTEM_IDS": (
        ".artifacts",
        "MODERN_POINT_ARTIFACT_REQUIRED_POLICY_SYSTEM_IDS",
    ),
    "MODERN_POINT_ARTIFACT_REQUIRED_SYSTEM_IDS": (
        ".artifacts",
        "MODERN_POINT_ARTIFACT_REQUIRED_SYSTEM_IDS",
    ),
    "MODERN_POINT_ARTIFACT_SCHEMA_ID": (".artifacts", "MODERN_POINT_ARTIFACT_SCHEMA_ID"),
    "MODERN_POINT_ARTIFACT_SCHEMA_VERSION": (
        ".artifacts",
        "MODERN_POINT_ARTIFACT_SCHEMA_VERSION",
    ),
    "MODERN_POINT_ARTIFACT_VERDICT_SCHEMA_ID": (
        ".artifacts",
        "MODERN_POINT_ARTIFACT_VERDICT_SCHEMA_ID",
    ),
    "MODERN_POINT_BRIDGE_ARTIFACT_SCHEMA_ID": (
        ".artifacts",
        "MODERN_POINT_BRIDGE_ARTIFACT_SCHEMA_ID",
    ),
    "MODERN_POINT_BRIDGE_ARTIFACT_SCHEMA_VERSION": (
        ".artifacts",
        "MODERN_POINT_BRIDGE_ARTIFACT_SCHEMA_VERSION",
    ),
    "MODERN_STRICT_PAPER_RESOLUTION_POLICY_ID": (
        ".inputs",
        "MODERN_STRICT_PAPER_RESOLUTION_POLICY_ID",
    ),
    "MODERN_STRICT_PAPER_SOURCE_SNAPSHOT_ID": (
        ".inputs",
        "MODERN_STRICT_PAPER_SOURCE_SNAPSHOT_ID",
    ),
    "MODERN_STRICT_PAPER_SOURCE_SNAPSHOT_NAMES": (
        ".inputs",
        "MODERN_STRICT_PAPER_SOURCE_SNAPSHOT_NAMES",
    ),
    "MODERN_STRICT_PAPER_SOURCE_SCHEMA_IDS": (
        ".inputs",
        "MODERN_STRICT_PAPER_SOURCE_SCHEMA_IDS",
    ),
    "MODERN_DEFAULT_BUNDLE_SLUG": (".conventions", "MODERN_DEFAULT_BUNDLE_SLUG"),
    "MODERN_INPUT_NAMESPACE": (".conventions", "MODERN_INPUT_NAMESPACE"),
    "MODERN_INPUT_REGISTRY_SCHEMA_ID": (".conventions", "MODERN_INPUT_REGISTRY_SCHEMA_ID"),
    "MODERN_LANE_ID": (".conventions", "MODERN_LANE_ID"),
    "MODERN_POINT_EVALUATION_BACKEND_KEYS": (".evaluation", "MODERN_POINT_EVALUATION_BACKEND_KEYS"),
    "MODERN_POINT_EVALUATION_BACKEND_SYSTEM_IDS": (
        ".evaluation",
        "MODERN_POINT_EVALUATION_BACKEND_SYSTEM_IDS",
    ),
    "MODERN_POINT_EVALUATION_OBSERVABLE_IDS": (
        ".evaluation",
        "MODERN_POINT_EVALUATION_OBSERVABLE_IDS",
    ),
    "MODERN_POINT_EVALUATION_POLICY_ID": (".evaluation", "MODERN_POINT_EVALUATION_POLICY_ID"),
    "MODERN_POINT_EVALUATION_POLICY_SYSTEM_IDS": (
        ".evaluation",
        "MODERN_POINT_EVALUATION_POLICY_SYSTEM_IDS",
    ),
    "MODERN_POINT_EVALUATION_SCHEMA_ID": (".evaluation", "MODERN_POINT_EVALUATION_SCHEMA_ID"),
    "MODERN_POINT_EVALUATION_SYSTEM_IDS": (".evaluation", "MODERN_POINT_EVALUATION_SYSTEM_IDS"),
    "MODERN_POINT_EVALUATION_VERDICT_SCHEMA_ID": (
        ".evaluation",
        "MODERN_POINT_EVALUATION_VERDICT_SCHEMA_ID",
    ),
    "MODERN_POINT_COUPLINGS_SCHEMA_ID": (".couplings", "MODERN_POINT_COUPLINGS_SCHEMA_ID"),
    "MODERN_POINT_COUPLING_CONVENTION_ID": (
        ".couplings",
        "MODERN_POINT_COUPLING_CONVENTION_ID",
    ),
    "MODERN_POINT_MATCHING_SCHEMA_ID": (".matching", "MODERN_POINT_MATCHING_SCHEMA_ID"),
    "MODERN_POINT_MATCHING_SYSTEM_SCHEMA_ID": (
        ".matching",
        "MODERN_POINT_MATCHING_SYSTEM_SCHEMA_ID",
    ),
    "MODERN_POINT_MATCHING_BACKEND_ID": (
        ".matching",
        "MODERN_POINT_MATCHING_BACKEND_ID",
    ),
    "MODERN_SCAN_CLAIM_LEVEL": (".scan", "MODERN_SCAN_CLAIM_LEVEL"),
    "MODERN_SCAN_CONFIG_SCHEMA_ID": (".scan", "MODERN_SCAN_CONFIG_SCHEMA_ID"),
    "MODERN_SCAN_MERGE_MANIFEST_SCHEMA_ID": (
        ".scan",
        "MODERN_SCAN_MERGE_MANIFEST_SCHEMA_ID",
    ),
    "MODERN_SCAN_NOTES": (".scan", "MODERN_SCAN_NOTES"),
    "MODERN_SCAN_POINT_SCHEMA_ID": (".scan", "MODERN_SCAN_POINT_SCHEMA_ID"),
    "MODERN_SCAN_RESULT_SCHEMA_ID": (".scan", "MODERN_SCAN_RESULT_SCHEMA_ID"),
    "MODERN_SCAN_RUN_CONFIG_SCHEMA_ID": (".scan", "MODERN_SCAN_RUN_CONFIG_SCHEMA_ID"),
    "MODERN_SCAN_SHARD_MANIFEST_SCHEMA_ID": (
        ".scan",
        "MODERN_SCAN_SHARD_MANIFEST_SCHEMA_ID",
    ),
    "MODERN_SCAN_VERIFICATION_SCHEMA_ID": (".scan", "MODERN_SCAN_VERIFICATION_SCHEMA_ID"),
    "MODERN_PROVENANCE_ID_TEMPLATE": (".conventions", "MODERN_PROVENANCE_ID_TEMPLATE"),
    "MODERN_PROVENANCE_NAMESPACE": (".conventions", "MODERN_PROVENANCE_NAMESPACE"),
    "MODERN_STRICT_PAPER_BUNDLE_FAMILY_ID": (
        ".conventions",
        "MODERN_STRICT_PAPER_BUNDLE_FAMILY_ID",
    ),
    "MODERN_STRICT_PAPER_BUNDLE_SLUG": (".conventions", "MODERN_STRICT_PAPER_BUNDLE_SLUG"),
    "MODERN_STRICT_PAPER_INPUT_BUNDLE_ID": (".inputs", "MODERN_STRICT_PAPER_INPUT_BUNDLE_ID"),
    "MODERN_STRICT_PAPER_INPUT_PROVENANCE_ID": (
        ".inputs",
        "MODERN_STRICT_PAPER_INPUT_PROVENANCE_ID",
    ),
    "MODERN_STRICT_PAPER_INPUTS_SCHEMA_ID": (".inputs", "MODERN_STRICT_PAPER_INPUTS_SCHEMA_ID"),
    "ModernDefaultInputs": (".inputs", "ModernDefaultInputs"),
    "ModernPointArtifactHeader": (".artifacts", "ModernPointArtifactHeader"),
    "ModernPointArtifactV1": (".artifacts", "ModernPointArtifactV1"),
    "ModernPointArtifactVerdict": (".artifacts", "ModernPointArtifactVerdict"),
    "ModernPointBridgeArtifactV1": (
        ".artifacts",
        "ModernPointBridgeArtifactV1",
    ),
    "ModernInputRegistry": (".inputs", "ModernInputRegistry"),
    "ModernLaneConventions": (".conventions", "ModernLaneConventions"),
    "ModernMergedScanManifest": (".scan", "ModernMergedScanManifest"),
    "ModernPointCouplings": (".couplings", "ModernPointCouplings"),
    "ModernPointEvaluation": (".evaluation", "ModernPointEvaluation"),
    "ModernPointMatching": (".matching", "ModernPointMatching"),
    "ModernPointMatchingSystem": (".matching", "ModernPointMatchingSystem"),
    "ModernPointVerdict": (".evaluation", "ModernPointVerdict"),
    "ModernScanConfig": (".scan", "ModernScanConfig"),
    "ModernScanPoint": (".scan", "ModernScanPoint"),
    "ModernScanResultRecord": (".scan", "ModernScanResultRecord"),
    "ModernScanShardManifest": (".scan", "ModernScanShardManifest"),
    "ModernPhenomenologyPolicy": (".phenomenology", "ModernPhenomenologyPolicy"),
    "ModernPhenomenologySystemResult": (
        ".phenomenology",
        "ModernPhenomenologySystemResult",
    ),
    "ModernPhenomenologySystemPolicy": (
        ".phenomenology",
        "ModernPhenomenologySystemPolicy",
    ),
    "ModernPointPhenomenologyArtifactV1": (
        ".phenomenology",
        "ModernPointPhenomenologyArtifactV1",
    ),
    "ModernStrictPaperInputs": (".inputs", "ModernStrictPaperInputs"),
    "VerifierSubprocessResult": (".scan", "VerifierSubprocessResult"),
    "default_modern_default_inputs": (".inputs", "default_modern_default_inputs"),
    "build_modern_point_artifact": (".artifacts", "build_modern_point_artifact"),
    "build_modern_point_bridge_artifact": (
        ".artifacts",
        "build_modern_point_bridge_artifact",
    ),
    "build_modern_point_couplings": (".couplings", "build_modern_point_couplings"),
    "build_modern_point_matching": (".matching", "build_modern_point_matching"),
    "default_modern_point_artifact": (".artifacts", "default_modern_point_artifact"),
    "default_modern_point_bridge_artifact": (
        ".artifacts",
        "default_modern_point_bridge_artifact",
    ),
    "default_modern_input_registry": (".inputs", "default_modern_input_registry"),
    "default_modern_point_couplings": (".couplings", "default_modern_point_couplings"),
    "default_modern_point_evaluation": (".evaluation", "default_modern_point_evaluation"),
    "default_modern_point_matching": (".matching", "default_modern_point_matching"),
    "default_modern_lane_conventions": (".conventions", "default_modern_lane_conventions"),
    "default_modern_phenomenology_policy": (
        ".phenomenology",
        "default_modern_phenomenology_policy",
    ),
    "default_modern_phenomenology_system_policies": (
        ".phenomenology",
        "default_modern_phenomenology_system_policies",
    ),
    "default_modern_strict_paper_inputs": (
        ".inputs",
        "default_modern_strict_paper_inputs",
    ),
    "build_point_id": (".scan", "build_point_id"),
    "enumerate_modern_scan_points": (".scan", "enumerate_modern_scan_points"),
    "enumerate_modern_scan_shard_points": (".scan", "enumerate_modern_scan_shard_points"),
    "enumerate_scan_points": (".scan", "enumerate_scan_points"),
    "merge_modern_scan_shards": (".scan", "merge_modern_scan_shards"),
    "merge_scan_shards": (".scan", "merge_scan_shards"),
    "pilot_scan_config": (".scan", "pilot_scan_config"),
    "point_in_shard": (".scan", "point_in_shard"),
    "read_modern_point_artifact": (".artifacts", "read_modern_point_artifact"),
    "read_modern_point_bridge_artifact": (
        ".artifacts",
        "read_modern_point_bridge_artifact",
    ),
    "read_modern_point_phenomenology_artifact": (
        ".phenomenology",
        "read_modern_point_phenomenology_artifact",
    ),
    "read_scan_config": (".scan", "read_scan_config"),
    "run_scan_shard": (".scan", "run_scan_shard"),
    "run_modern_scan_shard": (".scan", "run_modern_scan_shard"),
    "scan_points_for_shard": (".scan", "scan_points_for_shard"),
    "seed_from_dict": (".scan", "seed_from_dict"),
    "seed_to_dict": (".scan", "seed_to_dict"),
    "smoke_scan_config": (".scan", "smoke_scan_config"),
    "verify_merged_scan": (".scan", "verify_merged_scan"),
    "verify_artifact": (".verifier", "verify_artifact"),
    "verify_artifact_path": (".verifier", "verify_artifact_path"),
    "verify_bridge_artifact": (".verifier", "verify_bridge_artifact"),
    "verify_bridge_artifact_path": (".verifier", "verify_bridge_artifact_path"),
    "verify_phenomenology_artifact": (
        ".verifier",
        "verify_phenomenology_artifact",
    ),
    "verify_phenomenology_artifact_path": (
        ".verifier",
        "verify_phenomenology_artifact_path",
    ),
    "write_scan_config": (".scan", "write_scan_config"),
    "write_modern_point_artifact": (".artifacts", "write_modern_point_artifact"),
    "write_modern_point_bridge_artifact": (
        ".artifacts",
        "write_modern_point_bridge_artifact",
    ),
    "write_modern_point_phenomenology_artifact": (
        ".phenomenology",
        "write_modern_point_phenomenology_artifact",
    ),
    "build_modern_point_phenomenology_artifact": (
        ".phenomenology",
        "build_modern_point_phenomenology_artifact",
    ),
    "evaluate_modern_point": (".evaluation", "evaluate_modern_point"),
}

__all__ = list(_EXPORTS)


def __getattr__(name: str):
    try:
        module_name, attr_name = _EXPORTS[name]
    except KeyError as exc:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from exc
    module = import_module(module_name, __name__)
    value = getattr(module, attr_name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
