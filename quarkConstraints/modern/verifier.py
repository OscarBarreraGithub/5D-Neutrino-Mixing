"""Artifact-only verification helpers for the modern quark lane.

HONESTY BANNER: these verifiers check schema and exported-field
self-consistency only. They do not recompute Wilson matching, QCD running,
hadronic matrix elements, or physics bounds from independent oracles; a green
modern verifier result is not physics validation. Independent physics checks
live in the Delta-F=2/GGMS and QCD-running oracle tests, for example
``tests/test_epsilon_k_physics.py`` and ``tests/test_qcd_running.py``.
"""

from __future__ import annotations

import ast
import math
import sys
from dataclasses import dataclass
from pathlib import Path

from .artifacts import (
    MODERN_POINT_ARTIFACT_BACKEND_KEYS,
    MODERN_POINT_ARTIFACT_BACKEND_SYSTEM_IDS,
    MODERN_POINT_ARTIFACT_EVALUATION_SCHEMA_ID,
    MODERN_POINT_ARTIFACT_HEADER_SCHEMA_ID,
    MODERN_POINT_ARTIFACT_REQUIRED_OBSERVABLE_IDS,
    MODERN_POINT_ARTIFACT_REQUIRED_POLICY_SYSTEM_IDS,
    MODERN_POINT_ARTIFACT_REQUIRED_SYSTEM_IDS,
    MODERN_POINT_ARTIFACT_SCHEMA_ID,
    MODERN_POINT_ARTIFACT_SCHEMA_VERSION,
    MODERN_POINT_ARTIFACT_VERDICT_SCHEMA_ID,
    MODERN_POINT_BRIDGE_ARTIFACT_SCHEMA_ID,
    MODERN_POINT_BRIDGE_ARTIFACT_SCHEMA_VERSION,
    MODERN_POINT_BRIDGE_COUPLINGS_SCHEMA_ID as MODERN_BRIDGE_COUPLINGS_SCHEMA_ID,
    MODERN_POINT_BRIDGE_MATCHING_SCHEMA_ID as MODERN_BRIDGE_MATCHING_SCHEMA_ID,
    ModernPointArtifactHeader,
    ModernPointArtifactV1,
    ModernPointBridgeArtifactV1,
    read_modern_point_artifact,
    read_modern_point_bridge_artifact,
)
from .phenomenology import (
    MODERN_POINT_PHENOMENOLOGY_ARTIFACT_SCHEMA_ID,
    MODERN_POINT_PHENOMENOLOGY_ARTIFACT_SCHEMA_VERSION,
    MODERN_POINT_PHENOMENOLOGY_NON_CP_ACCEPTANCE_SYSTEM_IDS,
    MODERN_POINT_PHENOMENOLOGY_RELEASE_SCOPE_ID,
    MODERN_POINT_PHENOMENOLOGY_SYSTEM_RESULT_SCHEMA_ID,
    MODERN_POINT_PHENOMENOLOGY_SYSTEM_TREATMENT_IDS,
    ModernPointPhenomenologyArtifactV1,
    read_modern_point_phenomenology_artifact,
)

PACKAGE_ROOT = "quarkConstraints.modern"
ALLOWED_RELATIVE_IMPORTS = ("artifacts", "phenomenology")
ALLOWED_ABSOLUTE_LOCAL_IMPORTS = (
    f"{PACKAGE_ROOT}.artifacts",
    f"{PACKAGE_ROOT}.phenomenology",
)
MODERN_BRIDGE_INPUTS_SCHEMA_ID = "quarkConstraints.modern.inputs.modern_default_inputs.v2"
FORBIDDEN_EXTERNAL_MODULE_PREFIXES = (
    "quarkConstraints.benchmarks",
    "quarkConstraints.couplings",
    "quarkConstraints.deltaf2",
    "quarkConstraints.fit",
    "quarkConstraints.model",
    "quarkConstraints.paper_0710_1869",
    "quarkConstraints.scan",
    "quarkConstraints.modern.evaluation",
    "quarkConstraints.modern.scan",
)


class ArtifactSchemaError(ValueError):
    """Raised when a modern artifact violates the frozen JSON contract."""


def _require_text(name: str, value: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ArtifactSchemaError(f"{name} must be a non-empty string")
    return value.strip()


def _require_positive_float(name: str, value: float) -> float:
    numeric = float(value)
    if numeric <= 0.0:
        raise ArtifactSchemaError(f"{name} must be positive")
    return numeric


def _append_issue(
    issues: list[VerificationIssue],
    *,
    code: str,
    message: str,
    subject: str | None = None,
) -> None:
    issues.append(VerificationIssue(code=code, message=message, subject=subject))


@dataclass(frozen=True)
class VerificationIssue:
    """One verification failure."""

    code: str
    message: str
    subject: str | None = None


@dataclass(frozen=True)
class ImportIsolationReport:
    """Static/runtime import-isolation status for the standalone verifier."""

    verifier_module: str
    verifier_source_path: str
    allowed_relative_imports: tuple[str, ...]
    allowed_absolute_imports: tuple[str, ...]
    static_ok: bool
    runtime_ok: bool
    unexpected_import_targets: tuple[str, ...]
    forbidden_loaded_modules: tuple[str, ...]

    @property
    def ok(self) -> bool:
        return self.static_ok and self.runtime_ok


@dataclass(frozen=True)
class FrozenVerifierContract:
    """Frozen per-point modern artifact contract enforced by the verifier."""

    schema_id: str = MODERN_POINT_ARTIFACT_SCHEMA_ID
    schema_version: int = MODERN_POINT_ARTIFACT_SCHEMA_VERSION
    header_schema_id: str = MODERN_POINT_ARTIFACT_HEADER_SCHEMA_ID
    policy_schema_id: str = "quarkConstraints.modern.phenomenology.v1"
    policy_id: str = "quarkConstraints.modern.phenomenology.policy.v1"
    evaluation_schema_id: str = MODERN_POINT_ARTIFACT_EVALUATION_SCHEMA_ID
    verdict_schema_id: str = MODERN_POINT_ARTIFACT_VERDICT_SCHEMA_ID
    policy_system_ids: tuple[str, ...] = MODERN_POINT_ARTIFACT_REQUIRED_POLICY_SYSTEM_IDS
    verdict_system_ids: tuple[str, ...] = MODERN_POINT_ARTIFACT_REQUIRED_SYSTEM_IDS
    verdict_observable_ids: tuple[str, ...] = MODERN_POINT_ARTIFACT_REQUIRED_OBSERVABLE_IDS
    backend_system_ids: tuple[str, ...] = MODERN_POINT_ARTIFACT_BACKEND_SYSTEM_IDS
    backend_keys: tuple[str, ...] = MODERN_POINT_ARTIFACT_BACKEND_KEYS
    lane_id: str = "modern"


@dataclass(frozen=True)
class ModernPointArtifactVerificationReport:
    """Self-describing verifier result for one exported modern point artifact."""

    point_id: str
    point_label: str
    schema_id: str
    schema_version: int
    verdict_count: int
    verdict_system_ids: tuple[str, ...]
    verdict_observable_ids: tuple[str, ...]
    policy_system_ids: tuple[str, ...]
    header: ModernPointArtifactHeader
    contract: FrozenVerifierContract
    import_isolation: ImportIsolationReport
    issues: tuple[VerificationIssue, ...]

    @property
    def ok(self) -> bool:
        return not self.issues

    @property
    def issue_codes(self) -> tuple[str, ...]:
        return tuple(issue.code for issue in self.issues)

    def require_ok(self) -> "ModernPointArtifactVerificationReport":
        if self.issues:
            joined = "; ".join(f"{issue.code}: {issue.message}" for issue in self.issues)
            raise ArtifactSchemaError(joined)
        return self


@dataclass(frozen=True)
class FrozenBridgeVerifierContract:
    """Frozen bridge-sidecar contract enforced by the standalone verifier."""

    schema_id: str = MODERN_POINT_BRIDGE_ARTIFACT_SCHEMA_ID
    schema_version: int = MODERN_POINT_BRIDGE_ARTIFACT_SCHEMA_VERSION
    lane_id: str = "modern"
    input_bundle_schema_id: str = MODERN_BRIDGE_INPUTS_SCHEMA_ID
    coupling_schema_id: str = MODERN_BRIDGE_COUPLINGS_SCHEMA_ID
    matching_schema_id: str = MODERN_BRIDGE_MATCHING_SCHEMA_ID


@dataclass(frozen=True)
class ModernPointBridgeArtifactVerificationReport:
    """Self-describing verifier result for one exported modern bridge sidecar."""

    point_id: str
    point_label: str
    schema_id: str
    schema_version: int
    coupling_schema_id: str
    matching_schema_id: str
    system_ids: tuple[str, ...]
    contract: FrozenBridgeVerifierContract
    import_isolation: ImportIsolationReport
    issues: tuple[VerificationIssue, ...]

    @property
    def ok(self) -> bool:
        return not self.issues

    @property
    def issue_codes(self) -> tuple[str, ...]:
        return tuple(issue.code for issue in self.issues)

    def require_ok(self) -> "ModernPointBridgeArtifactVerificationReport":
        if self.issues:
            joined = "; ".join(f"{issue.code}: {issue.message}" for issue in self.issues)
            raise ArtifactSchemaError(joined)
        return self


@dataclass(frozen=True)
class FrozenPhenomenologyVerifierContract:
    """Frozen modern QS5 sidecar contract enforced by the standalone verifier."""

    schema_id: str = MODERN_POINT_PHENOMENOLOGY_ARTIFACT_SCHEMA_ID
    schema_version: int = MODERN_POINT_PHENOMENOLOGY_ARTIFACT_SCHEMA_VERSION
    lane_id: str = "modern"
    release_scope_id: str = MODERN_POINT_PHENOMENOLOGY_RELEASE_SCOPE_ID
    policy_schema_id: str = "quarkConstraints.modern.phenomenology.v1"
    policy_id: str = "quarkConstraints.modern.phenomenology.policy.v1"
    policy_system_ids: tuple[str, ...] = ("epsilon_K", "K", "B_d", "B_s", "D0")
    non_cp_acceptance_system_ids: tuple[str, ...] = (
        MODERN_POINT_PHENOMENOLOGY_NON_CP_ACCEPTANCE_SYSTEM_IDS
    )
    system_result_schema_id: str = MODERN_POINT_PHENOMENOLOGY_SYSTEM_RESULT_SCHEMA_ID


@dataclass(frozen=True)
class ModernPointPhenomenologyArtifactVerificationReport:
    """Self-describing verifier result for one exported modern QS5 sidecar."""

    point_id: str
    point_label: str
    schema_id: str
    schema_version: int
    system_ids: tuple[str, ...]
    non_cp_acceptance_system_ids: tuple[str, ...]
    contract: FrozenPhenomenologyVerifierContract
    import_isolation: ImportIsolationReport
    issues: tuple[VerificationIssue, ...]

    @property
    def ok(self) -> bool:
        return not self.issues

    @property
    def issue_codes(self) -> tuple[str, ...]:
        return tuple(issue.code for issue in self.issues)

    def require_ok(self) -> "ModernPointPhenomenologyArtifactVerificationReport":
        if self.issues:
            joined = "; ".join(f"{issue.code}: {issue.message}" for issue in self.issues)
            raise ArtifactSchemaError(joined)
        return self


def _static_import_violations() -> tuple[str, ...]:
    source_path = Path(__file__)
    source = source_path.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(source_path))
    violations: list[str] = []

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                module_name = alias.name
                if module_name == "quarkConstraints" or module_name.startswith("quarkConstraints."):
                    if module_name not in ALLOWED_ABSOLUTE_LOCAL_IMPORTS:
                        violations.append(f"import {module_name}")
        elif isinstance(node, ast.ImportFrom):
            module_name = node.module or ""
            if node.level == 1:
                if module_name not in ALLOWED_RELATIVE_IMPORTS:
                    violations.append(f"from .{module_name} import ...")
            elif node.level > 1:
                prefix = "." * node.level
                violations.append(f"from {prefix}{module_name} import ...")
            elif module_name == "quarkConstraints" or module_name.startswith("quarkConstraints."):
                if module_name not in ALLOWED_ABSOLUTE_LOCAL_IMPORTS:
                    violations.append(f"from {module_name} import ...")

    return tuple(sorted(set(violations)))


def _runtime_import_violations() -> tuple[str, ...]:
    allowed_modules = {
        __name__,
        __package__ or "",
        "quarkConstraints",
        *ALLOWED_ABSOLUTE_LOCAL_IMPORTS,
        "quarkConstraints.modern.conventions",
        "quarkConstraints.modern.inputs",
        "quarkConstraints.modern.phenomenology",
    }
    violations: set[str] = set()

    for module_name in sys.modules:
        if not module_name:
            continue
        if any(
            module_name == prefix or module_name.startswith(f"{prefix}.")
            for prefix in FORBIDDEN_EXTERNAL_MODULE_PREFIXES
        ):
            violations.add(module_name)
            continue
        if module_name.startswith("quarkConstraints.modern."):
            if module_name in allowed_modules:
                continue
            if module_name == __name__ or module_name.startswith(f"{__name__}."):
                continue
            violations.add(module_name)

    return tuple(sorted(violations))


def _build_import_isolation_report() -> ImportIsolationReport:
    unexpected_import_targets = _static_import_violations()
    forbidden_loaded_modules = _runtime_import_violations()
    return ImportIsolationReport(
        verifier_module=__name__,
        verifier_source_path=str(Path(__file__).resolve()),
        allowed_relative_imports=ALLOWED_RELATIVE_IMPORTS,
        allowed_absolute_imports=ALLOWED_ABSOLUTE_LOCAL_IMPORTS,
        static_ok=not unexpected_import_targets,
        runtime_ok=not forbidden_loaded_modules,
        unexpected_import_targets=unexpected_import_targets,
        forbidden_loaded_modules=forbidden_loaded_modules,
    )


def _verify_import_isolation(issues: list[VerificationIssue], report: ImportIsolationReport) -> None:
    if not report.static_ok:
        _append_issue(
            issues,
            code="import_isolation_static_violation",
            message=(
                "verifier source imports modules outside the allowed artifact-only "
                f"boundary: {', '.join(report.unexpected_import_targets)}"
            ),
        )
    if not report.runtime_ok:
        _append_issue(
            issues,
            code="import_isolation_runtime_violation",
            message=(
                "runtime loaded forbidden modules while executing the standalone "
                f"modern verifier: {', '.join(report.forbidden_loaded_modules)}"
            ),
        )
    if not report.static_ok or not report.runtime_ok:
        _append_issue(
            issues,
            code="import_isolation_failed",
            message="standalone modern verifier import isolation failed",
        )


def verify_artifact(artifact: ModernPointArtifactV1) -> ModernPointArtifactVerificationReport:
    """Validate exported-field self-consistency, not physics correctness."""

    contract = FrozenVerifierContract()
    import_isolation = _build_import_isolation_report()
    issues: list[VerificationIssue] = []

    if not isinstance(artifact, ModernPointArtifactV1):
        raise ArtifactSchemaError("artifact must be a ModernPointArtifactV1 instance")

    _verify_import_isolation(issues, import_isolation)

    if artifact.schema_id != contract.schema_id:
        _append_issue(
            issues,
            code="artifact_schema_id_mismatch",
            message=f"artifact schema_id must be exactly {contract.schema_id!r}",
            subject="schema_id",
        )
    if artifact.schema_version != contract.schema_version:
        _append_issue(
            issues,
            code="artifact_schema_version_mismatch",
            message=f"artifact schema_version must be exactly {contract.schema_version}",
            subject="schema_version",
        )

    header = artifact.header
    if header.schema_id != contract.header_schema_id:
        _append_issue(
            issues,
            code="header_schema_id_mismatch",
            message=f"header.schema_id must be exactly {contract.header_schema_id!r}",
            subject="header.schema_id",
        )
    if header.lane_id != contract.lane_id:
        _append_issue(
            issues,
            code="header_lane_id_mismatch",
            message=f"header.lane_id must be exactly {contract.lane_id!r}",
            subject="header.lane_id",
        )
    if header.policy_schema_id != contract.policy_schema_id:
        _append_issue(
            issues,
            code="header_policy_schema_id_mismatch",
            message=f"header.policy_schema_id must be exactly {contract.policy_schema_id!r}",
            subject="header.policy_schema_id",
        )
    if header.policy_id != contract.policy_id:
        _append_issue(
            issues,
            code="header_policy_id_mismatch",
            message=f"header.policy_id must be exactly {contract.policy_id!r}",
            subject="header.policy_id",
        )
    if header.evaluation_schema_id != contract.evaluation_schema_id:
        _append_issue(
            issues,
            code="header_evaluation_schema_id_mismatch",
            message=f"header.evaluation_schema_id must be exactly {contract.evaluation_schema_id!r}",
            subject="header.evaluation_schema_id",
        )
    if header.verdict_schema_id != contract.verdict_schema_id:
        _append_issue(
            issues,
            code="header_verdict_schema_id_mismatch",
            message=f"header.verdict_schema_id must be exactly {contract.verdict_schema_id!r}",
            subject="header.verdict_schema_id",
        )
    if header.policy_system_ids != contract.policy_system_ids:
        _append_issue(
            issues,
            code="header_policy_coverage_mismatch",
            message="header.policy_system_ids must explicitly cover epsilon_K, K, B_d, B_s, and D0",
            subject="header.policy_system_ids",
        )
    if header.verdict_system_ids != contract.verdict_system_ids:
        _append_issue(
            issues,
            code="header_verdict_coverage_mismatch",
            message="header.verdict_system_ids must explicitly cover K, B_d, B_s, and D0",
            subject="header.verdict_system_ids",
        )
    if header.verdict_observable_ids != contract.verdict_observable_ids:
        _append_issue(
            issues,
            code="header_observable_coverage_mismatch",
            message=(
                "header.verdict_observable_ids must explicitly cover epsilon_K, B_d, B_s, and D0"
            ),
            subject="header.verdict_observable_ids",
        )
    if header.verdict_count != len(contract.verdict_system_ids):
        _append_issue(
            issues,
            code="header_verdict_count_mismatch",
            message="header.verdict_count must be 4 for the frozen modern point artifact",
            subject="header.verdict_count",
        )
    if not math.isfinite(header.M_KK) or header.M_KK <= 0.0:
        _append_issue(
            issues,
            code="header_mkk_invalid",
            message="header.M_KK must be a positive finite float",
            subject="header.M_KK",
        )
    if not math.isfinite(header.xi_KK) or header.xi_KK <= 0.0:
        _append_issue(
            issues,
            code="header_xikk_invalid",
            message="header.xi_KK must be a positive finite float",
            subject="header.xi_KK",
        )

    if artifact.policy.schema_id != contract.policy_schema_id:
        _append_issue(
            issues,
            code="policy_schema_id_mismatch",
            message=f"policy.schema_id must be exactly {contract.policy_schema_id!r}",
            subject="policy.schema_id",
        )
    if artifact.policy.policy_id != contract.policy_id:
        _append_issue(
            issues,
            code="policy_id_mismatch",
            message=f"policy.policy_id must be exactly {contract.policy_id!r}",
            subject="policy.policy_id",
        )
    if artifact.policy.system_ids != contract.policy_system_ids:
        _append_issue(
            issues,
            code="policy_coverage_mismatch",
            message="policy.system_ids must explicitly cover epsilon_K, K, B_d, B_s, and D0",
            subject="policy.system_ids",
        )

    verdicts = artifact.verdicts
    if len(verdicts) != len(contract.verdict_system_ids):
        _append_issue(
            issues,
            code="verdict_count_mismatch",
            message="artifact must contain exactly four verdicts",
            subject="verdicts",
        )

    observed_system_ids = tuple(verdict.system_id for verdict in verdicts)
    observed_observable_ids = tuple(verdict.observable_id for verdict in verdicts)
    expected_index_by_system = {system_id: index for index, system_id in enumerate(contract.verdict_system_ids)}
    if observed_system_ids != contract.verdict_system_ids:
        _append_issue(
            issues,
            code="verdict_system_coverage_mismatch",
            message="verdicts must preserve the frozen K/B_d/B_s/D0 order",
            subject="verdicts.system_id",
        )
    if observed_observable_ids != contract.verdict_observable_ids:
        _append_issue(
            issues,
            code="verdict_observable_coverage_mismatch",
            message="verdicts must preserve the frozen epsilon_K/B_d/B_s/D0 observable order",
            subject="verdicts.observable_id",
        )

    for verdict in verdicts:
        if verdict.schema_id != contract.verdict_schema_id:
            _append_issue(
                issues,
                code="verdict_schema_id_mismatch",
                message=f"verdict.schema_id must be exactly {contract.verdict_schema_id!r}",
                subject=verdict.system_id,
            )
            continue
        try:
            policy_system = artifact.policy.system_policy(verdict.policy_system_id)
        except ValueError as exc:
            _append_issue(
                issues,
                code="verdict_policy_system_id_mismatch",
                message=str(exc),
                subject=verdict.system_id,
            )
            continue
        system_index = expected_index_by_system.get(verdict.system_id)
        if system_index is None:
            _append_issue(
                issues,
                code="verdict_system_id_mismatch",
                message="verdict.system_id must be one of the frozen K/B_d/B_s/D0 systems",
                subject=verdict.system_id,
            )
            continue
        if verdict.policy_id != policy_system.policy_id:
            _append_issue(
                issues,
                code="verdict_policy_id_mismatch",
                message="verdict.policy_id must match the exported policy entry",
                subject=verdict.system_id,
            )
        if verdict.policy_display_name != policy_system.display_name:
            _append_issue(
                issues,
                code="verdict_policy_display_name_mismatch",
                message="verdict.policy_display_name must match the exported policy entry",
                subject=verdict.system_id,
            )
        if verdict.policy_notes != policy_system.notes:
            _append_issue(
                issues,
                code="verdict_policy_notes_mismatch",
                message="verdict.policy_notes must match the exported policy entry",
                subject=verdict.system_id,
            )
        if verdict.backend_system_id != contract.backend_system_ids[
            system_index
        ]:
            _append_issue(
                issues,
                code="verdict_backend_system_id_mismatch",
                message="verdict.backend_system_id must match the frozen backend mapping",
                subject=verdict.system_id,
            )
        if verdict.backend_key != contract.backend_keys[
            system_index
        ]:
            _append_issue(
                issues,
                code="verdict_backend_key_mismatch",
                message="verdict.backend_key must match the frozen backend mapping",
                subject=verdict.system_id,
            )
        if verdict.bound <= 0.0 or not math.isfinite(verdict.bound):
            _append_issue(
                issues,
                code="verdict_bound_invalid",
                message="verdict.bound must be positive and finite",
                subject=verdict.system_id,
            )
        if verdict.ratio_to_bound < 0.0 or not math.isfinite(verdict.ratio_to_bound):
            _append_issue(
                issues,
                code="verdict_ratio_invalid",
                message="verdict.ratio_to_bound must be non-negative and finite",
                subject=verdict.system_id,
            )
        if not math.isfinite(verdict.effective_amplitude):
            _append_issue(
                issues,
                code="verdict_effective_amplitude_invalid",
                message="verdict.effective_amplitude must be finite",
                subject=verdict.system_id,
            )
        if not math.isfinite(verdict.coherent_amplitude):
            _append_issue(
                issues,
                code="verdict_coherent_amplitude_invalid",
                message="verdict.coherent_amplitude must be finite",
                subject=verdict.system_id,
            )
        if verdict.dominant_operator_size < 0.0 or not math.isfinite(verdict.dominant_operator_size):
            _append_issue(
                issues,
                code="verdict_dominant_operator_size_invalid",
                message="verdict.dominant_operator_size must be non-negative and finite",
                subject=verdict.system_id,
            )
        if not hasattr(verdict.weighted_operator_sizes, "items"):
            _append_issue(
                issues,
                code="verdict_weighted_operator_sizes_invalid",
                message="verdict.weighted_operator_sizes must be a mapping",
                subject=verdict.system_id,
            )
            continue
        invalid_weighted_values = False
        weighted_operator_sizes = {
            str(name): float(value)
            for name, value in verdict.weighted_operator_sizes.items()
        }
        if not weighted_operator_sizes:
            _append_issue(
                issues,
                code="verdict_weighted_operator_sizes_empty",
                message="verdict.weighted_operator_sizes must contain at least one operator",
                subject=verdict.system_id,
            )
            continue
        for operator_name, operator_size in weighted_operator_sizes.items():
            if not math.isfinite(operator_size) or operator_size < 0.0:
                _append_issue(
                    issues,
                    code="verdict_weighted_operator_size_invalid",
                    message=(
                        f"verdict.weighted_operator_sizes[{operator_name!r}] must be finite "
                        "and non-negative"
                    ),
                    subject=verdict.system_id,
                )
                invalid_weighted_values = True
        if invalid_weighted_values:
            continue
        if verdict.dominant_operator not in weighted_operator_sizes:
            _append_issue(
                issues,
                code="verdict_dominant_operator_missing",
                message="verdict.dominant_operator must be present in weighted_operator_sizes",
                subject=verdict.system_id,
            )
            continue
        dominant_size = weighted_operator_sizes[verdict.dominant_operator]
        max_size = max(weighted_operator_sizes.values())
        if not math.isclose(
            verdict.dominant_operator_size,
            dominant_size,
            rel_tol=0.0,
            abs_tol=0.0,
        ):
            _append_issue(
                issues,
                code="verdict_dominant_operator_size_mismatch",
                message="verdict.dominant_operator_size must match weighted_operator_sizes[dominant_operator]",
                subject=verdict.system_id,
            )
        if not math.isclose(dominant_size, max_size, rel_tol=0.0, abs_tol=0.0):
            _append_issue(
                issues,
                code="verdict_dominant_operator_not_maximum",
                message="verdict.dominant_operator must name an operator with the maximum size",
                subject=verdict.system_id,
            )
        if not math.isclose(
            verdict.dominant_operator_size,
            max_size,
            rel_tol=0.0,
            abs_tol=0.0,
        ):
            _append_issue(
                issues,
                code="verdict_dominant_operator_size_not_maximum",
                message="verdict.dominant_operator_size must equal the maximum operator size",
                subject=verdict.system_id,
            )
        if verdict.passes != (verdict.ratio_to_bound <= 1.0):
            _append_issue(
                issues,
                code="verdict_pass_ratio_mismatch",
                message="verdict.passes must agree with verdict.ratio_to_bound <= 1.0",
                subject=verdict.system_id,
            )

    return ModernPointArtifactVerificationReport(
        point_id=artifact.header.point_id,
        point_label=artifact.header.point_label,
        schema_id=artifact.schema_id,
        schema_version=artifact.schema_version,
        verdict_count=len(verdicts),
        verdict_system_ids=observed_system_ids,
        verdict_observable_ids=observed_observable_ids,
        policy_system_ids=artifact.policy.system_ids,
        header=header,
        contract=contract,
        import_isolation=import_isolation,
        issues=tuple(issues),
    )


def verify_artifact_path(path: str | Path) -> ModernPointArtifactVerificationReport:
    """Read and verify one exported modern point artifact from disk."""

    try:
        artifact = read_modern_point_artifact(path)
    except Exception as exc:  # pragma: no cover - defensive wrapper
        raise ArtifactSchemaError(str(exc)) from exc
    return verify_artifact(artifact)


def verify_bridge_artifact(
    artifact: ModernPointBridgeArtifactV1,
) -> ModernPointBridgeArtifactVerificationReport:
    """Validate that one exported modern bridge sidecar is self-consistent."""

    contract = FrozenBridgeVerifierContract()
    import_isolation = _build_import_isolation_report()
    issues: list[VerificationIssue] = []

    if not isinstance(artifact, ModernPointBridgeArtifactV1):
        raise ArtifactSchemaError("artifact must be a ModernPointBridgeArtifactV1 instance")

    _verify_import_isolation(issues, import_isolation)

    if artifact.schema_id != contract.schema_id:
        _append_issue(
            issues,
            code="bridge_schema_id_mismatch",
            message=f"artifact schema_id must be exactly {contract.schema_id!r}",
            subject="schema_id",
        )
    if artifact.schema_version != contract.schema_version:
        _append_issue(
            issues,
            code="bridge_schema_version_mismatch",
            message=f"artifact schema_version must be exactly {contract.schema_version}",
            subject="schema_version",
        )
    if artifact.lane_id != contract.lane_id:
        _append_issue(
            issues,
            code="bridge_lane_id_mismatch",
            message=f"artifact lane_id must be exactly {contract.lane_id!r}",
            subject="lane_id",
        )
    if artifact.input_bundle_schema_id != contract.input_bundle_schema_id:
        _append_issue(
            issues,
            code="bridge_input_bundle_schema_id_mismatch",
            message=(
                "artifact input_bundle_schema_id must be exactly "
                f"{contract.input_bundle_schema_id!r}"
            ),
            subject="input_bundle_schema_id",
        )
    if artifact.coupling_schema_id != contract.coupling_schema_id:
        _append_issue(
            issues,
            code="bridge_coupling_schema_id_mismatch",
            message=f"artifact coupling_schema_id must be exactly {contract.coupling_schema_id!r}",
            subject="coupling_schema_id",
        )
    if artifact.matching_schema_id != contract.matching_schema_id:
        _append_issue(
            issues,
            code="bridge_matching_schema_id_mismatch",
            message=f"artifact matching_schema_id must be exactly {contract.matching_schema_id!r}",
            subject="matching_schema_id",
        )

    couplings = artifact.couplings
    matching = artifact.matching
    if couplings.schema_id != artifact.coupling_schema_id:
        _append_issue(
            issues,
            code="bridge_nested_couplings_schema_id_mismatch",
            message="nested couplings.schema_id must match coupling_schema_id",
            subject="couplings.schema_id",
        )
    if matching.schema_id != artifact.matching_schema_id:
        _append_issue(
            issues,
            code="bridge_nested_matching_schema_id_mismatch",
            message="nested matching.schema_id must match matching_schema_id",
            subject="matching.schema_id",
        )
    if matching.couplings_schema_id != artifact.coupling_schema_id:
        _append_issue(
            issues,
            code="bridge_matching_couplings_schema_id_mismatch",
            message="nested matching.couplings_schema_id must match coupling_schema_id",
            subject="matching.couplings_schema_id",
        )
    if couplings.point_id != artifact.point_id or matching.point_id != artifact.point_id:
        _append_issue(
            issues,
            code="bridge_point_id_mismatch",
            message="artifact point_id must match nested couplings and matching point_id",
            subject="point_id",
        )
    if couplings.point_label != artifact.point_label or matching.point_label != artifact.point_label:
        _append_issue(
            issues,
            code="bridge_point_label_mismatch",
            message="artifact point_label must match nested couplings and matching point_label",
            subject="point_label",
        )
    if couplings.input_bundle_id != artifact.input_bundle_id or matching.input_bundle_id != artifact.input_bundle_id:
        _append_issue(
            issues,
            code="bridge_input_bundle_id_mismatch",
            message="artifact input_bundle_id must match nested couplings and matching",
            subject="input_bundle_id",
        )
    if couplings.input_provenance_id != artifact.input_provenance_id or matching.input_provenance_id != artifact.input_provenance_id:
        _append_issue(
            issues,
            code="bridge_input_provenance_id_mismatch",
            message="artifact input_provenance_id must match nested couplings and matching",
            subject="input_provenance_id",
        )
    if couplings.input_resolution_policy_id != artifact.input_resolution_policy_id or matching.input_resolution_policy_id != artifact.input_resolution_policy_id:
        _append_issue(
            issues,
            code="bridge_input_resolution_policy_id_mismatch",
            message="artifact input_resolution_policy_id must match nested couplings and matching",
            subject="input_resolution_policy_id",
        )
    if couplings.qcd_metadata_id != artifact.qcd_metadata_id or matching.qcd_metadata_id != artifact.qcd_metadata_id:
        _append_issue(
            issues,
            code="bridge_qcd_metadata_id_mismatch",
            message="artifact qcd_metadata_id must match nested couplings and matching",
            subject="qcd_metadata_id",
        )
    if couplings.alpha_s_policy_id != artifact.alpha_s_policy_id or matching.alpha_s_policy_id != artifact.alpha_s_policy_id:
        _append_issue(
            issues,
            code="bridge_alpha_s_policy_id_mismatch",
            message="artifact alpha_s_policy_id must match nested couplings and matching",
            subject="alpha_s_policy_id",
        )
    if not math.isclose(
        float(couplings.M_KK),
        float(matching.M_KK),
        rel_tol=0.0,
        abs_tol=0.0,
    ):
        _append_issue(
            issues,
            code="bridge_mkk_mismatch",
            message="nested couplings.M_KK and matching.M_KK must agree",
            subject="M_KK",
        )
    if not math.isclose(
        float(couplings.xi_KK),
        float(matching.xi_KK),
        rel_tol=0.0,
        abs_tol=0.0,
    ):
        _append_issue(
            issues,
            code="bridge_xi_kk_mismatch",
            message="nested couplings.xi_KK and matching.xi_KK must agree",
            subject="xi_KK",
        )

    systems = matching.system_ids
    if systems != ("K", "B_d", "B_s", "D0"):
        _append_issue(
            issues,
            code="bridge_matching_system_coverage_mismatch",
            message="nested matching.system_matches must preserve the frozen K/B_d/B_s/D0 order",
            subject="matching.system_matches",
        )
    return ModernPointBridgeArtifactVerificationReport(
        point_id=artifact.point_id,
        point_label=artifact.point_label,
        schema_id=artifact.schema_id,
        schema_version=artifact.schema_version,
        coupling_schema_id=artifact.coupling_schema_id,
        matching_schema_id=artifact.matching_schema_id,
        system_ids=systems,
        contract=contract,
        import_isolation=import_isolation,
        issues=tuple(issues),
    )


def verify_bridge_artifact_path(path: str | Path) -> ModernPointBridgeArtifactVerificationReport:
    """Read and verify one exported modern bridge sidecar artifact from disk."""

    try:
        artifact = read_modern_point_bridge_artifact(path)
    except Exception as exc:  # pragma: no cover - defensive wrapper
        raise ArtifactSchemaError(str(exc)) from exc
    return verify_bridge_artifact(artifact)


def verify_phenomenology_artifact(
    artifact: ModernPointPhenomenologyArtifactV1,
) -> ModernPointPhenomenologyArtifactVerificationReport:
    """Validate QS5 sidecar self-consistency, not physics correctness."""

    contract = FrozenPhenomenologyVerifierContract()
    import_isolation = _build_import_isolation_report()
    issues: list[VerificationIssue] = []

    if not isinstance(artifact, ModernPointPhenomenologyArtifactV1):
        raise ArtifactSchemaError(
            "artifact must be a ModernPointPhenomenologyArtifactV1 instance"
        )

    _verify_import_isolation(issues, import_isolation)

    if artifact.schema_id != contract.schema_id:
        _append_issue(
            issues,
            code="phenomenology_schema_id_mismatch",
            message=f"artifact schema_id must be exactly {contract.schema_id!r}",
            subject="schema_id",
        )
    if artifact.schema_version != contract.schema_version:
        _append_issue(
            issues,
            code="phenomenology_schema_version_mismatch",
            message=f"artifact schema_version must be exactly {contract.schema_version}",
            subject="schema_version",
        )
    if artifact.lane_id != contract.lane_id:
        _append_issue(
            issues,
            code="phenomenology_lane_id_mismatch",
            message=f"artifact lane_id must be exactly {contract.lane_id!r}",
            subject="lane_id",
        )
    if artifact.release_scope_id != contract.release_scope_id:
        _append_issue(
            issues,
            code="phenomenology_release_scope_id_mismatch",
            message=(
                "artifact release_scope_id must be exactly "
                f"{contract.release_scope_id!r}"
            ),
            subject="release_scope_id",
        )
    if artifact.policy_schema_id != contract.policy_schema_id:
        _append_issue(
            issues,
            code="phenomenology_policy_schema_id_mismatch",
            message=(
                "artifact policy_schema_id must be exactly "
                f"{contract.policy_schema_id!r}"
            ),
            subject="policy_schema_id",
        )
    if artifact.policy_id != contract.policy_id:
        _append_issue(
            issues,
            code="phenomenology_policy_id_mismatch",
            message=f"artifact policy_id must be exactly {contract.policy_id!r}",
            subject="policy_id",
        )
    if artifact.policy_system_ids != contract.policy_system_ids:
        _append_issue(
            issues,
            code="phenomenology_policy_system_ids_mismatch",
            message="artifact must preserve the frozen epsilon_K/K/B_d/B_s/D0 policy order",
            subject="policy_system_ids",
        )
    if artifact.non_cp_acceptance_system_ids != contract.non_cp_acceptance_system_ids:
        _append_issue(
            issues,
            code="phenomenology_non_cp_acceptance_system_ids_mismatch",
            message=(
                "artifact non_cp_acceptance_system_ids must be exactly "
                f"{contract.non_cp_acceptance_system_ids!r}"
            ),
            subject="non_cp_acceptance_system_ids",
        )
    if not artifact.kaon_viability_claimed:
        _append_issue(
            issues,
            code="phenomenology_kaon_viability_claim_missing",
            message="artifact must claim kaon viability in the full Delta F = 2 release scope",
            subject="kaon_viability_claimed",
        )
    if not math.isfinite(artifact.M_KK) or artifact.M_KK <= 0.0:
        _append_issue(
            issues,
            code="phenomenology_mkk_invalid",
            message="artifact M_KK must be a positive finite float",
            subject="M_KK",
        )
    if not math.isfinite(artifact.xi_KK) or artifact.xi_KK <= 0.0:
        _append_issue(
            issues,
            code="phenomenology_xi_kk_invalid",
            message="artifact xi_KK must be a positive finite float",
            subject="xi_KK",
        )

    observed_system_ids = tuple(result.system_id for result in artifact.system_results)
    if observed_system_ids != contract.policy_system_ids:
        _append_issue(
            issues,
            code="phenomenology_system_result_order_mismatch",
            message="system_results must preserve the frozen epsilon_K/K/B_d/B_s/D0 order",
            subject="system_results",
        )

    failing_non_cp: list[str] = []
    for result in artifact.system_results:
        if result.schema_id != contract.system_result_schema_id:
            _append_issue(
                issues,
                code="phenomenology_system_result_schema_id_mismatch",
                message=(
                    "system_result.schema_id must be exactly "
                    f"{contract.system_result_schema_id!r}"
                ),
                subject=result.system_id,
            )
            continue
        if result.policy_system_id != result.system_id:
            _append_issue(
                issues,
                code="phenomenology_policy_system_id_mismatch",
                message="system_result.policy_system_id must equal system_result.system_id",
                subject=result.system_id,
            )
        expected_treatment = MODERN_POINT_PHENOMENOLOGY_SYSTEM_TREATMENT_IDS[result.system_id]
        if result.treatment_id != expected_treatment:
            _append_issue(
                issues,
                code="phenomenology_treatment_id_mismatch",
                message=(
                    "system_result.treatment_id must match the frozen modern "
                    f"treatment for {result.system_id}"
                ),
                subject=result.system_id,
            )
        if not result.evaluated_from_bridge:
            _append_issue(
                issues,
                code="phenomenology_bridge_evaluation_missing",
                message="all systems must be sourced from the exported bridge artifact",
                subject=result.system_id,
            )
            continue
        if not result.included_in_non_cp_acceptance:
            _append_issue(
                issues,
                code="phenomenology_acceptance_system_not_included",
                message="all systems must be included in acceptance",
                subject=result.system_id,
            )
        if result.passes is False and result.included_in_non_cp_acceptance:
            failing_non_cp.append(result.system_id)

    if tuple(failing_non_cp) != artifact.failing_non_cp_system_ids:
        _append_issue(
            issues,
            code="phenomenology_failing_non_cp_system_ids_mismatch",
            message="artifact failing_non_cp_system_ids must match the failing acceptance systems",
            subject="failing_non_cp_system_ids",
        )
    if artifact.non_cp_passes != (not failing_non_cp):
        _append_issue(
            issues,
            code="phenomenology_non_cp_passes_mismatch",
            message="artifact non_cp_passes must agree with the acceptance-bearing system results",
            subject="non_cp_passes",
        )

    return ModernPointPhenomenologyArtifactVerificationReport(
        point_id=artifact.point_id,
        point_label=artifact.point_label,
        schema_id=artifact.schema_id,
        schema_version=artifact.schema_version,
        system_ids=observed_system_ids,
        non_cp_acceptance_system_ids=artifact.non_cp_acceptance_system_ids,
        contract=contract,
        import_isolation=import_isolation,
        issues=tuple(issues),
    )


def verify_phenomenology_artifact_path(
    path: str | Path,
) -> ModernPointPhenomenologyArtifactVerificationReport:
    """Read and verify one exported modern QS5 sidecar artifact from disk."""

    try:
        artifact = read_modern_point_phenomenology_artifact(path)
    except Exception as exc:  # pragma: no cover - defensive wrapper
        raise ArtifactSchemaError(str(exc)) from exc
    return verify_phenomenology_artifact(artifact)


__all__ = [
    "ArtifactSchemaError",
    "FrozenBridgeVerifierContract",
    "FrozenPhenomenologyVerifierContract",
    "FrozenVerifierContract",
    "ImportIsolationReport",
    "ModernPointBridgeArtifactVerificationReport",
    "ModernPointPhenomenologyArtifactVerificationReport",
    "ModernPointArtifactVerificationReport",
    "VerificationIssue",
    "verify_artifact",
    "verify_artifact_path",
    "verify_bridge_artifact",
    "verify_bridge_artifact_path",
    "verify_phenomenology_artifact",
    "verify_phenomenology_artifact_path",
]
