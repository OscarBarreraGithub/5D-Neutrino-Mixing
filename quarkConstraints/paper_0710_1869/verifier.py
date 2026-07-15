"""Artifact-only verification helpers for ``paper_0710_1869`` exports."""

from __future__ import annotations

import ast
import math
import sys
from dataclasses import dataclass
from pathlib import Path

from .artifacts import (
    ArtifactSchemaError,
    ComplexValue,
    HadronicArtifactBundleV1,
    ObservableArtifactBundleV1,
    ObservableRecord,
    ProvenanceBundleV1,
    WilsonArtifactBundleV1,
    read_artifact,
)
from .conventions import PAPER_0710_1869_MODE_ID, PAPER_0710_1869_PAPER_ID

PACKAGE_ROOT = "quarkConstraints.paper_0710_1869"
ALLOWED_RELATIVE_IMPORTS = ("artifacts", "conventions")
ALLOWED_ABSOLUTE_LOCAL_IMPORTS = (
    f"{PACKAGE_ROOT}.artifacts",
    f"{PACKAGE_ROOT}.conventions",
)
FORBIDDEN_EXTERNAL_MODULE_PREFIXES = (
    "quarkConstraints.deltaf2",
    "quarkConstraints.couplings",
    "quarkConstraints.model",
)

EXPECTED_OPERATOR_BASIS = "kk_gluon_tree_np_only.v1"
EXPECTED_RENORMALIZATION_SCHEME = "bmu.hep-ph-0005183.ndr-ms.lo.v1"
EXPECTED_INTERPRETATION = "kaon.np_only.v1"
EXPECTED_SYSTEM = "kaon"
EXPECTED_SECTOR = "down"
EXPECTED_REQUIRED_OBSERVABLES = ("M12_K_NP.re", "M12_K_NP.im", "Delta_m_K_NP")
EXPECTED_SUPPORTED_OPERATORS = ("Q1_VLL", "Q1_VRR", "Q4_LR", "Q5_LR")
EXPECTED_UNSUPPORTED_OPERATORS: tuple[str, ...] = ()
EXPECTED_HAMILTONIAN_CONVENTION_ID = "heff.sum_ci_qi.no_hc_factor.v1"
EXPECTED_MATRIX_ELEMENT_FORMULA_ID = (
    "kaon.q1_vll_vrr.2over3_fk2_mk2_bk_mu.plpr_projectors.v2"
)
EXPECTED_PARITY_RELATION_ID = "kaon.q1_vll_equals_q1_vrr.by_parity.v1"
EXPECTED_OBSERVABLE_UNITS = "GeV"
EXPECTED_BAG_PARAMETER_SOURCE_SCHEME_ID = "rgi.flag21.average.v1"
EXPECTED_BAG_PARAMETER_TRANSFORMATION_ID = "hat_bk_to_bk_mu_had.lo_inverse_q1_running.v1"
EXPECTED_ALPHA_S_POLICY_ID = "qcd.alpha_s.lo.beta0.thresholds_continuous.v1"
EXPECTED_PROVENANCE_MODE_ID = "paper_sourced_defaults.v1"
DEFAULT_THRESHOLD_MASSES_GEV = (
    (1.27, 4),
    (4.18, 5),
    (163.5, 6),
)
FROZEN_PUBLIC_TOLERANCES = {"abs": 1.0e-18, "rel": 0.0}


def frozen_public_verifier_tolerances() -> dict[str, float]:
    """Return the benchmark-facing observable tolerance contract."""

    return dict(FROZEN_PUBLIC_TOLERANCES)


@dataclass(frozen=True)
class VerificationIssue:
    """One verification failure."""

    code: str
    message: str
    subject: str | None = None


@dataclass(frozen=True)
class VerifierInputSet:
    """Verifier input restricted to exported artifact structures."""

    wilson_bundle: WilsonArtifactBundleV1
    hadronic_bundle: HadronicArtifactBundleV1
    observable_bundle: ObservableArtifactBundleV1
    provenance_bundle: ProvenanceBundleV1


@dataclass(frozen=True)
class TolerancePolicy:
    """Explicit numeric tolerance policy for artifact-only checks."""

    scale_abs_tol: float = 0.0
    matrix_element_abs_tol: float = 1e-18
    observable_abs_tol: float = FROZEN_PUBLIC_TOLERANCES["abs"]
    observable_rel_tol: float = FROZEN_PUBLIC_TOLERANCES["rel"]
    unsupported_lr_abs_tol: float = 1e-30
    description: str = (
        "Scales use exact agreement; reconstructed matrix elements and observables "
        "use the frozen public absolute-only tolerance contract; LR Wilson "
        "coefficients are part of the supported Wilson surface."
    )


@dataclass(frozen=True)
class FrozenVerifierContract:
    """Frozen theory contract enforced by the standalone verifier."""

    operator_basis: str = EXPECTED_OPERATOR_BASIS
    renormalization_scheme: str = EXPECTED_RENORMALIZATION_SCHEME
    interpretation: str = EXPECTED_INTERPRETATION
    system: str = EXPECTED_SYSTEM
    sector: str = EXPECTED_SECTOR
    supported_operators: tuple[str, ...] = EXPECTED_SUPPORTED_OPERATORS
    unsupported_operators: tuple[str, ...] = EXPECTED_UNSUPPORTED_OPERATORS
    required_observables: tuple[str, ...] = EXPECTED_REQUIRED_OBSERVABLES
    observable_units: str = EXPECTED_OBSERVABLE_UNITS
    hamiltonian_convention_id: str = EXPECTED_HAMILTONIAN_CONVENTION_ID
    matrix_element_formula_id: str = EXPECTED_MATRIX_ELEMENT_FORMULA_ID
    parity_relation_id: str = EXPECTED_PARITY_RELATION_ID
    bag_parameter_source_scheme_id: str = EXPECTED_BAG_PARAMETER_SOURCE_SCHEME_ID
    bag_parameter_transformation_id: str = EXPECTED_BAG_PARAMETER_TRANSFORMATION_ID
    alpha_s_policy_id: str = EXPECTED_ALPHA_S_POLICY_ID
    provenance_mode_id: str = EXPECTED_PROVENANCE_MODE_ID
    mode_id: str = PAPER_0710_1869_MODE_ID
    paper_id: str = PAPER_0710_1869_PAPER_ID
    h_eff_formula: str = "H_eff = sum_i C_i Q_i"
    h_eff_notes: str = "No additional hermitian-conjugate factor is applied."
    m12_formula: str = "M12 = <H_eff> / (2 m_K)"
    delta_m_formula: str = "Delta_m = 2 Re(M12)"
    q1_matrix_element_formula: str = "<Q1> = (2/3) f_K^2 m_K^2 B_K(mu_had)"
    bag_parameter_formula: str = (
        "B_K(mu_had) = hat_B_K * alpha_s(mu_had)^(4 / (2 beta_0(n_f(mu_had))))"
    )


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
class NumericCheckResult:
    """One numeric comparison performed by the verifier."""

    name: str
    issue_code: str
    expected: float
    actual: float | None
    present: bool
    expected_units: str
    actual_units: str | None
    abs_diff: float | None
    rel_diff: float | None
    abs_tol: float
    rel_tol: float
    ok: bool


@dataclass(frozen=True)
class ArtifactVerificationReport:
    """Self-describing verifier result for one exported paper benchmark."""

    point_id: str
    observed_point_ids: tuple[str, ...]
    wilson_bundle_id: str
    hadronic_bundle_id: str
    observable_bundle_id: str
    provenance_bundle_id: str
    coefficient_count: int
    observable_count: int
    wilson_operator_names: tuple[str, ...]
    hadronic_supported_operator_names: tuple[str, ...]
    hadronic_unsupported_operator_names: tuple[str, ...]
    observable_names: tuple[str, ...]
    contract: FrozenVerifierContract
    tolerance_policy: TolerancePolicy
    import_isolation: ImportIsolationReport
    recomputed_q1_matrix_element_GeV4: float
    recomputed_b_k_mu_had_from_hat_bk: float | None
    recomputed_m12_GeV: ComplexValue
    recomputed_delta_m_GeV: float
    numeric_checks: tuple[NumericCheckResult, ...]
    issues: tuple[VerificationIssue, ...]

    @property
    def ok(self) -> bool:
        return not self.issues

    @property
    def issue_codes(self) -> tuple[str, ...]:
        return tuple(issue.code for issue in self.issues)

    def require_ok(self) -> "ArtifactVerificationReport":
        if self.issues:
            joined = "; ".join(
                f"{issue.code}: {issue.message}" for issue in self.issues
            )
            raise ArtifactSchemaError(joined)
        return self


def _load_typed_bundle(path: str | Path, *, expected_type: type[object]) -> object:
    artifact = read_artifact(path)
    if not isinstance(artifact, expected_type):
        raise ArtifactSchemaError(
            f"{Path(path).name} is {type(artifact).__name__}, expected "
            f"{expected_type.__name__}"
        )
    return artifact


def read_verifier_inputs(
    wilson_path: str | Path,
    hadronic_path: str | Path,
    observable_path: str | Path,
    provenance_path: str | Path,
) -> VerifierInputSet:
    """Load the verifier inputs from exported paper artifacts only."""

    return VerifierInputSet(
        wilson_bundle=_load_typed_bundle(wilson_path, expected_type=WilsonArtifactBundleV1),
        hadronic_bundle=_load_typed_bundle(hadronic_path, expected_type=HadronicArtifactBundleV1),
        observable_bundle=_load_typed_bundle(
            observable_path, expected_type=ObservableArtifactBundleV1
        ),
        provenance_bundle=_load_typed_bundle(
            provenance_path, expected_type=ProvenanceBundleV1
        ),
    )


def _scale_map(bundle_scales) -> dict[str, float]:
    return {item.name: float(item.value_gev) for item in bundle_scales}


def _hadronic_scale_map(bundle: HadronicArtifactBundleV1) -> dict[str, float]:
    bundle_scales = getattr(bundle, "scales", None)
    if bundle_scales is not None:
        return _scale_map(bundle_scales)
    scale_map: dict[str, float] = {}
    mu_had_value = getattr(bundle, "mu_had_GeV", None)
    if mu_had_value is not None:
        scale_map["mu_had"] = float(mu_had_value)
    return scale_map


def _coefficient_map(bundle: WilsonArtifactBundleV1) -> dict[str, complex]:
    return {record.operator: record.value.to_complex() for record in bundle.coefficients}


def _observable_record_map(bundle: ObservableArtifactBundleV1) -> dict[str, ObservableRecord]:
    return {record.name: record for record in bundle.observables}


def _hadronic_source_wilson_bundle_id(bundle: HadronicArtifactBundleV1) -> str | None:
    value = getattr(bundle, "source_wilson_bundle_id", None)
    if value is None:
        return None
    return str(value)


def _hadronic_system(bundle: HadronicArtifactBundleV1) -> str:
    try:
        return str(bundle.system)
    except AttributeError:
        return str(bundle.system_id)


def _relative_diff(expected: float, actual: float) -> float:
    scale = max(abs(expected), abs(actual))
    if scale == 0.0:
        return 0.0
    return abs(expected - actual) / scale


def _beta_0(n_f: int) -> float:
    return 11.0 - (2.0 / 3.0) * float(n_f)


def _infer_n_f_from_scale(mu_had_gev: float) -> int:
    n_f = 3
    for threshold_mass_gev, n_f_above in DEFAULT_THRESHOLD_MASSES_GEV:
        if mu_had_gev >= threshold_mass_gev:
            n_f = n_f_above
        else:
            break
    return n_f


def _derive_implied_alpha_s_mu_had(
    *,
    hat_b_k_rgi: float,
    b_k_mu_had: float,
    mu_had_gev: float,
) -> tuple[float, float, int]:
    n_f = _infer_n_f_from_scale(mu_had_gev)
    exponent = 4.0 / (2.0 * _beta_0(n_f))
    if exponent <= 0.0:
        raise ArtifactSchemaError("frozen bag-parameter exponent must be positive")
    implied_alpha_s = (b_k_mu_had / hat_b_k_rgi) ** (1.0 / exponent)
    return implied_alpha_s, exponent, n_f


def _numeric_check(
    *,
    name: str,
    issue_code: str,
    expected: float,
    record: ObservableRecord | None,
    tolerance_policy: TolerancePolicy,
) -> NumericCheckResult:
    if record is None:
        return NumericCheckResult(
            name=name,
            issue_code=issue_code,
            expected=expected,
            actual=None,
            present=False,
            expected_units=EXPECTED_OBSERVABLE_UNITS,
            actual_units=None,
            abs_diff=None,
            rel_diff=None,
            abs_tol=tolerance_policy.observable_abs_tol,
            rel_tol=tolerance_policy.observable_rel_tol,
            ok=False,
        )
    actual = float(record.value)
    abs_diff = abs(expected - actual)
    rel_diff = _relative_diff(expected, actual)
    return NumericCheckResult(
        name=name,
        issue_code=issue_code,
        expected=expected,
        actual=actual,
        present=True,
        expected_units=EXPECTED_OBSERVABLE_UNITS,
        actual_units=record.units,
        abs_diff=abs_diff,
        rel_diff=rel_diff,
        abs_tol=tolerance_policy.observable_abs_tol,
        rel_tol=tolerance_policy.observable_rel_tol,
        ok=math.isclose(
            expected,
            actual,
            rel_tol=tolerance_policy.observable_rel_tol,
            abs_tol=tolerance_policy.observable_abs_tol,
        ),
    )


def _append_issue(
    issues: list[VerificationIssue],
    *,
    code: str,
    message: str,
    subject: str | None = None,
) -> None:
    issues.append(VerificationIssue(code=code, message=message, subject=subject))


def _duplicates(values: tuple[str, ...]) -> tuple[str, ...]:
    seen: set[str] = set()
    duplicates: list[str] = []
    for value in values:
        if value in seen and value not in duplicates:
            duplicates.append(value)
        seen.add(value)
    return tuple(duplicates)


def _static_import_violations() -> tuple[str, ...]:
    source_path = Path(__file__)
    source = source_path.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(source_path))
    violations: list[str] = []

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                module_name = alias.name
                if module_name == "quarkConstraints" or module_name.startswith(
                    "quarkConstraints."
                ):
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
            elif module_name == "quarkConstraints" or module_name.startswith(
                "quarkConstraints."
            ):
                if module_name not in ALLOWED_ABSOLUTE_LOCAL_IMPORTS:
                    violations.append(f"from {module_name} import ...")

    return tuple(sorted(set(violations)))


def _runtime_import_violations() -> tuple[str, ...]:
    allowed_modules = {
        __name__,
        __package__ or "",
        "quarkConstraints",
        *ALLOWED_ABSOLUTE_LOCAL_IMPORTS,
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
        if module_name.startswith(f"{PACKAGE_ROOT}."):
            if module_name in allowed_modules:
                continue
            if any(
                module_name == allowed or module_name.startswith(f"{allowed}.")
                for allowed in ALLOWED_ABSOLUTE_LOCAL_IMPORTS
            ):
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


def verify_inputs(inputs: VerifierInputSet) -> ArtifactVerificationReport:
    """Validate that the exported bundles are self-consistent and numerically sound."""

    wilson = inputs.wilson_bundle
    hadronic = inputs.hadronic_bundle
    observable = inputs.observable_bundle
    provenance = inputs.provenance_bundle
    contract = FrozenVerifierContract()
    tolerance_policy = TolerancePolicy()
    import_isolation = _build_import_isolation_report()
    issues: list[VerificationIssue] = []

    if not import_isolation.static_ok:
        _append_issue(
            issues,
            code="import_isolation_static_violation",
            message=(
                "verifier source imports modules outside the allowed artifacts/conventions "
                f"boundary: {', '.join(import_isolation.unexpected_import_targets)}"
            ),
        )
    if not import_isolation.runtime_ok:
        _append_issue(
            issues,
            code="import_isolation_runtime_violation",
            message=(
                "runtime imported forbidden modules while loading or executing the "
                f"standalone verifier: {', '.join(import_isolation.forbidden_loaded_modules)}"
            ),
        )
    if not import_isolation.ok:
        _append_issue(
            issues,
            code="import_isolation_failed",
            message="standalone verifier import isolation failed",
        )

    observed_point_ids = (
        wilson.metadata.point_id,
        hadronic.metadata.point_id,
        observable.metadata.point_id,
        provenance.metadata.point_id,
    )
    if len(set(observed_point_ids)) != 1:
        _append_issue(
            issues,
            code="point_id_mismatch",
            message="all artifact bundles must carry the same point_id",
        )

    if observable.source_wilson_bundle_id != wilson.metadata.bundle_id:
        _append_issue(
            issues,
            code="source_wilson_bundle_id_mismatch",
            message="observable bundle does not reference the exported Wilson bundle",
            subject="observable.source_wilson_bundle_id",
        )

    hadronic_source_wilson_bundle_id = _hadronic_source_wilson_bundle_id(hadronic)
    if (
        hadronic_source_wilson_bundle_id is not None
        and hadronic_source_wilson_bundle_id != wilson.metadata.bundle_id
    ):
        _append_issue(
            issues,
            code="hadronic_source_wilson_bundle_id_mismatch",
            message="hadronic bundle does not reference the exported Wilson bundle",
            subject="hadronic.source_wilson_bundle_id",
        )

    if observable.source_hadronic_bundle_id != hadronic.metadata.bundle_id:
        _append_issue(
            issues,
            code="source_hadronic_bundle_id_mismatch",
            message="observable bundle does not reference the exported hadronic bundle",
            subject="observable.source_hadronic_bundle_id",
        )

    provenance_ids = {record.record_id for record in provenance.records}
    referenced_ids = (
        set(wilson.provenance_ids)
        | set(hadronic.provenance_ids)
        | set(observable.provenance_ids)
    )
    for record_id in sorted(referenced_ids - provenance_ids):
        _append_issue(
            issues,
            code="missing_provenance_record",
            message=f"referenced provenance record {record_id!r} is missing",
            subject=record_id,
        )

    if hadronic.operator_basis != wilson.operator_basis:
        _append_issue(
            issues,
            code="hadronic_operator_basis_mismatch",
            message="hadronic and Wilson bundles must agree on operator_basis",
        )
    if observable.operator_basis != wilson.operator_basis:
        _append_issue(
            issues,
            code="operator_basis_mismatch",
            message="observable and Wilson bundles must agree on operator_basis",
        )
    if hadronic.renormalization_scheme != wilson.renormalization_scheme:
        _append_issue(
            issues,
            code="hadronic_renormalization_scheme_mismatch",
            message="hadronic and Wilson bundles must agree on renormalization_scheme",
        )
    if observable.renormalization_scheme != wilson.renormalization_scheme:
        _append_issue(
            issues,
            code="renormalization_scheme_mismatch",
            message="observable and Wilson bundles must agree on renormalization_scheme",
        )

    if wilson.operator_basis != contract.operator_basis:
        _append_issue(
            issues,
            code="convention_operator_basis_unexpected",
            message=(
                "wilson bundle operator_basis does not match the frozen verifier contract "
                f"{contract.operator_basis!r}"
            ),
            subject="wilson.operator_basis",
        )
    if hadronic.operator_basis != contract.operator_basis:
        _append_issue(
            issues,
            code="convention_operator_basis_unexpected",
            message=(
                "hadronic bundle operator_basis does not match the frozen verifier contract "
                f"{contract.operator_basis!r}"
            ),
            subject="hadronic.operator_basis",
        )
    if observable.operator_basis != contract.operator_basis:
        _append_issue(
            issues,
            code="convention_operator_basis_unexpected",
            message=(
                "observable bundle operator_basis does not match the frozen verifier "
                f"contract {contract.operator_basis!r}"
            ),
            subject="observable.operator_basis",
        )

    if wilson.renormalization_scheme != contract.renormalization_scheme:
        _append_issue(
            issues,
            code="convention_renormalization_scheme_unexpected",
            message=(
                "wilson bundle renormalization_scheme does not match the frozen verifier "
                f"contract {contract.renormalization_scheme!r}"
            ),
            subject="wilson.renormalization_scheme",
        )
    if hadronic.renormalization_scheme != contract.renormalization_scheme:
        _append_issue(
            issues,
            code="convention_renormalization_scheme_unexpected",
            message=(
                "hadronic bundle renormalization_scheme does not match the frozen "
                f"verifier contract {contract.renormalization_scheme!r}"
            ),
            subject="hadronic.renormalization_scheme",
        )
    if observable.renormalization_scheme != contract.renormalization_scheme:
        _append_issue(
            issues,
            code="convention_renormalization_scheme_unexpected",
            message=(
                "observable bundle renormalization_scheme does not match the frozen "
                f"verifier contract {contract.renormalization_scheme!r}"
            ),
            subject="observable.renormalization_scheme",
        )

    if observable.interpretation != contract.interpretation:
        _append_issue(
            issues,
            code="convention_interpretation_unexpected",
            message=(
                "observable bundle interpretation does not match the frozen verifier "
                f"contract {contract.interpretation!r}"
            ),
            subject="observable.interpretation",
        )

    hadronic_system = _hadronic_system(hadronic)
    if hadronic_system != contract.system:
        _append_issue(
            issues,
            code="convention_system_unexpected",
            message=(
                "hadronic bundle system does not match the frozen verifier contract "
                f"{contract.system!r}"
            ),
            subject="hadronic.system",
        )

    coefficient_names = tuple(record.operator for record in wilson.coefficients)
    observable_names = tuple(record.name for record in observable.observables)
    for duplicate_name in _duplicates(coefficient_names):
        _append_issue(
            issues,
            code="duplicate_wilson_operator",
            message=f"Wilson bundle repeats operator {duplicate_name!r}",
            subject=duplicate_name,
        )
    for duplicate_name in _duplicates(observable_names):
        _append_issue(
            issues,
            code="duplicate_observable_row",
            message=f"observable bundle repeats row {duplicate_name!r}",
            subject=duplicate_name,
        )

    coefficient_sectors = {record.sector for record in wilson.coefficients}
    coefficient_systems = {record.system for record in wilson.coefficients}
    for sector_name in sorted(coefficient_sectors):
        if sector_name != contract.sector:
            _append_issue(
                issues,
                code="convention_wilson_sector_unexpected",
                message=(
                    "wilson bundle contains sector values outside the frozen verifier "
                    f"contract {contract.sector!r}"
                ),
                subject=sector_name,
            )
    for system_name in sorted(coefficient_systems):
        if system_name != contract.system:
            _append_issue(
                issues,
                code="convention_wilson_system_unexpected",
                message=(
                    "wilson bundle contains system values outside the frozen verifier "
                    f"contract {contract.system!r}"
                ),
                subject=system_name,
            )

    if hadronic.hamiltonian_convention_id != contract.hamiltonian_convention_id:
        _append_issue(
            issues,
            code="convention_hamiltonian_mismatch",
            message=(
                "hadronic bundle hamiltonian_convention_id does not match the frozen "
                f"contract {contract.hamiltonian_convention_id!r}"
            ),
            subject="hadronic.hamiltonian_convention_id",
        )
    if hadronic.matrix_element_formula_id != contract.matrix_element_formula_id:
        _append_issue(
            issues,
            code="convention_matrix_element_formula_mismatch",
            message=(
                "hadronic bundle matrix_element_formula_id does not match the frozen "
                f"contract {contract.matrix_element_formula_id!r}"
            ),
            subject="hadronic.matrix_element_formula_id",
        )
    if hadronic.parity_relation_id != contract.parity_relation_id:
        _append_issue(
            issues,
            code="convention_parity_relation_mismatch",
            message=(
                "hadronic bundle parity_relation_id does not match the frozen verifier "
                f"contract {contract.parity_relation_id!r}"
            ),
            subject="hadronic.parity_relation_id",
        )
    if (
        hadronic.bag_parameter_source_scheme_id is not None
        and hadronic.bag_parameter_source_scheme_id
        != contract.bag_parameter_source_scheme_id
    ):
        _append_issue(
            issues,
            code="bag_parameter_source_scheme_mismatch",
            message=(
                "hadronic bundle bag_parameter_source_scheme_id does not match the frozen "
                f"contract {contract.bag_parameter_source_scheme_id!r}"
            ),
            subject="hadronic.bag_parameter_source_scheme_id",
        )
    if (
        hadronic.bag_parameter_transformation_id is not None
        and hadronic.bag_parameter_transformation_id
        != contract.bag_parameter_transformation_id
    ):
        _append_issue(
            issues,
            code="bag_parameter_transformation_mismatch",
            message=(
                "hadronic bundle bag_parameter_transformation_id does not match the frozen "
                f"contract {contract.bag_parameter_transformation_id!r}"
            ),
            subject="hadronic.bag_parameter_transformation_id",
        )
    if (
        hadronic.alpha_s_policy_id is not None
        and hadronic.alpha_s_policy_id != contract.alpha_s_policy_id
    ):
        _append_issue(
            issues,
            code="bag_parameter_alpha_s_policy_mismatch",
            message=(
                "hadronic bundle alpha_s_policy_id does not match the frozen contract "
                f"{contract.alpha_s_policy_id!r}"
            ),
            subject="hadronic.alpha_s_policy_id",
        )
    if (
        hadronic.input_provenance_mode_id is not None
        and hadronic.input_provenance_mode_id != contract.provenance_mode_id
    ):
        _append_issue(
            issues,
            code="hadronic_provenance_mode_mismatch",
            message=(
                "hadronic bundle input_provenance_mode_id does not match the frozen "
                f"contract {contract.provenance_mode_id!r}"
            ),
            subject="hadronic.input_provenance_mode_id",
        )

    if tuple(hadronic.supported_operator_names) != contract.supported_operators:
        _append_issue(
            issues,
            code="supported_operator_subset_mismatch",
            message=(
                "hadronic bundle must advertise the supported Q1/LR subset "
                "in the frozen order"
            ),
            subject="hadronic.supported_operator_names",
        )
    if tuple(hadronic.unsupported_operator_names) != contract.unsupported_operators:
        _append_issue(
            issues,
            code="unsupported_operator_subset_mismatch",
            message=(
                "hadronic bundle must not advertise any guarded unsupported operators"
            ),
            subject="hadronic.unsupported_operator_names",
        )

    wilson_scales = _scale_map(wilson.scales)
    hadronic_scales = _hadronic_scale_map(hadronic)
    observable_scales = _scale_map(observable.scales)
    if not set(hadronic_scales).issubset(set(wilson_scales)):
        _append_issue(
            issues,
            code="hadronic_scale_name_set_mismatch",
            message="hadronic scales must be a named subset of the Wilson scales",
        )
    if set(wilson_scales) != set(observable_scales):
        _append_issue(
            issues,
            code="scale_name_set_mismatch",
            message="observable and Wilson bundles must carry the same named scales",
        )

    shared_scale_names = (
        set(wilson_scales) & set(hadronic_scales) & set(observable_scales)
    )
    for scale_name in sorted(shared_scale_names):
        if not math.isclose(
            wilson_scales[scale_name],
            hadronic_scales[scale_name],
            rel_tol=0.0,
            abs_tol=tolerance_policy.scale_abs_tol,
        ):
            _append_issue(
                issues,
                code="hadronic_scale_value_mismatch",
                message=(
                    f"shared scale {scale_name!r} differs across Wilson and hadronic bundles"
                ),
                subject=scale_name,
            )
        if not math.isclose(
            wilson_scales[scale_name],
            observable_scales[scale_name],
            rel_tol=0.0,
            abs_tol=tolerance_policy.scale_abs_tol,
        ):
            _append_issue(
                issues,
                code="scale_value_mismatch",
                message=f"shared scale {scale_name!r} differs across bundles",
                subject=scale_name,
            )

    if wilson.coefficient_scale_name not in wilson_scales:
        _append_issue(
            issues,
            code="wilson_scale_reference_mismatch",
            message=(
                "wilson coefficient_scale_name does not resolve to a named entry in "
                "wilson.scales"
            ),
            subject="wilson.coefficient_scale_name",
        )
    if wilson.matching_scale_name not in wilson_scales:
        _append_issue(
            issues,
            code="wilson_scale_reference_mismatch",
            message=(
                "wilson matching_scale_name does not resolve to a named entry in "
                "wilson.scales"
            ),
            subject="wilson.matching_scale_name",
        )

    coefficient_map = _coefficient_map(wilson)
    expected_operator_names = set(contract.supported_operators) | set(
        contract.unsupported_operators
    )
    for operator_name in sorted(set(coefficient_map) - expected_operator_names):
        _append_issue(
            issues,
            code="unexpected_operator_name",
            message=(
                "standalone verifier only supports the frozen Q1/Q4/Q5 export boundary; "
                f"encountered unexpected Wilson operator {operator_name!r}"
            ),
            subject=operator_name,
        )

    q1_vll = complex(coefficient_map.get("Q1_VLL", 0.0 + 0.0j))
    q1_vrr = complex(coefficient_map.get("Q1_VRR", 0.0 + 0.0j))
    q4_lr = complex(coefficient_map.get("Q4_LR", 0.0 + 0.0j))
    q5_lr = complex(coefficient_map.get("Q5_LR", 0.0 + 0.0j))
    # RESIDUAL(C-2): default RH-down alignment model choice pending paper 0710.1869.

    recomputed_q1_matrix_element = (
        (2.0 / 3.0)
        * (hadronic.f_K_GeV**2)
        * (hadronic.m_K0_GeV**2)
        * hadronic.B_K_mu_had
    )
    recomputed_b_k_mu_had_from_hat_bk: float | None = None
    if hadronic.hat_B_K_rgi_source_value is not None:
        if hadronic.bag_parameter_source_scheme_id is None:
            _append_issue(
                issues,
                code="bag_parameter_source_scheme_missing",
                message=(
                    "hadronic bundle exports hat_B_K_rgi_source_value but omits "
                    "bag_parameter_source_scheme_id"
                ),
                subject="hadronic.bag_parameter_source_scheme_id",
            )
        if hadronic.bag_parameter_transformation_id is None:
            _append_issue(
                issues,
                code="bag_parameter_transformation_missing",
                message=(
                    "hadronic bundle exports hat_B_K_rgi_source_value but omits "
                    "bag_parameter_transformation_id"
                ),
                subject="hadronic.bag_parameter_transformation_id",
            )
        if hadronic.alpha_s_policy_id is None:
            _append_issue(
                issues,
                code="bag_parameter_alpha_s_policy_missing",
                message=(
                    "hadronic bundle exports hat_B_K_rgi_source_value but omits "
                    "alpha_s_policy_id"
                ),
                subject="hadronic.alpha_s_policy_id",
            )
        implied_alpha_s_mu_had, exponent, n_f = _derive_implied_alpha_s_mu_had(
            hat_b_k_rgi=hadronic.hat_B_K_rgi_source_value,
            b_k_mu_had=hadronic.B_K_mu_had,
            mu_had_gev=hadronic.mu_had_GeV,
        )
        recomputed_b_k_mu_had_from_hat_bk = (
            hadronic.hat_B_K_rgi_source_value * (implied_alpha_s_mu_had**exponent)
        )
        if not (0.0 < implied_alpha_s_mu_had < 1.0):
            _append_issue(
                issues,
                code="bag_parameter_alpha_s_implied_out_of_range",
                message=(
                    "implied alpha_s(mu_had) from the exported hat_B_K -> B_K(mu_had) "
                    f"conversion is outside (0, 1): {implied_alpha_s_mu_had:.16e} "
                    f"for n_f={n_f}"
                ),
                subject="hadronic.B_K_mu_had",
            )
        if hadronic.B_K_mu_had >= hadronic.hat_B_K_rgi_source_value:
            _append_issue(
                issues,
                code="bag_parameter_conversion_semantics_mismatch",
                message=(
                    "inverse Q1 running from hat_B_K to B_K(mu_had) must reduce the "
                    "RGI bag parameter at the exported hadronic scale"
                ),
                subject="hadronic.B_K_mu_had",
            )
        if not math.isclose(
            recomputed_b_k_mu_had_from_hat_bk,
            hadronic.B_K_mu_had,
            rel_tol=tolerance_policy.observable_rel_tol,
            abs_tol=tolerance_policy.matrix_element_abs_tol,
        ):
            _append_issue(
                issues,
                code="bag_parameter_conversion_reconstruction_mismatch",
                message=(
                    "reconstructed B_K(mu_had) from exported hat_B_K, mu_had, and the "
                    "frozen inverse-Q1-running semantics does not match the exported "
                    "B_K_mu_had"
                ),
                subject="hadronic.B_K_mu_had",
            )
        if (
            hadronic.bag_parameter_source is not None
            and hadronic.bag_parameter_source.transformation_id
            != contract.bag_parameter_transformation_id
        ):
            _append_issue(
                issues,
                code="bag_parameter_source_transformation_mismatch",
                message=(
                    "bag_parameter_source.transformation_id does not match the frozen "
                    "inverse-Q1-running conversion id"
                ),
                subject="hadronic.bag_parameter_source.transformation_id",
            )
        if (
            hadronic.bag_parameter_source is not None
            and hadronic.bag_parameter_source.renormalization_scheme_id
            != contract.bag_parameter_source_scheme_id
        ):
            _append_issue(
                issues,
                code="bag_parameter_source_scheme_semantics_mismatch",
                message=(
                    "bag_parameter_source.renormalization_scheme_id does not match "
                    "the frozen RGI source scheme id"
                ),
                subject="hadronic.bag_parameter_source.renormalization_scheme_id",
            )
    if not math.isclose(
        recomputed_q1_matrix_element,
        hadronic.q1_matrix_element_GeV4,
        rel_tol=tolerance_policy.observable_rel_tol,
        abs_tol=tolerance_policy.matrix_element_abs_tol,
    ):
        _append_issue(
            issues,
            code="q1_matrix_element_mismatch",
            message=(
                "hadronic q1_matrix_element_GeV4 does not match the frozen "
                "8/3 f_K^2 m_K^2 B_K(mu_had) identity"
            ),
            subject="hadronic.q1_matrix_element_GeV4",
        )

    observable_records = _observable_record_map(observable)
    missing_observables = tuple(
        name for name in contract.required_observables if name not in observable_records
    )
    for name in missing_observables:
        _append_issue(
            issues,
            code="observable_row_missing",
            message=f"observable bundle is missing required row {name!r}",
            subject=name,
        )
    for name in sorted(set(observable_records) - set(contract.required_observables)):
        _append_issue(
            issues,
            code="unexpected_observable_row",
            message=(
                "standalone verifier only supports the frozen kaon NP observable rows; "
                f"encountered unexpected row {name!r}"
            ),
            subject=name,
        )

    for row_name in contract.required_observables:
        record = observable_records.get(row_name)
        if record is None:
            continue
        if record.system != contract.system:
            _append_issue(
                issues,
                code="observable_system_mismatch",
                message=(
                    f"observable row {row_name!r} must use system {contract.system!r}"
                ),
                subject=row_name,
            )
        if record.units != contract.observable_units:
            _append_issue(
                issues,
                code="observable_units_mismatch",
                message=(
                    f"observable row {row_name!r} must use units "
                    f"{contract.observable_units!r}"
                ),
                subject=row_name,
            )

    recomputed_m12 = (
        (q1_vll + q1_vrr) * recomputed_q1_matrix_element / (2.0 * hadronic.m_K0_GeV)
    )
    recomputed_delta_m = 2.0 * float(recomputed_m12.real)

    numeric_checks = (
        _numeric_check(
            name="M12_K_NP.re",
            issue_code="m12_reconstruction_mismatch",
            expected=float(recomputed_m12.real),
            record=observable_records.get("M12_K_NP.re"),
            tolerance_policy=tolerance_policy,
        ),
        _numeric_check(
            name="M12_K_NP.im",
            issue_code="m12_imag_reconstruction_mismatch",
            expected=float(recomputed_m12.imag),
            record=observable_records.get("M12_K_NP.im"),
            tolerance_policy=tolerance_policy,
        ),
        _numeric_check(
            name="Delta_m_K_NP",
            issue_code="delta_m_reconstruction_mismatch",
            expected=recomputed_delta_m,
            record=observable_records.get("Delta_m_K_NP"),
            tolerance_policy=tolerance_policy,
        ),
    )
    for check in numeric_checks:
        if check.present and not check.ok:
            _append_issue(
                issues,
                code=check.issue_code,
                message=(
                    f"recomputed {check.name}={check.expected:.16e} does not match the "
                    f"observable bundle value {check.actual:.16e} within "
                    f"abs_tol={check.abs_tol:.1e}, rel_tol={check.rel_tol:.1e}"
                ),
                subject=check.name,
            )

    return ArtifactVerificationReport(
        point_id=wilson.metadata.point_id,
        observed_point_ids=observed_point_ids,
        wilson_bundle_id=wilson.metadata.bundle_id,
        hadronic_bundle_id=hadronic.metadata.bundle_id,
        observable_bundle_id=observable.metadata.bundle_id,
        provenance_bundle_id=provenance.metadata.bundle_id,
        coefficient_count=len(wilson.coefficients),
        observable_count=len(observable.observables),
        wilson_operator_names=tuple(sorted(coefficient_names)),
        hadronic_supported_operator_names=tuple(hadronic.supported_operator_names),
        hadronic_unsupported_operator_names=tuple(hadronic.unsupported_operator_names),
        observable_names=tuple(sorted(observable_names)),
        contract=contract,
        tolerance_policy=tolerance_policy,
        import_isolation=import_isolation,
        recomputed_q1_matrix_element_GeV4=recomputed_q1_matrix_element,
        recomputed_b_k_mu_had_from_hat_bk=recomputed_b_k_mu_had_from_hat_bk,
        recomputed_m12_GeV=ComplexValue.from_complex(recomputed_m12),
        recomputed_delta_m_GeV=recomputed_delta_m,
        numeric_checks=numeric_checks,
        issues=tuple(issues),
    )


def verify_artifact_paths(
    wilson_path: str | Path,
    hadronic_path: str | Path,
    observable_path: str | Path,
    provenance_path: str | Path,
) -> ArtifactVerificationReport:
    """Load and verify one paper artifact quartet from disk."""

    return verify_inputs(
        read_verifier_inputs(wilson_path, hadronic_path, observable_path, provenance_path)
    )


def verify_artifact_paths_public(
    wilson_path: str | Path,
    hadronic_path: str | Path,
    observable_path: str | Path,
    provenance_path: str | Path,
) -> dict[str, object]:
    """Return a JSON-safe benchmark-facing verifier payload."""

    report = verify_artifact_paths(
        wilson_path,
        hadronic_path,
        observable_path,
        provenance_path,
    )
    return {
        "ok": report.ok,
        "issue_codes": list(report.issue_codes),
        "tolerances": frozen_public_verifier_tolerances(),
        "import_isolation_ok": report.import_isolation.ok,
        "loaded_forbidden_modules": list(report.import_isolation.forbidden_loaded_modules),
        "point_id": report.point_id,
        "bundle_ids": {
            "wilsons": report.wilson_bundle_id,
            "hadronic": report.hadronic_bundle_id,
            "observables": report.observable_bundle_id,
            "provenance": report.provenance_bundle_id,
        },
        "numeric_checks": {
            check.name: {
                "ok": check.ok,
                "expected": check.expected,
                "actual": check.actual,
                "abs_diff": check.abs_diff,
                "rel_diff": check.rel_diff,
                "abs_tol": check.abs_tol,
                "rel_tol": check.rel_tol,
            }
            for check in report.numeric_checks
        },
    }


__all__ = [
    "ArtifactVerificationReport",
    "FrozenVerifierContract",
    "ImportIsolationReport",
    "NumericCheckResult",
    "TolerancePolicy",
    "VerifierInputSet",
    "VerificationIssue",
    "frozen_public_verifier_tolerances",
    "read_verifier_inputs",
    "verify_artifact_paths",
    "verify_artifact_paths_public",
    "verify_inputs",
]
