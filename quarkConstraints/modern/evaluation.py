"""Deterministic pointwise evaluation contract for the modern quark lane."""

from __future__ import annotations

from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Mapping

from ..deltaf2 import (
    DeltaF2ConstraintSummary,
    DeltaF2Input,
    DeltaF2ObservableSummary,
    DeltaF2WilsonCoefficients,
)
from ..fit import QuarkFitResult, QuarkFitSolution
from ..scales import DEFAULT_QUARK_XI_KK, default_quark_m_kk_from_lambda_ir
from .couplings import (
    MODERN_POINT_COUPLINGS_SCHEMA_ID,
    ModernPointCouplings,
    build_modern_point_couplings,
)
from .inputs import (
    MODERN_DEFAULT_INPUT_BUNDLE_ID,
    MODERN_DEFAULT_INPUT_PROVENANCE_ID,
    MODERN_DEFAULT_INPUTS_SCHEMA_ID,
    MODERN_DEFAULT_RESOLUTION_POLICY_ID,
    ModernDefaultInputs,
    default_modern_default_inputs,
)
from .matching import (
    MODERN_POINT_MATCHING_SCHEMA_ID,
    ModernPointMatching,
    build_modern_point_matching,
)
from .phenomenology import (
    MODERN_PHENOMENOLOGY_SYSTEM_IDS,
    ModernPhenomenologyPolicy,
    ModernPhenomenologySystemPolicy,
    default_modern_phenomenology_policy,
)

MODERN_POINT_EVALUATION_SCHEMA_ID = "quarkConstraints.modern.evaluation.v1"
MODERN_POINT_EVALUATION_POLICY_ID = "quarkConstraints.modern.evaluation.policy.v1"
MODERN_POINT_EVALUATION_VERDICT_SCHEMA_ID = "quarkConstraints.modern.evaluation.verdict.v1"
MODERN_POINT_EVALUATION_SYSTEM_IDS = ("K", "B_d", "B_s", "D0")
MODERN_POINT_EVALUATION_OBSERVABLE_IDS = (
    "epsilon_K",
    "B_d",
    "B_s",
    "D0",
)
MODERN_POINT_EVALUATION_BACKEND_SYSTEM_IDS = ("K", "B_d", "B_s", "D")
MODERN_POINT_EVALUATION_BACKEND_KEYS = ("epsilon_k", "b_d", "b_s", "d")
MODERN_POINT_EVALUATION_POLICY_SYSTEM_IDS = (
    "epsilon_K",
    "B_d",
    "B_s",
    "D0",
)
MODERN_POINT_EVALUATION_POLICY_ID = "quarkConstraints.modern.evaluation.policy.v1"


def _require_text(name: str, value: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty string")
    return value.strip()


def _require_exact(name: str, value: str, *, expected: str) -> str:
    if value != expected:
        raise ValueError(f"{name} must be exactly {expected!r}")
    return value


def _require_positive_float(name: str, value: float) -> float:
    numeric = float(value)
    if numeric <= 0.0:
        raise ValueError(f"{name} must be positive")
    return numeric


def _require_nonnegative_float(name: str, value: float) -> float:
    numeric = float(value)
    if numeric < 0.0:
        raise ValueError(f"{name} must be non-negative")
    return numeric


def _index_for_system(system_id: str) -> int:
    try:
        return MODERN_POINT_EVALUATION_SYSTEM_IDS.index(system_id)
    except ValueError as exc:
        raise ValueError(
            f"system_id must be one of {MODERN_POINT_EVALUATION_SYSTEM_IDS!r}"
        ) from exc


def _policy_system_id_for_system(system_id: str) -> str:
    if system_id == "K":
        return "epsilon_K"
    if system_id in {"B_d", "B_s", "D0"}:
        return system_id
    raise ValueError(f"system_id must be one of {MODERN_POINT_EVALUATION_SYSTEM_IDS!r}")


def _backend_system_id_for_system(system_id: str) -> str:
    return MODERN_POINT_EVALUATION_BACKEND_SYSTEM_IDS[_index_for_system(system_id)]


def _backend_key_for_system(system_id: str) -> str:
    return MODERN_POINT_EVALUATION_BACKEND_KEYS[_index_for_system(system_id)]


def _observable_id_for_system(system_id: str) -> str:
    return MODERN_POINT_EVALUATION_OBSERVABLE_IDS[_index_for_system(system_id)]


def _policy_system_for_system(
    policy: ModernPhenomenologyPolicy,
    system_id: str,
) -> ModernPhenomenologySystemPolicy:
    return policy.system_policy(_policy_system_id_for_system(system_id))


def _coerce_fit_result(
    source: QuarkFitResult | QuarkFitSolution,
) -> QuarkFitResult:
    if isinstance(source, QuarkFitResult):
        return source
    if isinstance(source, QuarkFitSolution):
        return source.result
    raise TypeError(
        "source must be a QuarkFitResult or QuarkFitSolution instance"
    )


def _deltaf2_inputs_from_modern_inputs(
    bundle: ModernDefaultInputs,
) -> tuple[DeltaF2Input, ...]:
    weights = bundle.operator_weight_policy
    return tuple(
        DeltaF2Input(
            key=system.backend_key,
            display_name=system.display_name,
            column_name=system.column_name,
            reject_reason=system.reject_reason,
            sector=system.sector_id,
            generations=system.generations,
            bound=system.bound,
            ll_weight=weights.ll_weight,
            rr_weight=weights.rr_weight,
            lr1_weight=weights.lr1_weight,
            lr2_weight=weights.lr2_weight,
            reference_scale=weights.reference_scale_GeV,
            note=(
                f"{system.note} Uses modern bundle {bundle.bundle_id} with "
                f"provenance {bundle.provenance_id}."
            ),
        )
        for system in bundle.neutral_meson_inputs
    )


def _deltaf2_summary_from_matching(
    matching: ModernPointMatching,
    *,
    bundle: ModernDefaultInputs,
) -> DeltaF2ConstraintSummary:
    backend_inputs = {
        item.key: item for item in _deltaf2_inputs_from_modern_inputs(bundle)
    }
    observables: list[DeltaF2ObservableSummary] = []
    for system in matching.systems:
        item = backend_inputs[system.backend_key]
        wilsons = DeltaF2WilsonCoefficients(
            input=item,
            M_KK=float(matching.M_KK),
            matching_scale=float(system.matching_scale_GeV),
            left_coupling=complex(system.left_coupling),
            right_coupling=complex(system.right_coupling),
            c1_vll=complex(system.c1_vll),
            c1_vrr=complex(system.c1_vrr),
            c4_lr=complex(system.c4_lr),
            c5_lr=complex(system.c5_lr),
        )
        weighted = {
            "C1_VLL": item.ll_weight * wilsons.c1_vll,
            "C1_VRR": item.rr_weight * wilsons.c1_vrr,
            "C4_LR": item.lr1_weight * wilsons.c4_lr,
            "C5_LR": item.lr2_weight * wilsons.c5_lr,
        }
        operator_sizes = {
            name: float(item.reference_scale**2 * abs(value))
            for name, value in weighted.items()
        }
        dominant_operator = max(operator_sizes, key=operator_sizes.get)
        coherent_amplitude = float(item.reference_scale**2 * abs(sum(weighted.values())))
        effective_amplitude = float(operator_sizes[dominant_operator])
        ratio_to_bound = float(effective_amplitude / item.bound)
        observables.append(
            DeltaF2ObservableSummary(
                input=item,
                wilsons=wilsons,
                effective_amplitude=effective_amplitude,
                coherent_amplitude=coherent_amplitude,
                ratio_to_bound=ratio_to_bound,
                passes=ratio_to_bound <= 1.0,
                dominant_operator=dominant_operator,
                dominant_operator_size=operator_sizes[dominant_operator],
                weighted_operator_sizes=operator_sizes,
            )
        )
    return DeltaF2ConstraintSummary(
        model_label=matching.operator_basis_id,
        input_bundle_label=bundle.bundle_id,
        M_KK=float(matching.M_KK),
        xi_KK=float(matching.xi_KK),
        observables=tuple(observables),
    )


@dataclass(frozen=True)
class ModernPointVerdict:
    """Frozen verdict for one explicit modern neutral-meson system."""

    schema_id: str = MODERN_POINT_EVALUATION_VERDICT_SCHEMA_ID
    system_id: str = ""
    observable_id: str = ""
    policy_system_id: str = ""
    backend_system_id: str = ""
    backend_key: str = ""
    policy_id: str = ""
    policy_display_name: str = ""
    policy_notes: str = ""
    passes: bool = False
    ratio_to_bound: float = 0.0
    effective_amplitude: float = 0.0
    coherent_amplitude: float = 0.0
    bound: float = 0.0
    dominant_operator: str = ""
    dominant_operator_size: float = 0.0
    weighted_operator_sizes: Mapping[str, float] = field(default_factory=dict)
    note: str = ""

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            _require_exact(
                "schema_id",
                self.schema_id,
                expected=MODERN_POINT_EVALUATION_VERDICT_SCHEMA_ID,
            ),
        )
        object.__setattr__(self, "system_id", _require_text("system_id", self.system_id))
        object.__setattr__(
            self,
            "observable_id",
            _require_text("observable_id", self.observable_id),
        )
        object.__setattr__(
            self,
            "policy_system_id",
            _require_text("policy_system_id", self.policy_system_id),
        )
        object.__setattr__(
            self,
            "backend_system_id",
            _require_text("backend_system_id", self.backend_system_id),
        )
        object.__setattr__(self, "backend_key", _require_text("backend_key", self.backend_key))
        object.__setattr__(self, "policy_id", _require_text("policy_id", self.policy_id))
        object.__setattr__(
            self,
            "policy_display_name",
            _require_text("policy_display_name", self.policy_display_name),
        )
        object.__setattr__(self, "policy_notes", _require_text("policy_notes", self.policy_notes))
        object.__setattr__(self, "dominant_operator", _require_text("dominant_operator", self.dominant_operator))
        object.__setattr__(self, "note", _require_text("note", self.note))
        object.__setattr__(self, "ratio_to_bound", float(self.ratio_to_bound))
        object.__setattr__(self, "effective_amplitude", float(self.effective_amplitude))
        object.__setattr__(self, "coherent_amplitude", float(self.coherent_amplitude))
        object.__setattr__(self, "bound", _require_positive_float("bound", self.bound))
        object.__setattr__(
            self,
            "dominant_operator_size",
            _require_nonnegative_float("dominant_operator_size", self.dominant_operator_size),
        )
        if self.system_id not in MODERN_POINT_EVALUATION_SYSTEM_IDS:
            raise ValueError(
                f"system_id must be one of {MODERN_POINT_EVALUATION_SYSTEM_IDS!r}"
            )
        expected_observable = _observable_id_for_system(self.system_id)
        expected_policy_system = _policy_system_id_for_system(self.system_id)
        expected_backend_system = _backend_system_id_for_system(self.system_id)
        expected_backend_key = _backend_key_for_system(self.system_id)
        if self.observable_id != expected_observable:
            raise ValueError(
                f"observable_id for {self.system_id!r} must be exactly {expected_observable!r}"
            )
        if self.policy_system_id != expected_policy_system:
            raise ValueError(
                f"policy_system_id for {self.system_id!r} must be exactly {expected_policy_system!r}"
            )
        if self.backend_system_id != expected_backend_system:
            raise ValueError(
                f"backend_system_id for {self.system_id!r} must be exactly {expected_backend_system!r}"
            )
        if self.backend_key != expected_backend_key:
            raise ValueError(
                f"backend_key for {self.system_id!r} must be exactly {expected_backend_key!r}"
            )
        if not isinstance(self.weighted_operator_sizes, Mapping):
            raise ValueError("weighted_operator_sizes must be a mapping")
        canonical_operator_sizes = dict(
            sorted((str(key), float(value)) for key, value in self.weighted_operator_sizes.items())
        )
        object.__setattr__(
            self,
            "weighted_operator_sizes",
            MappingProxyType(canonical_operator_sizes),
        )

    @property
    def system(self) -> str:
        return self.system_id

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "system_id": self.system_id,
            "observable_id": self.observable_id,
            "policy_system_id": self.policy_system_id,
            "backend_system_id": self.backend_system_id,
            "backend_key": self.backend_key,
            "policy_id": self.policy_id,
            "policy_display_name": self.policy_display_name,
            "policy_notes": self.policy_notes,
            "passes": self.passes,
            "ratio_to_bound": self.ratio_to_bound,
            "effective_amplitude": self.effective_amplitude,
            "coherent_amplitude": self.coherent_amplitude,
            "bound": self.bound,
            "dominant_operator": self.dominant_operator,
            "dominant_operator_size": self.dominant_operator_size,
            "weighted_operator_sizes": dict(self.weighted_operator_sizes),
            "note": self.note,
        }


@dataclass(frozen=True)
class ModernPointEvaluation:
    """Frozen per-point modern evaluation record."""

    schema_id: str = MODERN_POINT_EVALUATION_SCHEMA_ID
    policy: ModernPhenomenologyPolicy = field(default_factory=default_modern_phenomenology_policy)
    input_bundle_schema_id: str = MODERN_DEFAULT_INPUTS_SCHEMA_ID
    input_bundle_id: str = MODERN_DEFAULT_INPUT_BUNDLE_ID
    input_provenance_id: str = MODERN_DEFAULT_INPUT_PROVENANCE_ID
    input_resolution_policy_id: str = MODERN_DEFAULT_RESOLUTION_POLICY_ID
    point_id: str = ""
    point_label: str = ""
    M_KK: float = 0.0
    xi_KK: float = 0.0
    couplings: ModernPointCouplings | None = None
    matching: ModernPointMatching | None = None
    deltaf2: DeltaF2ConstraintSummary | None = None
    verdicts: tuple[ModernPointVerdict, ...] = ()
    notes: str = (
        "Deterministic per-point modern evaluation. It carries the explicit "
        "modern QS2 bridge plus verdicts for K, B_d, B_s, and D0 only."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            _require_exact("schema_id", self.schema_id, expected=MODERN_POINT_EVALUATION_SCHEMA_ID),
        )
        if not isinstance(self.policy, ModernPhenomenologyPolicy):
            raise ValueError("policy must be a ModernPhenomenologyPolicy instance")
        self.policy.require_system_ids(MODERN_PHENOMENOLOGY_SYSTEM_IDS)
        object.__setattr__(
            self,
            "input_bundle_schema_id",
            _require_exact(
                "input_bundle_schema_id",
                self.input_bundle_schema_id,
                expected=MODERN_DEFAULT_INPUTS_SCHEMA_ID,
            ),
        )
        object.__setattr__(
            self,
            "input_bundle_id",
            _require_exact(
                "input_bundle_id",
                self.input_bundle_id,
                expected=MODERN_DEFAULT_INPUT_BUNDLE_ID,
            ),
        )
        object.__setattr__(
            self,
            "input_provenance_id",
            _require_exact(
                "input_provenance_id",
                self.input_provenance_id,
                expected=MODERN_DEFAULT_INPUT_PROVENANCE_ID,
            ),
        )
        object.__setattr__(
            self,
            "input_resolution_policy_id",
            _require_exact(
                "input_resolution_policy_id",
                self.input_resolution_policy_id,
                expected=MODERN_DEFAULT_RESOLUTION_POLICY_ID,
            ),
        )
        object.__setattr__(self, "point_id", _require_text("point_id", self.point_id))
        object.__setattr__(self, "point_label", _require_text("point_label", self.point_label))
        object.__setattr__(self, "M_KK", _require_positive_float("M_KK", self.M_KK))
        object.__setattr__(self, "xi_KK", _require_positive_float("xi_KK", self.xi_KK))
        if not isinstance(self.couplings, ModernPointCouplings):
            raise ValueError("couplings must be a ModernPointCouplings instance")
        if not isinstance(self.matching, ModernPointMatching):
            raise ValueError("matching must be a ModernPointMatching instance")
        if self.couplings.schema_id != MODERN_POINT_COUPLINGS_SCHEMA_ID:
            raise ValueError("couplings must use the frozen modern couplings schema")
        if self.matching.schema_id != MODERN_POINT_MATCHING_SCHEMA_ID:
            raise ValueError("matching must use the frozen modern matching schema")
        if self.couplings.point_id != self.point_id:
            raise ValueError("couplings.point_id must match point_id")
        if self.matching.point_id != self.point_id:
            raise ValueError("matching.point_id must match point_id")
        if self.couplings.point_label != self.point_label:
            raise ValueError("couplings.point_label must match point_label")
        if self.matching.point_label != self.point_label:
            raise ValueError("matching.point_label must match point_label")
        if self.couplings.input_bundle_id != self.input_bundle_id:
            raise ValueError("couplings.input_bundle_id must match input_bundle_id")
        if self.matching.input_bundle_id != self.input_bundle_id:
            raise ValueError("matching.input_bundle_id must match input_bundle_id")
        if self.couplings.input_provenance_id != self.input_provenance_id:
            raise ValueError("couplings.input_provenance_id must match input_provenance_id")
        if self.matching.input_provenance_id != self.input_provenance_id:
            raise ValueError("matching.input_provenance_id must match input_provenance_id")
        if self.couplings.input_resolution_policy_id != self.input_resolution_policy_id:
            raise ValueError("couplings.input_resolution_policy_id must match input_resolution_policy_id")
        if self.matching.input_resolution_policy_id != self.input_resolution_policy_id:
            raise ValueError("matching.input_resolution_policy_id must match input_resolution_policy_id")
        if self.couplings.M_KK != self.M_KK or self.matching.M_KK != self.M_KK:
            raise ValueError("bridge M_KK must match evaluation M_KK")
        if self.couplings.xi_KK != self.xi_KK or self.matching.xi_KK != self.xi_KK:
            raise ValueError("bridge xi_KK must match evaluation xi_KK")
        if self.deltaf2 is None:
            raise ValueError("deltaf2 must be provided")
        if not isinstance(self.deltaf2, DeltaF2ConstraintSummary):
            raise ValueError("deltaf2 must be a DeltaF2ConstraintSummary instance")
        if not isinstance(self.verdicts, tuple):
            raise ValueError("verdicts must be a tuple of ModernPointVerdict instances")
        normalized_verdicts = tuple(self.verdicts)
        if len(normalized_verdicts) != len(MODERN_POINT_EVALUATION_SYSTEM_IDS):
            raise ValueError(
                "verdicts must contain exactly K, B_d, B_s, and D0"
            )
        for verdict in normalized_verdicts:
            if not isinstance(verdict, ModernPointVerdict):
                raise ValueError("verdicts must contain only ModernPointVerdict instances")
        normalized_system_ids = tuple(verdict.system_id for verdict in normalized_verdicts)
        if normalized_system_ids != MODERN_POINT_EVALUATION_SYSTEM_IDS:
            unknown = tuple(
                system_id
                for system_id in normalized_system_ids
                if system_id not in MODERN_POINT_EVALUATION_SYSTEM_IDS
            )
            missing = tuple(
                system_id
                for system_id in MODERN_POINT_EVALUATION_SYSTEM_IDS
                if system_id not in normalized_system_ids
            )
            if unknown:
                raise ValueError(f"verdicts contains unknown systems: {unknown!r}")
            if missing:
                raise ValueError(f"verdicts is missing systems: {missing!r}")
            raise ValueError("verdicts must preserve the frozen modern evaluation order")
        if tuple(verdict.observable_id for verdict in normalized_verdicts) != MODERN_POINT_EVALUATION_OBSERVABLE_IDS:
            raise ValueError("verdicts must preserve the frozen modern observable order")
        for verdict in normalized_verdicts:
            policy_system = self.policy.system_policy(verdict.policy_system_id)
            if verdict.policy_id != policy_system.policy_id:
                raise ValueError(
                    f"verdict policy_id for {verdict.system_id!r} must match the supplied policy entry"
                )
            if verdict.policy_display_name != policy_system.display_name:
                raise ValueError(
                    f"verdict policy_display_name for {verdict.system_id!r} must match the supplied policy entry"
                )
            if verdict.policy_notes != policy_system.notes:
                raise ValueError(
                    f"verdict policy_notes for {verdict.system_id!r} must match the supplied policy entry"
                )
        object.__setattr__(self, "verdicts", normalized_verdicts)
        verdict_passes = all(verdict.passes for verdict in self.verdicts)
        if self.deltaf2.passes_all != verdict_passes:
            raise ValueError("deltaf2 and verdicts must agree on the pass/fail outcome")

    @property
    def system_ids(self) -> tuple[str, ...]:
        return tuple(verdict.system_id for verdict in self.verdicts)

    @property
    def verdict_by_system(self) -> dict[str, ModernPointVerdict]:
        return {verdict.system_id: verdict for verdict in self.verdicts}

    def verdict_for(self, system_id: str) -> ModernPointVerdict:
        requested = _require_text("system_id", system_id)
        for verdict in self.verdicts:
            if verdict.system_id == requested:
                return verdict
        raise ValueError(
            f"system_id must be one of {MODERN_POINT_EVALUATION_SYSTEM_IDS!r}"
        )

    def require_system_ids(self, system_ids: tuple[str, ...]) -> tuple[str, ...]:
        provided = tuple(_require_text("system_id", system_id) for system_id in system_ids)
        unknown = tuple(
            system_id for system_id in provided if system_id not in MODERN_POINT_EVALUATION_SYSTEM_IDS
        )
        missing = tuple(
            system_id for system_id in MODERN_POINT_EVALUATION_SYSTEM_IDS if system_id not in provided
        )
        if unknown:
            raise ValueError(f"system_ids contains unknown systems: {unknown!r}")
        if missing:
            raise ValueError(f"system_ids is missing systems: {missing!r}")
        if provided != MODERN_POINT_EVALUATION_SYSTEM_IDS:
            raise ValueError("system_ids must preserve the frozen modern evaluation order")
        return provided

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "policy": self.policy.as_dict(),
            "input_bundle_schema_id": self.input_bundle_schema_id,
            "input_bundle_id": self.input_bundle_id,
            "input_provenance_id": self.input_provenance_id,
            "input_resolution_policy_id": self.input_resolution_policy_id,
            "point_id": self.point_id,
            "point_label": self.point_label,
            "M_KK": self.M_KK,
            "xi_KK": self.xi_KK,
            "couplings": self.couplings.as_dict(),
            "matching": self.matching.as_dict(),
            "deltaf2": {
                "model_label": self.deltaf2.model_label,
                "input_bundle_label": self.deltaf2.input_bundle_label,
                "M_KK": self.deltaf2.M_KK,
                "xi_KK": self.deltaf2.xi_KK,
                "passes_all": self.deltaf2.passes_all,
                "max_ratio_to_bound": self.deltaf2.max_ratio_to_bound,
                "failing_reasons": list(self.deltaf2.failing_reasons),
                "ratios": self.deltaf2.as_ratio_dict(),
            },
            "verdicts": [verdict.as_dict() for verdict in self.verdicts],
            "notes": self.notes,
        }


def evaluate_modern_point(
    source: QuarkFitResult | QuarkFitSolution,
    *,
    policy: ModernPhenomenologyPolicy | None = None,
    inputs: ModernDefaultInputs | None = None,
    point_id: str | None = None,
    point_label: str | None = None,
    M_KK: float | None = None,
    xi_KK: float = DEFAULT_QUARK_XI_KK,
) -> ModernPointEvaluation:
    """Evaluate one fitted modern point against the explicit per-point contract."""
    fit_result = _coerce_fit_result(source)
    resolved_policy = default_modern_phenomenology_policy() if policy is None else policy
    resolved_inputs = default_modern_default_inputs() if inputs is None else inputs
    resolved_policy.require_system_ids(MODERN_PHENOMENOLOGY_SYSTEM_IDS)
    if not isinstance(resolved_inputs, ModernDefaultInputs):
        raise TypeError("inputs must be a ModernDefaultInputs instance")
    if abs(float(xi_KK) - resolved_inputs.qcd_metadata.xi_KK) > 1e-12:
        raise ValueError("xi_KK must match modern default input bundle qcd_metadata.xi_KK")
    resolved_m_kk = (
        default_quark_m_kk_from_lambda_ir(fit_result.point.Lambda_IR, xi_KK=xi_KK)
        if M_KK is None
        else float(M_KK)
    )
    resolved_point_label = fit_result.point.label if point_label is None else str(point_label)
    resolved_point_id = resolved_point_label if point_id is None else str(point_id)
    couplings = build_modern_point_couplings(
        fit_result,
        inputs=resolved_inputs,
        point_id=resolved_point_id,
        point_label=resolved_point_label,
        M_KK=resolved_m_kk,
        xi_KK=xi_KK,
    )
    matching = build_modern_point_matching(
        couplings,
        inputs=resolved_inputs,
    )
    deltaf2 = _deltaf2_summary_from_matching(matching, bundle=resolved_inputs)
    verdicts = []
    for system_id in MODERN_POINT_EVALUATION_SYSTEM_IDS:
        backend_summary: DeltaF2ObservableSummary = deltaf2.by_system[
            _backend_system_id_for_system(system_id)
        ]
        policy_system = _policy_system_for_system(resolved_policy, system_id)
        verdicts.append(
            ModernPointVerdict(
                system_id=system_id,
                observable_id=_observable_id_for_system(system_id),
                policy_system_id=policy_system.system_id,
                backend_system_id=_backend_system_id_for_system(system_id),
                backend_key=backend_summary.input.key,
                policy_id=policy_system.policy_id,
                policy_display_name=policy_system.display_name,
                policy_notes=policy_system.notes,
                passes=backend_summary.passes,
                ratio_to_bound=backend_summary.ratio_to_bound,
                effective_amplitude=backend_summary.effective_amplitude,
                coherent_amplitude=backend_summary.coherent_amplitude,
                bound=backend_summary.bound,
                dominant_operator=backend_summary.dominant_operator,
                dominant_operator_size=backend_summary.dominant_operator_size,
                weighted_operator_sizes=dict(backend_summary.weighted_operator_sizes),
                note=(
                    f"{system_id} verdict uses observable {_observable_id_for_system(system_id)} "
                    f"and backend key {backend_summary.input.key}"
                ),
            )
        )
    return ModernPointEvaluation(
        policy=resolved_policy,
        input_bundle_schema_id=resolved_inputs.schema_id,
        input_bundle_id=resolved_inputs.bundle_id,
        input_provenance_id=resolved_inputs.provenance_id,
        input_resolution_policy_id=resolved_inputs.source_resolution_policy_id,
        point_id=resolved_point_id,
        point_label=resolved_point_label,
        M_KK=resolved_m_kk,
        xi_KK=xi_KK,
        couplings=couplings,
        matching=matching,
        deltaf2=deltaf2,
        verdicts=tuple(verdicts),
    )


def default_modern_point_evaluation(
    source: QuarkFitResult | QuarkFitSolution,
    *,
    policy: ModernPhenomenologyPolicy | None = None,
    inputs: ModernDefaultInputs | None = None,
    point_id: str | None = None,
    point_label: str | None = None,
    M_KK: float | None = None,
    xi_KK: float = DEFAULT_QUARK_XI_KK,
) -> ModernPointEvaluation:
    """Return the deterministic modern point evaluation for one source point."""
    return evaluate_modern_point(
        source,
        policy=policy,
        inputs=inputs,
        point_id=point_id,
        point_label=point_label,
        M_KK=M_KK,
        xi_KK=xi_KK,
    )
