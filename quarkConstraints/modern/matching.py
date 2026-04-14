"""Schemaed tree-level matching bridge for the modern quark lane."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from .conventions import MODERN_LANE_ID
from .couplings import (
    MODERN_POINT_COUPLINGS_SCHEMA_ID,
    ModernPointCouplings,
    build_modern_point_couplings,
)
from .inputs import (
    MODERN_DEFAULT_ALPHA_S_POLICY_ID,
    MODERN_DEFAULT_INPUT_BUNDLE_ID,
    MODERN_DEFAULT_INPUT_PROVENANCE_ID,
    MODERN_DEFAULT_INPUTS_SCHEMA_ID,
    MODERN_DEFAULT_RESOLUTION_POLICY_ID,
    ModernDefaultInputs,
    default_modern_default_inputs,
)

MODERN_POINT_MATCHING_SCHEMA_ID = "quarkConstraints.modern.matching.point.v1"
MODERN_MATCHING_SYSTEM_SCHEMA_ID = "quarkConstraints.modern.matching.system.v1"
MODERN_MATCHING_SYSTEM_IDS = ("K", "B_d", "B_s", "D0")
MODERN_POINT_MATCHING_SYSTEM_SCHEMA_ID = MODERN_MATCHING_SYSTEM_SCHEMA_ID
MODERN_POINT_MATCHING_BACKEND_ID = "quarkConstraints.modern.matching.repo_deltaf2_tree_v1_wrapper.v1"


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
    if not np.isfinite(numeric) or numeric <= 0.0:
        raise ValueError(f"{name} must be a positive finite float")
    return numeric


def _require_nonnegative_float(name: str, value: float) -> float:
    numeric = float(value)
    if not np.isfinite(numeric) or numeric < 0.0:
        raise ValueError(f"{name} must be a non-negative finite float")
    return numeric


def _require_complex(name: str, value: Any) -> complex:
    number = complex(value)
    if not np.isfinite(number.real) or not np.isfinite(number.imag):
        raise ValueError(f"{name} must be finite")
    return number


def _complex_payload(value: complex) -> dict[str, float]:
    return {"real": float(value.real), "imag": float(value.imag)}


def _complex_from_payload(name: str, payload: Any) -> complex:
    if not isinstance(payload, dict):
        raise ValueError(f"{name} must be a payload dictionary")
    return _require_complex(name, float(payload["real"]) + 1j * float(payload["imag"]))


def _system_id_for_input_id(input_id: str) -> str:
    if input_id == "epsilon_K":
        return "K"
    if input_id in {"B_d", "B_s", "D0"}:
        return input_id
    raise ValueError(f"unknown modern default input system_id {input_id!r}")


def _observable_id_for_system(system_id: str) -> str:
    if system_id == "K":
        return "epsilon_K"
    if system_id in {"B_d", "B_s", "D0"}:
        return system_id
    raise ValueError(f"system_id must be one of {MODERN_MATCHING_SYSTEM_IDS!r}")


def _backend_system_id_for_system(system_id: str) -> str:
    if system_id in {"K", "B_d", "B_s"}:
        return system_id
    if system_id == "D0":
        return "D"
    raise ValueError(f"system_id must be one of {MODERN_MATCHING_SYSTEM_IDS!r}")


def _matrix_entry(matrix: tuple[tuple[complex, ...], ...], i: int, j: int) -> complex:
    return complex(matrix[i][j])


@dataclass(frozen=True)
class ModernSystemMatching:
    """Raw tree-level Wilsons for one explicit neutral-meson system."""

    schema_id: str = MODERN_MATCHING_SYSTEM_SCHEMA_ID
    system_id: str = ""
    observable_id: str = ""
    backend_system_id: str = ""
    backend_key: str = ""
    display_name: str = ""
    sector_id: str = ""
    generations: tuple[int, int] = (0, 1)
    operator_basis_id: str = ""
    matching_scale_GeV: float = 0.0
    weight_policy_id: str = ""
    left_coupling: complex = 0.0j
    right_coupling: complex = 0.0j
    c1_vll: complex = 0.0j
    c1_vrr: complex = 0.0j
    c4_lr: complex = 0.0j
    c5_lr: complex = 0.0j
    note: str = ""

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            _require_exact("schema_id", self.schema_id, expected=MODERN_MATCHING_SYSTEM_SCHEMA_ID),
        )
        object.__setattr__(self, "system_id", _require_text("system_id", self.system_id))
        object.__setattr__(self, "observable_id", _require_text("observable_id", self.observable_id))
        object.__setattr__(self, "backend_system_id", _require_text("backend_system_id", self.backend_system_id))
        object.__setattr__(self, "backend_key", _require_text("backend_key", self.backend_key))
        object.__setattr__(self, "display_name", _require_text("display_name", self.display_name))
        object.__setattr__(self, "sector_id", _require_text("sector_id", self.sector_id))
        if self.system_id not in MODERN_MATCHING_SYSTEM_IDS:
            raise ValueError(f"system_id must be one of {MODERN_MATCHING_SYSTEM_IDS!r}")
        expected_observable = _observable_id_for_system(self.system_id)
        if self.observable_id != expected_observable:
            raise ValueError(
                f"observable_id for {self.system_id!r} must be exactly {expected_observable!r}"
            )
        expected_backend_system = _backend_system_id_for_system(self.system_id)
        if self.backend_system_id != expected_backend_system:
            raise ValueError(
                f"backend_system_id for {self.system_id!r} must be exactly {expected_backend_system!r}"
            )
        if not isinstance(self.generations, tuple) or len(self.generations) != 2:
            raise ValueError("generations must be a tuple of two indices")
        i, j = self.generations
        if not isinstance(i, int) or not isinstance(j, int) or i == j or i not in {0, 1, 2} or j not in {0, 1, 2}:
            raise ValueError("generations must be distinct integers in {0, 1, 2}")
        object.__setattr__(
            self,
            "matching_scale_GeV",
            _require_positive_float("matching_scale_GeV", self.matching_scale_GeV),
        )
        object.__setattr__(self, "operator_basis_id", _require_text("operator_basis_id", self.operator_basis_id))
        object.__setattr__(self, "weight_policy_id", _require_text("weight_policy_id", self.weight_policy_id))
        object.__setattr__(self, "left_coupling", _require_complex("left_coupling", self.left_coupling))
        object.__setattr__(self, "right_coupling", _require_complex("right_coupling", self.right_coupling))
        object.__setattr__(self, "c1_vll", _require_complex("c1_vll", self.c1_vll))
        object.__setattr__(self, "c1_vrr", _require_complex("c1_vrr", self.c1_vrr))
        object.__setattr__(self, "c4_lr", _require_complex("c4_lr", self.c4_lr))
        object.__setattr__(self, "c5_lr", _require_complex("c5_lr", self.c5_lr))
        object.__setattr__(self, "note", _require_text("note", self.note))

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "system_id": self.system_id,
            "observable_id": self.observable_id,
            "backend_system_id": self.backend_system_id,
            "backend_key": self.backend_key,
            "display_name": self.display_name,
            "sector_id": self.sector_id,
            "generations": list(self.generations),
            "operator_basis_id": self.operator_basis_id,
            "matching_scale_GeV": self.matching_scale_GeV,
            "weight_policy_id": self.weight_policy_id,
            "left_coupling": _complex_payload(self.left_coupling),
            "right_coupling": _complex_payload(self.right_coupling),
            "c1_vll": _complex_payload(self.c1_vll),
            "c1_vrr": _complex_payload(self.c1_vrr),
            "c4_lr": _complex_payload(self.c4_lr),
            "c5_lr": _complex_payload(self.c5_lr),
            "note": self.note,
        }

    @classmethod
    def from_dict(cls, payload: dict[str, object]) -> "ModernSystemMatching":
        return cls(
            schema_id=str(payload["schema_id"]),
            system_id=str(payload["system_id"]),
            observable_id=str(payload["observable_id"]),
            backend_system_id=str(payload["backend_system_id"]),
            backend_key=str(payload["backend_key"]),
            display_name=str(payload["display_name"]),
            sector_id=str(payload["sector_id"]),
            generations=tuple(int(item) for item in payload["generations"]),
            operator_basis_id=str(payload["operator_basis_id"]),
            matching_scale_GeV=float(payload["matching_scale_GeV"]),
            weight_policy_id=str(payload["weight_policy_id"]),
            left_coupling=_complex_from_payload("left_coupling", payload["left_coupling"]),
            right_coupling=_complex_from_payload("right_coupling", payload["right_coupling"]),
            c1_vll=_complex_from_payload("c1_vll", payload["c1_vll"]),
            c1_vrr=_complex_from_payload("c1_vrr", payload["c1_vrr"]),
            c4_lr=_complex_from_payload("c4_lr", payload["c4_lr"]),
            c5_lr=_complex_from_payload("c5_lr", payload["c5_lr"]),
            note=str(payload["note"]),
        )


@dataclass(frozen=True)
class ModernPointMatching:
    """Frozen modern pointwise matching bridge for one fitted quark point."""

    schema_id: str = MODERN_POINT_MATCHING_SCHEMA_ID
    lane_id: str = MODERN_LANE_ID
    couplings_schema_id: str = MODERN_POINT_COUPLINGS_SCHEMA_ID
    input_bundle_schema_id: str = MODERN_DEFAULT_INPUTS_SCHEMA_ID
    input_bundle_id: str = MODERN_DEFAULT_INPUT_BUNDLE_ID
    input_provenance_id: str = MODERN_DEFAULT_INPUT_PROVENANCE_ID
    input_resolution_policy_id: str = MODERN_DEFAULT_RESOLUTION_POLICY_ID
    qcd_metadata_id: str = ""
    alpha_s_policy_id: str = MODERN_DEFAULT_ALPHA_S_POLICY_ID
    operator_basis_id: str = ""
    weight_policy_id: str = ""
    point_id: str = ""
    point_label: str = ""
    Lambda_IR: float = 0.0
    M_KK: float = 0.0
    xi_KK: float = 0.0
    alpha_s: float = 0.0
    g_s: float = 0.0
    system_matches: tuple[ModernSystemMatching, ...] = field(default_factory=tuple)
    notes: str = (
        "Schemaed modern tree-level matching bridge. This is a versioned "
        "wrapper over the current compact repo KK-gluon formulas, not a full "
        "operator-complete EFT or RG package."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            _require_exact("schema_id", self.schema_id, expected=MODERN_POINT_MATCHING_SCHEMA_ID),
        )
        object.__setattr__(self, "lane_id", _require_exact("lane_id", self.lane_id, expected=MODERN_LANE_ID))
        object.__setattr__(
            self,
            "couplings_schema_id",
            _require_exact(
                "couplings_schema_id",
                self.couplings_schema_id,
                expected=MODERN_POINT_COUPLINGS_SCHEMA_ID,
            ),
        )
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
            _require_exact("input_bundle_id", self.input_bundle_id, expected=MODERN_DEFAULT_INPUT_BUNDLE_ID),
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
        object.__setattr__(self, "qcd_metadata_id", _require_text("qcd_metadata_id", self.qcd_metadata_id))
        object.__setattr__(self, "alpha_s_policy_id", _require_text("alpha_s_policy_id", self.alpha_s_policy_id))
        object.__setattr__(self, "operator_basis_id", _require_text("operator_basis_id", self.operator_basis_id))
        object.__setattr__(self, "weight_policy_id", _require_text("weight_policy_id", self.weight_policy_id))
        object.__setattr__(self, "point_id", _require_text("point_id", self.point_id))
        object.__setattr__(self, "point_label", _require_text("point_label", self.point_label))
        object.__setattr__(self, "Lambda_IR", _require_positive_float("Lambda_IR", self.Lambda_IR))
        object.__setattr__(self, "M_KK", _require_positive_float("M_KK", self.M_KK))
        object.__setattr__(self, "xi_KK", _require_positive_float("xi_KK", self.xi_KK))
        object.__setattr__(self, "alpha_s", _require_positive_float("alpha_s", self.alpha_s))
        object.__setattr__(self, "g_s", _require_positive_float("g_s", self.g_s))
        if not isinstance(self.system_matches, tuple):
            raise ValueError("system_matches must be a tuple")
        normalized = tuple(self.system_matches)
        if len(normalized) != len(MODERN_MATCHING_SYSTEM_IDS):
            raise ValueError("system_matches must contain the frozen modern system set")
        for item in normalized:
            if not isinstance(item, ModernSystemMatching):
                raise ValueError("system_matches must contain only ModernSystemMatching instances")
        if tuple(item.system_id for item in normalized) != MODERN_MATCHING_SYSTEM_IDS:
            raise ValueError("system_matches must preserve the frozen modern system order")
        object.__setattr__(self, "system_matches", normalized)
        object.__setattr__(self, "notes", _require_text("notes", self.notes))

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "lane_id": self.lane_id,
            "couplings_schema_id": self.couplings_schema_id,
            "input_bundle_schema_id": self.input_bundle_schema_id,
            "input_bundle_id": self.input_bundle_id,
            "input_provenance_id": self.input_provenance_id,
            "input_resolution_policy_id": self.input_resolution_policy_id,
            "qcd_metadata_id": self.qcd_metadata_id,
            "alpha_s_policy_id": self.alpha_s_policy_id,
            "operator_basis_id": self.operator_basis_id,
            "weight_policy_id": self.weight_policy_id,
            "point_id": self.point_id,
            "point_label": self.point_label,
            "Lambda_IR": self.Lambda_IR,
            "M_KK": self.M_KK,
            "xi_KK": self.xi_KK,
            "alpha_s": self.alpha_s,
            "g_s": self.g_s,
            "system_matches": [item.as_dict() for item in self.system_matches],
            "notes": self.notes,
        }

    @property
    def systems(self) -> tuple[ModernSystemMatching, ...]:
        return self.system_matches

    def system(self, system_id: str) -> ModernSystemMatching:
        requested = _require_text("system_id", system_id)
        for system in self.system_matches:
            if system.system_id == requested:
                return system
        raise ValueError(f"system_id must be one of {MODERN_MATCHING_SYSTEM_IDS!r}")

    @classmethod
    def from_dict(cls, payload: dict[str, object]) -> "ModernPointMatching":
        return cls(
            schema_id=str(payload["schema_id"]),
            lane_id=str(payload["lane_id"]),
            couplings_schema_id=str(payload["couplings_schema_id"]),
            input_bundle_schema_id=str(payload["input_bundle_schema_id"]),
            input_bundle_id=str(payload["input_bundle_id"]),
            input_provenance_id=str(payload["input_provenance_id"]),
            input_resolution_policy_id=str(payload["input_resolution_policy_id"]),
            qcd_metadata_id=str(payload["qcd_metadata_id"]),
            alpha_s_policy_id=str(payload["alpha_s_policy_id"]),
            operator_basis_id=str(payload["operator_basis_id"]),
            weight_policy_id=str(payload["weight_policy_id"]),
            point_id=str(payload["point_id"]),
            point_label=str(payload["point_label"]),
            Lambda_IR=float(payload["Lambda_IR"]),
            M_KK=float(payload["M_KK"]),
            xi_KK=float(payload["xi_KK"]),
            alpha_s=float(payload["alpha_s"]),
            g_s=float(payload["g_s"]),
            system_matches=tuple(
                ModernSystemMatching.from_dict(item) for item in payload["system_matches"]
            ),
            notes=str(payload["notes"]),
        )


def build_modern_point_matching(
    source: Any,
    *,
    inputs: ModernDefaultInputs | None = None,
    M_KK: float | None = None,
    xi_KK: float | None = None,
) -> ModernPointMatching:
    """Build the schemaed modern matching bridge for one fitted point."""

    bundle = default_modern_default_inputs() if inputs is None else inputs
    if not isinstance(bundle, ModernDefaultInputs):
        raise TypeError("inputs must be a ModernDefaultInputs instance")
    couplings = (
        source
        if isinstance(source, ModernPointCouplings)
        else build_modern_point_couplings(source, inputs=bundle, M_KK=M_KK, xi_KK=xi_KK)
    )
    if not isinstance(couplings, ModernPointCouplings):
        raise TypeError("source must be a ModernPointCouplings, QuarkFitResult, or QuarkFitSolution instance")
    prefactor = 1.0 / (couplings.M_KK**2)
    system_matches: list[ModernSystemMatching] = []
    for system_input in bundle.neutral_meson_inputs:
        i, j = system_input.generations
        if system_input.sector_id == "down":
            left = _matrix_entry(couplings.left_down, i, j)
            right = _matrix_entry(couplings.right_down, i, j)
        else:
            left = _matrix_entry(couplings.left_up, i, j)
            right = _matrix_entry(couplings.right_up, i, j)
        system_matches.append(
            ModernSystemMatching(
                system_id=_system_id_for_input_id(system_input.system_id),
                observable_id=system_input.system_id,
                backend_system_id=_backend_system_id_for_system(
                    _system_id_for_input_id(system_input.system_id)
                ),
                backend_key=system_input.backend_key,
                display_name=system_input.display_name,
                sector_id=system_input.sector_id,
                generations=system_input.generations,
                operator_basis_id=bundle.operator_weight_policy.operator_convention_id,
                matching_scale_GeV=couplings.M_KK,
                weight_policy_id=bundle.operator_weight_policy.policy_id,
                left_coupling=left,
                right_coupling=right,
                c1_vll=left * left * prefactor / 6.0,
                c1_vrr=right * right * prefactor / 6.0,
                c4_lr=-(left * right) * prefactor,
                c5_lr=(left * right) * prefactor / 3.0,
                note=f"{system_input.note} Observable {system_input.system_id}.",
            )
        )
    return ModernPointMatching(
        qcd_metadata_id=bundle.qcd_metadata.metadata_id,
        alpha_s_policy_id=bundle.qcd_metadata.alpha_s_policy_id,
        operator_basis_id=bundle.operator_weight_policy.operator_convention_id,
        weight_policy_id=bundle.operator_weight_policy.policy_id,
        point_id=couplings.point_id,
        point_label=couplings.point_label,
        Lambda_IR=couplings.Lambda_IR,
        M_KK=couplings.M_KK,
        xi_KK=couplings.xi_KK,
        alpha_s=couplings.alpha_s,
        g_s=couplings.g_s,
        system_matches=tuple(system_matches),
    )


def default_modern_point_matching(
    source: Any,
    *,
    inputs: ModernDefaultInputs | None = None,
    M_KK: float | None = None,
    xi_KK: float | None = None,
) -> ModernPointMatching:
    """Return the canonical modern matching bridge for one fitted point."""

    return build_modern_point_matching(source, inputs=inputs, M_KK=M_KK, xi_KK=xi_KK)


def matching_to_deltaf2_summary(
    matching: ModernPointMatching,
    *,
    inputs: ModernDefaultInputs | None = None,
) -> Any:
    """Project raw modern matching data into the current exclusion surrogate."""

    bundle = default_modern_default_inputs() if inputs is None else inputs
    if not isinstance(bundle, ModernDefaultInputs):
        raise TypeError("inputs must be a ModernDefaultInputs instance")
    weights = bundle.operator_weight_policy
    input_by_system = {
        _system_id_for_input_id(system_input.system_id): system_input
        for system_input in bundle.neutral_meson_inputs
    }
    from ..deltaf2 import (
        DeltaF2ConstraintSummary,
        DeltaF2Input,
        DeltaF2ObservableSummary,
        DeltaF2WilsonCoefficients,
    )

    observables: list[DeltaF2ObservableSummary] = []
    for system in matching.system_matches:
        system_input = input_by_system[system.system_id]
        deltaf2_input = DeltaF2Input(
            key=system.backend_key,
            display_name=system.display_name,
            column_name=system_input.column_name,
            reject_reason=system_input.reject_reason,
            sector=system.sector_id,
            generations=system.generations,
            bound=system_input.bound,
            ll_weight=weights.ll_weight,
            rr_weight=weights.rr_weight,
            lr1_weight=weights.lr1_weight,
            lr2_weight=weights.lr2_weight,
            reference_scale=weights.reference_scale_GeV,
            note=system.note,
        )
        wilsons = DeltaF2WilsonCoefficients(
            input=deltaf2_input,
            M_KK=matching.M_KK,
            matching_scale=system.matching_scale_GeV,
            left_coupling=system.left_coupling,
            right_coupling=system.right_coupling,
            c1_vll=system.c1_vll,
            c1_vrr=system.c1_vrr,
            c4_lr=system.c4_lr,
            c5_lr=system.c5_lr,
        )
        weighted = {
            "C1_VLL": weights.ll_weight * system.c1_vll,
            "C1_VRR": weights.rr_weight * system.c1_vrr,
            "C4_LR": weights.lr1_weight * system.c4_lr,
            "C5_LR": weights.lr2_weight * system.c5_lr,
        }
        operator_sizes = {
            name: float(weights.reference_scale_GeV**2 * abs(value))
            for name, value in weighted.items()
        }
        dominant_operator = max(operator_sizes, key=operator_sizes.get)
        coherent_amplitude = float(weights.reference_scale_GeV**2 * abs(sum(weighted.values())))
        effective_amplitude = float(operator_sizes[dominant_operator])
        ratio_to_bound = float(effective_amplitude / system_input.bound)
        observables.append(
            DeltaF2ObservableSummary(
                input=deltaf2_input,
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
        M_KK=matching.M_KK,
        xi_KK=matching.xi_KK,
        observables=tuple(observables),
    )


__all__ = [
    "MODERN_MATCHING_SYSTEM_IDS",
    "MODERN_MATCHING_SYSTEM_SCHEMA_ID",
    "MODERN_POINT_MATCHING_BACKEND_ID",
    "MODERN_POINT_MATCHING_SCHEMA_ID",
    "MODERN_POINT_MATCHING_SYSTEM_SCHEMA_ID",
    "ModernPointMatching",
    "ModernPointMatchingSystem",
    "ModernSystemMatching",
    "build_modern_point_matching",
    "default_modern_point_matching",
    "matching_to_deltaf2_summary",
]

ModernPointMatchingSystem = ModernSystemMatching
