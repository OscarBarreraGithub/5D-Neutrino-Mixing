"""Schemaed pointwise KK-gluon couplings for the modern quark lane."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from ..couplings import QuarkMassBasisCouplings, compute_quark_kk_gluon_couplings

_SENTINEL = object()  # used to distinguish "not passed" from explicit None
from ..fit import QuarkFitResult, QuarkFitSolution
from .conventions import MODERN_LANE_ID
from .inputs import (
    MODERN_DEFAULT_ALPHA_S_POLICY_ID,
    MODERN_DEFAULT_COUPLING_POLICY_ID,
    MODERN_DEFAULT_INPUT_BUNDLE_ID,
    MODERN_DEFAULT_INPUT_PROVENANCE_ID,
    MODERN_DEFAULT_INPUTS_SCHEMA_ID,
    MODERN_DEFAULT_OPERATOR_CONVENTION_ID,
    MODERN_DEFAULT_RESOLUTION_POLICY_ID,
    ModernDefaultInputs,
    default_modern_default_inputs,
)

MODERN_POINT_COUPLINGS_SCHEMA_ID = "quarkConstraints.modern.couplings.point.v1"
MODERN_POINT_COUPLING_CONVENTION_ID = MODERN_DEFAULT_OPERATOR_CONVENTION_ID


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


def _require_complex(name: str, value: complex) -> complex:
    number = complex(value)
    if not np.isfinite(number.real) or not np.isfinite(number.imag):
        raise ValueError(f"{name} must be finite")
    return number


def _canonical_complex_matrix(name: str, value: Any) -> tuple[tuple[complex, ...], ...]:
    matrix = np.asarray(value, dtype=np.complex128)
    if matrix.shape != (3, 3):
        raise ValueError(f"{name} must have shape (3, 3)")
    if not np.all(np.isfinite(np.real(matrix))) or not np.all(np.isfinite(np.imag(matrix))):
        raise ValueError(f"{name} must contain only finite complex entries")
    return tuple(tuple(complex(entry) for entry in row) for row in matrix)


def _matrix_payload(matrix: tuple[tuple[complex, ...], ...]) -> dict[str, object]:
    return {
        "real": [[float(entry.real) for entry in row] for row in matrix],
        "imag": [[float(entry.imag) for entry in row] for row in matrix],
    }


def _matrix_from_payload(name: str, payload: Any) -> tuple[tuple[complex, ...], ...]:
    if not isinstance(payload, dict):
        raise ValueError(f"{name} must be a payload dictionary")
    real = np.asarray(payload["real"], dtype=float)
    imag = np.asarray(payload["imag"], dtype=float)
    if real.shape != (3, 3) or imag.shape != (3, 3):
        raise ValueError(f"{name} payload must contain 3x3 real and imag arrays")
    return _canonical_complex_matrix(name, real + 1j * imag)


def _coerce_fit_result(source: Any) -> Any:
    if hasattr(source, "point") and hasattr(source, "state"):
        return source
    if hasattr(source, "result") and hasattr(source.result, "point") and hasattr(source.result, "state"):
        return source.result
    raise TypeError("source must be a QuarkFitResult or QuarkFitSolution instance")


@dataclass(frozen=True)
class ModernPointCouplings:
    """Frozen modern pointwise coupling bridge for one fitted quark point."""

    schema_id: str = MODERN_POINT_COUPLINGS_SCHEMA_ID
    lane_id: str = MODERN_LANE_ID
    input_bundle_schema_id: str = MODERN_DEFAULT_INPUTS_SCHEMA_ID
    input_bundle_id: str = MODERN_DEFAULT_INPUT_BUNDLE_ID
    input_provenance_id: str = MODERN_DEFAULT_INPUT_PROVENANCE_ID
    input_resolution_policy_id: str = MODERN_DEFAULT_RESOLUTION_POLICY_ID
    qcd_metadata_id: str = ""
    alpha_s_policy_id: str = ""
    operator_convention_id: str = MODERN_POINT_COUPLING_CONVENTION_ID
    ckm_target_id: str = ""
    quark_mass_target_id: str = ""
    point_id: str = ""
    point_label: str = ""
    Lambda_IR: float = 0.0
    M_KK: float = 0.0
    xi_KK: float = 0.0
    alpha_s: float = 0.0
    g_s: float = 0.0
    g_s_4d: float = 0.0
    g_eff: float = 0.0
    g_s_multiplier: float = 0.0
    coupling_policy_id: str = MODERN_DEFAULT_COUPLING_POLICY_ID
    left_overlap: tuple[tuple[complex, ...], ...] = ()
    right_up_overlap: tuple[tuple[complex, ...], ...] = ()
    right_down_overlap: tuple[tuple[complex, ...], ...] = ()
    left_up: tuple[tuple[complex, ...], ...] = ()
    left_down: tuple[tuple[complex, ...], ...] = ()
    right_up: tuple[tuple[complex, ...], ...] = ()
    right_down: tuple[tuple[complex, ...], ...] = ()
    notes: str = (
        "Schemaed modern pointwise KK-gluon couplings. This is a versioned "
        "wrapper over the current repo mass-basis coupling formula, not a full "
        "EFT or RG package."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            _require_exact("schema_id", self.schema_id, expected=MODERN_POINT_COUPLINGS_SCHEMA_ID),
        )
        object.__setattr__(self, "lane_id", _require_exact("lane_id", self.lane_id, expected=MODERN_LANE_ID))
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
        object.__setattr__(
            self,
            "operator_convention_id",
            _require_text("operator_convention_id", self.operator_convention_id),
        )
        object.__setattr__(self, "ckm_target_id", _require_text("ckm_target_id", self.ckm_target_id))
        object.__setattr__(
            self,
            "quark_mass_target_id",
            _require_text("quark_mass_target_id", self.quark_mass_target_id),
        )
        object.__setattr__(self, "point_id", _require_text("point_id", self.point_id))
        object.__setattr__(self, "point_label", _require_text("point_label", self.point_label))
        object.__setattr__(self, "Lambda_IR", _require_positive_float("Lambda_IR", self.Lambda_IR))
        object.__setattr__(self, "M_KK", _require_positive_float("M_KK", self.M_KK))
        object.__setattr__(self, "xi_KK", _require_positive_float("xi_KK", self.xi_KK))
        object.__setattr__(self, "alpha_s", _require_positive_float("alpha_s", self.alpha_s))
        object.__setattr__(self, "g_s", _require_positive_float("g_s", self.g_s))
        object.__setattr__(self, "g_s_4d", _require_positive_float("g_s_4d", self.g_s_4d))
        object.__setattr__(self, "g_eff", _require_positive_float("g_eff", self.g_eff))
        object.__setattr__(
            self,
            "g_s_multiplier",
            _require_positive_float("g_s_multiplier", self.g_s_multiplier),
        )
        object.__setattr__(
            self,
            "coupling_policy_id",
            _require_text("coupling_policy_id", self.coupling_policy_id),
        )
        object.__setattr__(self, "left_overlap", _canonical_complex_matrix("left_overlap", self.left_overlap))
        object.__setattr__(
            self,
            "right_up_overlap",
            _canonical_complex_matrix("right_up_overlap", self.right_up_overlap),
        )
        object.__setattr__(
            self,
            "right_down_overlap",
            _canonical_complex_matrix("right_down_overlap", self.right_down_overlap),
        )
        object.__setattr__(self, "left_up", _canonical_complex_matrix("left_up", self.left_up))
        object.__setattr__(self, "left_down", _canonical_complex_matrix("left_down", self.left_down))
        object.__setattr__(self, "right_up", _canonical_complex_matrix("right_up", self.right_up))
        object.__setattr__(self, "right_down", _canonical_complex_matrix("right_down", self.right_down))
        object.__setattr__(self, "notes", _require_text("notes", self.notes))

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "lane_id": self.lane_id,
            "input_bundle_schema_id": self.input_bundle_schema_id,
            "input_bundle_id": self.input_bundle_id,
            "input_provenance_id": self.input_provenance_id,
            "input_resolution_policy_id": self.input_resolution_policy_id,
            "qcd_metadata_id": self.qcd_metadata_id,
            "alpha_s_policy_id": self.alpha_s_policy_id,
            "operator_convention_id": self.operator_convention_id,
            "ckm_target_id": self.ckm_target_id,
            "quark_mass_target_id": self.quark_mass_target_id,
            "point_id": self.point_id,
            "point_label": self.point_label,
            "Lambda_IR": self.Lambda_IR,
            "M_KK": self.M_KK,
            "xi_KK": self.xi_KK,
            "alpha_s": self.alpha_s,
            "g_s": self.g_s,
            "g_s_4d": self.g_s_4d,
            "g_eff": self.g_eff,
            "g_s_multiplier": self.g_s_multiplier,
            "coupling_policy_id": self.coupling_policy_id,
            "left_overlap": _matrix_payload(self.left_overlap),
            "right_up_overlap": _matrix_payload(self.right_up_overlap),
            "right_down_overlap": _matrix_payload(self.right_down_overlap),
            "left_up": _matrix_payload(self.left_up),
            "left_down": _matrix_payload(self.left_down),
            "right_up": _matrix_payload(self.right_up),
            "right_down": _matrix_payload(self.right_down),
            "notes": self.notes,
        }

    @classmethod
    def from_dict(cls, payload: dict[str, object]) -> "ModernPointCouplings":
        return cls(
            schema_id=str(payload["schema_id"]),
            lane_id=str(payload["lane_id"]),
            input_bundle_schema_id=str(payload["input_bundle_schema_id"]),
            input_bundle_id=str(payload["input_bundle_id"]),
            input_provenance_id=str(payload["input_provenance_id"]),
            input_resolution_policy_id=str(payload["input_resolution_policy_id"]),
            qcd_metadata_id=str(payload["qcd_metadata_id"]),
            alpha_s_policy_id=str(payload["alpha_s_policy_id"]),
            operator_convention_id=str(payload["operator_convention_id"]),
            ckm_target_id=str(payload["ckm_target_id"]),
            quark_mass_target_id=str(payload["quark_mass_target_id"]),
            point_id=str(payload["point_id"]),
            point_label=str(payload["point_label"]),
            Lambda_IR=float(payload["Lambda_IR"]),
            M_KK=float(payload["M_KK"]),
            xi_KK=float(payload["xi_KK"]),
            alpha_s=float(payload["alpha_s"]),
            g_s=float(payload["g_s"]),
            g_s_4d=float(payload["g_s_4d"]),
            g_eff=float(payload["g_eff"]),
            g_s_multiplier=float(payload["g_s_multiplier"]),
            coupling_policy_id=str(payload["coupling_policy_id"]),
            left_overlap=_matrix_from_payload("left_overlap", payload["left_overlap"]),
            right_up_overlap=_matrix_from_payload("right_up_overlap", payload["right_up_overlap"]),
            right_down_overlap=_matrix_from_payload("right_down_overlap", payload["right_down_overlap"]),
            left_up=_matrix_from_payload("left_up", payload["left_up"]),
            left_down=_matrix_from_payload("left_down", payload["left_down"]),
            right_up=_matrix_from_payload("right_up", payload["right_up"]),
            right_down=_matrix_from_payload("right_down", payload["right_down"]),
            notes=str(payload["notes"]),
        )


def build_modern_point_couplings(
    source: Any,
    *,
    inputs: ModernDefaultInputs | None = None,
    point_id: str | None = None,
    point_label: str | None = None,
    M_KK: float | None = None,
    xi_KK: float | None = None,
    g_s_star: float | None = _SENTINEL,
) -> ModernPointCouplings:
    """Build the schemaed modern coupling bridge for one fitted point.

    Parameters
    ----------
    g_s_star
        Enhanced 5D KK-gluon coupling. When omitted, the value is read from
        the input bundle's QCD metadata (default 3.0). When explicitly set
        (including ``None``), that value is forwarded directly to
        :func:`compute_quark_kk_gluon_couplings`; passing ``None`` falls back
        to the perturbative ``g_s(M_KK)``.
    """

    fit_result = _coerce_fit_result(source)
    bundle = default_modern_default_inputs() if inputs is None else inputs
    if not isinstance(bundle, ModernDefaultInputs):
        raise TypeError("inputs must be a ModernDefaultInputs instance")
    resolved_xi_kk = bundle.qcd_metadata.xi_KK if xi_KK is None else float(xi_KK)
    if abs(resolved_xi_kk - bundle.qcd_metadata.xi_KK) > 1e-12:
        raise ValueError("xi_KK must match modern input bundle qcd_metadata.xi_KK")
    using_bundle_coupling_policy = g_s_star is _SENTINEL
    resolved_g_s_star = bundle.qcd_metadata.g_s_star if using_bundle_coupling_policy else g_s_star
    from ..couplings import compute_quark_kk_gluon_couplings

    couplings = compute_quark_kk_gluon_couplings(
        fit_result,
        M_KK=M_KK,
        xi_KK=resolved_xi_kk,
        g_s_star=resolved_g_s_star,
        coupling_policy_id=(
            bundle.qcd_metadata.coupling_policy_id
            if using_bundle_coupling_policy
            else None
        ),
    )
    state = fit_result.state
    resolved_point_label = fit_result.point.label if point_label is None else str(point_label)
    resolved_point_id = resolved_point_label if point_id is None else str(point_id)
    return ModernPointCouplings(
        qcd_metadata_id=bundle.qcd_metadata.metadata_id,
        alpha_s_policy_id=bundle.qcd_metadata.alpha_s_policy_id,
        operator_convention_id=couplings.operator_convention_id,
        ckm_target_id=bundle.ckm_target.target_id,
        quark_mass_target_id=bundle.quark_mass_target.target_id,
        point_id=resolved_point_id,
        point_label=resolved_point_label,
        Lambda_IR=float(state.point.Lambda_IR),
        M_KK=float(couplings.M_KK),
        xi_KK=float(couplings.xi_KK),
        alpha_s=float(couplings.alpha_s),
        g_s=float(couplings.g_s),
        g_s_4d=float(couplings.g_s_4d),
        g_eff=float(couplings.g_eff),
        g_s_multiplier=float(couplings.g_s_multiplier),
        coupling_policy_id=couplings.coupling_policy_id,
        left_overlap=couplings.left_overlap,
        right_up_overlap=couplings.right_up_overlap,
        right_down_overlap=couplings.right_down_overlap,
        left_up=couplings.left_up,
        left_down=couplings.left_down,
        right_up=couplings.right_up,
        right_down=couplings.right_down,
    )


def default_modern_point_couplings(
    source: Any,
    *,
    inputs: ModernDefaultInputs | None = None,
    point_id: str | None = None,
    point_label: str | None = None,
    M_KK: float | None = None,
    xi_KK: float | None = None,
    g_s_star: float | None = _SENTINEL,
) -> ModernPointCouplings:
    """Return the canonical modern coupling bridge for one fitted point."""

    return build_modern_point_couplings(
        source,
        inputs=inputs,
        point_id=point_id,
        point_label=point_label,
        M_KK=M_KK,
        xi_KK=xi_KK,
        g_s_star=g_s_star,
    )


__all__ = [
    "MODERN_POINT_COUPLINGS_SCHEMA_ID",
    "MODERN_POINT_COUPLING_CONVENTION_ID",
    "ModernPointCouplings",
    "build_modern_point_couplings",
    "default_modern_point_couplings",
]
