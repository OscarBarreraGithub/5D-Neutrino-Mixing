"""Paper-owned LO QCD running for the 0710.1869 ``Delta F = 2`` EFT slice.

This module intentionally implements only the subset that is frozen by the
current paper-mode contracts:

- full LO thresholded running for the paper basis
- a frozen paper ``O4/O5`` to BMU ``Q1^LR/Q2^LR`` map recorded under the
  BMU NDR-MS scheme contract
- BMU LO LR evolution retained as an audit object in the BMU LR basis
- active LR running on the public paper path through
  ``paper -> BMU -> BMU LO RG -> paper``
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from qcd import alpha_s
from qcd.beta_function import beta_0

from ..conventions import PAPER_0710_1869_MODE_ID, PAPER_0710_1869_PAPER_ID
from ..validation import (
    require_known_schema_id,
    require_member,
    require_nonempty_identifier,
    require_positive_finite,
)
from .matching_kkgluon import (
    default_paper_0710_1869_kaon_matching,
    default_paper_0710_1869_kaon_matching_summary,
)
from .operators import (
    PAPER_0710_1869_DELTAF2_MATCHING_ID,
    PAPER_0710_1869_DELTAF2_OPERATOR_BASIS_ID,
    PAPER_0710_1869_DELTAF2_OPERATOR_DEFINITION_Q4_LR_ID,
    PAPER_0710_1869_DELTAF2_OPERATOR_DEFINITION_Q5_LR_ID,
    PAPER_0710_1869_DELTAF2_OPERATOR_METADATA_COMPATIBILITY_NOTE,
    PAPER_0710_1869_DELTAF2_OPERATOR_NORMALIZATION_ID,
    PAPER_0710_1869_DELTAF2_OPERATOR_ORDER,
    PAPER_0710_1869_DELTAF2_PROJECTOR_NORMALIZATION_ID,
    PAPER_0710_1869_DELTAF2_Q1_VLL,
    PAPER_0710_1869_DELTAF2_Q1_VRR,
    PAPER_0710_1869_DELTAF2_Q4_LR,
    PAPER_0710_1869_DELTAF2_Q5_LR,
    PAPER_0710_1869_DELTAF2_RENORMALIZATION_SCHEME_ID,
    PAPER_0710_1869_DELTAF2_WILSON_SCHEMA_ID,
    Paper07101869DeltaF2WilsonCoefficients,
    Paper07101869DeltaF2WilsonContract,
)
from .rg_inputs import (
    PAPER_0710_1869_DELTAF2_DEFAULT_MU_HAD_GEV,
    PAPER_0710_1869_DELTAF2_RG_ADM_SOURCE_ID,
    PAPER_0710_1869_DELTAF2_RG_ALPHA_S_POLICY_ID,
    PAPER_0710_1869_DELTAF2_RG_LR_BASIS_CONTRACT_ID,
    PAPER_0710_1869_DELTAF2_RG_LR_BASIS_CONTRACT_SCHEMA_ID,
    PAPER_0710_1869_DELTAF2_RG_LR_BASIS_DIRECTION_ID,
    PAPER_0710_1869_DELTAF2_RG_LR_BASIS_STATUS_ID,
    PAPER_0710_1869_DELTAF2_RG_LR_BMU_BASIS_ID,
    PAPER_0710_1869_DELTAF2_RG_LR_BMU_OPERATOR_ORDER,
    PAPER_0710_1869_DELTAF2_RG_LR_FIERZ_VALIDITY_STATEMENT_ID,
    PAPER_0710_1869_DELTAF2_RG_LR_RUNNING_BRIDGE_ID,
    PAPER_0710_1869_DELTAF2_RG_ORDER_ID,
    PAPER_0710_1869_DELTAF2_RG_SCALE_SEMANTICS_ID,
    PAPER_0710_1869_DELTAF2_RG_SCHEME_ID,
    PAPER_0710_1869_DELTAF2_RG_SUPPORTED_OPERATOR_SUBSET_ID,
    PAPER_0710_1869_DELTAF2_RG_THRESHOLD_POLICY_ID,
    default_paper_0710_1869_rg_thresholds,
)

PAPER_0710_1869_DELTAF2_RG_CONTRACT_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.rg_contract.v1"
)
PAPER_0710_1869_DELTAF2_LO_ADM_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.lo_adm.v1"
)
PAPER_0710_1869_DELTAF2_RG_SEGMENT_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.rg_segment.v1"
)
PAPER_0710_1869_DELTAF2_RG_EVOLUTION_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.rg_evolution.v1"
)
PAPER_0710_1869_DELTAF2_RG_SUMMARY_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.rg_summary.v1"
)
PAPER_0710_1869_DELTAF2_RG_WILSON_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.rg_wilsons.v1"
)
PAPER_0710_1869_DELTAF2_RG_RESULT_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.rg_result.v1"
)

_SUPPORTED_OPERATORS = PAPER_0710_1869_DELTAF2_OPERATOR_ORDER
_UNSUPPORTED_OPERATORS: tuple[str, ...] = ()
_ZERO_TOLERANCE = 1e-30
_GAMMA0_VLL = 4.0
_GAMMA0_VRR = 4.0
_GAMMA0_BMU_LR = (
    (2.0, 12.0),
    (0.0, -16.0),
)
_PAPER_LR_OPERATOR_ORDER = (
    PAPER_0710_1869_DELTAF2_Q4_LR,
    PAPER_0710_1869_DELTAF2_Q5_LR,
)
_PAPER_LR_OPERATOR_DEFINITION_IDS = (
    PAPER_0710_1869_DELTAF2_OPERATOR_DEFINITION_Q4_LR_ID,
    PAPER_0710_1869_DELTAF2_OPERATOR_DEFINITION_Q5_LR_ID,
)
_PAPER_TO_BMU_LR_OPERATOR_MAP = (
    (0.0, -2.0),
    (1.0, 0.0),
)
_BMU_TO_PAPER_LR_OPERATOR_MAP = (
    (0.0, 1.0),
    (-0.5, 0.0),
)
_PAPER_TO_BMU_LR_WILSON_MAP = (
    (0.0, -0.5),
    (1.0, 0.0),
)
_BMU_TO_PAPER_LR_WILSON_MAP = (
    (0.0, 1.0),
    (-2.0, 0.0),
)
_LR_BASIS_SOURCE_ID = (
    "becirevic-villadoro.hep-lat-0408029.plus.ciuchini.hep-ph-9711402.plus."
    "bmu.hep-ph-0005183.lr_map.ndr_ms.frozen.v3"
)
_LR_BASIS_SOURCE_CITATION = (
    "Becirevic and Villadoro, hep-lat/0408029 (O4/O5 scalar LR basis used in the "
    "DeltaF=2 literature and scheme caveats) + Ciuchini et al., hep-ph/9711402 "
    "(notation cross-check for the O4/O5 scalar LR basis) + Buras, Misiak, Urban, "
    "hep-ph/0005183 (BMU Q1^LR/Q2^LR basis, LO LR ADM, and NDR-MS LR-sector "
    "Fierz-validity caveats)"
)
_LR_BASIS_SOURCE_LOCATOR = (
    "Freeze the paper O4/O5 scalar LR basis used in the DeltaF=2 literature onto the "
    "BMU Q1^LR/Q2^LR basis only under the BMU NDR-MS LR-sector Fierz-validity "
    "statement. Becirevic-Villadoro anchors the paper O4/O5 convention and notes "
    "scheme dependence of Dirac-basis relations; Ciuchini is retained as a notation "
    "cross-check; BMU anchors the LR basis, LO ADM block, and Fierz-validity caveat."
)
_LR_OPERATOR_MAP_FORMULA = "Q_BMU = M_op Q_paper"
_LR_WILSON_MAP_FORMULA = "C_BMU = (M_op^{-1})^T C_paper"
_LR_FIERZ_VALIDITY_NOTE = (
    "This map is frozen only under the BMU NDR-MS LR-sector Fierz-validity statement "
    "used to identify the paper O4/O5 scalar LR basis with the BMU Q1^LR/Q2^LR "
    "basis. It is not a general scheme-independent D-dimensional operator identity."
)


def _require_finite_complex(name: str, value: complex) -> complex:
    complex_value = complex(value)
    if not math.isfinite(complex_value.real) or not math.isfinite(complex_value.imag):
        raise ValueError(f"{name} must be finite")
    real_part = 0.0 if complex_value.real == 0.0 else float(complex_value.real)
    imag_part = 0.0 if complex_value.imag == 0.0 else float(complex_value.imag)
    return complex(real_part, imag_part)


def _complex_as_dict(value: complex) -> dict[str, float]:
    complex_value = _require_finite_complex("value", value)
    return {"real": float(complex_value.real), "imag": float(complex_value.imag)}


def _normalize_thresholds(
    thresholds: tuple[tuple[float, int, int], ...] | list[tuple[float, int, int]] | None,
) -> tuple[tuple[float, int, int], ...]:
    if thresholds is None:
        return default_paper_0710_1869_rg_thresholds()
    normalized: list[tuple[float, int, int]] = []
    for mass, nf_below, nf_above in thresholds:
        numeric_mass = require_positive_finite("threshold_mass_GeV", mass)
        if not isinstance(nf_below, int) or not isinstance(nf_above, int):
            raise ValueError("threshold flavor counts must be integers")
        if nf_above != nf_below + 1:
            raise ValueError("threshold flavor counts must increase by one across each threshold")
        normalized.append((numeric_mass, nf_below, nf_above))
    return tuple(sorted(normalized, key=lambda item: item[0]))


def _n_f_at_scale(
    mu_GeV: float,
    thresholds: tuple[tuple[float, int, int], ...],
) -> int:
    if not thresholds:
        return 5
    n_f = thresholds[0][1]
    for mass, _, n_f_above in thresholds:
        if mu_GeV >= mass:
            n_f = n_f_above
        else:
            break
    return n_f


def _matrix_multiply(
    left: tuple[tuple[complex, ...], ...],
    right: tuple[tuple[complex, ...], ...],
) -> tuple[tuple[complex, ...], ...]:
    if len(left[0]) != len(right):
        raise ValueError("matrix dimensions are incompatible")
    rows = len(left)
    cols = len(right[0])
    inner = len(right)
    result: list[tuple[complex, ...]] = []
    for row_index in range(rows):
        row: list[complex] = []
        for column_index in range(cols):
            entry = 0.0 + 0.0j
            for inner_index in range(inner):
                entry += left[row_index][inner_index] * right[inner_index][column_index]
            row.append(_require_finite_complex("matrix_entry", entry))
        result.append(tuple(row))
    return tuple(result)


def _matrix_vector_multiply(
    matrix: tuple[tuple[complex, ...], ...],
    vector: tuple[complex, ...],
) -> tuple[complex, ...]:
    if len(matrix[0]) != len(vector):
        raise ValueError("matrix and vector dimensions are incompatible")
    result: list[complex] = []
    for row in matrix:
        entry = 0.0 + 0.0j
        for row_entry, vector_entry in zip(row, vector, strict=True):
            entry += row_entry * vector_entry
        result.append(_require_finite_complex("vector_entry", entry))
    return tuple(result)


def _identity_matrix(size: int) -> tuple[tuple[complex, ...], ...]:
    return tuple(
        tuple(
            1.0 + 0.0j if row_index == column_index else 0.0 + 0.0j
            for column_index in range(size)
        )
        for row_index in range(size)
    )


def _matrix_as_dict(
    matrix: tuple[tuple[complex, ...], ...],
) -> list[list[dict[str, float]]]:
    return [[_complex_as_dict(value) for value in row] for row in matrix]


def _real_matrix_as_dict(
    matrix: tuple[tuple[float, float], tuple[float, float]],
) -> list[list[float]]:
    return [[float(value) for value in row] for row in matrix]


def _normalize_real_matrix_2x2(
    name: str,
    matrix: tuple[tuple[float, float], tuple[float, float]]
    | tuple[tuple[int, int], tuple[int, int]],
) -> tuple[tuple[float, float], tuple[float, float]]:
    if len(matrix) != 2 or any(len(row) != 2 for row in matrix):
        raise ValueError(f"{name} must be 2x2")
    normalized: list[tuple[float, float]] = []
    for row in matrix:
        normalized.append(tuple(float(value) for value in row))
    return (normalized[0], normalized[1])


def _complex_matrix_from_real(
    matrix: tuple[tuple[float, float], tuple[float, float]],
) -> tuple[tuple[complex, complex], tuple[complex, complex]]:
    return tuple(
        tuple(_require_finite_complex("real_matrix_entry", complex(value, 0.0)) for value in row)
        for row in matrix
    )


def _vector_as_dict(
    vector: tuple[complex, ...],
    operator_order: tuple[str, ...],
) -> dict[str, dict[str, float]]:
    return {
        operator_name: _complex_as_dict(value)
        for operator_name, value in zip(operator_order, vector, strict=True)
    }


def _closed_form_bmu_lr_segment_matrix(
    *,
    alpha_start: float,
    alpha_end: float,
    n_f: int,
) -> tuple[tuple[complex, complex], tuple[complex, complex]]:
    beta0_value = float(beta_0(n_f))
    gamma11 = _GAMMA0_BMU_LR[0][0]
    gamma12 = _GAMMA0_BMU_LR[0][1]
    gamma22 = _GAMMA0_BMU_LR[1][1]
    exponent_1 = gamma11 / (2.0 * beta0_value)
    exponent_2 = gamma22 / (2.0 * beta0_value)
    eta = alpha_start / alpha_end
    eta_a1 = eta**exponent_1
    eta_a2 = eta**exponent_2
    off_diagonal = (gamma12 / (gamma11 - gamma22)) * (eta_a1 - eta_a2)
    return (
        (
            _require_finite_complex("u11", eta_a1 + 0.0j),
            0.0 + 0.0j,
        ),
        (
            _require_finite_complex("u21", off_diagonal + 0.0j),
            _require_finite_complex("u22", eta_a2 + 0.0j),
        ),
    )


def _build_segment_boundaries(
    *,
    mu_high_GeV: float,
    mu_low_GeV: float,
    thresholds: tuple[tuple[float, int, int], ...],
) -> tuple[float, ...]:
    if mu_high_GeV == mu_low_GeV:
        return (mu_high_GeV, mu_low_GeV)
    direction_descending = mu_high_GeV > mu_low_GeV
    masses = [
        mass
        for mass, _, _ in thresholds
        if min(mu_high_GeV, mu_low_GeV) < mass < max(mu_high_GeV, mu_low_GeV)
    ]
    masses.sort(reverse=direction_descending)
    return (mu_high_GeV, *masses, mu_low_GeV)


@dataclass(frozen=True)
class Paper07101869DeltaF2LRBasisContract:
    """Frozen LR-map contract between the paper and BMU LR operator orderings."""

    contract_id: str = PAPER_0710_1869_DELTAF2_RG_LR_BASIS_CONTRACT_ID
    paper_operator_basis_id: str = PAPER_0710_1869_DELTAF2_OPERATOR_BASIS_ID
    paper_operator_order: tuple[str, str] = _PAPER_LR_OPERATOR_ORDER
    paper_operator_definition_ids: tuple[str, str] = _PAPER_LR_OPERATOR_DEFINITION_IDS
    bmu_lr_basis_id: str = PAPER_0710_1869_DELTAF2_RG_LR_BMU_BASIS_ID
    bmu_lr_operator_order: tuple[str, str] = PAPER_0710_1869_DELTAF2_RG_LR_BMU_OPERATOR_ORDER
    map_direction_id: str = PAPER_0710_1869_DELTAF2_RG_LR_BASIS_DIRECTION_ID
    operator_map_formula: str = _LR_OPERATOR_MAP_FORMULA
    paper_to_bmu_operator_map_matrix: tuple[tuple[float, float], tuple[float, float]] = (
        _PAPER_TO_BMU_LR_OPERATOR_MAP
    )
    bmu_to_paper_operator_map_matrix: tuple[tuple[float, float], tuple[float, float]] = (
        _BMU_TO_PAPER_LR_OPERATOR_MAP
    )
    wilson_map_formula: str = _LR_WILSON_MAP_FORMULA
    paper_to_bmu_wilson_map_matrix: tuple[tuple[float, float], tuple[float, float]] = (
        _PAPER_TO_BMU_LR_WILSON_MAP
    )
    bmu_to_paper_wilson_map_matrix: tuple[tuple[float, float], tuple[float, float]] = (
        _BMU_TO_PAPER_LR_WILSON_MAP
    )
    renormalization_scheme_id: str = PAPER_0710_1869_DELTAF2_RG_SCHEME_ID
    operator_normalization_id: str = PAPER_0710_1869_DELTAF2_OPERATOR_NORMALIZATION_ID
    projector_normalization_id: str = PAPER_0710_1869_DELTAF2_PROJECTOR_NORMALIZATION_ID
    fierz_validity_statement_id: str = (
        PAPER_0710_1869_DELTAF2_RG_LR_FIERZ_VALIDITY_STATEMENT_ID
    )
    fierz_validity_note: str = _LR_FIERZ_VALIDITY_NOTE
    source_id: str = _LR_BASIS_SOURCE_ID
    source_citation: str = _LR_BASIS_SOURCE_CITATION
    source_locator: str = _LR_BASIS_SOURCE_LOCATOR
    status_id: str = PAPER_0710_1869_DELTAF2_RG_LR_BASIS_STATUS_ID
    lr_definitions_frozen: bool = True
    mapping_matrix_frozen: bool = True
    lr_running_activated: bool = True
    schema_id: str = PAPER_0710_1869_DELTAF2_RG_LR_BASIS_CONTRACT_SCHEMA_ID
    metadata_compatibility_note: str = (
        PAPER_0710_1869_DELTAF2_OPERATOR_METADATA_COMPATIBILITY_NOTE
    )
    notes: str = (
        "Paper-mode Q4_LR/Q5_LR are frozen as the O4/O5 scalar LR basis used in the "
        "DeltaF=2 literature. The exact paper-to-BMU map is now frozen under the BMU "
        "NDR-MS LR-sector Fierz-validity statement, and the canonical LR status now "
        "means LR RG is active on the public paper path through the frozen BMU "
        "conjugation bridge, custom LR hadronic, custom LR-only observable, and "
        "custom combined Q1+LR observable surfaces are active under exact "
        "scheme/scale alignment, and the default/exported observable path "
        "intentionally remains Q1-only."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_DELTAF2_RG_LR_BASIS_CONTRACT_SCHEMA_ID,
            ),
        )
        for field_name in (
            "contract_id",
            "paper_operator_basis_id",
            "bmu_lr_basis_id",
            "map_direction_id",
            "operator_map_formula",
            "wilson_map_formula",
            "renormalization_scheme_id",
            "operator_normalization_id",
            "projector_normalization_id",
            "fierz_validity_statement_id",
            "fierz_validity_note",
            "source_id",
            "source_citation",
            "source_locator",
            "status_id",
            "metadata_compatibility_note",
            "notes",
        ):
            object.__setattr__(
                self,
                field_name,
                require_nonempty_identifier(field_name, getattr(self, field_name)),
            )
        if tuple(self.paper_operator_order) != _PAPER_LR_OPERATOR_ORDER:
            raise ValueError("paper_operator_order must match the frozen paper LR ordering")
        if tuple(self.paper_operator_definition_ids) != _PAPER_LR_OPERATOR_DEFINITION_IDS:
            raise ValueError(
                "paper_operator_definition_ids must match the frozen paper O4/O5 scalar "
                "LR operator definitions"
            )
        if tuple(self.bmu_lr_operator_order) != PAPER_0710_1869_DELTAF2_RG_LR_BMU_OPERATOR_ORDER:
            raise ValueError("bmu_lr_operator_order must match the frozen BMU LR ordering")
        object.__setattr__(
            self,
            "paper_to_bmu_operator_map_matrix",
            _normalize_real_matrix_2x2(
                "paper_to_bmu_operator_map_matrix",
                self.paper_to_bmu_operator_map_matrix,
            ),
        )
        object.__setattr__(
            self,
            "bmu_to_paper_operator_map_matrix",
            _normalize_real_matrix_2x2(
                "bmu_to_paper_operator_map_matrix",
                self.bmu_to_paper_operator_map_matrix,
            ),
        )
        object.__setattr__(
            self,
            "paper_to_bmu_wilson_map_matrix",
            _normalize_real_matrix_2x2(
                "paper_to_bmu_wilson_map_matrix",
                self.paper_to_bmu_wilson_map_matrix,
            ),
        )
        object.__setattr__(
            self,
            "bmu_to_paper_wilson_map_matrix",
            _normalize_real_matrix_2x2(
                "bmu_to_paper_wilson_map_matrix",
                self.bmu_to_paper_wilson_map_matrix,
            ),
        )
        if self.paper_to_bmu_operator_map_matrix != _PAPER_TO_BMU_LR_OPERATOR_MAP:
            raise ValueError(
                "paper_to_bmu_operator_map_matrix must match the frozen exact operator-space "
                "map Q_BMU = M_op Q_paper"
            )
        if self.bmu_to_paper_operator_map_matrix != _BMU_TO_PAPER_LR_OPERATOR_MAP:
            raise ValueError(
                "bmu_to_paper_operator_map_matrix must match the frozen exact inverse "
                "operator-space map"
            )
        if self.paper_to_bmu_wilson_map_matrix != _PAPER_TO_BMU_LR_WILSON_MAP:
            raise ValueError(
                "paper_to_bmu_wilson_map_matrix must match the frozen exact Wilson-space "
                "map C_BMU = (M_op^{-1})^T C_paper"
            )
        if self.bmu_to_paper_wilson_map_matrix != _BMU_TO_PAPER_LR_WILSON_MAP:
            raise ValueError(
                "bmu_to_paper_wilson_map_matrix must match the frozen exact inverse "
                "Wilson-space map"
            )
        if not self.lr_definitions_frozen:
            raise ValueError("lr_definitions_frozen must remain True in the LR-freeze slice")
        if not self.mapping_matrix_frozen:
            raise ValueError(
                "mapping_matrix_frozen must remain True once the exact numeric "
                "paper-to-BMU LR map is frozen"
            )
        if not self.lr_running_activated:
            raise ValueError(
                "lr_running_activated must remain True in the LR-RG-1 slice"
            )

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "contract_id": self.contract_id,
            "paper_operator_basis_id": self.paper_operator_basis_id,
            "paper_operator_order": list(self.paper_operator_order),
            "paper_operator_definition_ids": list(self.paper_operator_definition_ids),
            "bmu_lr_basis_id": self.bmu_lr_basis_id,
            "bmu_lr_operator_order": list(self.bmu_lr_operator_order),
            "map_direction_id": self.map_direction_id,
            "operator_map_formula": self.operator_map_formula,
            "paper_to_bmu_operator_map_matrix": _real_matrix_as_dict(
                self.paper_to_bmu_operator_map_matrix
            ),
            "bmu_to_paper_operator_map_matrix": _real_matrix_as_dict(
                self.bmu_to_paper_operator_map_matrix
            ),
            "wilson_map_formula": self.wilson_map_formula,
            "paper_to_bmu_wilson_map_matrix": _real_matrix_as_dict(
                self.paper_to_bmu_wilson_map_matrix
            ),
            "bmu_to_paper_wilson_map_matrix": _real_matrix_as_dict(
                self.bmu_to_paper_wilson_map_matrix
            ),
            "renormalization_scheme_id": self.renormalization_scheme_id,
            "operator_normalization_id": self.operator_normalization_id,
            "projector_normalization_id": self.projector_normalization_id,
            "fierz_validity_statement_id": self.fierz_validity_statement_id,
            "fierz_validity_note": self.fierz_validity_note,
            "source_id": self.source_id,
            "source_citation": self.source_citation,
            "source_locator": self.source_locator,
            "status_id": self.status_id,
            "lr_definitions_frozen": self.lr_definitions_frozen,
            "mapping_matrix_frozen": self.mapping_matrix_frozen,
            "lr_running_activated": self.lr_running_activated,
            "metadata_compatibility_note": self.metadata_compatibility_note,
            "notes": self.notes,
        }


@dataclass(frozen=True)
class Paper07101869DeltaF2RGContract:
    """Frozen LO RG contract for the paper-owned ``Delta F = 2`` EFT slice."""

    operator_basis_id: str = PAPER_0710_1869_DELTAF2_OPERATOR_BASIS_ID
    renormalization_scheme_id: str = PAPER_0710_1869_DELTAF2_RG_SCHEME_ID
    operator_normalization_id: str = PAPER_0710_1869_DELTAF2_OPERATOR_NORMALIZATION_ID
    operator_order: tuple[str, ...] = PAPER_0710_1869_DELTAF2_OPERATOR_ORDER
    supported_operator_names: tuple[str, ...] = _SUPPORTED_OPERATORS
    unsupported_operator_names: tuple[str, ...] = _UNSUPPORTED_OPERATORS
    matching_input_scheme_id: str = PAPER_0710_1869_DELTAF2_RENORMALIZATION_SCHEME_ID
    matching_input_matching_id: str = PAPER_0710_1869_DELTAF2_MATCHING_ID
    matching_input_wilson_schema_id: str = PAPER_0710_1869_DELTAF2_WILSON_SCHEMA_ID
    rg_order_id: str = PAPER_0710_1869_DELTAF2_RG_ORDER_ID
    alpha_s_policy_id: str = PAPER_0710_1869_DELTAF2_RG_ALPHA_S_POLICY_ID
    threshold_policy_id: str = PAPER_0710_1869_DELTAF2_RG_THRESHOLD_POLICY_ID
    adm_source_id: str = PAPER_0710_1869_DELTAF2_RG_ADM_SOURCE_ID
    lr_basis_contract_id: str = PAPER_0710_1869_DELTAF2_RG_LR_BASIS_CONTRACT_ID
    lr_basis_status_id: str = PAPER_0710_1869_DELTAF2_RG_LR_BASIS_STATUS_ID
    lr_running_bridge_id: str = PAPER_0710_1869_DELTAF2_RG_LR_RUNNING_BRIDGE_ID
    supported_operator_subset_id: str = (
        PAPER_0710_1869_DELTAF2_RG_SUPPORTED_OPERATOR_SUBSET_ID
    )
    evolved_scale_semantics_id: str = PAPER_0710_1869_DELTAF2_RG_SCALE_SEMANTICS_ID
    schema_id: str = PAPER_0710_1869_DELTAF2_RG_CONTRACT_SCHEMA_ID
    mode_id: str = PAPER_0710_1869_MODE_ID
    paper_id: str = PAPER_0710_1869_PAPER_ID
    notes: str = (
        "LO RG uses BMU hep-ph/0005183 for Q1_VLL/Q1_VRR directly and evolves "
        "Q4_LR/Q5_LR through the frozen BMU conjugation bridge. "
        "The explicit Q4_LR/Q5_LR definitions and paper-to-BMU LR map are frozen "
        "under the BMU NDR-MS LR-basis contract. The canonical LR status now covers "
        "RG support plus the custom LR hadronic, custom LR-only observable, and "
        "custom combined Q1+LR observable surfaces, while the default/exported "
        "observable path remains Q1-only."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_DELTAF2_RG_CONTRACT_SCHEMA_ID,
            ),
        )
        require_member("mode_id", self.mode_id, (PAPER_0710_1869_MODE_ID,))
        require_member("paper_id", self.paper_id, (PAPER_0710_1869_PAPER_ID,))
        for field_name in (
            "operator_basis_id",
            "renormalization_scheme_id",
            "operator_normalization_id",
            "matching_input_scheme_id",
            "matching_input_matching_id",
            "matching_input_wilson_schema_id",
            "rg_order_id",
            "alpha_s_policy_id",
            "threshold_policy_id",
            "adm_source_id",
            "lr_basis_contract_id",
            "lr_basis_status_id",
            "lr_running_bridge_id",
            "supported_operator_subset_id",
            "evolved_scale_semantics_id",
            "notes",
        ):
            object.__setattr__(
                self,
                field_name,
                require_nonempty_identifier(field_name, getattr(self, field_name)),
            )
        if tuple(self.operator_order) != PAPER_0710_1869_DELTAF2_OPERATOR_ORDER:
            raise ValueError(
                "operator_order must match the frozen paper Delta F = 2 basis order"
            )
        if tuple(self.supported_operator_names) != _SUPPORTED_OPERATORS:
            raise ValueError("supported_operator_names must match the RG-supported paper basis")
        if tuple(self.unsupported_operator_names) != _UNSUPPORTED_OPERATORS:
            raise ValueError(
                "unsupported_operator_names must match the current RG-unsupported subset"
            )

    @property
    def tags(self) -> dict[str, object]:
        return {
            "mode_id": self.mode_id,
            "paper_id": self.paper_id,
            "operator_basis_id": self.operator_basis_id,
            "renormalization_scheme_id": self.renormalization_scheme_id,
            "rg_scheme_id": self.rg_scheme_id,
            "operator_normalization_id": self.operator_normalization_id,
            "operator_order": list(self.operator_order),
            "matching_input_scheme_id": self.matching_input_scheme_id,
            "matching_input_matching_id": self.matching_input_matching_id,
            "matching_input_wilson_schema_id": self.matching_input_wilson_schema_id,
            "rg_order_id": self.rg_order_id,
            "alpha_s_policy_id": self.alpha_s_policy_id,
            "threshold_policy_id": self.threshold_policy_id,
            "adm_source_id": self.adm_source_id,
            "lr_basis_contract_id": self.lr_basis_contract_id,
            "lr_basis_status_id": self.lr_basis_status_id,
            "lr_running_bridge_id": self.lr_running_bridge_id,
            "supported_operator_subset_id": self.supported_operator_subset_id,
            "evolved_scale_semantics_id": self.evolved_scale_semantics_id,
        }

    @property
    def rg_scheme_id(self) -> str:
        return self.renormalization_scheme_id

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "mode_id": self.mode_id,
            "paper_id": self.paper_id,
            "operator_basis_id": self.operator_basis_id,
            "renormalization_scheme_id": self.renormalization_scheme_id,
            "rg_scheme_id": self.rg_scheme_id,
            "operator_normalization_id": self.operator_normalization_id,
            "operator_order": list(self.operator_order),
            "supported_operator_names": list(self.supported_operator_names),
            "unsupported_operator_names": list(self.unsupported_operator_names),
            "matching_input_scheme_id": self.matching_input_scheme_id,
            "matching_input_matching_id": self.matching_input_matching_id,
            "matching_input_wilson_schema_id": self.matching_input_wilson_schema_id,
            "rg_order_id": self.rg_order_id,
            "alpha_s_policy_id": self.alpha_s_policy_id,
            "threshold_policy_id": self.threshold_policy_id,
            "adm_source_id": self.adm_source_id,
            "lr_basis_contract_id": self.lr_basis_contract_id,
            "lr_basis_status_id": self.lr_basis_status_id,
            "lr_running_bridge_id": self.lr_running_bridge_id,
            "supported_operator_subset_id": self.supported_operator_subset_id,
            "evolved_scale_semantics_id": self.evolved_scale_semantics_id,
            "tags": self.tags,
            "notes": self.notes,
        }


@dataclass(frozen=True)
class Paper07101869DeltaF2LOADM:
    """Frozen LO ADM summary anchored to BMU for the RG-supported paper basis."""

    contract: Paper07101869DeltaF2RGContract
    gamma0_vll: float = _GAMMA0_VLL
    gamma0_vrr: float = _GAMMA0_VRR
    gamma0_bmu_lr_basis: tuple[tuple[float, float], tuple[float, float]] = _GAMMA0_BMU_LR
    paper_basis_lr_map_frozen: bool = True
    schema_id: str = PAPER_0710_1869_DELTAF2_LO_ADM_SCHEMA_ID

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_DELTAF2_LO_ADM_SCHEMA_ID,
            ),
        )
        if not isinstance(self.contract, Paper07101869DeltaF2RGContract):
            raise ValueError("contract must be a Paper07101869DeltaF2RGContract")
        for field_name in ("gamma0_vll", "gamma0_vrr"):
            numeric = float(getattr(self, field_name))
            if not math.isfinite(numeric):
                raise ValueError(f"{field_name} must be finite")
            object.__setattr__(self, field_name, numeric)
        if (
            tuple(tuple(float(value) for value in row) for row in self.gamma0_bmu_lr_basis)
            != _GAMMA0_BMU_LR
        ):
            raise ValueError("gamma0_bmu_lr_basis must match the frozen BMU LO LR block")
        if not self.paper_basis_lr_map_frozen:
            raise ValueError(
                "paper_basis_lr_map_frozen must remain True once the exact LR basis map "
                "is frozen"
            )

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "contract": self.contract.as_dict(),
            "gamma0_vll": self.gamma0_vll,
            "gamma0_vrr": self.gamma0_vrr,
            "gamma0_bmu_lr_basis": [list(row) for row in self.gamma0_bmu_lr_basis],
            "paper_basis_lr_map_frozen": self.paper_basis_lr_map_frozen,
        }


@dataclass(frozen=True)
class Paper07101869DeltaF2RGSegment:
    """One fixed-``n_f`` LO evolution segment."""

    mu_start_GeV: float
    mu_end_GeV: float
    n_f: int
    alpha_s_start: float
    alpha_s_end: float
    vll_factor: float
    vrr_factor: float
    bmu_lr_matrix: tuple[tuple[complex, complex], tuple[complex, complex]]
    schema_id: str = PAPER_0710_1869_DELTAF2_RG_SEGMENT_SCHEMA_ID

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_DELTAF2_RG_SEGMENT_SCHEMA_ID,
            ),
        )
        for field_name in (
            "mu_start_GeV",
            "mu_end_GeV",
            "alpha_s_start",
            "alpha_s_end",
            "vll_factor",
            "vrr_factor",
        ):
            object.__setattr__(
                self,
                field_name,
                require_positive_finite(field_name, getattr(self, field_name)),
            )
        if not isinstance(self.n_f, int) or self.n_f < 3 or self.n_f > 6:
            raise ValueError("n_f must be an integer in the range 3..6")
        if len(self.bmu_lr_matrix) != 2 or any(len(row) != 2 for row in self.bmu_lr_matrix):
            raise ValueError("bmu_lr_matrix must be 2x2")
        normalized_matrix: list[tuple[complex, complex]] = []
        for row in self.bmu_lr_matrix:
            normalized_matrix.append(
                tuple(_require_finite_complex("bmu_lr_entry", value) for value in row)
            )
        object.__setattr__(self, "bmu_lr_matrix", tuple(normalized_matrix))

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "mu_start_GeV": self.mu_start_GeV,
            "mu_end_GeV": self.mu_end_GeV,
            "n_f": self.n_f,
            "alpha_s_start": self.alpha_s_start,
            "alpha_s_end": self.alpha_s_end,
            "vll_factor": self.vll_factor,
            "vrr_factor": self.vrr_factor,
            "bmu_lr_matrix": _matrix_as_dict(self.bmu_lr_matrix),
        }


@dataclass(frozen=True)
class Paper07101869DeltaF2RGEvolution:
    """Deterministic LO evolution object for the paper-owned RG slice."""

    contract: Paper07101869DeltaF2RGContract
    adm: Paper07101869DeltaF2LOADM
    mu_high_GeV: float
    mu_low_GeV: float
    segments: tuple[Paper07101869DeltaF2RGSegment, ...]
    paper_basis_evolution_matrix: tuple[tuple[complex, ...], ...]
    bmu_lr_evolution_matrix: tuple[tuple[complex, complex], tuple[complex, complex]]
    lr_basis_map_supported: bool = True
    schema_id: str = PAPER_0710_1869_DELTAF2_RG_EVOLUTION_SCHEMA_ID

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_DELTAF2_RG_EVOLUTION_SCHEMA_ID,
            ),
        )
        if not isinstance(self.contract, Paper07101869DeltaF2RGContract):
            raise ValueError("contract must be a Paper07101869DeltaF2RGContract")
        if not isinstance(self.adm, Paper07101869DeltaF2LOADM):
            raise ValueError("adm must be a Paper07101869DeltaF2LOADM")
        for field_name in ("mu_high_GeV", "mu_low_GeV"):
            object.__setattr__(
                self,
                field_name,
                require_positive_finite(field_name, getattr(self, field_name)),
            )
        if not all(isinstance(segment, Paper07101869DeltaF2RGSegment) for segment in self.segments):
            raise ValueError("segments must contain only Paper07101869DeltaF2RGSegment entries")
        if len(self.paper_basis_evolution_matrix) != 4 or any(
            len(row) != 4 for row in self.paper_basis_evolution_matrix
        ):
            raise ValueError("paper_basis_evolution_matrix must be 4x4")
        if len(self.bmu_lr_evolution_matrix) != 2 or any(
            len(row) != 2 for row in self.bmu_lr_evolution_matrix
        ):
            raise ValueError("bmu_lr_evolution_matrix must be 2x2")
        object.__setattr__(
            self,
            "paper_basis_evolution_matrix",
            tuple(
                tuple(_require_finite_complex("paper_basis_entry", value) for value in row)
                for row in self.paper_basis_evolution_matrix
            ),
        )
        object.__setattr__(
            self,
            "bmu_lr_evolution_matrix",
            tuple(
                tuple(_require_finite_complex("bmu_lr_entry", value) for value in row)
                for row in self.bmu_lr_evolution_matrix
            ),
        )
        if not self.lr_basis_map_supported:
            raise ValueError(
                "lr_basis_map_supported must remain True once LR running is activated "
                "through the frozen paper-to-BMU bridge"
            )

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "contract": self.contract.as_dict(),
            "adm": self.adm.as_dict(),
            "mu_high_GeV": self.mu_high_GeV,
            "mu_low_GeV": self.mu_low_GeV,
            "segments": [segment.as_dict() for segment in self.segments],
            "paper_basis_evolution_matrix": _matrix_as_dict(self.paper_basis_evolution_matrix),
            "bmu_lr_evolution_matrix": _matrix_as_dict(self.bmu_lr_evolution_matrix),
            "lr_basis_map_supported": self.lr_basis_map_supported,
        }

    def summary(self) -> dict[str, object]:
        return {
            "schema_id": PAPER_0710_1869_DELTAF2_RG_SUMMARY_SCHEMA_ID,
            "mode_id": self.contract.mode_id,
            "paper_id": self.contract.paper_id,
            "mu_high_GeV": self.mu_high_GeV,
            "mu_low_GeV": self.mu_low_GeV,
            "operator_order": list(self.contract.operator_order),
            "supported_operator_names": list(self.contract.supported_operator_names),
            "unsupported_operator_names": list(self.contract.unsupported_operator_names),
            "renormalization_scheme_id": self.contract.renormalization_scheme_id,
            "scheme_id": self.contract.renormalization_scheme_id,
            "rg_scheme_id": self.contract.rg_scheme_id,
            "rg_order_id": self.contract.rg_order_id,
            "alpha_s_policy_id": self.contract.alpha_s_policy_id,
            "threshold_policy_id": self.contract.threshold_policy_id,
            "adm_source_id": self.contract.adm_source_id,
            "lr_basis_status_id": self.contract.lr_basis_status_id,
            "lr_running_bridge_id": self.contract.lr_running_bridge_id,
            "lr_basis_map_supported": self.lr_basis_map_supported,
            "segments": [segment.as_dict() for segment in self.segments],
            "paper_basis_evolution_matrix": _matrix_as_dict(self.paper_basis_evolution_matrix),
            "bmu_lr_evolution_matrix": _matrix_as_dict(self.bmu_lr_evolution_matrix),
        }


@dataclass(frozen=True)
class Paper07101869DeltaF2RGWilsonSnapshot:
    """Self-describing RG-evolved Wilson snapshot tagged with RG provenance."""

    contract: Paper07101869DeltaF2RGContract
    benchmark_id: str
    scale_label: str
    system_id: str
    sector_id: str
    generations: tuple[int, int]
    matching_scale_GeV: float
    source_matching_scale_GeV: float
    propagator_mass_GeV: float
    left_coupling: complex
    right_coupling: complex
    source_renormalization_scheme_id: str
    source_wilson_schema_id: str
    q1_vll: complex
    q1_vrr: complex
    q4_lr: complex
    q5_lr: complex
    schema_id: str = PAPER_0710_1869_DELTAF2_RG_WILSON_SCHEMA_ID

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_DELTAF2_RG_WILSON_SCHEMA_ID,
            ),
        )
        if not isinstance(self.contract, Paper07101869DeltaF2RGContract):
            raise ValueError("contract must be a Paper07101869DeltaF2RGContract")
        for field_name in ("benchmark_id", "scale_label", "system_id"):
            object.__setattr__(
                self,
                field_name,
                require_nonempty_identifier(field_name, getattr(self, field_name)),
            )
        require_member("sector_id", self.sector_id, ("down", "up"))
        if tuple(self.generations) not in ((0, 1), (0, 2), (1, 2), (1, 0), (2, 0), (2, 1)):
            raise ValueError("generations must contain two distinct flavor indices in {0,1,2}")
        object.__setattr__(self, "generations", tuple(int(item) for item in self.generations))
        for field_name in (
            "matching_scale_GeV",
            "source_matching_scale_GeV",
            "propagator_mass_GeV",
        ):
            object.__setattr__(
                self,
                field_name,
                require_positive_finite(field_name, getattr(self, field_name)),
            )
        for field_name in (
            "source_renormalization_scheme_id",
            "source_wilson_schema_id",
        ):
            object.__setattr__(
                self,
                field_name,
                require_nonempty_identifier(field_name, getattr(self, field_name)),
            )
        for field_name in (
            "left_coupling",
            "right_coupling",
            "q1_vll",
            "q1_vrr",
            "q4_lr",
            "q5_lr",
        ):
            object.__setattr__(
                self,
                field_name,
                _require_finite_complex(field_name, getattr(self, field_name)),
            )

    @property
    def coefficients(self) -> dict[str, complex]:
        return {
            PAPER_0710_1869_DELTAF2_Q1_VLL: self.q1_vll,
            PAPER_0710_1869_DELTAF2_Q1_VRR: self.q1_vrr,
            PAPER_0710_1869_DELTAF2_Q4_LR: self.q4_lr,
            PAPER_0710_1869_DELTAF2_Q5_LR: self.q5_lr,
        }

    @property
    def operator_basis_id(self) -> str:
        return self.contract.operator_basis_id

    @property
    def renormalization_scheme_id(self) -> str:
        return self.contract.renormalization_scheme_id

    @property
    def operator_normalization_id(self) -> str:
        return self.contract.operator_normalization_id

    @property
    def tags(self) -> dict[str, object]:
        return {
            **self.contract.tags,
            "benchmark_id": self.benchmark_id,
            "scale_label": self.scale_label,
            "system_id": self.system_id,
            "sector_id": self.sector_id,
            "generations": list(self.generations),
            "matching_scale_GeV": self.matching_scale_GeV,
            "evaluation_scale_GeV": self.matching_scale_GeV,
            "source_matching_scale_GeV": self.source_matching_scale_GeV,
            "propagator_mass_GeV": self.propagator_mass_GeV,
            "source_renormalization_scheme_id": self.source_renormalization_scheme_id,
            "source_wilson_schema_id": self.source_wilson_schema_id,
        }

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "benchmark_id": self.benchmark_id,
            "scale_label": self.scale_label,
            "system_id": self.system_id,
            "sector_id": self.sector_id,
            "generations": list(self.generations),
            "basis_id": self.operator_basis_id,
            "operator_basis_id": self.operator_basis_id,
            "scheme_id": self.renormalization_scheme_id,
            "renormalization_scheme_id": self.renormalization_scheme_id,
            "rg_scheme_id": self.contract.rg_scheme_id,
            "operator_normalization_id": self.operator_normalization_id,
            "rg_order_id": self.contract.rg_order_id,
            "matching_scale_GeV": self.matching_scale_GeV,
            "evaluation_scale_GeV": self.matching_scale_GeV,
            "source_matching_scale_GeV": self.source_matching_scale_GeV,
            "propagator_mass_GeV": self.propagator_mass_GeV,
            "left_coupling": _complex_as_dict(self.left_coupling),
            "right_coupling": _complex_as_dict(self.right_coupling),
            "source_renormalization_scheme_id": self.source_renormalization_scheme_id,
            "source_wilson_schema_id": self.source_wilson_schema_id,
            "coefficients": {
                name: _complex_as_dict(value) for name, value in self.coefficients.items()
            },
            "tags": self.tags,
        }


@dataclass(frozen=True)
class Paper07101869DeltaF2RGResult:
    """Public RG result carrying the evolved Wilsons and full RG provenance."""

    contract: Paper07101869DeltaF2RGContract
    evolution: Paper07101869DeltaF2RGEvolution
    wilsons: Paper07101869DeltaF2RGWilsonSnapshot
    schema_id: str = PAPER_0710_1869_DELTAF2_RG_RESULT_SCHEMA_ID

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_DELTAF2_RG_RESULT_SCHEMA_ID,
            ),
        )
        if not isinstance(self.contract, Paper07101869DeltaF2RGContract):
            raise ValueError("contract must be a Paper07101869DeltaF2RGContract")
        if not isinstance(self.evolution, Paper07101869DeltaF2RGEvolution):
            raise ValueError("evolution must be a Paper07101869DeltaF2RGEvolution")
        if not isinstance(self.wilsons, Paper07101869DeltaF2RGWilsonSnapshot):
            raise ValueError("wilsons must be a Paper07101869DeltaF2RGWilsonSnapshot")
        if self.evolution.contract != self.contract:
            raise ValueError("evolution.contract must match contract")
        if self.wilsons.contract != self.contract:
            raise ValueError("wilsons.contract must match contract")
        if self.wilsons.matching_scale_GeV != self.evolution.mu_low_GeV:
            raise ValueError("wilsons.matching_scale_GeV must match evolution.mu_low_GeV")
        if self.wilsons.source_matching_scale_GeV != self.evolution.mu_high_GeV:
            raise ValueError(
                "wilsons.source_matching_scale_GeV must match evolution.mu_high_GeV"
            )

    @property
    def coefficients(self) -> dict[str, complex]:
        return self.wilsons.coefficients

    @property
    def matching_scale_GeV(self) -> float:
        return self.wilsons.matching_scale_GeV

    @property
    def tags(self) -> dict[str, object]:
        return self.wilsons.tags

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "basis_id": self.contract.operator_basis_id,
            "operator_basis_id": self.contract.operator_basis_id,
            "scheme_id": self.contract.renormalization_scheme_id,
            "renormalization_scheme_id": self.contract.renormalization_scheme_id,
            "rg_scheme_id": self.contract.rg_scheme_id,
            "operator_normalization_id": self.contract.operator_normalization_id,
            "matching_scale_GeV": self.wilsons.matching_scale_GeV,
            "evaluation_scale_GeV": self.wilsons.matching_scale_GeV,
            "mu_low_GeV": self.evolution.mu_low_GeV,
            "mu_high_GeV": self.evolution.mu_high_GeV,
            "source_matching_scale_GeV": self.wilsons.source_matching_scale_GeV,
            "source_renormalization_scheme_id": self.wilsons.source_renormalization_scheme_id,
            "source_wilson_schema_id": self.wilsons.source_wilson_schema_id,
            "rg_order_id": self.contract.rg_order_id,
            "alpha_s_policy_id": self.contract.alpha_s_policy_id,
            "threshold_policy_id": self.contract.threshold_policy_id,
            "adm_source_id": self.contract.adm_source_id,
            "lr_basis_status_id": self.contract.lr_basis_status_id,
            "lr_running_bridge_id": self.contract.lr_running_bridge_id,
            "contract": self.contract.as_dict(),
            "evolution": self.evolution.summary(),
            "wilsons": self.wilsons.as_dict(),
            "coefficients": {
                name: _complex_as_dict(value) for name, value in self.coefficients.items()
            },
            "evolved_coefficients": {
                name: _complex_as_dict(value) for name, value in self.coefficients.items()
            },
            "tags": self.tags,
        }

    def summary(self) -> dict[str, object]:
        return self.as_dict()


def default_paper_0710_1869_deltaf2_rg_contract() -> Paper07101869DeltaF2RGContract:
    """Return the frozen LO RG contract for the paper-owned Delta F = 2 slice."""
    return Paper07101869DeltaF2RGContract()


def default_paper_0710_1869_deltaf2_lr_basis_contract() -> Paper07101869DeltaF2LRBasisContract:
    """Return the frozen LR-map contract for the paper-owned EFT slice."""

    return Paper07101869DeltaF2LRBasisContract()


def paper_lr_to_bmu_lr_operator_map_matrix() -> tuple[
    tuple[complex, complex], tuple[complex, complex]
]:
    """Return the frozen operator-space map from paper ``(Q4_LR,Q5_LR)`` to BMU LR."""

    return _complex_matrix_from_real(_PAPER_TO_BMU_LR_OPERATOR_MAP)


def bmu_lr_to_paper_lr_operator_map_matrix() -> tuple[
    tuple[complex, complex], tuple[complex, complex]
]:
    """Return the frozen inverse operator-space map from BMU LR to paper LR."""

    return _complex_matrix_from_real(_BMU_TO_PAPER_LR_OPERATOR_MAP)


def paper_lr_to_bmu_lr_wilson_map_matrix() -> tuple[
    tuple[complex, complex], tuple[complex, complex]
]:
    """Return the frozen Wilson-space map from paper LR coefficients to BMU LR."""

    return _complex_matrix_from_real(_PAPER_TO_BMU_LR_WILSON_MAP)


def bmu_lr_to_paper_lr_wilson_map_matrix() -> tuple[
    tuple[complex, complex], tuple[complex, complex]
]:
    """Return the frozen inverse Wilson-space map from BMU LR coefficients to paper LR."""

    return _complex_matrix_from_real(_BMU_TO_PAPER_LR_WILSON_MAP)


def map_paper_lr_wilsons_to_bmu_lr(
    paper_lr_wilsons: tuple[complex, complex] | list[complex],
) -> tuple[complex, complex]:
    """Map paper-basis LR Wilsons into the frozen BMU LR Wilson ordering."""

    values = tuple(_require_finite_complex("paper_lr_wilson", value) for value in paper_lr_wilsons)
    if len(values) != 2:
        raise ValueError("paper_lr_wilsons must contain exactly two LR coefficients")
    return _matrix_vector_multiply(paper_lr_to_bmu_lr_wilson_map_matrix(), values)


def map_bmu_lr_wilsons_to_paper_lr(
    bmu_lr_wilsons: tuple[complex, complex] | list[complex],
) -> tuple[complex, complex]:
    """Map BMU-basis LR Wilsons into the frozen paper LR Wilson ordering."""

    values = tuple(_require_finite_complex("bmu_lr_wilson", value) for value in bmu_lr_wilsons)
    if len(values) != 2:
        raise ValueError("bmu_lr_wilsons must contain exactly two LR coefficients")
    return _matrix_vector_multiply(bmu_lr_to_paper_lr_wilson_map_matrix(), values)


def build_paper_0710_1869_deltaf2_lo_adm(
    *,
    contract: Paper07101869DeltaF2RGContract | None = None,
) -> Paper07101869DeltaF2LOADM:
    """Return the frozen BMU-anchored LO ADM bundle for the RG-supported basis."""
    resolved_contract = (
        default_paper_0710_1869_deltaf2_rg_contract() if contract is None else contract
    )
    if not isinstance(resolved_contract, Paper07101869DeltaF2RGContract):
        raise ValueError("contract must be a Paper07101869DeltaF2RGContract")
    return Paper07101869DeltaF2LOADM(contract=resolved_contract)


def _require_wilson_contract_compatibility(
    wilson_contract: Paper07101869DeltaF2WilsonContract,
    rg_contract: Paper07101869DeltaF2RGContract,
) -> None:
    if wilson_contract.operator_basis_id != rg_contract.operator_basis_id:
        raise ValueError("wilson contract operator_basis_id is incompatible with the RG contract")
    if wilson_contract.renormalization_scheme_id != rg_contract.matching_input_scheme_id:
        raise ValueError(
            "wilson contract renormalization_scheme_id is incompatible with the RG contract"
        )
    if wilson_contract.matching_id != rg_contract.matching_input_matching_id:
        raise ValueError(
            "wilson contract matching_id is incompatible with the RG contract; "
            "LO RG requires the frozen PR3 input matching_id "
            f"{rg_contract.matching_input_matching_id!r}"
        )
    if wilson_contract.operator_normalization_id != rg_contract.operator_normalization_id:
        raise ValueError(
            "wilson contract operator_normalization_id is incompatible with the RG contract"
        )
    if tuple(wilson_contract.operator_order) != tuple(rg_contract.operator_order):
        raise ValueError("wilson contract operator order is incompatible with the RG contract")


def _vector_from_wilsons(
    wilsons: Paper07101869DeltaF2WilsonCoefficients | Paper07101869DeltaF2RGWilsonSnapshot,
) -> tuple[complex, complex, complex, complex]:
    return (
        _require_finite_complex("q1_vll", wilsons.q1_vll),
        _require_finite_complex("q1_vrr", wilsons.q1_vrr),
        _require_finite_complex("q4_lr", wilsons.q4_lr),
        _require_finite_complex("q5_lr", wilsons.q5_lr),
    )


def _paper_lr_wilson_evolution_matrix(
    bmu_lr_matrix: tuple[tuple[complex, complex], tuple[complex, complex]],
) -> tuple[tuple[complex, complex], tuple[complex, complex]]:
    return _matrix_multiply(
        bmu_lr_to_paper_lr_wilson_map_matrix(),
        _matrix_multiply(bmu_lr_matrix, paper_lr_to_bmu_lr_wilson_map_matrix()),
    )


def _paper_basis_segment_matrix(
    *,
    vll_factor: float,
    vrr_factor: float,
    bmu_lr_matrix: tuple[tuple[complex, complex], tuple[complex, complex]],
) -> tuple[tuple[complex, ...], ...]:
    paper_lr_matrix = _paper_lr_wilson_evolution_matrix(bmu_lr_matrix)
    return (
        (
            _require_finite_complex("u00", complex(vll_factor, 0.0)),
            0.0 + 0.0j,
            0.0 + 0.0j,
            0.0 + 0.0j,
        ),
        (
            0.0 + 0.0j,
            _require_finite_complex("u11", complex(vrr_factor, 0.0)),
            0.0 + 0.0j,
            0.0 + 0.0j,
        ),
        (
            0.0 + 0.0j,
            0.0 + 0.0j,
            paper_lr_matrix[0][0],
            paper_lr_matrix[0][1],
        ),
        (
            0.0 + 0.0j,
            0.0 + 0.0j,
            paper_lr_matrix[1][0],
            paper_lr_matrix[1][1],
        ),
    )


def compute_deltaf2_lo_evolution_matrix(
    mu_high_GeV: float,
    mu_low_GeV: float,
    *,
    contract: Paper07101869DeltaF2RGContract | None = None,
    adm: Paper07101869DeltaF2LOADM | None = None,
    thresholds: tuple[tuple[float, int, int], ...] | list[tuple[float, int, int]] | None = None,
) -> Paper07101869DeltaF2RGEvolution:
    """Compute the LO QCD evolution matrix between two scales.

    The returned 4x4 matrix evolves the full paper operator basis. The LR block
    is obtained by conjugating the BMU LO LR evolution matrix with the frozen
    Wilson-space map between the paper ``(Q4_LR,Q5_LR)`` ordering and the BMU
    ``(Q1^LR,Q2^LR)`` ordering.
    """

    resolved_contract = (
        default_paper_0710_1869_deltaf2_rg_contract() if contract is None else contract
    )
    if not isinstance(resolved_contract, Paper07101869DeltaF2RGContract):
        raise ValueError("contract must be a Paper07101869DeltaF2RGContract")
    resolved_adm = (
        build_paper_0710_1869_deltaf2_lo_adm(contract=resolved_contract) if adm is None else adm
    )
    if not isinstance(resolved_adm, Paper07101869DeltaF2LOADM):
        raise ValueError("adm must be a Paper07101869DeltaF2LOADM")
    if resolved_adm.contract != resolved_contract:
        raise ValueError("adm.contract must match the RG contract")
    high_scale = require_positive_finite("mu_high_GeV", mu_high_GeV)
    low_scale = require_positive_finite("mu_low_GeV", mu_low_GeV)
    resolved_thresholds = _normalize_thresholds(thresholds)

    if math.isclose(high_scale, low_scale, rel_tol=0.0, abs_tol=0.0):
        identity4 = _identity_matrix(4)
        identity2 = (
            (1.0 + 0.0j, 0.0 + 0.0j),
            (0.0 + 0.0j, 1.0 + 0.0j),
        )
        return Paper07101869DeltaF2RGEvolution(
            contract=resolved_contract,
            adm=resolved_adm,
            mu_high_GeV=high_scale,
            mu_low_GeV=low_scale,
            segments=tuple(),
            paper_basis_evolution_matrix=identity4,
            bmu_lr_evolution_matrix=identity2,
        )

    boundaries = _build_segment_boundaries(
        mu_high_GeV=high_scale,
        mu_low_GeV=low_scale,
        thresholds=resolved_thresholds,
    )
    paper_total = _identity_matrix(4)
    bmu_lr_total = (
        (1.0 + 0.0j, 0.0 + 0.0j),
        (0.0 + 0.0j, 1.0 + 0.0j),
    )
    segments: list[Paper07101869DeltaF2RGSegment] = []

    for start_scale, end_scale in zip(boundaries[:-1], boundaries[1:], strict=True):
        midpoint = math.sqrt(start_scale * end_scale)
        n_f = _n_f_at_scale(midpoint, resolved_thresholds)
        alpha_start = float(
            alpha_s(
                start_scale,
                n_loops=1,
                matching_loops=0,
                thresholds=list(resolved_thresholds),
            )
        )
        alpha_end = float(
            alpha_s(
                end_scale,
                n_loops=1,
                matching_loops=0,
                thresholds=list(resolved_thresholds),
            )
        )
        exponent_vll = resolved_adm.gamma0_vll / (2.0 * float(beta_0(n_f)))
        exponent_vrr = resolved_adm.gamma0_vrr / (2.0 * float(beta_0(n_f)))
        eta = alpha_start / alpha_end
        vll_factor = float(eta**exponent_vll)
        vrr_factor = float(eta**exponent_vrr)
        bmu_lr_segment = _closed_form_bmu_lr_segment_matrix(
            alpha_start=alpha_start,
            alpha_end=alpha_end,
            n_f=n_f,
        )
        segment = Paper07101869DeltaF2RGSegment(
            mu_start_GeV=start_scale,
            mu_end_GeV=end_scale,
            n_f=n_f,
            alpha_s_start=alpha_start,
            alpha_s_end=alpha_end,
            vll_factor=vll_factor,
            vrr_factor=vrr_factor,
            bmu_lr_matrix=bmu_lr_segment,
        )
        segments.append(segment)
        paper_total = _matrix_multiply(
            _paper_basis_segment_matrix(
                vll_factor=vll_factor,
                vrr_factor=vrr_factor,
                bmu_lr_matrix=bmu_lr_segment,
            ),
            paper_total,
        )
        bmu_lr_total = _matrix_multiply(bmu_lr_segment, bmu_lr_total)

    return Paper07101869DeltaF2RGEvolution(
        contract=resolved_contract,
        adm=resolved_adm,
        mu_high_GeV=high_scale,
        mu_low_GeV=low_scale,
        segments=tuple(segments),
        paper_basis_evolution_matrix=paper_total,
        bmu_lr_evolution_matrix=bmu_lr_total,
    )


def _require_supported_wilson_subset(
    wilsons: Paper07101869DeltaF2WilsonCoefficients | Paper07101869DeltaF2RGWilsonSnapshot,
    rg_contract: Paper07101869DeltaF2RGContract,
) -> None:
    unsupported_values = {
        operator_name: coefficient
        for operator_name, coefficient in (
            (PAPER_0710_1869_DELTAF2_Q1_VLL, complex(wilsons.q1_vll)),
            (PAPER_0710_1869_DELTAF2_Q1_VRR, complex(wilsons.q1_vrr)),
            (PAPER_0710_1869_DELTAF2_Q4_LR, complex(wilsons.q4_lr)),
            (PAPER_0710_1869_DELTAF2_Q5_LR, complex(wilsons.q5_lr)),
        )
        if operator_name in rg_contract.unsupported_operator_names
    }
    if not unsupported_values:
        return
    nonzero = [
        name
        for name, value in unsupported_values.items()
        if abs(value.real) > _ZERO_TOLERANCE or abs(value.imag) > _ZERO_TOLERANCE
    ]
    if nonzero:
        joined = ", ".join(nonzero)
        raise NotImplementedError(
            "LO RG currently supports only the "
            f"{rg_contract.supported_operator_names!r} subset. "
            f"Non-zero unsupported operators {joined} remain blocked by LR contract "
            f"{rg_contract.lr_basis_contract_id!r} under status "
            f"{rg_contract.lr_basis_status_id!r}."
        )


def _coerce_rg_input_wilsons(
    wilsons: object,
) -> Paper07101869DeltaF2WilsonCoefficients | Paper07101869DeltaF2RGWilsonSnapshot:
    if isinstance(wilsons, Paper07101869DeltaF2RGResult):
        return wilsons.wilsons
    if isinstance(
        wilsons,
        (Paper07101869DeltaF2WilsonCoefficients, Paper07101869DeltaF2RGWilsonSnapshot),
    ):
        return wilsons
    raise ValueError(
        "wilsons must be a Paper07101869DeltaF2WilsonCoefficients, "
        "Paper07101869DeltaF2RGWilsonSnapshot, or Paper07101869DeltaF2RGResult"
    )


def _require_rg_input_compatibility(
    wilsons: Paper07101869DeltaF2WilsonCoefficients | Paper07101869DeltaF2RGWilsonSnapshot,
    rg_contract: Paper07101869DeltaF2RGContract,
) -> None:
    if isinstance(wilsons, Paper07101869DeltaF2WilsonCoefficients):
        _require_wilson_contract_compatibility(wilsons.contract, rg_contract)
        return
    if wilsons.contract != rg_contract:
        raise ValueError("RG-tagged wilsons must carry the same RG contract as the run request")
    if wilsons.contract.renormalization_scheme_id != rg_contract.renormalization_scheme_id:
        raise ValueError("RG-tagged wilsons carry an incompatible renormalization scheme")
    if tuple(wilsons.contract.operator_order) != tuple(rg_contract.operator_order):
        raise ValueError("RG-tagged wilsons carry an incompatible operator order")


def _source_renormalization_scheme_id(
    wilsons: Paper07101869DeltaF2WilsonCoefficients | Paper07101869DeltaF2RGWilsonSnapshot,
) -> str:
    if isinstance(wilsons, Paper07101869DeltaF2WilsonCoefficients):
        return wilsons.contract.renormalization_scheme_id
    return wilsons.contract.renormalization_scheme_id


def _source_wilson_schema_id(
    wilsons: Paper07101869DeltaF2WilsonCoefficients | Paper07101869DeltaF2RGWilsonSnapshot,
) -> str:
    return wilsons.schema_id


def evolve_deltaf2_wilsons_lo(
    wilsons: (
        Paper07101869DeltaF2WilsonCoefficients
        | Paper07101869DeltaF2RGWilsonSnapshot
        | Paper07101869DeltaF2RGResult
    ),
    *,
    mu_low_GeV: float = PAPER_0710_1869_DELTAF2_DEFAULT_MU_HAD_GEV,
    contract: Paper07101869DeltaF2RGContract | None = None,
    adm: Paper07101869DeltaF2LOADM | None = None,
    thresholds: tuple[tuple[float, int, int], ...] | list[tuple[float, int, int]] | None = None,
) -> Paper07101869DeltaF2RGResult:
    """Run paper-owned Wilson coefficients to a lower scale at LO.

    The evolved Wilson object reuses ``matching_scale_GeV`` to carry the current
    evaluation scale. This is frozen explicitly in the RG contract via
    ``evolved_scale_semantics_id`` until the paper artifact schema grows a
    dedicated evaluation-scale field.
    """

    resolved_input = _coerce_rg_input_wilsons(wilsons)
    resolved_contract = (
        default_paper_0710_1869_deltaf2_rg_contract() if contract is None else contract
    )
    if not isinstance(resolved_contract, Paper07101869DeltaF2RGContract):
        raise ValueError("contract must be a Paper07101869DeltaF2RGContract")
    _require_rg_input_compatibility(resolved_input, resolved_contract)
    _require_supported_wilson_subset(resolved_input, resolved_contract)
    evolution = compute_deltaf2_lo_evolution_matrix(
        resolved_input.matching_scale_GeV,
        mu_low_GeV,
        contract=resolved_contract,
        adm=adm,
        thresholds=thresholds,
    )
    coefficients = _vector_from_wilsons(resolved_input)
    evolved_coefficients = _matrix_vector_multiply(
        evolution.paper_basis_evolution_matrix,
        coefficients,
    )
    evolved_wilsons = Paper07101869DeltaF2RGWilsonSnapshot(
        contract=resolved_contract,
        benchmark_id=resolved_input.benchmark_id,
        scale_label=resolved_input.scale_label,
        system_id=resolved_input.system_id,
        sector_id=resolved_input.sector_id,
        generations=resolved_input.generations,
        matching_scale_GeV=evolution.mu_low_GeV,
        source_matching_scale_GeV=evolution.mu_high_GeV,
        propagator_mass_GeV=resolved_input.propagator_mass_GeV,
        left_coupling=resolved_input.left_coupling,
        right_coupling=resolved_input.right_coupling,
        source_renormalization_scheme_id=_source_renormalization_scheme_id(resolved_input),
        source_wilson_schema_id=_source_wilson_schema_id(resolved_input),
        q1_vll=evolved_coefficients[0],
        q1_vrr=evolved_coefficients[1],
        q4_lr=evolved_coefficients[2],
        q5_lr=evolved_coefficients[3],
    )
    return Paper07101869DeltaF2RGResult(
        contract=resolved_contract,
        evolution=evolution,
        wilsons=evolved_wilsons,
    )


def run_deltaf2_wilsons_lo(
    wilsons: (
        Paper07101869DeltaF2WilsonCoefficients
        | Paper07101869DeltaF2RGWilsonSnapshot
        | Paper07101869DeltaF2RGResult
    ),
    *,
    mu_low_GeV: float = PAPER_0710_1869_DELTAF2_DEFAULT_MU_HAD_GEV,
    contract: Paper07101869DeltaF2RGContract | None = None,
    adm: Paper07101869DeltaF2LOADM | None = None,
    thresholds: tuple[tuple[float, int, int], ...] | list[tuple[float, int, int]] | None = None,
) -> Paper07101869DeltaF2RGResult:
    """Compatibility alias for LO paper-owned Wilson running."""
    return evolve_deltaf2_wilsons_lo(
        wilsons,
        mu_low_GeV=mu_low_GeV,
        contract=contract,
        adm=adm,
        thresholds=thresholds,
    )


def build_paper_0710_1869_default_kaon_rg_summary(
    *,
    mu_low_GeV: float = PAPER_0710_1869_DELTAF2_DEFAULT_MU_HAD_GEV,
    contract: Paper07101869DeltaF2RGContract | None = None,
    adm: Paper07101869DeltaF2LOADM | None = None,
    thresholds: tuple[tuple[float, int, int], ...] | list[tuple[float, int, int]] | None = None,
) -> dict[str, object]:
    """Return a deterministic summary for the default kaon LO RG path."""
    default_match = default_paper_0710_1869_kaon_matching()
    result = evolve_deltaf2_wilsons_lo(
        default_match.wilsons,
        mu_low_GeV=mu_low_GeV,
        contract=contract,
        adm=adm,
        thresholds=thresholds,
    )
    payload = result.summary()
    payload["input_matching_summary"] = default_paper_0710_1869_kaon_matching_summary()
    payload["mu_match_GeV"] = default_match.wilsons.matching_scale_GeV
    payload["mu_low_GeV"] = mu_low_GeV
    return payload


def default_paper_0710_1869_default_kaon_rg_summary() -> dict[str, object]:
    """Return the default deterministic kaon LO RG summary."""
    return build_paper_0710_1869_default_kaon_rg_summary()


__all__ = [
    "PAPER_0710_1869_DELTAF2_LO_ADM_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_RG_CONTRACT_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_RG_EVOLUTION_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_RG_RESULT_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_RG_SEGMENT_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_RG_SUMMARY_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_RG_WILSON_SCHEMA_ID",
    "Paper07101869DeltaF2LRBasisContract",
    "Paper07101869DeltaF2LOADM",
    "Paper07101869DeltaF2RGContract",
    "Paper07101869DeltaF2RGEvolution",
    "Paper07101869DeltaF2RGResult",
    "Paper07101869DeltaF2RGSegment",
    "Paper07101869DeltaF2RGWilsonSnapshot",
    "bmu_lr_to_paper_lr_operator_map_matrix",
    "bmu_lr_to_paper_lr_wilson_map_matrix",
    "build_paper_0710_1869_default_kaon_rg_summary",
    "build_paper_0710_1869_deltaf2_lo_adm",
    "compute_deltaf2_lo_evolution_matrix",
    "default_paper_0710_1869_default_kaon_rg_summary",
    "default_paper_0710_1869_deltaf2_lr_basis_contract",
    "default_paper_0710_1869_deltaf2_rg_contract",
    "evolve_deltaf2_wilsons_lo",
    "map_bmu_lr_wilsons_to_paper_lr",
    "map_paper_lr_wilsons_to_bmu_lr",
    "paper_lr_to_bmu_lr_operator_map_matrix",
    "paper_lr_to_bmu_lr_wilson_map_matrix",
    "run_deltaf2_wilsons_lo",
]
