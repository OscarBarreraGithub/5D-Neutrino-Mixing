"""Tree-level KK-gluon matching for the paper-owned ``Delta F = 2`` slice."""

from __future__ import annotations

import math
from dataclasses import dataclass

from ..benchmarks import Paper07101869Benchmark, default_paper_0710_1869_pr1_benchmark
from ..conventions import PAPER_0710_1869_MODE_ID, PAPER_0710_1869_PAPER_ID
from ..couplings import (
    Paper07101869CouplingContract,
    default_paper_0710_1869_coupling_contract,
)
from ..kkgluon import (
    Paper07101869KKGluonCouplings,
    Paper07101869KKGluonFlavorMatrix,
    build_paper_0710_1869_kk_gluon_couplings,
)
from ..scales import Paper07101869ScalePoint, default_paper_0710_1869_scales
from ..validation import (
    require_known_schema_id,
    require_member,
    require_nonempty_identifier,
    require_positive_finite,
)
from .operators import (
    PAPER_0710_1869_DELTAF2_Q1_VLL,
    PAPER_0710_1869_DELTAF2_Q1_VRR,
    PAPER_0710_1869_DELTAF2_Q4_LR,
    PAPER_0710_1869_DELTAF2_Q5_LR,
    PAPER_0710_1869_DELTAF2_WILSON_CONTRACT_SCHEMA_ID,
    Paper07101869DeltaF2WilsonCoefficients,
    Paper07101869DeltaF2WilsonContract,
    default_paper_0710_1869_deltaf2_wilson_contract,
)
from .rg_inputs import (
    PAPER_0710_1869_DELTAF2_RG_LR_BASIS_CONTRACT_ID,
    PAPER_0710_1869_DELTAF2_RG_LR_BASIS_STATUS_ID,
)

PAPER_0710_1869_DELTAF2_SYSTEM_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.system.v1"
)
PAPER_0710_1869_DELTAF2_KKGLUON_MATCH_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.matching_kkgluon.v1"
)
PAPER_0710_1869_DELTAF2_MATCHING_SUMMARY_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.matching_summary.v1"
)
PAPER_0710_1869_DELTAF2_MATCHING_OBSERVABLE_SUPPORT_STATUS_ID = (
    "paper_0710_1869.deltaf2.matching.lr_coefficients_exported.observable_surface_lr_capable.default_rh_down_alignment_residual.v5"
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
    return {
        "real": float(complex_value.real),
        "imag": float(complex_value.imag),
    }


def _require_flavor_pair(
    flavor_indices: tuple[int, int] | list[int],
) -> tuple[int, int]:
    values = tuple(flavor_indices)
    if len(values) != 2:
        raise ValueError("flavor_indices must contain exactly two indices")
    if any(not isinstance(item, int) or isinstance(item, bool) for item in values):
        raise ValueError("flavor_indices must contain integers")
    i, j = values
    if i == j or i not in (0, 1, 2) or j not in (0, 1, 2):
        raise ValueError("flavor_indices must be distinct indices in {0, 1, 2}")
    return (i, j)


def _require_matching_contract_compatibility(
    contract: Paper07101869DeltaF2WilsonContract,
    kk_gluon_couplings: Paper07101869KKGluonCouplings,
) -> None:
    if (
        contract.kk_gluon_normalization_id
        != kk_gluon_couplings.normalization.kk_gluon_normalization_id
    ):
        raise ValueError("Wilson contract kk_gluon_normalization_id is incompatible")
    if (
        contract.dimensionless_matrix_policy_id
        != kk_gluon_couplings.contract.dimensionless_matrix_policy_id
    ):
        raise ValueError("Wilson contract dimensionless_matrix_policy_id is incompatible")
    if (
        contract.universal_subtraction_policy_id
        != kk_gluon_couplings.contract.universal_subtraction_policy_id
    ):
        raise ValueError("Wilson contract universal_subtraction_policy_id is incompatible")
    if contract.propagator_mass_rule_id != kk_gluon_couplings.contract.propagator_mass_rule_id:
        raise ValueError("Wilson contract propagator_mass_rule_id is incompatible")


def _select_chiral_matrices(
    kk_gluon_couplings: Paper07101869KKGluonCouplings,
    system: "Paper07101869DeltaF2System",
) -> tuple[Paper07101869KKGluonFlavorMatrix, Paper07101869KKGluonFlavorMatrix]:
    if system.sector_id == "down":
        return kk_gluon_couplings.left_down_aligned, kk_gluon_couplings.right_down
    return kk_gluon_couplings.left_up_aligned, kk_gluon_couplings.right_up


@dataclass(frozen=True)
class Paper07101869DeltaF2System:
    """One neutral-meson system definition for the paper-owned matching layer."""

    system_id: str = "kaon"
    sector_id: str = "down"
    flavor_indices: tuple[int, int] = (0, 1)
    label: str = "K0-K0bar"
    schema_id: str = PAPER_0710_1869_DELTAF2_SYSTEM_SCHEMA_ID
    mode_id: str = PAPER_0710_1869_MODE_ID
    paper_id: str = PAPER_0710_1869_PAPER_ID
    notes: str = (
        "Default PR3 benchmark system: s-d mixing with left-down-aligned and "
        "right-diagonal KK-gluon matrices."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_DELTAF2_SYSTEM_SCHEMA_ID,
            ),
        )
        require_member("mode_id", self.mode_id, (PAPER_0710_1869_MODE_ID,))
        require_member("paper_id", self.paper_id, (PAPER_0710_1869_PAPER_ID,))
        require_member("sector_id", self.sector_id, ("down", "up"))
        object.__setattr__(
            self,
            "system_id",
            require_nonempty_identifier("system_id", self.system_id),
        )
        object.__setattr__(self, "label", require_nonempty_identifier("label", self.label))
        object.__setattr__(self, "notes", require_nonempty_identifier("notes", self.notes))
        object.__setattr__(
            self,
            "flavor_indices",
            _require_flavor_pair(self.flavor_indices),
        )

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "mode_id": self.mode_id,
            "paper_id": self.paper_id,
            "system_id": self.system_id,
            "sector_id": self.sector_id,
            "flavor_indices": list(self.flavor_indices),
            "label": self.label,
            "notes": self.notes,
        }


@dataclass(frozen=True)
class Paper07101869KKGluonDeltaF2Match:
    """Tree-level KK-gluon matching result for one paper-owned benchmark system."""

    system: Paper07101869DeltaF2System
    benchmark_id: str
    benchmark_label: str
    scale_label: str
    left_matrix_label: str
    left_basis_id: str
    right_matrix_label: str
    right_basis_id: str
    left_coupling: complex
    right_coupling: complex
    prefactor_GeV_minus2: float
    matching_scale_GeV: float
    propagator_mass_GeV: float
    wilsons: Paper07101869DeltaF2WilsonCoefficients
    supported_observable_operator_names: tuple[str, ...] = (
        PAPER_0710_1869_DELTAF2_Q1_VLL,
        PAPER_0710_1869_DELTAF2_Q1_VRR,
        PAPER_0710_1869_DELTAF2_Q4_LR,
        PAPER_0710_1869_DELTAF2_Q5_LR,
    )
    unsupported_observable_operator_names: tuple[str, ...] = ()
    lr_observable_support_contract_id: str = PAPER_0710_1869_DELTAF2_RG_LR_BASIS_CONTRACT_ID
    lr_observable_support_status_id: str = (
        PAPER_0710_1869_DELTAF2_MATCHING_OBSERVABLE_SUPPORT_STATUS_ID
    )
    lr_observable_support_note: str = (
        "Q4_LR and Q5_LR are matched and exported in the paper O4/O5 scalar LR basis "
        "used in the DeltaF=2 literature. The canonical LR status now means LO LR RG "
        "is available through the frozen BMU NDR-MS bridge, custom LR hadronic "
        "inputs are active under exact scheme/scale alignment, and separate custom "
        "LR-only plus custom combined Q1+LR observable surfaces are available. "
        "RESIDUAL(C-2): default RH-down alignment model choice pending paper 0710.1869."
    )
    contract_schema_id: str = PAPER_0710_1869_DELTAF2_WILSON_CONTRACT_SCHEMA_ID
    schema_id: str = PAPER_0710_1869_DELTAF2_KKGLUON_MATCH_SCHEMA_ID

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_DELTAF2_KKGLUON_MATCH_SCHEMA_ID,
            ),
        )
        object.__setattr__(
            self,
            "contract_schema_id",
            require_known_schema_id(
                "contract_schema_id",
                self.contract_schema_id,
                expected=PAPER_0710_1869_DELTAF2_WILSON_CONTRACT_SCHEMA_ID,
            ),
        )
        if not isinstance(self.system, Paper07101869DeltaF2System):
            raise ValueError("system must be a Paper07101869DeltaF2System")
        if not isinstance(self.wilsons, Paper07101869DeltaF2WilsonCoefficients):
            raise ValueError("wilsons must be a Paper07101869DeltaF2WilsonCoefficients")
        for field_name in (
            "benchmark_id",
            "benchmark_label",
            "scale_label",
            "left_matrix_label",
            "left_basis_id",
            "right_matrix_label",
            "right_basis_id",
            "lr_observable_support_contract_id",
            "lr_observable_support_status_id",
            "lr_observable_support_note",
        ):
            object.__setattr__(
                self,
                field_name,
                require_nonempty_identifier(field_name, getattr(self, field_name)),
            )
        for field_name in ("left_coupling", "right_coupling"):
            object.__setattr__(
                self,
                field_name,
                _require_finite_complex(field_name, getattr(self, field_name)),
            )
        for field_name in (
            "prefactor_GeV_minus2",
            "matching_scale_GeV",
            "propagator_mass_GeV",
        ):
            object.__setattr__(
                self,
                field_name,
                require_positive_finite(field_name, getattr(self, field_name)),
            )
        if self.wilsons.benchmark_id != self.benchmark_id:
            raise ValueError("wilsons.benchmark_id must match benchmark_id")
        if self.wilsons.scale_label != self.scale_label:
            raise ValueError("wilsons.scale_label must match scale_label")
        if self.wilsons.system_id != self.system.system_id:
            raise ValueError("wilsons.system_id must match system.system_id")
        if self.wilsons.matching_scale_GeV != self.matching_scale_GeV:
            raise ValueError("wilsons.matching_scale_GeV must match matching_scale_GeV")
        if self.wilsons.propagator_mass_GeV != self.propagator_mass_GeV:
            raise ValueError("wilsons.propagator_mass_GeV must match propagator_mass_GeV")
        if self.wilsons.left_coupling != self.left_coupling:
            raise ValueError("wilsons.left_coupling must match left_coupling")
        if self.wilsons.right_coupling != self.right_coupling:
            raise ValueError("wilsons.right_coupling must match right_coupling")
        if tuple(self.supported_observable_operator_names) != (
            PAPER_0710_1869_DELTAF2_Q1_VLL,
            PAPER_0710_1869_DELTAF2_Q1_VRR,
            PAPER_0710_1869_DELTAF2_Q4_LR,
            PAPER_0710_1869_DELTAF2_Q5_LR,
        ):
            raise ValueError(
                "supported_observable_operator_names must include the Q1 and LR observable subset"
            )
        if tuple(self.unsupported_observable_operator_names) != ():
            raise ValueError("unsupported_observable_operator_names must be empty")

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "contract_schema_id": self.contract_schema_id,
            "system": self.system.as_dict(),
            "benchmark_id": self.benchmark_id,
            "benchmark_label": self.benchmark_label,
            "scale_label": self.scale_label,
            "left_matrix_label": self.left_matrix_label,
            "left_basis_id": self.left_basis_id,
            "right_matrix_label": self.right_matrix_label,
            "right_basis_id": self.right_basis_id,
            "left_coupling": _complex_as_dict(self.left_coupling),
            "right_coupling": _complex_as_dict(self.right_coupling),
            "prefactor_GeV_minus2": self.prefactor_GeV_minus2,
            "matching_scale_GeV": self.matching_scale_GeV,
            "propagator_mass_GeV": self.propagator_mass_GeV,
            "supported_observable_operator_names": list(self.supported_observable_operator_names),
            "unsupported_observable_operator_names": list(
                self.unsupported_observable_operator_names
            ),
            "lr_observable_support_contract_id": self.lr_observable_support_contract_id,
            "lr_observable_support_status_id": self.lr_observable_support_status_id,
            "lr_observable_support_note": self.lr_observable_support_note,
            "wilsons": self.wilsons.as_dict(),
        }

    def summary(self) -> dict[str, object]:
        """Return a flattened deterministic summary for acceptance/export discovery."""
        contract = self.wilsons.contract
        return {
            "schema_id": PAPER_0710_1869_DELTAF2_MATCHING_SUMMARY_SCHEMA_ID,
            "match_schema_id": self.schema_id,
            "wilson_schema_id": self.wilsons.schema_id,
            "mode_id": contract.mode_id,
            "paper_id": contract.paper_id,
            "benchmark_id": self.benchmark_id,
            "benchmark_label": self.benchmark_label,
            "system_id": self.system.system_id,
            "sector_id": self.system.sector_id,
            "generations": list(self.system.flavor_indices),
            "scale_label": self.scale_label,
            "basis_id": contract.operator_basis_id,
            "operator_basis_id": contract.operator_basis_id,
            "scheme_id": contract.renormalization_scheme_id,
            "renormalization_scheme_id": contract.renormalization_scheme_id,
            "operator_normalization_id": contract.operator_normalization_id,
            "matching_id": contract.matching_id,
            "kk_gluon_normalization_id": contract.kk_gluon_normalization_id,
            "dimensionless_matrix_policy_id": contract.dimensionless_matrix_policy_id,
            "universal_subtraction_policy_id": contract.universal_subtraction_policy_id,
            "propagator_mass_rule_id": contract.propagator_mass_rule_id,
            "matching_scale_name": contract.matching_scale_name,
            "propagator_mass_name": contract.propagator_mass_name,
            "mu_match_GeV": self.matching_scale_GeV,
            "matching_scale_GeV": self.matching_scale_GeV,
            "propagator_mass_GeV": self.propagator_mass_GeV,
            "prefactor_GeV_minus2": self.prefactor_GeV_minus2,
            "left_matrix_label": self.left_matrix_label,
            "left_basis_id": self.left_basis_id,
            "right_matrix_label": self.right_matrix_label,
            "right_basis_id": self.right_basis_id,
            "left_coupling": _complex_as_dict(self.left_coupling),
            "right_coupling": _complex_as_dict(self.right_coupling),
            "coefficients": {
                name: _complex_as_dict(value)
                for name, value in self.wilsons.coefficients.items()
            },
            "supported_observable_operator_names": list(self.supported_observable_operator_names),
            "unsupported_observable_operator_names": list(
                self.unsupported_observable_operator_names
            ),
            "lr_observable_support_contract_id": self.lr_observable_support_contract_id,
            "lr_observable_support_status_id": self.lr_observable_support_status_id,
            "lr_basis_status_id": PAPER_0710_1869_DELTAF2_RG_LR_BASIS_STATUS_ID,
            "lr_observable_support_note": self.lr_observable_support_note,
            "tags": self.wilsons.tags,
        }


def _resolve_alias_kk_gluon_couplings(
    *,
    kk_couplings: Paper07101869KKGluonCouplings | None = None,
    kk_gluon_couplings: Paper07101869KKGluonCouplings | None = None,
    couplings: Paper07101869KKGluonCouplings | None = None,
) -> Paper07101869KKGluonCouplings | None:
    provided = [
        value for value in (kk_couplings, kk_gluon_couplings, couplings) if value is not None
    ]
    if len(provided) > 1:
        raise ValueError(
            "pass at most one of kk_couplings, kk_gluon_couplings, or couplings"
        )
    return provided[0] if provided else None


def default_paper_0710_1869_kaon_deltaf2_system() -> Paper07101869DeltaF2System:
    """Return the default kaon-system benchmark used for PR3 verification."""
    return Paper07101869DeltaF2System()


def default_paper_0710_1869_bd_deltaf2_system() -> Paper07101869DeltaF2System:
    """Return the custom B_d-system definition for paper-mode DeltaF=2 matching."""

    return Paper07101869DeltaF2System(
        system_id="B_d",
        sector_id="down",
        flavor_indices=(0, 2),
        label="Bd-Bdbar",
        notes=(
            "Custom B_d benchmark system: b-d mixing with left-down-aligned and "
            "right-diagonal KK-gluon matrices."
        ),
    )


def default_paper_0710_1869_bs_deltaf2_system() -> Paper07101869DeltaF2System:
    """Return the custom B_s-system definition for paper-mode DeltaF=2 matching."""

    return Paper07101869DeltaF2System(
        system_id="B_s",
        sector_id="down",
        flavor_indices=(1, 2),
        label="Bs-Bsbar",
        notes=(
            "Custom B_s benchmark system: b-s mixing with left-down-aligned and "
            "right-diagonal KK-gluon matrices."
        ),
    )


def default_paper_0710_1869_d0_deltaf2_system() -> Paper07101869DeltaF2System:
    """Return the custom D0-system definition for paper-mode DeltaF=2 matching."""

    return Paper07101869DeltaF2System(
        system_id="D0",
        sector_id="up",
        flavor_indices=(0, 1),
        label="D0-D0bar",
        notes=(
            "Custom D0 benchmark system: c-u mixing with left-up-aligned and "
            "right-up-diagonal KK-gluon matrices."
        ),
    )


def match_paper_0710_1869_kk_gluon_deltaf2(
    system: Paper07101869DeltaF2System | None = None,
    *,
    benchmark: Paper07101869Benchmark | None = None,
    scale_point: Paper07101869ScalePoint | None = None,
    coupling_contract: Paper07101869CouplingContract | None = None,
    wilson_contract: Paper07101869DeltaF2WilsonContract | None = None,
    kk_gluon_couplings: Paper07101869KKGluonCouplings | None = None,
) -> Paper07101869KKGluonDeltaF2Match:
    """Match the minimal frozen four-operator basis at tree level."""
    resolved_system = (
        default_paper_0710_1869_kaon_deltaf2_system() if system is None else system
    )
    if not isinstance(resolved_system, Paper07101869DeltaF2System):
        raise ValueError("system must be a Paper07101869DeltaF2System")

    if kk_gluon_couplings is not None:
        if benchmark is not None or scale_point is not None or coupling_contract is not None:
            raise ValueError(
                "benchmark, scale_point, and coupling_contract cannot be passed with "
                "kk_gluon_couplings"
            )
        resolved_kk_gluon = kk_gluon_couplings
    else:
        resolved_benchmark = (
            default_paper_0710_1869_pr1_benchmark() if benchmark is None else benchmark
        )
        if not isinstance(resolved_benchmark, Paper07101869Benchmark):
            raise ValueError("benchmark must be a Paper07101869Benchmark")
        resolved_scale = default_paper_0710_1869_scales() if scale_point is None else scale_point
        if not isinstance(resolved_scale, Paper07101869ScalePoint):
            raise ValueError("scale_point must be a Paper07101869ScalePoint")
        resolved_coupling_contract = (
            default_paper_0710_1869_coupling_contract()
            if coupling_contract is None
            else coupling_contract
        )
        if not isinstance(resolved_coupling_contract, Paper07101869CouplingContract):
            raise ValueError("coupling_contract must be a Paper07101869CouplingContract")
        resolved_kk_gluon = build_paper_0710_1869_kk_gluon_couplings(
            resolved_benchmark,
            scale_point=resolved_scale,
            contract=resolved_coupling_contract,
        )

    if not isinstance(resolved_kk_gluon, Paper07101869KKGluonCouplings):
        raise ValueError("kk_gluon_couplings must be a Paper07101869KKGluonCouplings")

    resolved_wilson_contract = (
        default_paper_0710_1869_deltaf2_wilson_contract()
        if wilson_contract is None
        else wilson_contract
    )
    if not isinstance(resolved_wilson_contract, Paper07101869DeltaF2WilsonContract):
        raise ValueError("wilson_contract must be a Paper07101869DeltaF2WilsonContract")

    _require_matching_contract_compatibility(resolved_wilson_contract, resolved_kk_gluon)

    left_matrix, right_matrix = _select_chiral_matrices(resolved_kk_gluon, resolved_system)
    i, j = resolved_system.flavor_indices
    left_coupling = complex(left_matrix.universal_subtracted_gs_normalized[i, j])
    right_coupling = complex(right_matrix.universal_subtracted_gs_normalized[i, j])
    propagator_mass = resolved_kk_gluon.normalization.propagator_mass_GeV
    prefactor = 1.0 / (propagator_mass**2)
    matching_scale = resolved_kk_gluon.normalization.mu_match_GeV

    wilsons = Paper07101869DeltaF2WilsonCoefficients(
        contract=resolved_wilson_contract,
        benchmark_id=resolved_kk_gluon.benchmark_id,
        scale_label=resolved_kk_gluon.normalization.scale_label,
        system_id=resolved_system.system_id,
        sector_id=resolved_system.sector_id,
        generations=resolved_system.flavor_indices,
        matching_scale_GeV=matching_scale,
        propagator_mass_GeV=propagator_mass,
        left_coupling=left_coupling,
        right_coupling=right_coupling,
        q1_vll=(left_coupling * left_coupling) * (prefactor / 6.0),
        q1_vrr=(right_coupling * right_coupling) * (prefactor / 6.0),
        q4_lr=-(left_coupling * right_coupling) * prefactor,
        q5_lr=(left_coupling * right_coupling) * (prefactor / 3.0),
    )

    return Paper07101869KKGluonDeltaF2Match(
        system=resolved_system,
        benchmark_id=resolved_kk_gluon.benchmark_id,
        benchmark_label=resolved_kk_gluon.benchmark_label,
        scale_label=resolved_kk_gluon.normalization.scale_label,
        left_matrix_label=left_matrix.label,
        left_basis_id=left_matrix.basis_id,
        right_matrix_label=right_matrix.label,
        right_basis_id=right_matrix.basis_id,
        left_coupling=left_coupling,
        right_coupling=right_coupling,
        prefactor_GeV_minus2=prefactor,
        matching_scale_GeV=matching_scale,
        propagator_mass_GeV=propagator_mass,
        wilsons=wilsons,
    )


def match_default_paper_0710_1869_kaon_kk_gluon_deltaf2() -> Paper07101869KKGluonDeltaF2Match:
    """Return the default kaon matching point for the frozen PR3 slice."""
    return match_paper_0710_1869_kk_gluon_deltaf2()


def match_default_paper_0710_1869_bd_kk_gluon_deltaf2() -> Paper07101869KKGluonDeltaF2Match:
    """Return the custom B_d matching point for the paper-owned Q1 slice."""

    return match_paper_0710_1869_kk_gluon_deltaf2(
        system=default_paper_0710_1869_bd_deltaf2_system()
    )


def match_default_paper_0710_1869_bs_kk_gluon_deltaf2() -> Paper07101869KKGluonDeltaF2Match:
    """Return the custom B_s matching point for the paper-owned Q1 slice."""

    return match_paper_0710_1869_kk_gluon_deltaf2(
        system=default_paper_0710_1869_bs_deltaf2_system()
    )


def match_default_paper_0710_1869_d0_kk_gluon_deltaf2() -> Paper07101869KKGluonDeltaF2Match:
    """Return the custom D0 matching point for the paper-owned Q1 slice."""

    return match_paper_0710_1869_kk_gluon_deltaf2(
        system=default_paper_0710_1869_d0_deltaf2_system()
    )


def match_kkgluon_to_deltaf2(
    *,
    system: Paper07101869DeltaF2System | None = None,
    benchmark: Paper07101869Benchmark | None = None,
    scale_point: Paper07101869ScalePoint | None = None,
    coupling_contract: Paper07101869CouplingContract | None = None,
    wilson_contract: Paper07101869DeltaF2WilsonContract | None = None,
    kk_couplings: Paper07101869KKGluonCouplings | None = None,
    kk_gluon_couplings: Paper07101869KKGluonCouplings | None = None,
    couplings: Paper07101869KKGluonCouplings | None = None,
) -> Paper07101869KKGluonDeltaF2Match:
    """Compatibility alias for parameterized KK-gluon to Delta F=2 matching."""
    resolved_couplings = _resolve_alias_kk_gluon_couplings(
        kk_couplings=kk_couplings,
        kk_gluon_couplings=kk_gluon_couplings,
        couplings=couplings,
    )
    return match_paper_0710_1869_kk_gluon_deltaf2(
        system=system,
        benchmark=benchmark,
        scale_point=scale_point,
        coupling_contract=coupling_contract,
        wilson_contract=wilson_contract,
        kk_gluon_couplings=resolved_couplings,
    )


def match_paper_0710_1869_kkgluon_to_deltaf2(
    *,
    system: Paper07101869DeltaF2System | None = None,
    benchmark: Paper07101869Benchmark | None = None,
    scale_point: Paper07101869ScalePoint | None = None,
    coupling_contract: Paper07101869CouplingContract | None = None,
    wilson_contract: Paper07101869DeltaF2WilsonContract | None = None,
    kk_couplings: Paper07101869KKGluonCouplings | None = None,
    kk_gluon_couplings: Paper07101869KKGluonCouplings | None = None,
    couplings: Paper07101869KKGluonCouplings | None = None,
) -> Paper07101869KKGluonDeltaF2Match:
    """Paper-namespaced alias for parameterized KK-gluon to Delta F=2 matching."""
    return match_kkgluon_to_deltaf2(
        system=system,
        benchmark=benchmark,
        scale_point=scale_point,
        coupling_contract=coupling_contract,
        wilson_contract=wilson_contract,
        kk_couplings=kk_couplings,
        kk_gluon_couplings=kk_gluon_couplings,
        couplings=couplings,
    )


def build_paper_0710_1869_kaon_matching(
    *,
    kk_couplings: Paper07101869KKGluonCouplings | None = None,
    kk_gluon_couplings: Paper07101869KKGluonCouplings | None = None,
    couplings: Paper07101869KKGluonCouplings | None = None,
    benchmark: Paper07101869Benchmark | None = None,
    scale_point: Paper07101869ScalePoint | None = None,
    coupling_contract: Paper07101869CouplingContract | None = None,
    wilson_contract: Paper07101869DeltaF2WilsonContract | None = None,
) -> Paper07101869KKGluonDeltaF2Match:
    """Build the deterministic default kaon matching object."""
    return match_kkgluon_to_deltaf2(
        benchmark=benchmark,
        scale_point=scale_point,
        coupling_contract=coupling_contract,
        wilson_contract=wilson_contract,
        kk_couplings=kk_couplings,
        kk_gluon_couplings=kk_gluon_couplings,
        couplings=couplings,
    )


def default_paper_0710_1869_kaon_matching() -> Paper07101869KKGluonDeltaF2Match:
    """Return the deterministic default kaon matching object."""
    return build_paper_0710_1869_kaon_matching()


def build_paper_0710_1869_bd_matching(
    *,
    kk_couplings: Paper07101869KKGluonCouplings | None = None,
    kk_gluon_couplings: Paper07101869KKGluonCouplings | None = None,
    couplings: Paper07101869KKGluonCouplings | None = None,
    benchmark: Paper07101869Benchmark | None = None,
    scale_point: Paper07101869ScalePoint | None = None,
    coupling_contract: Paper07101869CouplingContract | None = None,
    wilson_contract: Paper07101869DeltaF2WilsonContract | None = None,
) -> Paper07101869KKGluonDeltaF2Match:
    """Build the deterministic custom B_d matching object."""

    return match_kkgluon_to_deltaf2(
        system=default_paper_0710_1869_bd_deltaf2_system(),
        benchmark=benchmark,
        scale_point=scale_point,
        coupling_contract=coupling_contract,
        wilson_contract=wilson_contract,
        kk_couplings=kk_couplings,
        kk_gluon_couplings=kk_gluon_couplings,
        couplings=couplings,
    )


def default_paper_0710_1869_bd_matching() -> Paper07101869KKGluonDeltaF2Match:
    """Return the deterministic custom B_d matching object."""

    return build_paper_0710_1869_bd_matching()


def build_paper_0710_1869_bs_matching(
    *,
    kk_couplings: Paper07101869KKGluonCouplings | None = None,
    kk_gluon_couplings: Paper07101869KKGluonCouplings | None = None,
    couplings: Paper07101869KKGluonCouplings | None = None,
    benchmark: Paper07101869Benchmark | None = None,
    scale_point: Paper07101869ScalePoint | None = None,
    coupling_contract: Paper07101869CouplingContract | None = None,
    wilson_contract: Paper07101869DeltaF2WilsonContract | None = None,
) -> Paper07101869KKGluonDeltaF2Match:
    """Build the deterministic custom B_s matching object."""

    return match_kkgluon_to_deltaf2(
        system=default_paper_0710_1869_bs_deltaf2_system(),
        benchmark=benchmark,
        scale_point=scale_point,
        coupling_contract=coupling_contract,
        wilson_contract=wilson_contract,
        kk_couplings=kk_couplings,
        kk_gluon_couplings=kk_gluon_couplings,
        couplings=couplings,
    )


def default_paper_0710_1869_bs_matching() -> Paper07101869KKGluonDeltaF2Match:
    """Return the deterministic custom B_s matching object."""

    return build_paper_0710_1869_bs_matching()


def build_paper_0710_1869_d0_matching(
    *,
    kk_couplings: Paper07101869KKGluonCouplings | None = None,
    kk_gluon_couplings: Paper07101869KKGluonCouplings | None = None,
    couplings: Paper07101869KKGluonCouplings | None = None,
    benchmark: Paper07101869Benchmark | None = None,
    scale_point: Paper07101869ScalePoint | None = None,
    coupling_contract: Paper07101869CouplingContract | None = None,
    wilson_contract: Paper07101869DeltaF2WilsonContract | None = None,
) -> Paper07101869KKGluonDeltaF2Match:
    """Build the deterministic custom D0 matching object."""

    return match_kkgluon_to_deltaf2(
        system=default_paper_0710_1869_d0_deltaf2_system(),
        benchmark=benchmark,
        scale_point=scale_point,
        coupling_contract=coupling_contract,
        wilson_contract=wilson_contract,
        kk_couplings=kk_couplings,
        kk_gluon_couplings=kk_gluon_couplings,
        couplings=couplings,
    )


def default_paper_0710_1869_d0_matching() -> Paper07101869KKGluonDeltaF2Match:
    """Return the deterministic custom D0 matching object."""

    return build_paper_0710_1869_d0_matching()


def build_paper_0710_1869_deltaf2_matching(
    *,
    kk_couplings: Paper07101869KKGluonCouplings | None = None,
    kk_gluon_couplings: Paper07101869KKGluonCouplings | None = None,
    couplings: Paper07101869KKGluonCouplings | None = None,
    benchmark: Paper07101869Benchmark | None = None,
    scale_point: Paper07101869ScalePoint | None = None,
    coupling_contract: Paper07101869CouplingContract | None = None,
    wilson_contract: Paper07101869DeltaF2WilsonContract | None = None,
) -> Paper07101869KKGluonDeltaF2Match:
    """Compatibility alias for the deterministic default Delta F=2 match object."""
    return build_paper_0710_1869_kaon_matching(
        kk_couplings=kk_couplings,
        kk_gluon_couplings=kk_gluon_couplings,
        couplings=couplings,
        benchmark=benchmark,
        scale_point=scale_point,
        coupling_contract=coupling_contract,
        wilson_contract=wilson_contract,
    )


def default_paper_0710_1869_deltaf2_matching() -> Paper07101869KKGluonDeltaF2Match:
    """Return the deterministic default Delta F=2 match object."""
    return build_paper_0710_1869_deltaf2_matching()


def build_paper_0710_1869_kkgluon_matching(
    *,
    kk_couplings: Paper07101869KKGluonCouplings | None = None,
    kk_gluon_couplings: Paper07101869KKGluonCouplings | None = None,
    couplings: Paper07101869KKGluonCouplings | None = None,
    benchmark: Paper07101869Benchmark | None = None,
    scale_point: Paper07101869ScalePoint | None = None,
    coupling_contract: Paper07101869CouplingContract | None = None,
    wilson_contract: Paper07101869DeltaF2WilsonContract | None = None,
) -> Paper07101869KKGluonDeltaF2Match:
    """Compatibility alias for the deterministic default KK-gluon match object."""
    return build_paper_0710_1869_kaon_matching(
        kk_couplings=kk_couplings,
        kk_gluon_couplings=kk_gluon_couplings,
        couplings=couplings,
        benchmark=benchmark,
        scale_point=scale_point,
        coupling_contract=coupling_contract,
        wilson_contract=wilson_contract,
    )


def default_paper_0710_1869_kkgluon_matching() -> Paper07101869KKGluonDeltaF2Match:
    """Return the deterministic default KK-gluon match object."""
    return build_paper_0710_1869_kkgluon_matching()


def build_paper_0710_1869_kaon_matching_summary(
    *,
    kk_couplings: Paper07101869KKGluonCouplings | None = None,
    kk_gluon_couplings: Paper07101869KKGluonCouplings | None = None,
    couplings: Paper07101869KKGluonCouplings | None = None,
    benchmark: Paper07101869Benchmark | None = None,
    scale_point: Paper07101869ScalePoint | None = None,
    coupling_contract: Paper07101869CouplingContract | None = None,
    wilson_contract: Paper07101869DeltaF2WilsonContract | None = None,
) -> dict[str, object]:
    """Build the deterministic default kaon matching summary mapping."""
    return build_paper_0710_1869_kaon_matching(
        kk_couplings=kk_couplings,
        kk_gluon_couplings=kk_gluon_couplings,
        couplings=couplings,
        benchmark=benchmark,
        scale_point=scale_point,
        coupling_contract=coupling_contract,
        wilson_contract=wilson_contract,
    ).summary()


def default_paper_0710_1869_kaon_matching_summary() -> dict[str, object]:
    """Return the deterministic default kaon matching summary mapping."""
    return build_paper_0710_1869_kaon_matching_summary()


def build_paper_0710_1869_bd_matching_summary(
    *,
    kk_couplings: Paper07101869KKGluonCouplings | None = None,
    kk_gluon_couplings: Paper07101869KKGluonCouplings | None = None,
    couplings: Paper07101869KKGluonCouplings | None = None,
    benchmark: Paper07101869Benchmark | None = None,
    scale_point: Paper07101869ScalePoint | None = None,
    coupling_contract: Paper07101869CouplingContract | None = None,
    wilson_contract: Paper07101869DeltaF2WilsonContract | None = None,
) -> dict[str, object]:
    """Build the deterministic custom B_d matching summary mapping."""

    return build_paper_0710_1869_bd_matching(
        kk_couplings=kk_couplings,
        kk_gluon_couplings=kk_gluon_couplings,
        couplings=couplings,
        benchmark=benchmark,
        scale_point=scale_point,
        coupling_contract=coupling_contract,
        wilson_contract=wilson_contract,
    ).summary()


def default_paper_0710_1869_bd_matching_summary() -> dict[str, object]:
    """Return the deterministic custom B_d matching summary mapping."""

    return build_paper_0710_1869_bd_matching_summary()


def build_paper_0710_1869_bs_matching_summary(
    *,
    kk_couplings: Paper07101869KKGluonCouplings | None = None,
    kk_gluon_couplings: Paper07101869KKGluonCouplings | None = None,
    couplings: Paper07101869KKGluonCouplings | None = None,
    benchmark: Paper07101869Benchmark | None = None,
    scale_point: Paper07101869ScalePoint | None = None,
    coupling_contract: Paper07101869CouplingContract | None = None,
    wilson_contract: Paper07101869DeltaF2WilsonContract | None = None,
) -> dict[str, object]:
    """Build the deterministic custom B_s matching summary mapping."""

    return build_paper_0710_1869_bs_matching(
        kk_couplings=kk_couplings,
        kk_gluon_couplings=kk_gluon_couplings,
        couplings=couplings,
        benchmark=benchmark,
        scale_point=scale_point,
        coupling_contract=coupling_contract,
        wilson_contract=wilson_contract,
    ).summary()


def default_paper_0710_1869_bs_matching_summary() -> dict[str, object]:
    """Return the deterministic custom B_s matching summary mapping."""

    return build_paper_0710_1869_bs_matching_summary()


def build_paper_0710_1869_d0_matching_summary(
    *,
    kk_couplings: Paper07101869KKGluonCouplings | None = None,
    kk_gluon_couplings: Paper07101869KKGluonCouplings | None = None,
    couplings: Paper07101869KKGluonCouplings | None = None,
    benchmark: Paper07101869Benchmark | None = None,
    scale_point: Paper07101869ScalePoint | None = None,
    coupling_contract: Paper07101869CouplingContract | None = None,
    wilson_contract: Paper07101869DeltaF2WilsonContract | None = None,
) -> dict[str, object]:
    """Build the deterministic custom D0 matching summary mapping."""

    return build_paper_0710_1869_d0_matching(
        kk_couplings=kk_couplings,
        kk_gluon_couplings=kk_gluon_couplings,
        couplings=couplings,
        benchmark=benchmark,
        scale_point=scale_point,
        coupling_contract=coupling_contract,
        wilson_contract=wilson_contract,
    ).summary()


def default_paper_0710_1869_d0_matching_summary() -> dict[str, object]:
    """Return the deterministic custom D0 matching summary mapping."""

    return build_paper_0710_1869_d0_matching_summary()


def build_paper_0710_1869_deltaf2_matching_summary(
    *,
    kk_couplings: Paper07101869KKGluonCouplings | None = None,
    kk_gluon_couplings: Paper07101869KKGluonCouplings | None = None,
    couplings: Paper07101869KKGluonCouplings | None = None,
    benchmark: Paper07101869Benchmark | None = None,
    scale_point: Paper07101869ScalePoint | None = None,
    coupling_contract: Paper07101869CouplingContract | None = None,
    wilson_contract: Paper07101869DeltaF2WilsonContract | None = None,
) -> dict[str, object]:
    """Compatibility alias for the deterministic default Delta F=2 summary mapping."""
    return build_paper_0710_1869_kaon_matching_summary(
        kk_couplings=kk_couplings,
        kk_gluon_couplings=kk_gluon_couplings,
        couplings=couplings,
        benchmark=benchmark,
        scale_point=scale_point,
        coupling_contract=coupling_contract,
        wilson_contract=wilson_contract,
    )


def default_paper_0710_1869_deltaf2_matching_summary() -> dict[str, object]:
    """Return the deterministic default Delta F=2 summary mapping."""
    return build_paper_0710_1869_deltaf2_matching_summary()


def build_paper_0710_1869_kkgluon_matching_summary(
    *,
    kk_couplings: Paper07101869KKGluonCouplings | None = None,
    kk_gluon_couplings: Paper07101869KKGluonCouplings | None = None,
    couplings: Paper07101869KKGluonCouplings | None = None,
    benchmark: Paper07101869Benchmark | None = None,
    scale_point: Paper07101869ScalePoint | None = None,
    coupling_contract: Paper07101869CouplingContract | None = None,
    wilson_contract: Paper07101869DeltaF2WilsonContract | None = None,
) -> dict[str, object]:
    """Compatibility alias for the deterministic default KK-gluon summary mapping."""
    return build_paper_0710_1869_kaon_matching_summary(
        kk_couplings=kk_couplings,
        kk_gluon_couplings=kk_gluon_couplings,
        couplings=couplings,
        benchmark=benchmark,
        scale_point=scale_point,
        coupling_contract=coupling_contract,
        wilson_contract=wilson_contract,
    )


def default_paper_0710_1869_kkgluon_matching_summary() -> dict[str, object]:
    """Return the deterministic default KK-gluon summary mapping."""
    return build_paper_0710_1869_kkgluon_matching_summary()


__all__ = [
    "PAPER_0710_1869_DELTAF2_KKGLUON_MATCH_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_MATCHING_SUMMARY_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_SYSTEM_SCHEMA_ID",
    "Paper07101869DeltaF2System",
    "Paper07101869KKGluonDeltaF2Match",
    "build_paper_0710_1869_deltaf2_matching",
    "build_paper_0710_1869_deltaf2_matching_summary",
    "build_paper_0710_1869_bd_matching",
    "build_paper_0710_1869_bd_matching_summary",
    "build_paper_0710_1869_bs_matching",
    "build_paper_0710_1869_bs_matching_summary",
    "build_paper_0710_1869_d0_matching",
    "build_paper_0710_1869_d0_matching_summary",
    "build_paper_0710_1869_kaon_matching",
    "build_paper_0710_1869_kaon_matching_summary",
    "build_paper_0710_1869_kkgluon_matching",
    "build_paper_0710_1869_kkgluon_matching_summary",
    "default_paper_0710_1869_bd_deltaf2_system",
    "default_paper_0710_1869_bd_matching",
    "default_paper_0710_1869_bd_matching_summary",
    "default_paper_0710_1869_bs_deltaf2_system",
    "default_paper_0710_1869_bs_matching",
    "default_paper_0710_1869_bs_matching_summary",
    "default_paper_0710_1869_d0_deltaf2_system",
    "default_paper_0710_1869_d0_matching",
    "default_paper_0710_1869_d0_matching_summary",
    "default_paper_0710_1869_deltaf2_matching",
    "default_paper_0710_1869_deltaf2_matching_summary",
    "default_paper_0710_1869_kaon_deltaf2_system",
    "default_paper_0710_1869_kaon_matching",
    "default_paper_0710_1869_kaon_matching_summary",
    "default_paper_0710_1869_kkgluon_matching",
    "default_paper_0710_1869_kkgluon_matching_summary",
    "match_kkgluon_to_deltaf2",
    "match_default_paper_0710_1869_bd_kk_gluon_deltaf2",
    "match_default_paper_0710_1869_bs_kk_gluon_deltaf2",
    "match_default_paper_0710_1869_d0_kk_gluon_deltaf2",
    "match_default_paper_0710_1869_kaon_kk_gluon_deltaf2",
    "match_paper_0710_1869_kkgluon_to_deltaf2",
    "match_paper_0710_1869_kk_gluon_deltaf2",
]
