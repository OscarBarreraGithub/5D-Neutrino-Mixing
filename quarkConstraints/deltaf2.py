"""Tree-level KK-gluon-inspired ``Delta F = 2`` benchmark observables.

This module intentionally implements a compact, repo-owned v1 matching layer:

- exact fitted quark mass-basis couplings come from ``quarkConstraints.couplings``
- KK-gluon exchange is represented by a small fixed operator basis
- hadronic/bound inputs live in a deterministic in-repo bundle

The goal is to turn the current MFV fit machinery into a first exclusion slice,
not to provide a general-purpose EFT/RG package.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from .couplings import QuarkMassBasisCouplings, compute_quark_kk_gluon_couplings
from .fit import QuarkFitResult
from .scales import DEFAULT_QUARK_TARGET_SCALE_GEV, DEFAULT_QUARK_XI_KK

DELTA_F2_MODEL_V1 = "kk_gluon_tree_v1"
DELTA_F2_INPUT_BUNDLE_V1 = "deltaf2_inputs_mu3tev_v1"
DELTA_F2_OPERATOR_CONVENTION = DELTA_F2_MODEL_V1
DELTA_F2_INPUT_BUNDLE = DELTA_F2_INPUT_BUNDLE_V1


@dataclass(frozen=True)
class DeltaF2Input:
    """Repo-owned v1 input for one neutral-meson mixing system."""

    key: str
    display_name: str
    column_name: str
    reject_reason: str
    sector: str
    generations: tuple[int, int]
    bound: float
    ll_weight: float = 1.0
    rr_weight: float = 1.0
    lr1_weight: float = 7.0
    lr2_weight: float = 2.0
    reference_scale: float = DEFAULT_QUARK_TARGET_SCALE_GEV
    note: str = ""

    def __post_init__(self) -> None:
        if self.sector not in {"down", "up"}:
            raise ValueError("sector must be 'down' or 'up'")
        if len(self.generations) != 2:
            raise ValueError("generations must be a pair of flavor indices")
        i, j = self.generations
        if i == j or not (0 <= i < 3) or not (0 <= j < 3):
            raise ValueError("generations must be distinct indices in {0,1,2}")
        if self.bound <= 0.0:
            raise ValueError("bound must be positive")
        if self.reference_scale <= 0.0:
            raise ValueError("reference_scale must be positive")

    @property
    def system(self) -> str:
        aliases = {
            "epsilon_k": "K",
            "b_d": "B_d",
            "b_s": "B_s",
            "d": "D",
        }
        return aliases.get(self.key, self.display_name)

    @property
    def flavor_indices(self) -> tuple[int, int]:
        return self.generations


@dataclass(frozen=True)
class DeltaF2WilsonCoefficients:
    """Wilson coefficients for one neutral-meson system."""

    input: DeltaF2Input
    M_KK: float
    matching_scale: float
    left_coupling: complex
    right_coupling: complex
    c1_vll: complex
    c1_vrr: complex
    c4_lr: complex
    c5_lr: complex

    @property
    def wilsons(self) -> Mapping[str, complex]:
        return {
            "C1_VLL": self.c1_vll,
            "C1_VRR": self.c1_vrr,
            "C4_LR": self.c4_lr,
            "C5_LR": self.c5_lr,
        }


@dataclass(frozen=True)
class DeltaF2ObservableSummary:
    """Compact per-system exclusion summary."""

    input: DeltaF2Input
    wilsons: DeltaF2WilsonCoefficients
    effective_amplitude: float
    ratio_to_bound: float
    passes: bool
    dominant_operator: str
    dominant_operator_size: float
    weighted_operator_sizes: Mapping[str, float]

    @property
    def key(self) -> str:
        return self.input.key

    @property
    def system(self) -> str:
        return self.input.system

    @property
    def display_name(self) -> str:
        return self.input.display_name

    @property
    def column_name(self) -> str:
        return self.input.column_name

    @property
    def reject_reason(self) -> str:
        return self.input.reject_reason

    @property
    def bound(self) -> float:
        return self.input.bound


@dataclass(frozen=True)
class DeltaF2ConstraintSummary:
    """Overall ``Delta F = 2`` evaluation for one fitted MFV point."""

    model_label: str
    input_bundle_label: str
    M_KK: float
    xi_KK: float
    observables: tuple[DeltaF2ObservableSummary, ...]

    @property
    def passes_all(self) -> bool:
        return all(item.passes for item in self.observables)

    @property
    def worst_ratio(self) -> float:
        if not self.observables:
            return 0.0
        return float(max(item.ratio_to_bound for item in self.observables))

    @property
    def max_ratio_to_bound(self) -> float:
        return self.worst_ratio

    @property
    def operator_convention(self) -> str:
        return self.model_label

    @property
    def input_bundle(self) -> str:
        return self.input_bundle_label

    @property
    def matching_scale(self) -> float:
        return self.M_KK

    @property
    def failing_reasons(self) -> tuple[str, ...]:
        return tuple(item.reject_reason for item in self.observables if not item.passes)

    @property
    def reject_reasons(self) -> tuple[str, ...]:
        return self.failing_reasons

    @property
    def by_system(self) -> dict[str, DeltaF2ObservableSummary]:
        aliases = {
            "epsilon_k": "K",
            "b_d": "B_d",
            "b_s": "B_s",
            "d": "D",
        }
        return {aliases.get(item.key, item.display_name): item for item in self.observables}

    def as_ratio_dict(self) -> dict[str, float]:
        return {item.column_name: float(item.ratio_to_bound) for item in self.observables}

    def get(self, key: str) -> DeltaF2ObservableSummary:
        for item in self.observables:
            if item.key == key or item.column_name == key or item.reject_reason == key:
                return item
        raise KeyError(key)


DEFAULT_DELTA_F2_INPUTS_V1: tuple[DeltaF2Input, ...] = (
    DeltaF2Input(
        key="epsilon_k",
        display_name="epsilon_K",
        column_name="epsilon_k_ratio",
        reject_reason="epsilon_k",
        sector="down",
        generations=(0, 1),
        bound=2.0e-8,
        ll_weight=1.0,
        rr_weight=1.0,
        lr1_weight=7.0,
        lr2_weight=2.0,
        note="Kaon mixing uses the strongest LR enhancement in the v1 input bundle.",
    ),
    DeltaF2Input(
        key="b_d",
        display_name="B_d mixing",
        column_name="b_d_ratio",
        reject_reason="b_d_mix",
        sector="down",
        generations=(0, 2),
        bound=4.0e-7,
        ll_weight=1.0,
        rr_weight=1.0,
        lr1_weight=7.0,
        lr2_weight=2.0,
        note="Repo-owned B_d benchmark normalization at the 3 TeV matching scale.",
    ),
    DeltaF2Input(
        key="b_s",
        display_name="B_s mixing",
        column_name="b_s_ratio",
        reject_reason="b_s_mix",
        sector="down",
        generations=(1, 2),
        bound=5.5e-6,
        ll_weight=1.0,
        rr_weight=1.0,
        lr1_weight=7.0,
        lr2_weight=2.0,
        note="Repo-owned B_s benchmark normalization at the 3 TeV matching scale.",
    ),
    DeltaF2Input(
        key="d",
        display_name="D mixing",
        column_name="d_ratio",
        reject_reason="d_mix",
        sector="up",
        generations=(0, 1),
        bound=8.5e-9,
        ll_weight=1.0,
        rr_weight=1.0,
        lr1_weight=7.0,
        lr2_weight=2.0,
        note="Repo-owned D-mixing normalization for the first exclusion slice.",
    ),
)


def default_delta_f2_inputs() -> tuple[DeltaF2Input, ...]:
    """Return the repo-owned v1 input bundle."""
    return DEFAULT_DELTA_F2_INPUTS_V1


def _coerce_couplings(
    source: QuarkFitResult | QuarkMassBasisCouplings,
    *,
    M_KK: float | None,
    xi_KK: float,
) -> QuarkMassBasisCouplings:
    if isinstance(source, QuarkMassBasisCouplings):
        return source
    required = ("M_KK", "left_down", "right_down", "left_up", "right_up")
    if all(hasattr(source, name) for name in required):
        return source  # type: ignore[return-value]
    return compute_quark_kk_gluon_couplings(source, M_KK=M_KK, xi_KK=xi_KK)


def _pair_couplings(
    couplings: QuarkMassBasisCouplings,
    item: DeltaF2Input,
) -> tuple[complex, complex]:
    i, j = item.generations
    if item.sector == "down":
        left = complex(couplings.left_down[i, j])
        right = complex(couplings.right_down[i, j])
    else:
        left = complex(couplings.left_up[i, j])
        right = complex(couplings.right_up[i, j])
    return left, right


def compute_delta_f2_wilsons(
    source: QuarkFitResult | QuarkMassBasisCouplings,
    *,
    M_KK: float | None = None,
    xi_KK: float = DEFAULT_QUARK_XI_KK,
    inputs: Sequence[DeltaF2Input] | None = None,
) -> tuple[DeltaF2WilsonCoefficients, ...]:
    """Compute tree-level KK-gluon-inspired ``Delta F = 2`` Wilsons.

    The returned coefficients follow a compact v1 convention with one vector
    left-left operator, one vector right-right operator, and two left-right
    operators. The coefficients are matched directly at ``mu = M_KK``.
    """
    if xi_KK <= 0.0:
        raise ValueError("xi_KK must be positive")
    couplings = _coerce_couplings(source, M_KK=M_KK, xi_KK=xi_KK)
    items = default_delta_f2_inputs() if inputs is None else tuple(inputs)
    out: list[DeltaF2WilsonCoefficients] = []
    prefactor = 1.0 / (couplings.M_KK**2)
    for item in items:
        left, right = _pair_couplings(couplings, item)
        out.append(
            DeltaF2WilsonCoefficients(
                input=item,
                M_KK=float(couplings.M_KK),
                matching_scale=float(couplings.M_KK),
                left_coupling=left,
                right_coupling=right,
                c1_vll=left * left * prefactor / 6.0,
                c1_vrr=right * right * prefactor / 6.0,
                c4_lr=-(left * right) * prefactor,
                c5_lr=(left * right) * prefactor / 3.0,
            )
        )
    return tuple(out)


def evaluate_delta_f2_constraints(
    source: QuarkFitResult | QuarkMassBasisCouplings,
    *,
    M_KK: float | None = None,
    xi_KK: float = DEFAULT_QUARK_XI_KK,
    inputs: Sequence[DeltaF2Input] | None = None,
) -> DeltaF2ConstraintSummary:
    """Evaluate the repo-owned ``Delta F = 2`` benchmark bundle."""
    couplings = _coerce_couplings(source, M_KK=M_KK, xi_KK=xi_KK)
    coefficients = compute_delta_f2_wilsons(
        couplings,
        M_KK=couplings.M_KK,
        xi_KK=couplings.xi_KK,
        inputs=inputs,
    )
    observables: list[DeltaF2ObservableSummary] = []
    for coeffs in coefficients:
        item = coeffs.input
        weighted = {
            "C1_VLL": item.ll_weight * coeffs.c1_vll,
            "C1_VRR": item.rr_weight * coeffs.c1_vrr,
            "C4_LR": item.lr1_weight * coeffs.c4_lr,
            "C5_LR": item.lr2_weight * coeffs.c5_lr,
        }
        operator_sizes = {
            name: float(item.reference_scale**2 * abs(value))
            for name, value in weighted.items()
        }
        dominant_operator = max(operator_sizes, key=operator_sizes.get)
        effective_amplitude = float(item.reference_scale**2 * abs(sum(weighted.values())))
        ratio_to_bound = float(effective_amplitude / item.bound)
        observables.append(
            DeltaF2ObservableSummary(
                input=item,
                wilsons=coeffs,
                effective_amplitude=effective_amplitude,
                ratio_to_bound=ratio_to_bound,
                passes=ratio_to_bound <= 1.0,
                dominant_operator=dominant_operator,
                dominant_operator_size=operator_sizes[dominant_operator],
                weighted_operator_sizes=operator_sizes,
            )
        )
    return DeltaF2ConstraintSummary(
        model_label=DELTA_F2_MODEL_V1,
        input_bundle_label=DELTA_F2_INPUT_BUNDLE_V1,
        M_KK=float(couplings.M_KK),
        xi_KK=float(getattr(couplings, "xi_KK", xi_KK)),
        observables=tuple(observables),
    )


DeltaF2WilsonSet = DeltaF2WilsonCoefficients
