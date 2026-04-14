"""Closed modern phenomenology policy surface for the modern quark lane."""

from __future__ import annotations

import json
import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from ..qcd_running import evolve_deltaf2_wilsons
from .conventions import MODERN_INPUT_REGISTRY_SCHEMA_ID, MODERN_LANE_ID
from .inputs import (
    MODERN_DEFAULT_INPUTS_SCHEMA_ID,
    ModernDefaultInputs,
    default_modern_default_inputs,
)

MODERN_PHENOMENOLOGY_SCHEMA_ID = "quarkConstraints.modern.phenomenology.v1"
MODERN_PHENOMENOLOGY_POLICY_ID = "quarkConstraints.modern.phenomenology.policy.v1"
MODERN_PHENOMENOLOGY_SYSTEM_SCHEMA_ID = "quarkConstraints.modern.phenomenology.system.v1"
MODERN_PHENOMENOLOGY_SYSTEM_IDS = ("epsilon_K", "K", "B_d", "B_s", "D0")
MODERN_PHENOMENOLOGY_SYSTEM_POLICY_ID_TEMPLATE = (
    "quarkConstraints.modern.phenomenology.system.{system_id}.v1"
)
MODERN_POINT_PHENOMENOLOGY_SYSTEM_RESULT_SCHEMA_ID = (
    "quarkConstraints.modern.phenomenology.system_result.v1"
)
MODERN_POINT_PHENOMENOLOGY_ARTIFACT_SCHEMA_ID = (
    "quarkConstraints.modern.phenomenology.artifact.v1"
)
MODERN_POINT_PHENOMENOLOGY_ARTIFACT_SCHEMA_VERSION = 1
MODERN_POINT_PHENOMENOLOGY_RELEASE_SCOPE_ID = (
    "quarkConstraints.modern.phenomenology.release.full_deltaf2.v1"
)
MODERN_POINT_PHENOMENOLOGY_BRIDGE_SCHEMA_ID = (
    "quarkConstraints.modern.artifacts.bridge.v1"
)
MODERN_POINT_PHENOMENOLOGY_NON_CP_ACCEPTANCE_SYSTEM_IDS = ("epsilon_K", "K", "B_d", "B_s", "D0")
MODERN_POINT_PHENOMENOLOGY_SYSTEM_TREATMENT_IDS = {
    "epsilon_K": "cp_violation_epsilon_k",
    "K": "non_cp_mixing_amplitude_kaon",
    "B_d": "non_cp_mixing_amplitude_surrogate",
    "B_s": "non_cp_mixing_amplitude_surrogate",
    "D0": "conservative_non_cp_mixing_amplitude_surrogate",
}


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
    if not math.isfinite(numeric) or numeric <= 0.0:
        raise ValueError(f"{name} must be a positive finite float")
    return numeric


def _system_policy_id(system_id: str) -> str:
    system = _require_text("system_id", system_id)
    if system not in MODERN_PHENOMENOLOGY_SYSTEM_IDS:
        raise ValueError(
            f"system_id must be one of {MODERN_PHENOMENOLOGY_SYSTEM_IDS!r}"
        )
    return MODERN_PHENOMENOLOGY_SYSTEM_POLICY_ID_TEMPLATE.format(system_id=system)


@dataclass(frozen=True)
class ModernPhenomenologySystemPolicy:
    """Frozen policy entry for one explicit modern phenomenology system."""

    schema_id: str = MODERN_PHENOMENOLOGY_SYSTEM_SCHEMA_ID
    system_id: str = ""
    policy_id: str = ""
    display_name: str = ""
    notes: str = ""

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            _require_exact(
                "schema_id",
                self.schema_id,
                expected=MODERN_PHENOMENOLOGY_SYSTEM_SCHEMA_ID,
            ),
        )
        object.__setattr__(self, "system_id", _require_text("system_id", self.system_id))
        object.__setattr__(
            self,
            "policy_id",
            _require_exact("policy_id", self.policy_id, expected=_system_policy_id(self.system_id)),
        )
        object.__setattr__(self, "display_name", _require_text("display_name", self.display_name))
        object.__setattr__(self, "notes", _require_text("notes", self.notes))

    def as_dict(self) -> dict[str, str]:
        return asdict(self)


def default_modern_phenomenology_system_policies() -> tuple[ModernPhenomenologySystemPolicy, ...]:
    """Return the frozen explicit system policy entries for the modern lane."""
    return DEFAULT_MODERN_PHENOMENOLOGY_SYSTEM_POLICIES


DEFAULT_MODERN_PHENOMENOLOGY_SYSTEM_POLICIES = (
    ModernPhenomenologySystemPolicy(
        system_id="epsilon_K",
        policy_id=_system_policy_id("epsilon_K"),
        display_name="epsilon_K CP-violating kaon-mixing policy",
        notes="Explicit CP-violating neutral-kaon mixing policy entry.",
    ),
    ModernPhenomenologySystemPolicy(
        system_id="K",
        policy_id=_system_policy_id("K"),
        display_name="K neutral-kaon mixing policy",
        notes="Explicit neutral-kaon mixing policy entry.",
    ),
    ModernPhenomenologySystemPolicy(
        system_id="B_d",
        policy_id=_system_policy_id("B_d"),
        display_name="B_d neutral-B mixing policy",
        notes="Explicit B_d neutral-meson mixing policy entry.",
    ),
    ModernPhenomenologySystemPolicy(
        system_id="B_s",
        policy_id=_system_policy_id("B_s"),
        display_name="B_s neutral-B mixing policy",
        notes="Explicit B_s neutral-meson mixing policy entry.",
    ),
    ModernPhenomenologySystemPolicy(
        system_id="D0",
        policy_id=_system_policy_id("D0"),
        display_name="D0 neutral-D mixing policy",
        notes="Explicit D0 neutral-meson mixing policy entry.",
    ),
)

MODERN_PHENOMENOLOGY_SYSTEM_POLICY_IDS = tuple(
    system.policy_id for system in DEFAULT_MODERN_PHENOMENOLOGY_SYSTEM_POLICIES
)


@dataclass(frozen=True)
class ModernPhenomenologyPolicy:
    """Closed versioned modern phenomenology policy surface."""

    schema_id: str = MODERN_PHENOMENOLOGY_SCHEMA_ID
    policy_id: str = MODERN_PHENOMENOLOGY_POLICY_ID
    lane_id: str = MODERN_LANE_ID
    registry_schema_id: str = MODERN_INPUT_REGISTRY_SCHEMA_ID
    systems: tuple[ModernPhenomenologySystemPolicy, ...] = field(
        default_factory=default_modern_phenomenology_system_policies
    )
    notes: str = (
        "Policy-only modern phenomenology fence. It enumerates the explicit "
        "systems epsilon_K, K, B_d, B_s, and D0, but it does not define any "
        "numeric EFT backend, scan engine, or phenomenology verdict model."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            _require_exact("schema_id", self.schema_id, expected=MODERN_PHENOMENOLOGY_SCHEMA_ID),
        )
        object.__setattr__(
            self,
            "policy_id",
            _require_exact("policy_id", self.policy_id, expected=MODERN_PHENOMENOLOGY_POLICY_ID),
        )
        object.__setattr__(self, "lane_id", _require_exact("lane_id", self.lane_id, expected=MODERN_LANE_ID))
        object.__setattr__(
            self,
            "registry_schema_id",
            _require_exact(
                "registry_schema_id",
                self.registry_schema_id,
                expected=MODERN_INPUT_REGISTRY_SCHEMA_ID,
            ),
        )
        object.__setattr__(self, "notes", _require_text("notes", self.notes))
        if not isinstance(self.systems, tuple):
            raise ValueError("systems must be a tuple of ModernPhenomenologySystemPolicy instances")
        normalized_systems = tuple(self.systems)
        if len(normalized_systems) != len(MODERN_PHENOMENOLOGY_SYSTEM_IDS):
            raise ValueError(
                "systems must contain exactly epsilon_K, K, B_d, B_s, and D0"
            )
        for system in normalized_systems:
            if not isinstance(system, ModernPhenomenologySystemPolicy):
                raise ValueError(
                    "systems must contain only ModernPhenomenologySystemPolicy instances"
                )
        normalized_ids = tuple(system.system_id for system in normalized_systems)
        if normalized_ids != MODERN_PHENOMENOLOGY_SYSTEM_IDS:
            unknown = tuple(system_id for system_id in normalized_ids if system_id not in MODERN_PHENOMENOLOGY_SYSTEM_IDS)
            missing = tuple(
                system_id
                for system_id in MODERN_PHENOMENOLOGY_SYSTEM_IDS
                if system_id not in normalized_ids
            )
            if unknown:
                raise ValueError(f"systems contains unknown systems: {unknown!r}")
            if missing:
                raise ValueError(f"systems is missing systems: {missing!r}")
            raise ValueError("systems must preserve the frozen modern phenomenology order")
        if tuple(system.policy_id for system in normalized_systems) != MODERN_PHENOMENOLOGY_SYSTEM_POLICY_IDS:
            raise ValueError("systems must preserve the frozen modern phenomenology policy ids")
        object.__setattr__(self, "systems", normalized_systems)

    @property
    def system_ids(self) -> tuple[str, ...]:
        return tuple(system.system_id for system in self.systems)

    def system_policy(self, system_id: str) -> ModernPhenomenologySystemPolicy:
        """Return the explicit policy entry for one known system."""
        requested = _require_text("system_id", system_id)
        for system in self.systems:
            if system.system_id == requested:
                return system
        raise ValueError(
            f"system_id must be one of {MODERN_PHENOMENOLOGY_SYSTEM_IDS!r}"
        )

    def require_system_ids(self, system_ids: Iterable[str]) -> tuple[str, ...]:
        """Require the exact frozen modern phenomenology system set."""
        provided = tuple(_require_text("system_id", system_id) for system_id in system_ids)
        unknown = tuple(system_id for system_id in provided if system_id not in MODERN_PHENOMENOLOGY_SYSTEM_IDS)
        missing = tuple(
            system_id
            for system_id in MODERN_PHENOMENOLOGY_SYSTEM_IDS
            if system_id not in provided
        )
        if unknown:
            raise ValueError(f"system_ids contains unknown systems: {unknown!r}")
        if missing:
            raise ValueError(f"system_ids is missing systems: {missing!r}")
        if provided != MODERN_PHENOMENOLOGY_SYSTEM_IDS:
            raise ValueError(
                "system_ids must preserve the frozen modern phenomenology order"
            )
        return provided

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "policy_id": self.policy_id,
            "lane_id": self.lane_id,
            "registry_schema_id": self.registry_schema_id,
            "systems": [system.as_dict() for system in self.systems],
            "notes": self.notes,
        }


def default_modern_phenomenology_policy() -> ModernPhenomenologyPolicy:
    """Return the frozen modern phenomenology policy surface."""
    return ModernPhenomenologyPolicy()


def _require_mapping(value: Any, *, context: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{context} must be a JSON object")
    return value


def _require_sequence(value: Any, *, context: str) -> Sequence[Any]:
    if isinstance(value, (str, bytes, bytearray)) or not isinstance(value, Sequence):
        raise ValueError(f"{context} must be a JSON array")
    return value


def _require_int(name: str, value: Any) -> int:
    if not isinstance(value, int) or isinstance(value, bool):
        raise ValueError(f"{name} must be an integer")
    return value


def _require_bool(name: str, value: Any) -> bool:
    if not isinstance(value, bool):
        raise ValueError(f"{name} must be a bool")
    return value


def _require_optional_bool(name: str, value: Any) -> bool | None:
    if value is None:
        return None
    return _require_bool(name, value)


def _require_nonnegative_float(name: str, value: float) -> float:
    numeric = float(value)
    if not math.isfinite(numeric) or numeric < 0.0:
        raise ValueError(f"{name} must be a non-negative finite float")
    return numeric


def _require_optional_nonnegative_float(name: str, value: Any) -> float | None:
    if value is None:
        return None
    return _require_nonnegative_float(name, float(value))


def _require_optional_text(name: str, value: Any) -> str | None:
    if value is None:
        return None
    return _require_text(name, str(value))


def _require_exact_keys(
    mapping: Mapping[str, Any],
    *,
    expected: tuple[str, ...],
    context: str,
) -> None:
    keys = tuple(mapping.keys())
    if set(keys) != set(expected):
        missing = tuple(key for key in expected if key not in mapping)
        unexpected = tuple(key for key in keys if key not in expected)
        parts: list[str] = []
        if missing:
            parts.append(f"missing keys {missing!r}")
        if unexpected:
            parts.append(f"unexpected keys {unexpected!r}")
        raise ValueError(f"{context} must contain exactly {expected!r} ({'; '.join(parts)})")


def _complex_from_payload(name: str, payload: Any) -> complex:
    mapping = _require_mapping(payload, context=name)
    _require_exact_keys(mapping, expected=("real", "imag"), context=name)
    real = float(mapping["real"])
    imag = float(mapping["imag"])
    if not math.isfinite(real) or not math.isfinite(imag):
        raise ValueError(f"{name} must contain finite real/imag parts")
    return real + 1j * imag


def _bridge_payload_from_source(source: Any) -> Mapping[str, Any]:
    if isinstance(source, Mapping):
        payload = source
    elif hasattr(source, "as_dict"):
        payload = source.as_dict()
    else:
        raise TypeError("source must be a mapping or provide as_dict()")
    mapping = _require_mapping(payload, context="source")
    if mapping.get("schema_id") != MODERN_POINT_PHENOMENOLOGY_BRIDGE_SCHEMA_ID:
        raise ValueError(
            "source.schema_id must be exactly "
            f"{MODERN_POINT_PHENOMENOLOGY_BRIDGE_SCHEMA_ID!r}"
        )
    return mapping


def _evolve_bridge_wilsons(
    system_match: Mapping[str, Any],
    *,
    mu_had: float = 2.0,
) -> Mapping[str, Any]:
    """Return a copy of *system_match* with Wilson coefficients evolved to *mu_had*.

    The bridge stores Wilson coefficients at the KK matching scale.  This helper
    applies leading-log QCD running via ``evolve_deltaf2_wilsons`` so that the
    phenomenology sidecar uses the same evolved coefficients as the evaluation
    layer.
    """
    c1_vll = _complex_from_payload("system_match.c1_vll", system_match["c1_vll"])
    c1_vrr = _complex_from_payload("system_match.c1_vrr", system_match["c1_vrr"])
    c4_lr = _complex_from_payload("system_match.c4_lr", system_match["c4_lr"])
    c5_lr = _complex_from_payload("system_match.c5_lr", system_match["c5_lr"])
    mu_high = float(system_match["matching_scale_GeV"])

    c1_vll_ev, c1_vrr_ev, c4_lr_ev, c5_lr_ev = evolve_deltaf2_wilsons(
        c1_vll, c1_vrr, c4_lr, c5_lr,
        mu_high=mu_high,
        mu_low=mu_had,
    )
    evolved = dict(system_match)
    evolved["c1_vll"] = {"real": c1_vll_ev.real, "imag": c1_vll_ev.imag}
    evolved["c1_vrr"] = {"real": c1_vrr_ev.real, "imag": c1_vrr_ev.imag}
    evolved["c4_lr"] = {"real": c4_lr_ev.real, "imag": c4_lr_ev.imag}
    evolved["c5_lr"] = {"real": c5_lr_ev.real, "imag": c5_lr_ev.imag}
    return evolved


def _operator_sizes_from_bridge_match(
    system_match: Mapping[str, Any],
    *,
    ll_weight: float,
    rr_weight: float,
    lr1_weight: float,
    lr2_weight: float,
    reference_scale: float,
) -> tuple[dict[str, float], str, float]:
    weighted = {
        "C1_VLL": ll_weight * _complex_from_payload("system_match.c1_vll", system_match["c1_vll"]),
        "C1_VRR": rr_weight * _complex_from_payload("system_match.c1_vrr", system_match["c1_vrr"]),
        "C4_LR": lr1_weight * _complex_from_payload("system_match.c4_lr", system_match["c4_lr"]),
        "C5_LR": lr2_weight * _complex_from_payload("system_match.c5_lr", system_match["c5_lr"]),
    }
    operator_sizes = {
        name: float(reference_scale**2 * abs(value)) for name, value in weighted.items()
    }
    dominant_operator = max(operator_sizes, key=operator_sizes.get)
    return operator_sizes, dominant_operator, float(operator_sizes[dominant_operator])


# --------------------------------------------------------------------------
# Kaon hadronic parameters (inline, no deltaf2 import)
# --------------------------------------------------------------------------

_KAON_F_K = 0.1557
_KAON_M_K = 0.49761
_KAON_DELTA_M_K = 3.484e-15
_KAON_M_S_2GEV = 0.0934
_KAON_M_D_2GEV = 0.00467
_KAON_B_1 = 0.717
_KAON_B_4 = 0.78
_KAON_B_5 = 0.57
_KAON_KAPPA_EPSILON = 0.94
_KAON_EPSILON_K_EXP = 2.228e-3
_KAON_EPSILON_K_SM = 1.81e-3


def _kaon_m12_np_from_bridge_match(
    system_match: Mapping[str, Any],
) -> complex:
    """Compute M_12^NP for kaon system from bridge Wilson coefficients."""
    c1_vll = _complex_from_payload("system_match.c1_vll", system_match["c1_vll"])
    c1_vrr = _complex_from_payload("system_match.c1_vrr", system_match["c1_vrr"])
    c4_lr = _complex_from_payload("system_match.c4_lr", system_match["c4_lr"])
    c5_lr = _complex_from_payload("system_match.c5_lr", system_match["c5_lr"])

    fk2_mk = _KAON_F_K**2 * _KAON_M_K
    m_ratio_sq = (_KAON_M_K / (_KAON_M_S_2GEV + _KAON_M_D_2GEV)) ** 2

    o1_vll = (2.0 / 3.0) * fk2_mk * _KAON_B_1
    o4_lr = (m_ratio_sq * (1.0 / 6.0) + 1.0 / 4.0) * fk2_mk * _KAON_B_4
    o5_lr = (m_ratio_sq * (1.0 / 2.0) + 1.0 / 12.0) * fk2_mk * _KAON_B_5

    return (
        c1_vll * o1_vll
        + c1_vrr * o1_vll  # VRR has same ME as VLL by parity
        + c4_lr * o4_lr
        + c5_lr * o5_lr
    )


def _evaluate_epsilon_k_from_bridge(
    system_match: Mapping[str, Any],
) -> tuple[float, dict[str, float], str, float]:
    """Evaluate epsilon_K^NP from bridge match and return (ratio_to_budget, operator_sizes, dominant, dominant_size)."""
    m12_np = _kaon_m12_np_from_bridge_match(system_match)
    im_m12_np = m12_np.imag
    epsilon_k_np = abs(_KAON_KAPPA_EPSILON / (math.sqrt(2.0) * _KAON_DELTA_M_K) * im_m12_np)
    budget = abs(_KAON_EPSILON_K_EXP - _KAON_EPSILON_K_SM)
    ratio = epsilon_k_np / budget if budget > 0.0 else float("inf")

    # Provide per-operator breakdown using imaginary parts
    c1_vll = _complex_from_payload("system_match.c1_vll", system_match["c1_vll"])
    c1_vrr = _complex_from_payload("system_match.c1_vrr", system_match["c1_vrr"])
    c4_lr = _complex_from_payload("system_match.c4_lr", system_match["c4_lr"])
    c5_lr = _complex_from_payload("system_match.c5_lr", system_match["c5_lr"])

    fk2_mk = _KAON_F_K**2 * _KAON_M_K
    m_ratio_sq = (_KAON_M_K / (_KAON_M_S_2GEV + _KAON_M_D_2GEV)) ** 2
    o1_vll = (2.0 / 3.0) * fk2_mk * _KAON_B_1
    o4_lr = (m_ratio_sq * (1.0 / 6.0) + 1.0 / 4.0) * fk2_mk * _KAON_B_4
    o5_lr = (m_ratio_sq * (1.0 / 2.0) + 1.0 / 12.0) * fk2_mk * _KAON_B_5
    prefactor = abs(_KAON_KAPPA_EPSILON / (math.sqrt(2.0) * _KAON_DELTA_M_K))

    operator_sizes = {
        "C1_VLL": float(prefactor * abs(c1_vll.imag * o1_vll)),
        "C1_VRR": float(prefactor * abs(c1_vrr.imag * o1_vll)),
        "C4_LR": float(prefactor * abs(c4_lr.imag * o4_lr)),
        "C5_LR": float(prefactor * abs(c5_lr.imag * o5_lr)),
    }
    dominant_operator = max(operator_sizes, key=operator_sizes.get)
    dominant_size = float(operator_sizes[dominant_operator])

    return ratio, operator_sizes, dominant_operator, dominant_size


def _evaluate_delta_mk_from_bridge(
    system_match: Mapping[str, Any],
) -> tuple[float, dict[str, float], str, float]:
    """Evaluate Delta m_K NP contribution from bridge match."""
    m12_np = _kaon_m12_np_from_bridge_match(system_match)
    abs_m12_np = abs(m12_np)
    half_dm = _KAON_DELTA_M_K / 2.0
    ratio = abs_m12_np / half_dm if half_dm > 0.0 else float("inf")

    # Per-operator breakdown using absolute values
    c1_vll = _complex_from_payload("system_match.c1_vll", system_match["c1_vll"])
    c1_vrr = _complex_from_payload("system_match.c1_vrr", system_match["c1_vrr"])
    c4_lr = _complex_from_payload("system_match.c4_lr", system_match["c4_lr"])
    c5_lr = _complex_from_payload("system_match.c5_lr", system_match["c5_lr"])

    fk2_mk = _KAON_F_K**2 * _KAON_M_K
    m_ratio_sq = (_KAON_M_K / (_KAON_M_S_2GEV + _KAON_M_D_2GEV)) ** 2
    o1_vll = (2.0 / 3.0) * fk2_mk * _KAON_B_1
    o4_lr = (m_ratio_sq * (1.0 / 6.0) + 1.0 / 4.0) * fk2_mk * _KAON_B_4
    o5_lr = (m_ratio_sq * (1.0 / 2.0) + 1.0 / 12.0) * fk2_mk * _KAON_B_5

    operator_sizes = {
        "C1_VLL": float(abs(c1_vll * o1_vll)),
        "C1_VRR": float(abs(c1_vrr * o1_vll)),
        "C4_LR": float(abs(c4_lr * o4_lr)),
        "C5_LR": float(abs(c5_lr * o5_lr)),
    }
    dominant_operator = max(operator_sizes, key=operator_sizes.get)
    dominant_size = float(operator_sizes[dominant_operator])

    return ratio, operator_sizes, dominant_operator, dominant_size


@dataclass(frozen=True)
class ModernPhenomenologySystemResult:
    """Explicit per-system modern phenomenology result for the QS5 sidecar."""

    schema_id: str = MODERN_POINT_PHENOMENOLOGY_SYSTEM_RESULT_SCHEMA_ID
    system_id: str = ""
    policy_system_id: str = ""
    policy_id: str = ""
    policy_display_name: str = ""
    observable_kind: str = ""
    treatment_id: str = ""
    bridge_system_id: str | None = None
    bridge_observable_id: str | None = None
    bridge_backend_system_id: str | None = None
    bridge_backend_key: str | None = None
    evaluated_from_bridge: bool = False
    included_in_non_cp_acceptance: bool = False
    passes: bool | None = None
    ratio_to_bound: float | None = None
    bound: float | None = None
    dominant_operator: str | None = None
    dominant_operator_size: float | None = None
    weighted_operator_sizes: Mapping[str, float] = field(default_factory=dict)
    note: str = ""

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            _require_exact(
                "schema_id",
                self.schema_id,
                expected=MODERN_POINT_PHENOMENOLOGY_SYSTEM_RESULT_SCHEMA_ID,
            ),
        )
        object.__setattr__(self, "system_id", _require_text("system_id", self.system_id))
        if self.system_id not in MODERN_PHENOMENOLOGY_SYSTEM_IDS:
            raise ValueError(
                f"system_id must be one of {MODERN_PHENOMENOLOGY_SYSTEM_IDS!r}"
            )
        expected_treatment = MODERN_POINT_PHENOMENOLOGY_SYSTEM_TREATMENT_IDS[self.system_id]
        object.__setattr__(
            self,
            "policy_system_id",
            _require_exact(
                "policy_system_id",
                _require_text("policy_system_id", self.policy_system_id),
                expected=self.system_id,
            ),
        )
        object.__setattr__(self, "policy_id", _require_text("policy_id", self.policy_id))
        object.__setattr__(
            self,
            "policy_display_name",
            _require_text("policy_display_name", self.policy_display_name),
        )
        object.__setattr__(
            self,
            "observable_kind",
            _require_text("observable_kind", self.observable_kind),
        )
        object.__setattr__(
            self,
            "treatment_id",
            _require_exact(
                "treatment_id",
                _require_text("treatment_id", self.treatment_id),
                expected=expected_treatment,
            ),
        )
        object.__setattr__(
            self,
            "bridge_system_id",
            _require_optional_text("bridge_system_id", self.bridge_system_id),
        )
        object.__setattr__(
            self,
            "bridge_observable_id",
            _require_optional_text("bridge_observable_id", self.bridge_observable_id),
        )
        object.__setattr__(
            self,
            "bridge_backend_system_id",
            _require_optional_text(
                "bridge_backend_system_id", self.bridge_backend_system_id
            ),
        )
        object.__setattr__(
            self,
            "bridge_backend_key",
            _require_optional_text("bridge_backend_key", self.bridge_backend_key),
        )
        object.__setattr__(
            self,
            "evaluated_from_bridge",
            _require_bool("evaluated_from_bridge", self.evaluated_from_bridge),
        )
        object.__setattr__(
            self,
            "included_in_non_cp_acceptance",
            _require_bool(
                "included_in_non_cp_acceptance",
                self.included_in_non_cp_acceptance,
            ),
        )
        object.__setattr__(self, "passes", _require_optional_bool("passes", self.passes))
        object.__setattr__(
            self,
            "ratio_to_bound",
            _require_optional_nonnegative_float("ratio_to_bound", self.ratio_to_bound),
        )
        object.__setattr__(self, "bound", _require_optional_nonnegative_float("bound", self.bound))
        object.__setattr__(
            self,
            "dominant_operator",
            _require_optional_text("dominant_operator", self.dominant_operator),
        )
        object.__setattr__(
            self,
            "dominant_operator_size",
            _require_optional_nonnegative_float(
                "dominant_operator_size", self.dominant_operator_size
            ),
        )
        if not isinstance(self.weighted_operator_sizes, Mapping):
            raise ValueError("weighted_operator_sizes must be a mapping")
        canonical_operator_sizes = {
            str(key): float(value)
            for key, value in sorted(self.weighted_operator_sizes.items())
        }
        for operator_size in canonical_operator_sizes.values():
            _require_nonnegative_float("weighted_operator_sizes", operator_size)
        object.__setattr__(self, "weighted_operator_sizes", canonical_operator_sizes)
        object.__setattr__(self, "note", _require_text("note", self.note))
        if self.evaluated_from_bridge:
            if self.bridge_system_id is None:
                raise ValueError("evaluated systems must define bridge_system_id")
            if self.bridge_backend_key is None:
                raise ValueError("evaluated systems must define bridge_backend_key")
            if self.passes is None:
                raise ValueError("evaluated systems must define passes")
            if self.ratio_to_bound is None or self.bound is None:
                raise ValueError("evaluated systems must define ratio_to_bound and bound")
            if self.dominant_operator is None or self.dominant_operator_size is None:
                raise ValueError(
                    "evaluated systems must define dominant_operator and dominant_operator_size"
                )
            if not self.weighted_operator_sizes:
                raise ValueError(
                    "evaluated systems must define weighted_operator_sizes"
                )
        else:
            if self.passes is not None or self.ratio_to_bound is not None or self.bound is not None:
                raise ValueError("blocked systems must not define pass/bound data")
            if self.dominant_operator is not None or self.dominant_operator_size is not None:
                raise ValueError(
                    "blocked systems must not define dominant operator data"
                )
            if self.weighted_operator_sizes:
                raise ValueError("blocked systems must not define operator sizes")

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema_id": self.schema_id,
            "system_id": self.system_id,
            "policy_system_id": self.policy_system_id,
            "policy_id": self.policy_id,
            "policy_display_name": self.policy_display_name,
            "observable_kind": self.observable_kind,
            "treatment_id": self.treatment_id,
            "bridge_system_id": self.bridge_system_id,
            "bridge_observable_id": self.bridge_observable_id,
            "bridge_backend_system_id": self.bridge_backend_system_id,
            "bridge_backend_key": self.bridge_backend_key,
            "evaluated_from_bridge": self.evaluated_from_bridge,
            "included_in_non_cp_acceptance": self.included_in_non_cp_acceptance,
            "passes": self.passes,
            "ratio_to_bound": self.ratio_to_bound,
            "bound": self.bound,
            "dominant_operator": self.dominant_operator,
            "dominant_operator_size": self.dominant_operator_size,
            "weighted_operator_sizes": dict(self.weighted_operator_sizes),
            "note": self.note,
        }

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "ModernPhenomenologySystemResult":
        mapping = _require_mapping(payload, context="system_result")
        _require_exact_keys(
            mapping,
            expected=(
                "schema_id",
                "system_id",
                "policy_system_id",
                "policy_id",
                "policy_display_name",
                "observable_kind",
                "treatment_id",
                "bridge_system_id",
                "bridge_observable_id",
                "bridge_backend_system_id",
                "bridge_backend_key",
                "evaluated_from_bridge",
                "included_in_non_cp_acceptance",
                "passes",
                "ratio_to_bound",
                "bound",
                "dominant_operator",
                "dominant_operator_size",
                "weighted_operator_sizes",
                "note",
            ),
            context="system_result",
        )
        return cls(
            schema_id=_require_text("system_result.schema_id", mapping["schema_id"]),
            system_id=_require_text("system_result.system_id", mapping["system_id"]),
            policy_system_id=_require_text(
                "system_result.policy_system_id", mapping["policy_system_id"]
            ),
            policy_id=_require_text("system_result.policy_id", mapping["policy_id"]),
            policy_display_name=_require_text(
                "system_result.policy_display_name", mapping["policy_display_name"]
            ),
            observable_kind=_require_text(
                "system_result.observable_kind", mapping["observable_kind"]
            ),
            treatment_id=_require_text(
                "system_result.treatment_id", mapping["treatment_id"]
            ),
            bridge_system_id=mapping["bridge_system_id"],
            bridge_observable_id=mapping["bridge_observable_id"],
            bridge_backend_system_id=mapping["bridge_backend_system_id"],
            bridge_backend_key=mapping["bridge_backend_key"],
            evaluated_from_bridge=bool(mapping["evaluated_from_bridge"]),
            included_in_non_cp_acceptance=bool(
                mapping["included_in_non_cp_acceptance"]
            ),
            passes=mapping["passes"],
            ratio_to_bound=mapping["ratio_to_bound"],
            bound=mapping["bound"],
            dominant_operator=mapping["dominant_operator"],
            dominant_operator_size=mapping["dominant_operator_size"],
            weighted_operator_sizes=_require_mapping(
                mapping["weighted_operator_sizes"],
                context="system_result.weighted_operator_sizes",
            ),
            note=_require_text("system_result.note", mapping["note"]),
        )


@dataclass(frozen=True)
class ModernPointPhenomenologyArtifactV1:
    """Explicit modern phenomenology sidecar for the first non-CP release."""

    schema_id: str = MODERN_POINT_PHENOMENOLOGY_ARTIFACT_SCHEMA_ID
    schema_version: int = MODERN_POINT_PHENOMENOLOGY_ARTIFACT_SCHEMA_VERSION
    lane_id: str = MODERN_LANE_ID
    release_scope_id: str = MODERN_POINT_PHENOMENOLOGY_RELEASE_SCOPE_ID
    input_bundle_schema_id: str = MODERN_DEFAULT_INPUTS_SCHEMA_ID
    input_bundle_id: str = ""
    input_provenance_id: str = ""
    input_resolution_policy_id: str = ""
    point_id: str = ""
    point_label: str = ""
    policy_schema_id: str = MODERN_PHENOMENOLOGY_SCHEMA_ID
    policy_id: str = MODERN_PHENOMENOLOGY_POLICY_ID
    policy_system_ids: tuple[str, ...] = MODERN_PHENOMENOLOGY_SYSTEM_IDS
    bridge_artifact_schema_id: str = MODERN_POINT_PHENOMENOLOGY_BRIDGE_SCHEMA_ID
    bridge_coupling_schema_id: str = ""
    bridge_matching_schema_id: str = ""
    M_KK: float = 0.0
    xi_KK: float = 0.0
    non_cp_acceptance_system_ids: tuple[str, ...] = (
        MODERN_POINT_PHENOMENOLOGY_NON_CP_ACCEPTANCE_SYSTEM_IDS
    )
    system_results: tuple[ModernPhenomenologySystemResult, ...] = field(
        default_factory=tuple
    )
    non_cp_passes: bool = False
    failing_non_cp_system_ids: tuple[str, ...] = ()
    kaon_viability_claimed: bool = True
    notes: str = (
        "Explicit modern QS5 sidecar for the full Delta F = 2 release. All five "
        "systems (epsilon_K, K, B_d, B_s, D0) are acceptance-bearing. epsilon_K "
        "uses proper hadronic matrix elements for CP violation, and K uses "
        "proper hadronic matrix elements for non-CP kaon mass difference."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            _require_exact(
                "schema_id",
                self.schema_id,
                expected=MODERN_POINT_PHENOMENOLOGY_ARTIFACT_SCHEMA_ID,
            ),
        )
        object.__setattr__(
            self,
            "schema_version",
            _require_int("schema_version", self.schema_version),
        )
        if self.schema_version != MODERN_POINT_PHENOMENOLOGY_ARTIFACT_SCHEMA_VERSION:
            raise ValueError(
                "schema_version must be exactly "
                f"{MODERN_POINT_PHENOMENOLOGY_ARTIFACT_SCHEMA_VERSION}"
            )
        object.__setattr__(
            self,
            "lane_id",
            _require_exact("lane_id", self.lane_id, expected=MODERN_LANE_ID),
        )
        object.__setattr__(
            self,
            "release_scope_id",
            _require_exact(
                "release_scope_id",
                self.release_scope_id,
                expected=MODERN_POINT_PHENOMENOLOGY_RELEASE_SCOPE_ID,
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
        object.__setattr__(self, "input_bundle_id", _require_text("input_bundle_id", self.input_bundle_id))
        object.__setattr__(
            self,
            "input_provenance_id",
            _require_text("input_provenance_id", self.input_provenance_id),
        )
        object.__setattr__(
            self,
            "input_resolution_policy_id",
            _require_text("input_resolution_policy_id", self.input_resolution_policy_id),
        )
        object.__setattr__(self, "point_id", _require_text("point_id", self.point_id))
        object.__setattr__(self, "point_label", _require_text("point_label", self.point_label))
        object.__setattr__(
            self,
            "policy_schema_id",
            _require_exact(
                "policy_schema_id",
                self.policy_schema_id,
                expected=MODERN_PHENOMENOLOGY_SCHEMA_ID,
            ),
        )
        object.__setattr__(
            self,
            "policy_id",
            _require_exact(
                "policy_id",
                self.policy_id,
                expected=MODERN_PHENOMENOLOGY_POLICY_ID,
            ),
        )
        if tuple(self.policy_system_ids) != MODERN_PHENOMENOLOGY_SYSTEM_IDS:
            raise ValueError("policy_system_ids must preserve the frozen modern policy order")
        object.__setattr__(
            self,
            "bridge_artifact_schema_id",
            _require_exact(
                "bridge_artifact_schema_id",
                self.bridge_artifact_schema_id,
                expected=MODERN_POINT_PHENOMENOLOGY_BRIDGE_SCHEMA_ID,
            ),
        )
        object.__setattr__(
            self,
            "bridge_coupling_schema_id",
            _require_text("bridge_coupling_schema_id", self.bridge_coupling_schema_id),
        )
        object.__setattr__(
            self,
            "bridge_matching_schema_id",
            _require_text("bridge_matching_schema_id", self.bridge_matching_schema_id),
        )
        object.__setattr__(self, "M_KK", _require_positive_float("M_KK", self.M_KK))
        object.__setattr__(self, "xi_KK", _require_positive_float("xi_KK", self.xi_KK))
        if tuple(self.non_cp_acceptance_system_ids) != MODERN_POINT_PHENOMENOLOGY_NON_CP_ACCEPTANCE_SYSTEM_IDS:
            raise ValueError(
                "non_cp_acceptance_system_ids must be exactly "
                f"{MODERN_POINT_PHENOMENOLOGY_NON_CP_ACCEPTANCE_SYSTEM_IDS!r}"
            )
        normalized_results = tuple(self.system_results)
        if tuple(result.system_id for result in normalized_results) != MODERN_PHENOMENOLOGY_SYSTEM_IDS:
            raise ValueError(
                "system_results must preserve the frozen epsilon_K/K/B_d/B_s/D0 order"
            )
        object.__setattr__(self, "system_results", normalized_results)
        object.__setattr__(
            self,
            "failing_non_cp_system_ids",
            tuple(str(system_id) for system_id in self.failing_non_cp_system_ids),
        )
        expected_failing = tuple(
            result.system_id
            for result in normalized_results
            if result.included_in_non_cp_acceptance and result.passes is False
        )
        if self.failing_non_cp_system_ids != expected_failing:
            raise ValueError(
                "failing_non_cp_system_ids must match the evaluated acceptance-bearing failures"
            )
        expected_non_cp_passes = not expected_failing
        if self.non_cp_passes != expected_non_cp_passes:
            raise ValueError(
                "non_cp_passes must agree with the acceptance-bearing system results"
            )
        object.__setattr__(
            self,
            "kaon_viability_claimed",
            _require_bool("kaon_viability_claimed", self.kaon_viability_claimed),
        )
        if not self.kaon_viability_claimed:
            raise ValueError("kaon_viability_claimed must be true for the full Delta F = 2 release scope")
        object.__setattr__(self, "notes", _require_text("notes", self.notes))

    @property
    def system_ids(self) -> tuple[str, ...]:
        return tuple(result.system_id for result in self.system_results)

    def system_result(self, system_id: str) -> ModernPhenomenologySystemResult:
        requested = _require_text("system_id", system_id)
        for result in self.system_results:
            if result.system_id == requested:
                return result
        raise ValueError(
            f"system_id must be one of {MODERN_PHENOMENOLOGY_SYSTEM_IDS!r}"
        )

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema_id": self.schema_id,
            "schema_version": self.schema_version,
            "lane_id": self.lane_id,
            "release_scope_id": self.release_scope_id,
            "input_bundle_schema_id": self.input_bundle_schema_id,
            "input_bundle_id": self.input_bundle_id,
            "input_provenance_id": self.input_provenance_id,
            "input_resolution_policy_id": self.input_resolution_policy_id,
            "point_id": self.point_id,
            "point_label": self.point_label,
            "policy_schema_id": self.policy_schema_id,
            "policy_id": self.policy_id,
            "policy_system_ids": list(self.policy_system_ids),
            "bridge_artifact_schema_id": self.bridge_artifact_schema_id,
            "bridge_coupling_schema_id": self.bridge_coupling_schema_id,
            "bridge_matching_schema_id": self.bridge_matching_schema_id,
            "M_KK": self.M_KK,
            "xi_KK": self.xi_KK,
            "non_cp_acceptance_system_ids": list(self.non_cp_acceptance_system_ids),
            "system_results": [result.as_dict() for result in self.system_results],
            "non_cp_passes": self.non_cp_passes,
            "failing_non_cp_system_ids": list(self.failing_non_cp_system_ids),
            "kaon_viability_claimed": self.kaon_viability_claimed,
            "notes": self.notes,
        }

    def to_json(self) -> str:
        return json.dumps(self.as_dict(), indent=2, sort_keys=True, allow_nan=False) + "\n"

    def write_json(self, path: str | Path) -> None:
        Path(path).write_text(self.to_json(), encoding="utf-8")

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "ModernPointPhenomenologyArtifactV1":
        mapping = _require_mapping(payload, context="phenomenology_artifact")
        _require_exact_keys(
            mapping,
            expected=(
                "schema_id",
                "schema_version",
                "lane_id",
                "release_scope_id",
                "input_bundle_schema_id",
                "input_bundle_id",
                "input_provenance_id",
                "input_resolution_policy_id",
                "point_id",
                "point_label",
                "policy_schema_id",
                "policy_id",
                "policy_system_ids",
                "bridge_artifact_schema_id",
                "bridge_coupling_schema_id",
                "bridge_matching_schema_id",
                "M_KK",
                "xi_KK",
                "non_cp_acceptance_system_ids",
                "system_results",
                "non_cp_passes",
                "failing_non_cp_system_ids",
                "kaon_viability_claimed",
                "notes",
            ),
            context="phenomenology_artifact",
        )
        return cls(
            schema_id=_require_text("phenomenology_artifact.schema_id", mapping["schema_id"]),
            schema_version=_require_int(
                "phenomenology_artifact.schema_version",
                mapping["schema_version"],
            ),
            lane_id=_require_text("phenomenology_artifact.lane_id", mapping["lane_id"]),
            release_scope_id=_require_text(
                "phenomenology_artifact.release_scope_id",
                mapping["release_scope_id"],
            ),
            input_bundle_schema_id=_require_text(
                "phenomenology_artifact.input_bundle_schema_id",
                mapping["input_bundle_schema_id"],
            ),
            input_bundle_id=_require_text(
                "phenomenology_artifact.input_bundle_id",
                mapping["input_bundle_id"],
            ),
            input_provenance_id=_require_text(
                "phenomenology_artifact.input_provenance_id",
                mapping["input_provenance_id"],
            ),
            input_resolution_policy_id=_require_text(
                "phenomenology_artifact.input_resolution_policy_id",
                mapping["input_resolution_policy_id"],
            ),
            point_id=_require_text("phenomenology_artifact.point_id", mapping["point_id"]),
            point_label=_require_text(
                "phenomenology_artifact.point_label", mapping["point_label"]
            ),
            policy_schema_id=_require_text(
                "phenomenology_artifact.policy_schema_id",
                mapping["policy_schema_id"],
            ),
            policy_id=_require_text(
                "phenomenology_artifact.policy_id", mapping["policy_id"]
            ),
            policy_system_ids=tuple(
                _require_text("phenomenology_artifact.policy_system_ids", item)
                for item in _require_sequence(
                    mapping["policy_system_ids"],
                    context="phenomenology_artifact.policy_system_ids",
                )
            ),
            bridge_artifact_schema_id=_require_text(
                "phenomenology_artifact.bridge_artifact_schema_id",
                mapping["bridge_artifact_schema_id"],
            ),
            bridge_coupling_schema_id=_require_text(
                "phenomenology_artifact.bridge_coupling_schema_id",
                mapping["bridge_coupling_schema_id"],
            ),
            bridge_matching_schema_id=_require_text(
                "phenomenology_artifact.bridge_matching_schema_id",
                mapping["bridge_matching_schema_id"],
            ),
            M_KK=float(mapping["M_KK"]),
            xi_KK=float(mapping["xi_KK"]),
            non_cp_acceptance_system_ids=tuple(
                _require_text(
                    "phenomenology_artifact.non_cp_acceptance_system_ids",
                    item,
                )
                for item in _require_sequence(
                    mapping["non_cp_acceptance_system_ids"],
                    context="phenomenology_artifact.non_cp_acceptance_system_ids",
                )
            ),
            system_results=tuple(
                ModernPhenomenologySystemResult.from_dict(item)
                for item in _require_sequence(
                    mapping["system_results"],
                    context="phenomenology_artifact.system_results",
                )
            ),
            non_cp_passes=bool(mapping["non_cp_passes"]),
            failing_non_cp_system_ids=tuple(
                _require_text("phenomenology_artifact.failing_non_cp_system_ids", item)
                for item in _require_sequence(
                    mapping["failing_non_cp_system_ids"],
                    context="phenomenology_artifact.failing_non_cp_system_ids",
                )
            ),
            kaon_viability_claimed=bool(mapping["kaon_viability_claimed"]),
            notes=_require_text("phenomenology_artifact.notes", mapping["notes"]),
        )

    @classmethod
    def from_source(
        cls,
        source: Any,
        *,
        inputs: ModernDefaultInputs | None = None,
        policy: ModernPhenomenologyPolicy | None = None,
    ) -> "ModernPointPhenomenologyArtifactV1":
        bridge = _bridge_payload_from_source(source)
        resolved_inputs = default_modern_default_inputs() if inputs is None else inputs
        resolved_policy = (
            default_modern_phenomenology_policy() if policy is None else policy
        )
        if bridge["input_bundle_schema_id"] != resolved_inputs.schema_id:
            raise ValueError("bridge input_bundle_schema_id must match the supplied modern inputs")
        if bridge["input_bundle_id"] != resolved_inputs.bundle_id:
            raise ValueError("bridge input_bundle_id must match the supplied modern inputs")
        if bridge["input_provenance_id"] != resolved_inputs.provenance_id:
            raise ValueError(
                "bridge input_provenance_id must match the supplied modern inputs"
            )
        matching = _require_mapping(bridge["matching"], context="source.matching")
        system_matches = {
            str(item["system_id"]): _require_mapping(item, context="source.matching.system_matches")
            for item in _require_sequence(
                matching["system_matches"],
                context="source.matching.system_matches",
            )
        }
        required_matching_systems = ("K", "B_d", "B_s", "D0")
        if tuple(system_matches) != required_matching_systems:
            raise ValueError(
                "bridge matching.system_matches must preserve the frozen K/B_d/B_s/D0 order"
            )
        bundle_inputs = {item.system_id: item for item in resolved_inputs.neutral_meson_inputs}
        weights = resolved_inputs.operator_weight_policy
        results: list[ModernPhenomenologySystemResult] = []
        for system_id in MODERN_PHENOMENOLOGY_SYSTEM_IDS:
            policy_system = resolved_policy.system_policy(system_id)
            treatment_id = MODERN_POINT_PHENOMENOLOGY_SYSTEM_TREATMENT_IDS[system_id]
            bridge_system_id = "K" if system_id in ("epsilon_K", "K") else system_id
            system_match = system_matches[bridge_system_id]
            # Evolve Wilson coefficients from matching scale to hadronic
            # scale so the sidecar is consistent with the evaluation layer
            # (which applies QCD running by default).
            evolved_match = _evolve_bridge_wilsons(system_match)
            input_item = bundle_inputs[system_id]

            if system_id == "epsilon_K":
                # CP-violating epsilon_K uses proper hadronic matrix elements
                ratio_to_bound, operator_sizes, dominant_operator, dominant_size = (
                    _evaluate_epsilon_k_from_bridge(evolved_match)
                )
                note = (
                    f"{input_item.note} Evaluated with proper hadronic matrix elements "
                    "for CP-violating kaon mixing (QCD-evolved). Included in acceptance."
                )
            elif system_id == "K":
                # Non-CP Delta m_K uses proper hadronic matrix elements
                ratio_to_bound, operator_sizes, dominant_operator, dominant_size = (
                    _evaluate_delta_mk_from_bridge(evolved_match)
                )
                note = (
                    f"{input_item.note} Evaluated with proper hadronic matrix elements "
                    "for non-CP kaon mass difference (QCD-evolved). Included in acceptance."
                )
            else:
                operator_sizes, dominant_operator, dominant_size = _operator_sizes_from_bridge_match(
                    evolved_match,
                    ll_weight=weights.ll_weight,
                    rr_weight=weights.rr_weight,
                    lr1_weight=weights.lr1_weight,
                    lr2_weight=weights.lr2_weight,
                    reference_scale=weights.reference_scale_GeV,
                )
                ratio_to_bound = float(dominant_size / input_item.bound)
                if system_id == "D0":
                    note = (
                        f"{input_item.note} Included in acceptance with conservative D0 treatment (QCD-evolved)."
                    )
                else:
                    note = f"{input_item.note} Included in acceptance (QCD-evolved)."

            results.append(
                ModernPhenomenologySystemResult(
                    system_id=system_id,
                    policy_system_id=policy_system.system_id,
                    policy_id=policy_system.policy_id,
                    policy_display_name=policy_system.display_name,
                    observable_kind=input_item.observable_kind,
                    treatment_id=treatment_id,
                    bridge_system_id=bridge_system_id,
                    bridge_observable_id=str(system_match["observable_id"]),
                    bridge_backend_system_id=str(system_match["backend_system_id"]),
                    bridge_backend_key=str(system_match["backend_key"]),
                    evaluated_from_bridge=True,
                    included_in_non_cp_acceptance=(
                        system_id in MODERN_POINT_PHENOMENOLOGY_NON_CP_ACCEPTANCE_SYSTEM_IDS
                    ),
                    passes=ratio_to_bound <= 1.0,
                    ratio_to_bound=ratio_to_bound,
                    bound=float(input_item.bound),
                    dominant_operator=dominant_operator,
                    dominant_operator_size=dominant_size,
                    weighted_operator_sizes=operator_sizes,
                    note=note,
                )
            )
        failing_non_cp = tuple(
            result.system_id
            for result in results
            if result.included_in_non_cp_acceptance and result.passes is False
        )
        couplings = _require_mapping(bridge["couplings"], context="source.couplings")
        return cls(
            input_bundle_schema_id=_require_text(
                "source.input_bundle_schema_id",
                bridge["input_bundle_schema_id"],
            ),
            input_bundle_id=_require_text("source.input_bundle_id", bridge["input_bundle_id"]),
            input_provenance_id=_require_text(
                "source.input_provenance_id", bridge["input_provenance_id"]
            ),
            input_resolution_policy_id=_require_text(
                "source.input_resolution_policy_id",
                bridge["input_resolution_policy_id"],
            ),
            point_id=_require_text("source.point_id", bridge["point_id"]),
            point_label=_require_text("source.point_label", bridge["point_label"]),
            bridge_coupling_schema_id=_require_text(
                "source.coupling_schema_id",
                bridge["coupling_schema_id"],
            ),
            bridge_matching_schema_id=_require_text(
                "source.matching_schema_id",
                bridge["matching_schema_id"],
            ),
            M_KK=float(couplings["M_KK"]),
            xi_KK=float(couplings["xi_KK"]),
            system_results=tuple(results),
            non_cp_passes=not failing_non_cp,
            failing_non_cp_system_ids=failing_non_cp,
        )


def build_modern_point_phenomenology_artifact(
    source: Any,
    *,
    inputs: ModernDefaultInputs | None = None,
    policy: ModernPhenomenologyPolicy | None = None,
) -> ModernPointPhenomenologyArtifactV1:
    """Build the explicit modern QS5 sidecar from the exported bridge artifact."""

    return ModernPointPhenomenologyArtifactV1.from_source(
        source,
        inputs=inputs,
        policy=policy,
    )


def read_modern_point_phenomenology_artifact(
    path: str | Path,
) -> ModernPointPhenomenologyArtifactV1:
    """Read one modern phenomenology sidecar from disk."""

    return ModernPointPhenomenologyArtifactV1.from_dict(
        json.loads(Path(path).read_text(encoding="utf-8"))
    )


def write_modern_point_phenomenology_artifact(
    artifact: ModernPointPhenomenologyArtifactV1,
    path: str | Path,
) -> Path:
    """Write one modern phenomenology sidecar to disk."""

    destination = Path(path)
    destination.write_text(artifact.to_json(), encoding="utf-8")
    return destination
