"""Z-pole pseudo-observables from effective ``Z f f`` couplings.

Physics convention
------------------
The effective couplings are normalized as

    L_Z = g_Z Z_mu fbar gamma^mu (g_L P_L + g_R P_R) f.

For a fermion species ``f``,

    A_f = (|g_L|^2 - |g_R|^2) / (|g_L|^2 + |g_R|^2),

and the quark ratio is evaluated from partial-width weights

    R_q = Gamma_q / sum_{u,d,s,c,b} Gamma_i,
    Gamma_i proportional to N_c * radiator_i * (|g_L^i|^2 + |g_R^i|^2).

This is a pseudo-observable layer, not a full electroweak two-loop Z-width
calculation.  The caller may provide effective couplings and radiator factors
already carrying SM corrections.

RS matching status
------------------
NEEDS-HUMAN-PHYSICS: the current ``ParameterPoint`` carries quark
``KK-gluon``-style mass-basis overlap matrices, but not the electroweak
KK/Z/Z' spectrum, custodial representations, brane kinetic terms, or fermion
embedding data needed for rigorous RS shifts in ``Z b_L b_L`` and
``Z b_R b_R``.  The ``zbb_coupling_shift_proxy`` helper is therefore only a
documented leading proxy for the classic RS ``Z -> b bbar`` tension: it maps
available bottom-vs-light overlap non-universality onto
``delta g_b ~ (m_Z / M_KK)^2 Delta overlap``.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
import math
from typing import Mapping

import numpy as np

from .couplings import QuarkMassBasisCouplings

ZPOLE_MODEL_V1 = "zpole_effective_couplings_width_ratios_v1"
ZPOLE_INPUT_BUNDLE_V1 = "zpole_tree_couplings_with_configurable_radiators_v1"
ZPOLE_RS_ZBB_PROXY_V1 = (
    "NEEDS-HUMAN-PHYSICS: available quark KK-gluon overlap matrices are used "
    "as bottom-vs-light neutral-current non-universality proxies; the shift "
    "delta g_b = proxy_strength * (m_Z/M_KK)^2 * (overlap_b - overlap_light) "
    "stands in for missing model-dependent RS EW KK/Z/Z' and custodial "
    "Zbb matching."
)

_QUANTUM_NUMBERS: Mapping[str, tuple[float, float, int]] = {
    "u": (0.5, 2.0 / 3.0, 3),
    "c": (0.5, 2.0 / 3.0, 3),
    "t": (0.5, 2.0 / 3.0, 3),
    "d": (-0.5, -1.0 / 3.0, 3),
    "s": (-0.5, -1.0 / 3.0, 3),
    "b": (-0.5, -1.0 / 3.0, 3),
    "e": (-0.5, -1.0, 1),
    "mu": (-0.5, -1.0, 1),
    "tau": (-0.5, -1.0, 1),
    "nu": (0.5, 0.0, 1),
}
_HADRONIC_FLAVORS: tuple[str, ...] = ("u", "d", "s", "c", "b")


@dataclass(frozen=True)
class ZPoleSMInputs:
    """Inputs for Z-pole pseudo-observable arithmetic."""

    input_bundle: str = ZPOLE_INPUT_BUNDLE_V1
    sin2_theta_eff: float = 0.2315
    m_z_gev: float = 91.1876
    proxy_strength: float = 1.0
    quark_radiators: Mapping[str, float] = field(
        default_factory=lambda: {
            "u": 1.0,
            "d": 1.0,
            "s": 1.0,
            "c": 1.0,
            "b": 1.0,
        }
    )
    constants_citation: str = (
        "Effective-coupling Z-pole pseudo-observables; sin2_theta_eff default "
        "is a LEP/SLC-scale effective weak mixing input.  Per-flavor "
        "radiators may be supplied by the catalog constraint to align the "
        "SM-limit width ratio with a quoted Z-pole SM fit."
    )

    def __post_init__(self) -> None:
        if not 0.0 < float(self.sin2_theta_eff) < 1.0:
            raise ValueError("sin2_theta_eff must lie between zero and one")
        for name in ("m_z_gev", "proxy_strength"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        for flavor in _HADRONIC_FLAVORS:
            radiator = float(self.quark_radiators.get(flavor, 1.0))
            if not math.isfinite(radiator) or radiator <= 0.0:
                raise ValueError(f"quark radiator for {flavor!r} must be positive")

    def radiator_for(self, flavor: str) -> float:
        """Return the positive width radiator for a quark flavor."""
        return float(self.quark_radiators.get(flavor, 1.0))


@dataclass(frozen=True)
class ZPoleCouplings:
    """Effective chiral Z couplings for one fermion species."""

    flavor: str
    g_left: complex
    g_right: complex

    @property
    def norm(self) -> float:
        return _coupling_norm(self)


@dataclass(frozen=True)
class ZPoleQuarkObservables:
    """Pseudo-observables for one quark flavor."""

    model_label: str
    input_bundle: str
    target_quark: str
    r_q: float
    a_q: float
    a_fb: float
    a_e: float
    couplings: ZPoleCouplings
    width_weights: Mapping[str, float]
    diagnostics: Mapping[str, float | complex | str]


@dataclass(frozen=True)
class ZbbCouplingShiftProxy:
    """Documented RS proxy for bottom Z-coupling shifts."""

    model_label: str
    matching_assumption: str
    M_KK: float
    matching_scale: float
    left_bottom_overlap: float
    right_bottom_overlap: float
    left_light_overlap: float
    right_light_overlap: float
    left_nonuniversality: float
    right_nonuniversality: float
    scale_factor: float
    delta_g_left_b: float
    delta_g_right_b: float
    diagnostics: Mapping[str, float | str]


def default_sm_inputs() -> ZPoleSMInputs:
    """Return the default Z-pole pseudo-observable input bundle."""
    return ZPoleSMInputs()


def inputs_with_bottom_radiator(
    target_r_b: float,
    inputs: ZPoleSMInputs | None = None,
) -> ZPoleSMInputs:
    """Return inputs whose bottom radiator reproduces ``target_r_b`` in the SM.

    This calibrates only the inclusive bottom width weight.  It is useful when
    a catalog entry quotes a reference SM ``R_b`` from a full electroweak fit,
    while this module keeps the reusable arithmetic at the pseudo-observable
    level.
    """
    p = default_sm_inputs() if inputs is None else inputs
    target = _bounded_probability("target_r_b", target_r_b)
    base_radiators = {flavor: p.radiator_for(flavor) for flavor in _HADRONIC_FLAVORS}
    base_radiators["b"] = 1.0
    sm_weights = {
        flavor: partial_width_weight(
            sm_couplings(flavor, p),
            radiator=base_radiators[flavor],
        )
        for flavor in _HADRONIC_FLAVORS
    }
    bottom_unit = sm_weights["b"]
    other_sum = sum(weight for flavor, weight in sm_weights.items() if flavor != "b")
    bottom_radiator = target * other_sum / (bottom_unit * (1.0 - target))
    if not math.isfinite(bottom_radiator) or bottom_radiator <= 0.0:
        raise ValueError("calibrated bottom radiator is not positive and finite")
    base_radiators["b"] = float(bottom_radiator)
    return replace(p, quark_radiators=base_radiators)


def sm_couplings(
    flavor: str,
    inputs: ZPoleSMInputs | None = None,
) -> ZPoleCouplings:
    """Return tree-level SM effective chiral couplings for ``flavor``."""
    p = default_sm_inputs() if inputs is None else inputs
    t3, charge, _ = _quantum_numbers(flavor)
    g_left = t3 - charge * p.sin2_theta_eff
    g_right = -charge * p.sin2_theta_eff
    return ZPoleCouplings(
        flavor=flavor,
        g_left=complex(g_left),
        g_right=complex(g_right),
    )


def shifted_couplings(
    base: ZPoleCouplings,
    *,
    delta_g_left: complex = 0.0j,
    delta_g_right: complex = 0.0j,
) -> ZPoleCouplings:
    """Return ``base`` shifted by additive chiral coupling corrections."""
    return ZPoleCouplings(
        flavor=base.flavor,
        g_left=complex(base.g_left + delta_g_left),
        g_right=complex(base.g_right + delta_g_right),
    )


def asymmetry_parameter(couplings: ZPoleCouplings) -> float:
    """Return ``A_f`` from effective chiral couplings."""
    denom = _coupling_norm(couplings)
    if denom <= 0.0:
        raise ValueError("coupling norm must be positive")
    left = abs(complex(couplings.g_left)) ** 2
    right = abs(complex(couplings.g_right)) ** 2
    return float((left - right) / denom)


def forward_backward_asymmetry(
    final_state: ZPoleCouplings,
    *,
    initial_state: ZPoleCouplings | None = None,
) -> float:
    """Return ``A_FB^0,f = 3/4 A_e A_f``."""
    electron = sm_couplings("e") if initial_state is None else initial_state
    return float(0.75 * asymmetry_parameter(electron) * asymmetry_parameter(final_state))


def partial_width_weight(
    couplings: ZPoleCouplings,
    *,
    radiator: float = 1.0,
    n_color: int | None = None,
) -> float:
    """Return a relative partial-width weight for one fermion species."""
    _, _, default_color = _quantum_numbers(couplings.flavor)
    color = default_color if n_color is None else int(n_color)
    if color <= 0:
        raise ValueError("n_color must be positive")
    rad = _positive_float("radiator", radiator)
    return float(color * rad * _coupling_norm(couplings))


def hadronic_width_weights(
    couplings_by_flavor: Mapping[str, ZPoleCouplings] | None = None,
    *,
    inputs: ZPoleSMInputs | None = None,
) -> dict[str, float]:
    """Return ``u,d,s,c,b`` partial-width weights."""
    p = default_sm_inputs() if inputs is None else inputs
    supplied = {} if couplings_by_flavor is None else dict(couplings_by_flavor)
    weights: dict[str, float] = {}
    for flavor in _HADRONIC_FLAVORS:
        couplings = supplied.get(flavor, sm_couplings(flavor, p))
        weights[flavor] = partial_width_weight(
            couplings,
            radiator=p.radiator_for(flavor),
        )
    return weights


def r_quark(
    target_quark: str,
    couplings_by_flavor: Mapping[str, ZPoleCouplings] | None = None,
    *,
    inputs: ZPoleSMInputs | None = None,
) -> float:
    """Return ``R_q = Gamma_q / Gamma_had`` for a light quark flavor."""
    flavor = _canonical_flavor(target_quark)
    if flavor not in _HADRONIC_FLAVORS:
        raise ValueError(f"{target_quark!r} is not a supported hadronic flavor")
    weights = hadronic_width_weights(couplings_by_flavor, inputs=inputs)
    total = sum(weights.values())
    if total <= 0.0:
        raise ValueError("total hadronic width weight must be positive")
    return float(weights[flavor] / total)


def evaluate_quark_pseudo_observables(
    target_quark: str,
    couplings_by_flavor: Mapping[str, ZPoleCouplings] | None = None,
    *,
    inputs: ZPoleSMInputs | None = None,
) -> ZPoleQuarkObservables:
    """Evaluate ``R_q``, ``A_q``, and ``A_FB^0,q`` for ``target_quark``."""
    p = default_sm_inputs() if inputs is None else inputs
    flavor = _canonical_flavor(target_quark)
    supplied = {} if couplings_by_flavor is None else dict(couplings_by_flavor)
    couplings = supplied.get(flavor, sm_couplings(flavor, p))
    weights = hadronic_width_weights(supplied, inputs=p)
    total = sum(weights.values())
    if total <= 0.0:
        raise ValueError("total hadronic width weight must be positive")
    a_q = asymmetry_parameter(couplings)
    electron = sm_couplings("e", p)
    a_e = asymmetry_parameter(electron)
    return ZPoleQuarkObservables(
        model_label=ZPOLE_MODEL_V1,
        input_bundle=p.input_bundle,
        target_quark=flavor,
        r_q=float(weights[flavor] / total),
        a_q=float(a_q),
        a_fb=float(0.75 * a_e * a_q),
        a_e=float(a_e),
        couplings=couplings,
        width_weights={key: float(value) for key, value in weights.items()},
        diagnostics={
            "sin2_theta_eff": float(p.sin2_theta_eff),
            "m_z_gev": float(p.m_z_gev),
            "target_width_weight": float(weights[flavor]),
            "total_hadronic_width_weight": float(total),
            "g_left": complex(couplings.g_left),
            "g_right": complex(couplings.g_right),
            "a_e": float(a_e),
        },
    )


def zbb_coupling_shift_proxy(
    source: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: ZPoleSMInputs | None = None,
) -> ZbbCouplingShiftProxy:
    """Map available mass-basis overlaps onto a documented ``Zbb`` proxy."""
    p = default_sm_inputs() if inputs is None else inputs
    resolved_m_kk = _positive_float(
        "m_kk_gev",
        getattr(source, "M_KK") if m_kk_gev is None else m_kk_gev,
    )
    left_overlap = _overlap_matrix(source, "left_overlap", "left_down")
    right_overlap = _overlap_matrix(source, "right_down_overlap", "right_down")
    left_b = _real_diagonal(left_overlap, 2, "left_bottom_overlap")
    right_b = _real_diagonal(right_overlap, 2, "right_bottom_overlap")
    left_light = 0.5 * (
        _real_diagonal(left_overlap, 0, "left_down_overlap")
        + _real_diagonal(left_overlap, 1, "strange_left_overlap")
    )
    right_light = 0.5 * (
        _real_diagonal(right_overlap, 0, "right_down_overlap")
        + _real_diagonal(right_overlap, 1, "right_strange_overlap")
    )
    left_nonuniversality = float(left_b - left_light)
    right_nonuniversality = float(right_b - right_light)
    scale = float((p.m_z_gev / resolved_m_kk) ** 2)
    delta_left = float(p.proxy_strength * scale * left_nonuniversality)
    delta_right = float(p.proxy_strength * scale * right_nonuniversality)
    return ZbbCouplingShiftProxy(
        model_label=ZPOLE_MODEL_V1,
        matching_assumption=ZPOLE_RS_ZBB_PROXY_V1,
        M_KK=resolved_m_kk,
        matching_scale=resolved_m_kk,
        left_bottom_overlap=float(left_b),
        right_bottom_overlap=float(right_b),
        left_light_overlap=float(left_light),
        right_light_overlap=float(right_light),
        left_nonuniversality=left_nonuniversality,
        right_nonuniversality=right_nonuniversality,
        scale_factor=scale,
        delta_g_left_b=delta_left,
        delta_g_right_b=delta_right,
        diagnostics={
            "m_kk_gev": float(resolved_m_kk),
            "matching_scale_gev": float(resolved_m_kk),
            "m_z_gev": float(p.m_z_gev),
            "proxy_strength": float(p.proxy_strength),
            "scale_factor": scale,
            "left_bottom_overlap": float(left_b),
            "right_bottom_overlap": float(right_b),
            "left_light_overlap": float(left_light),
            "right_light_overlap": float(right_light),
            "left_nonuniversality": left_nonuniversality,
            "right_nonuniversality": right_nonuniversality,
            "delta_g_left_b": delta_left,
            "delta_g_right_b": delta_right,
            "matching_assumption": ZPOLE_RS_ZBB_PROXY_V1,
        },
    )


def evaluate_zbb_with_proxy(
    source: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: ZPoleSMInputs | None = None,
) -> tuple[ZPoleQuarkObservables, ZbbCouplingShiftProxy]:
    """Evaluate ``Z -> b bbar`` pseudo-observables with the RS proxy shift."""
    p = default_sm_inputs() if inputs is None else inputs
    proxy = zbb_coupling_shift_proxy(source, m_kk_gev=m_kk_gev, inputs=p)
    shifted_bottom = shifted_couplings(
        sm_couplings("b", p),
        delta_g_left=proxy.delta_g_left_b,
        delta_g_right=proxy.delta_g_right_b,
    )
    observables = evaluate_quark_pseudo_observables(
        "b",
        {"b": shifted_bottom},
        inputs=p,
    )
    return observables, proxy


def _canonical_flavor(flavor: str) -> str:
    name = str(flavor)
    aliases = {
        "electron": "e",
        "muon": "mu",
        "bottom": "b",
        "beauty": "b",
        "charm": "c",
    }
    return aliases.get(name, name)


def _quantum_numbers(flavor: str) -> tuple[float, float, int]:
    name = _canonical_flavor(flavor)
    try:
        return _QUANTUM_NUMBERS[name]
    except KeyError as exc:
        raise ValueError(f"unsupported Z-pole fermion flavor {flavor!r}") from exc


def _coupling_norm(couplings: ZPoleCouplings) -> float:
    left = abs(complex(couplings.g_left)) ** 2
    right = abs(complex(couplings.g_right)) ** 2
    norm = float(left + right)
    if not math.isfinite(norm):
        raise ValueError("coupling norm must be finite")
    return norm


def _bounded_probability(name: str, value: float) -> float:
    number = float(value)
    if not math.isfinite(number) or not 0.0 < number < 1.0:
        raise ValueError(f"{name} must lie between zero and one")
    return number


def _positive_float(name: str, value: object) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return number


def _matrix(source: object, matrix_name: str) -> np.ndarray:
    matrix = np.asarray(getattr(source, matrix_name), dtype=np.complex128)
    if matrix.shape != (3, 3):
        raise ValueError(f"{matrix_name} must have shape (3, 3)")
    if not np.all(np.isfinite(matrix.real)) or not np.all(np.isfinite(matrix.imag)):
        raise ValueError(f"{matrix_name} entries must be finite")
    return matrix


def _overlap_matrix(
    source: QuarkMassBasisCouplings,
    overlap_name: str,
    coupling_name: str,
) -> np.ndarray:
    overlap = getattr(source, overlap_name, None)
    if overlap is not None:
        return _matrix(source, overlap_name)
    coupling = _matrix(source, coupling_name)
    g_s = _positive_float("g_s", getattr(source, "g_s", 1.0))
    return coupling / g_s


def _real_diagonal(matrix: np.ndarray, index: int, name: str) -> float:
    value = complex(matrix[index, index])
    if abs(value.imag) > 1.0e-12:
        raise ValueError(f"{name} diagonal overlap must be real")
    if not math.isfinite(value.real):
        raise ValueError(f"{name} diagonal overlap must be finite")
    return float(value.real)


ZPOLE_DOWN_FCNC_PROXY_V1 = (
    "NEEDS-HUMAN-PHYSICS: off-diagonal down-sector Z couplings are not "
    "available on ParameterPoint; available quark mass-basis overlap matrices "
    "are mapped to delta g_ij = proxy_strength * (m_Z/M_KK)^2 * overlap_ij, "
    "standing in for missing model-dependent RS EW KK/Z/Z' down-sector "
    "neutral-current matching."
)

_DOWN_FCNC_FLAVOR_INDEX: Mapping[str, int] = {"d": 0, "s": 1, "b": 2}
_DOWN_FCNC_CHARGED_LEPTONS: tuple[str, ...] = ("e", "mu", "tau")
_DOWN_FCNC_NEUTRINO_FLAVORS: tuple[str, ...] = ("nu_e", "nu_mu", "nu_tau")
_DOWN_FCNC_DEFAULT_CHARGE_STATE_FACTOR = 2.0


@dataclass(frozen=True)
class ZPoleDownFCNCCouplingProxy:
    """Documented proxy for an off-diagonal down-sector Z coupling."""

    model_label: str
    matching_assumption: str
    flavor_i: str
    flavor_j: str
    M_KK: float
    matching_scale: float
    left_overlap_ij: complex
    right_overlap_ij: complex
    scale_factor: float
    delta_g_left_ij: complex
    delta_g_right_ij: complex
    diagnostics: Mapping[str, object]


@dataclass(frozen=True)
class ZPoleDownFCNCBranchingResult:
    """Total-width-normalized result for ``Z -> q_i qbar_j + q_j qbar_i``."""

    model_label: str
    input_bundle: str
    flavor_i: str
    flavor_j: str
    branching_fraction: float
    sm_branching_fraction: float
    ratio_to_limit: float | None
    br_limit: float | None
    passes: bool | None
    delta_g_left: complex
    delta_g_right: complex
    coupling_norm: float
    fcnc_width_weight: float
    sm_hadronic_width_weight: float
    sm_total_width_weight: float
    total_hadronic_width_weight: float
    total_width_weight: float
    sm_hadronic_to_total_width_ratio: float
    charge_state_factor: float
    n_color: int
    radiator: float
    diagnostics: Mapping[str, object]


def down_fcnc_sm_hadronic_width_weight(
    inputs: ZPoleSMInputs | None = None,
) -> dict[str, float]:
    """Return SM hadronic Z-width weights from the shared Z-pole convention."""

    p = default_sm_inputs() if inputs is None else inputs
    return hadronic_width_weights(inputs=p)


def down_fcnc_sm_total_width_weight(
    inputs: ZPoleSMInputs | None = None,
) -> dict[str, float]:
    """Return SM total Z-width weights for quarks, charged leptons, and neutrinos."""

    p = default_sm_inputs() if inputs is None else inputs
    weights: dict[str, float] = {}
    weights.update(down_fcnc_sm_hadronic_width_weight(p))
    for flavor in _DOWN_FCNC_CHARGED_LEPTONS:
        weights[flavor] = partial_width_weight(sm_couplings(flavor, p))
    for flavor in _DOWN_FCNC_NEUTRINO_FLAVORS:
        weights[flavor] = partial_width_weight(sm_couplings("nu", p))
    return weights


def down_fcnc_branching_fraction_from_couplings(
    *,
    flavor_i: str,
    flavor_j: str,
    delta_g_left: complex,
    delta_g_right: complex,
    br_limit: float | None = None,
    inputs: ZPoleSMInputs | None = None,
    charge_state_factor: float = _DOWN_FCNC_DEFAULT_CHARGE_STATE_FACTOR,
) -> ZPoleDownFCNCBranchingResult:
    """Return total-Z-width ``BR(Z -> q_i qbar_j + q_j qbar_i)``."""

    p = default_sm_inputs() if inputs is None else inputs
    first = _canonical_down_fcnc_flavor(flavor_i)
    second = _canonical_down_fcnc_flavor(flavor_j)
    if first == second:
        raise ValueError("down-sector FCNC Z decay requires distinct flavors")

    charge_factor = _positive_float("charge_state_factor", charge_state_factor)
    left = complex(delta_g_left)
    right = complex(delta_g_right)
    coupling_norm = float(abs(left) ** 2 + abs(right) ** 2)
    if not math.isfinite(coupling_norm):
        raise ValueError("down-sector FCNC coupling norm must be finite")

    sm_hadronic_weights = down_fcnc_sm_hadronic_width_weight(p)
    sm_total_weights = down_fcnc_sm_total_width_weight(p)
    sm_hadronic = float(sum(sm_hadronic_weights.values()))
    sm_total = float(sum(sm_total_weights.values()))
    if sm_hadronic <= 0.0:
        raise ValueError("SM hadronic Z-width weight must be positive")
    if sm_total <= 0.0:
        raise ValueError("SM total Z-width weight must be positive")

    radiator = _down_fcnc_pair_radiator(first, second, p)
    fcnc_width = float(
        charge_factor
        * partial_width_weight(
            ZPoleCouplings(flavor=first, g_left=left, g_right=right),
            radiator=radiator,
            n_color=3,
        )
    )
    total_hadronic = float(sm_hadronic + fcnc_width)
    total_width = float(sm_total + fcnc_width)
    branching_fraction = float(fcnc_width / total_width)
    hadronic_to_total = float(sm_hadronic / sm_total)

    limit = None if br_limit is None else _bounded_probability("br_limit", br_limit)
    ratio = None if limit is None else float(branching_fraction / limit)
    passes = None if ratio is None else bool(ratio <= 1.0)

    return ZPoleDownFCNCBranchingResult(
        model_label=ZPOLE_MODEL_V1,
        input_bundle=p.input_bundle,
        flavor_i=first,
        flavor_j=second,
        branching_fraction=branching_fraction,
        sm_branching_fraction=0.0,
        ratio_to_limit=ratio,
        br_limit=limit,
        passes=passes,
        delta_g_left=left,
        delta_g_right=right,
        coupling_norm=coupling_norm,
        fcnc_width_weight=fcnc_width,
        sm_hadronic_width_weight=sm_hadronic,
        sm_total_width_weight=sm_total,
        total_hadronic_width_weight=total_hadronic,
        total_width_weight=total_width,
        sm_hadronic_to_total_width_ratio=hadronic_to_total,
        charge_state_factor=charge_factor,
        n_color=3,
        radiator=radiator,
        diagnostics={
            "base_model_label": ZPOLE_MODEL_V1,
            "branching_formula": (
                "BR(Z -> q_i qbar_j + q_j qbar_i) = "
                "charge_state_factor * N_c * radiator_ij * "
                "(|delta_g_L|^2 + |delta_g_R|^2) / "
                "(SM total Z width weight + FCNC width weight)"
            ),
            "sm_hadronic_width_weights": dict(sm_hadronic_weights),
            "sm_total_width_weights": dict(sm_total_weights),
            "sm_hadronic_width_weight": float(sm_hadronic),
            "sm_total_width_weight": float(sm_total),
            "sm_hadronic_to_total_width_ratio": float(hadronic_to_total),
            "sin2_theta_eff": float(p.sin2_theta_eff),
            "m_z_gev": float(p.m_z_gev),
        },
    )


def down_fcnc_effective_coupling_limit(
    br_limit: float,
    *,
    flavor_i: str,
    flavor_j: str,
    inputs: ZPoleSMInputs | None = None,
    charge_state_factor: float = _DOWN_FCNC_DEFAULT_CHARGE_STATE_FACTOR,
) -> float:
    """Return the limit on ``sqrt(|delta_g_L|^2 + |delta_g_R|^2)``."""

    p = default_sm_inputs() if inputs is None else inputs
    first = _canonical_down_fcnc_flavor(flavor_i)
    second = _canonical_down_fcnc_flavor(flavor_j)
    if first == second:
        raise ValueError("down-sector FCNC Z decay requires distinct flavors")
    limit = _bounded_probability("br_limit", br_limit)
    charge_factor = _positive_float("charge_state_factor", charge_state_factor)
    sm_total = float(sum(down_fcnc_sm_total_width_weight(p).values()))
    radiator = _down_fcnc_pair_radiator(first, second, p)
    norm_limit = limit * sm_total / (
        charge_factor * 3.0 * radiator * (1.0 - limit)
    )
    return float(math.sqrt(norm_limit))


def down_fcnc_coupling_proxy(
    source: QuarkMassBasisCouplings,
    *,
    flavor_i: str,
    flavor_j: str,
    m_kk_gev: float | None = None,
    inputs: ZPoleSMInputs | None = None,
) -> ZPoleDownFCNCCouplingProxy:
    """Map quark mass-basis overlaps onto an off-diagonal down-sector Z proxy."""

    p = default_sm_inputs() if inputs is None else inputs
    first = _canonical_down_fcnc_flavor(flavor_i)
    second = _canonical_down_fcnc_flavor(flavor_j)
    if first == second:
        raise ValueError("down-sector FCNC Z decay requires distinct flavors")
    first_index = _DOWN_FCNC_FLAVOR_INDEX[first]
    second_index = _DOWN_FCNC_FLAVOR_INDEX[second]
    resolved_m_kk = _positive_float(
        "m_kk_gev",
        getattr(source, "M_KK") if m_kk_gev is None else m_kk_gev,
    )
    left_overlap = _overlap_matrix(source, "left_overlap", "left_down")
    right_overlap = _overlap_matrix(source, "right_down_overlap", "right_down")
    left_ij = _down_fcnc_matrix_off_diagonal(
        left_overlap,
        first_index,
        second_index,
        "left_down_overlap",
    )
    right_ij = _down_fcnc_matrix_off_diagonal(
        right_overlap,
        first_index,
        second_index,
        "right_down_overlap",
    )
    scale = float((p.m_z_gev / resolved_m_kk) ** 2)
    delta_left = complex(p.proxy_strength * scale * left_ij)
    delta_right = complex(p.proxy_strength * scale * right_ij)
    return ZPoleDownFCNCCouplingProxy(
        model_label=ZPOLE_MODEL_V1,
        matching_assumption=ZPOLE_DOWN_FCNC_PROXY_V1,
        flavor_i=first,
        flavor_j=second,
        M_KK=resolved_m_kk,
        matching_scale=resolved_m_kk,
        left_overlap_ij=left_ij,
        right_overlap_ij=right_ij,
        scale_factor=scale,
        delta_g_left_ij=delta_left,
        delta_g_right_ij=delta_right,
        diagnostics={
            "m_kk_gev": float(resolved_m_kk),
            "matching_scale_gev": float(resolved_m_kk),
            "m_z_gev": float(p.m_z_gev),
            "proxy_strength": float(p.proxy_strength),
            "scale_factor": float(scale),
            "flavor_i": first,
            "flavor_j": second,
            "left_overlap_ij": left_ij,
            "right_overlap_ij": right_ij,
            "delta_g_left_ij": delta_left,
            "delta_g_right_ij": delta_right,
            "matching_assumption": ZPOLE_DOWN_FCNC_PROXY_V1,
        },
    )


def down_fcnc_branching_fraction_with_proxy(
    source: QuarkMassBasisCouplings,
    *,
    flavor_i: str,
    flavor_j: str,
    br_limit: float | None = None,
    m_kk_gev: float | None = None,
    inputs: ZPoleSMInputs | None = None,
    charge_state_factor: float = _DOWN_FCNC_DEFAULT_CHARGE_STATE_FACTOR,
) -> tuple[ZPoleDownFCNCBranchingResult, ZPoleDownFCNCCouplingProxy]:
    """Evaluate an off-diagonal down-sector Z decay with the documented proxy."""

    p = default_sm_inputs() if inputs is None else inputs
    proxy = down_fcnc_coupling_proxy(
        source,
        flavor_i=flavor_i,
        flavor_j=flavor_j,
        m_kk_gev=m_kk_gev,
        inputs=p,
    )
    result = down_fcnc_branching_fraction_from_couplings(
        flavor_i=flavor_i,
        flavor_j=flavor_j,
        delta_g_left=proxy.delta_g_left_ij,
        delta_g_right=proxy.delta_g_right_ij,
        br_limit=br_limit,
        inputs=p,
        charge_state_factor=charge_state_factor,
    )
    return result, proxy


def _canonical_down_fcnc_flavor(flavor: str) -> str:
    name = str(flavor).lower()
    aliases = {
        "down": "d",
        "strange": "s",
        "bottom": "b",
        "beauty": "b",
    }
    canonical = aliases.get(name, name)
    if canonical not in _DOWN_FCNC_FLAVOR_INDEX:
        raise ValueError(f"{flavor!r} is not a supported down-sector flavor")
    return canonical


def _down_fcnc_pair_radiator(
    flavor_i: str,
    flavor_j: str,
    inputs: ZPoleSMInputs,
) -> float:
    radiator = 0.5 * (inputs.radiator_for(flavor_i) + inputs.radiator_for(flavor_j))
    return _positive_float("down-pair radiator", radiator)


def _down_fcnc_matrix_off_diagonal(
    matrix: np.ndarray,
    row: int,
    column: int,
    name: str,
) -> complex:
    value = complex(matrix[row, column])
    if not math.isfinite(value.real) or not math.isfinite(value.imag):
        raise ValueError(f"{name}[{row},{column}] must be finite")
    return value
