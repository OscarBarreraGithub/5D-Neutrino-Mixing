"""Higgs lepton-flavor-violating branching fractions.

This module provides the process-independent part of the Higgs LFV
constraints ``h -> l_i l_j``.  It uses the effective off-diagonal Yukawa
convention

    L = -h (bar l_i Y_ij P_R l_j + bar l_i Y_ji^* P_L l_j) + h.c.

and the charge-summed width

    Gamma(h -> l_i l_j) = m_h (|Y_ij|^2 + |Y_ji|^2) / (8 pi),
    BR = Gamma / Gamma_h^total.

For charged-LFV Higgs decays the catalog SM prediction is zero.

RS matching status
------------------
NEEDS-HUMAN-PHYSICS: the current ``ParameterPoint`` does not carry the
charged-lepton Higgs-Yukawa misalignment matrix.  The proxy helper accepts
caller-supplied off-diagonal Higgs Yukawa entries directly, or extracts them
from a supplied charged-lepton Higgs-Yukawa matrix.  This is a diagnostic
effective-coupling proxy, not a complete RS calculation of Higgs localization,
KK-fermion mixing, brane kinetic terms, or total-width modifications.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any, Mapping

import numpy as np

HIGGS_LFV_MODEL_V1 = "higgs_lfv_effective_offdiagonal_yukawa_v1"
HIGGS_LFV_INPUT_BUNDLE_V1 = "higgs_lfv_pdg_like_sm_inputs_v1"
HIGGS_LFV_YUKAWA_CONVENTION = (
    "L=-h (bar l_i Y_ij P_R l_j + bar l_i Y_ji^* P_L l_j) + h.c.; "
    "Gamma(h -> l_i l_j)=m_h*(|Y_ij|^2+|Y_ji|^2)/(8*pi)"
)
HIGGS_LFV_RS_PROXY_V1 = (
    "NEEDS-HUMAN-PHYSICS: off-diagonal charged-lepton Higgs Yukawas are not "
    "available on ParameterPoint; caller-supplied Higgs-Yukawa spurions are "
    "used directly as Y_ij and Y_ji, standing in for missing model-dependent "
    "RS Higgs-Yukawa misalignment and total-width matching."
)

_CHARGED_LEPTONS: tuple[str, ...] = ("e", "mu", "tau")
_LEPTON_INDEX = {"e": 0, "mu": 1, "tau": 2}


@dataclass(frozen=True)
class HiggsLFVInputs:
    """Numerical inputs for effective Higgs LFV widths."""

    input_bundle: str = HIGGS_LFV_INPUT_BUNDLE_V1
    higgs_mass_gev: float = 125.25
    total_higgs_width_gev: float = 4.07e-3
    constants_citation: str = (
        "PDG-like 125 GeV Higgs mass and SM total-width inputs; width formula "
        "follows the effective LFV Higgs-Yukawa convention documented in this "
        "module. Callers may override these inputs for width-changing models."
    )

    def __post_init__(self) -> None:
        for name in ("higgs_mass_gev", "total_higgs_width_gev"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")


@dataclass(frozen=True)
class HiggsLFVYukawaProxyInput:
    """Explicit effective off-diagonal Higgs-Yukawa proxy input."""

    initial_flavor: str
    final_flavor: str
    yukawa_ij: complex
    yukawa_ji: complex
    source: str = "caller-supplied Higgs LFV Yukawa proxy"


@dataclass(frozen=True)
class HiggsLFVYukawaProxy:
    """Documented proxy view over off-diagonal Higgs-Yukawa entries."""

    model_label: str
    matching_assumption: str
    initial_flavor: str
    final_flavor: str
    yukawa_ij: complex
    yukawa_ji: complex
    source: str
    diagnostics: Mapping[str, Any]


@dataclass(frozen=True)
class HiggsLFVBranchingResult:
    """Branching-fraction result for one charged-LFV Higgs decay."""

    model_label: str
    input_bundle: str
    initial_flavor: str
    final_flavor: str
    branching_fraction: float
    sm_branching_fraction: float
    ratio_to_limit: float | None
    br_limit: float | None
    passes: bool | None
    partial_width_gev: float
    total_higgs_width_gev: float
    yukawa_ij: complex
    yukawa_ji: complex
    yukawa_norm: float
    yukawa_norm_squared: float
    proxy: HiggsLFVYukawaProxy | None = None
    diagnostics: Mapping[str, Any] = field(default_factory=dict)


def default_higgs_lfv_inputs() -> HiggsLFVInputs:
    """Return the repo-owned default Higgs LFV input bundle."""

    return HiggsLFVInputs()


def h_lfv_yukawa_proxy_input(
    initial_flavor: str,
    final_flavor: str,
    yukawa_ij: complex,
    yukawa_ji: complex = 0.0j,
    *,
    source: str = "caller-supplied Higgs LFV Yukawa proxy",
) -> HiggsLFVYukawaProxyInput:
    """Build a shape-checked explicit proxy input for ``h -> l_i l_j``."""

    initial = _canonical_charged_lepton(initial_flavor)
    final = _canonical_charged_lepton(final_flavor)
    if initial == final:
        raise ValueError("Higgs LFV decay requires distinct charged leptons")
    return HiggsLFVYukawaProxyInput(
        initial_flavor=initial,
        final_flavor=final,
        yukawa_ij=_finite_complex(yukawa_ij, "yukawa_ij"),
        yukawa_ji=_finite_complex(yukawa_ji, "yukawa_ji"),
        source=str(source),
    )


def h_lfv_partial_width(
    yukawa_ij: complex = 0.0j,
    yukawa_ji: complex = 0.0j,
    *,
    inputs: HiggsLFVInputs | None = None,
) -> float:
    """Return ``Gamma(h -> l_i l_j)`` from off-diagonal Yukawa entries."""

    p = default_higgs_lfv_inputs() if inputs is None else inputs
    y_ij = _finite_complex(yukawa_ij, "yukawa_ij")
    y_ji = _finite_complex(yukawa_ji, "yukawa_ji")
    return float(
        p.higgs_mass_gev
        / (8.0 * math.pi)
        * (abs(y_ij) ** 2 + abs(y_ji) ** 2)
    )


def h_lfv_branching_fraction_from_yukawas(
    *,
    initial_flavor: str,
    final_flavor: str,
    yukawa_ij: complex = 0.0j,
    yukawa_ji: complex = 0.0j,
    br_limit: float | None = None,
    inputs: HiggsLFVInputs | None = None,
    proxy: HiggsLFVYukawaProxy | None = None,
) -> HiggsLFVBranchingResult:
    """Return ``BR(h -> l_i l_j)`` from effective LFV Higgs Yukawas."""

    p = default_higgs_lfv_inputs() if inputs is None else inputs
    initial = _canonical_charged_lepton(initial_flavor)
    final = _canonical_charged_lepton(final_flavor)
    if initial == final:
        raise ValueError("Higgs LFV decay requires distinct charged leptons")
    y_ij = _finite_complex(yukawa_ij, "yukawa_ij")
    y_ji = _finite_complex(yukawa_ji, "yukawa_ji")
    norm_squared = float(abs(y_ij) ** 2 + abs(y_ji) ** 2)
    norm = float(math.sqrt(norm_squared))
    width = h_lfv_partial_width(y_ij, y_ji, inputs=p)
    branching_fraction = float(width / p.total_higgs_width_gev)

    limit = None if br_limit is None else _bounded_probability("br_limit", br_limit)
    ratio = None if limit is None else float(branching_fraction / limit)
    passes = None if ratio is None else bool(ratio <= 1.0)

    diagnostics: dict[str, Any] = {
        "operator_convention": HIGGS_LFV_YUKAWA_CONVENTION,
        "branching_formula": (
            "BR(h -> l_i l_j) = m_h*(|Y_ij|^2+|Y_ji|^2) / "
            "(8*pi*Gamma_h_total)"
        ),
        "m_h_gev": float(p.higgs_mass_gev),
        "total_higgs_width_gev": float(p.total_higgs_width_gev),
        "constants_citation": p.constants_citation,
        "sm_branching_fraction": 0.0,
        "yukawa_ij": y_ij,
        "yukawa_ji": y_ji,
        "yukawa_norm": norm,
        "yukawa_norm_squared": norm_squared,
        "partial_width_gev": float(width),
        "fixed_total_width_policy": True,
    }
    if proxy is not None:
        diagnostics["matching_assumption"] = proxy.matching_assumption
        diagnostics["proxy_source"] = proxy.source

    return HiggsLFVBranchingResult(
        model_label=HIGGS_LFV_MODEL_V1,
        input_bundle=p.input_bundle,
        initial_flavor=initial,
        final_flavor=final,
        branching_fraction=branching_fraction,
        sm_branching_fraction=0.0,
        ratio_to_limit=ratio,
        br_limit=limit,
        passes=passes,
        partial_width_gev=float(width),
        total_higgs_width_gev=float(p.total_higgs_width_gev),
        yukawa_ij=y_ij,
        yukawa_ji=y_ji,
        yukawa_norm=norm,
        yukawa_norm_squared=norm_squared,
        proxy=proxy,
        diagnostics=diagnostics,
    )


def h_lfv_effective_yukawa_limit(
    br_limit: float,
    *,
    inputs: HiggsLFVInputs | None = None,
) -> float:
    """Return the bound on ``sqrt(|Y_ij|^2 + |Y_ji|^2)``."""

    p = default_higgs_lfv_inputs() if inputs is None else inputs
    limit = _bounded_probability("br_limit", br_limit)
    norm_limit = limit * p.total_higgs_width_gev * 8.0 * math.pi / p.higgs_mass_gev
    return float(math.sqrt(norm_limit))


def h_lfv_yukawa_proxy(
    source: object,
    *,
    initial_flavor: str,
    final_flavor: str,
) -> HiggsLFVYukawaProxy:
    """Return the documented effective-Yukawa proxy for ``h -> l_i l_j``."""

    initial = _canonical_charged_lepton(initial_flavor)
    final = _canonical_charged_lepton(final_flavor)
    proxy_input = _coerce_proxy_input(source, initial=initial, final=final)
    y_ij = _finite_complex(proxy_input.yukawa_ij, "yukawa_ij")
    y_ji = _finite_complex(proxy_input.yukawa_ji, "yukawa_ji")
    return HiggsLFVYukawaProxy(
        model_label=HIGGS_LFV_MODEL_V1,
        matching_assumption=HIGGS_LFV_RS_PROXY_V1,
        initial_flavor=initial,
        final_flavor=final,
        yukawa_ij=y_ij,
        yukawa_ji=y_ji,
        source=proxy_input.source,
        diagnostics={
            "initial_flavor": initial,
            "final_flavor": final,
            "yukawa_ij": y_ij,
            "yukawa_ji": y_ji,
            "yukawa_norm": float(
                math.sqrt(abs(y_ij) ** 2 + abs(y_ji) ** 2)
            ),
            "yukawa_norm_squared": float(abs(y_ij) ** 2 + abs(y_ji) ** 2),
            "proxy_source": proxy_input.source,
            "matching_assumption": HIGGS_LFV_RS_PROXY_V1,
        },
    )


def h_lfv_branching_fraction_with_proxy(
    source: object,
    *,
    initial_flavor: str,
    final_flavor: str,
    br_limit: float | None = None,
    inputs: HiggsLFVInputs | None = None,
) -> tuple[HiggsLFVBranchingResult, HiggsLFVYukawaProxy]:
    """Evaluate ``h -> l_i l_j`` using the documented Yukawa proxy."""

    p = default_higgs_lfv_inputs() if inputs is None else inputs
    proxy = h_lfv_yukawa_proxy(
        source,
        initial_flavor=initial_flavor,
        final_flavor=final_flavor,
    )
    result = h_lfv_branching_fraction_from_yukawas(
        initial_flavor=proxy.initial_flavor,
        final_flavor=proxy.final_flavor,
        yukawa_ij=proxy.yukawa_ij,
        yukawa_ji=proxy.yukawa_ji,
        br_limit=br_limit,
        inputs=p,
        proxy=proxy,
    )
    return result, proxy


def _coerce_proxy_input(
    value: Any,
    *,
    initial: str,
    final: str,
) -> HiggsLFVYukawaProxyInput:
    if isinstance(value, HiggsLFVYukawaProxyInput):
        return _orient_proxy_input(value, initial=initial, final=final)
    if isinstance(value, Mapping):
        return _proxy_from_mapping(value, initial=initial, final=final)
    return _proxy_from_object(value, initial=initial, final=final)


def _orient_proxy_input(
    value: HiggsLFVYukawaProxyInput,
    *,
    initial: str,
    final: str,
) -> HiggsLFVYukawaProxyInput:
    source_initial = _canonical_charged_lepton(value.initial_flavor)
    source_final = _canonical_charged_lepton(value.final_flavor)
    if source_initial == initial and source_final == final:
        return value
    if source_initial == final and source_final == initial:
        return HiggsLFVYukawaProxyInput(
            initial_flavor=initial,
            final_flavor=final,
            yukawa_ij=complex(value.yukawa_ji),
            yukawa_ji=complex(value.yukawa_ij),
            source=value.source,
        )
    raise ValueError(
        f"proxy input is for {source_initial}-{source_final}, not {initial}-{final}"
    )


def _proxy_from_mapping(
    mapping: Mapping[str, Any],
    *,
    initial: str,
    final: str,
) -> HiggsLFVYukawaProxyInput:
    source = str(mapping.get("source", "mapping Higgs LFV Yukawa proxy"))
    direct = _yukawas_from_mapping(mapping, initial=initial, final=final)
    if direct is not None:
        y_ij, y_ji = direct
    else:
        matrix = _first_present_key(
            mapping,
            (
                "higgs_yukawa_matrix",
                "higgs_yukawa",
                "charged_lepton_higgs_yukawa",
                "yukawa_matrix",
                "Y_higgs",
                "Y_l",
            ),
        )
        if matrix is None:
            raise KeyError(
                f"mapping must provide {initial}-{final} Higgs-Yukawa entries "
                "or a charged-lepton Higgs-Yukawa matrix"
            )
        y_ij, y_ji = _matrix_yukawas(matrix, initial=initial, final=final)
    return h_lfv_yukawa_proxy_input(
        initial,
        final,
        y_ij,
        y_ji,
        source=source,
    )


def _proxy_from_object(
    value: Any,
    *,
    initial: str,
    final: str,
) -> HiggsLFVYukawaProxyInput:
    source = str(getattr(value, "source", "object Higgs LFV Yukawa proxy"))
    direct = _yukawas_from_object(value, initial=initial, final=final)
    if direct is not None:
        y_ij, y_ji = direct
    else:
        matrix = _first_present_attr(
            value,
            (
                "higgs_yukawa_matrix",
                "higgs_yukawa",
                "charged_lepton_higgs_yukawa",
                "yukawa_matrix",
                "Y_higgs",
                "Y_l",
            ),
        )
        if matrix is None:
            raise TypeError(
                f"lepton Higgs LFV proxy must provide {initial}-{final} "
                "Yukawa entries or a charged-lepton Higgs-Yukawa matrix"
            )
        y_ij, y_ji = _matrix_yukawas(matrix, initial=initial, final=final)
    return h_lfv_yukawa_proxy_input(
        initial,
        final,
        y_ij,
        y_ji,
        source=source,
    )


def _yukawas_from_mapping(
    mapping: Mapping[str, Any],
    *,
    initial: str,
    final: str,
) -> tuple[complex, complex] | None:
    forward = _first_present_key(mapping, _entry_keys(initial, final))
    reverse = _first_present_key(mapping, _entry_keys(final, initial))
    if forward is None and reverse is None:
        return None
    return complex(0.0j if forward is None else forward), complex(
        0.0j if reverse is None else reverse
    )


def _yukawas_from_object(
    value: Any,
    *,
    initial: str,
    final: str,
) -> tuple[complex, complex] | None:
    forward = _first_present_attr(value, _entry_keys(initial, final))
    reverse = _first_present_attr(value, _entry_keys(final, initial))
    if forward is None and reverse is None:
        return None
    return complex(0.0j if forward is None else forward), complex(
        0.0j if reverse is None else reverse
    )


def _entry_keys(initial: str, final: str) -> tuple[str, ...]:
    compact = f"{initial}{final}"
    underscored = f"{initial}_{final}"
    return (
        f"y_{compact}",
        f"Y_{compact}",
        f"y_{underscored}",
        f"Y_{underscored}",
        f"yukawa_{compact}",
        f"Yukawa_{compact}",
        f"yukawa_{underscored}",
        f"Yukawa_{underscored}",
        f"higgs_yukawa_{compact}",
        f"Higgs_yukawa_{compact}",
        f"higgs_yukawa_{underscored}",
        f"Higgs_yukawa_{underscored}",
    )


def _matrix_yukawas(
    value: Any,
    *,
    initial: str,
    final: str,
) -> tuple[complex, complex]:
    matrix = np.asarray(value, dtype=np.complex128)
    if matrix.shape != (3, 3):
        raise ValueError("charged-lepton Higgs-Yukawa matrix must have shape (3, 3)")
    if not np.all(np.isfinite(matrix.real)) or not np.all(np.isfinite(matrix.imag)):
        raise ValueError("charged-lepton Higgs-Yukawa matrix entries must be finite")
    i = _LEPTON_INDEX[initial]
    j = _LEPTON_INDEX[final]
    return complex(matrix[i, j]), complex(matrix[j, i])


def _canonical_charged_lepton(flavor: str) -> str:
    aliases = {
        "electron": "e",
        "e-": "e",
        "e+": "e",
        "muon": "mu",
        "mu-": "mu",
        "mu+": "mu",
        "tau-": "tau",
        "tau+": "tau",
    }
    name = aliases.get(str(flavor).lower(), str(flavor).lower())
    if name not in _CHARGED_LEPTONS:
        raise ValueError(f"unsupported charged-lepton flavor {flavor!r}")
    return name


def _finite_complex(value: complex, name: str) -> complex:
    number = complex(value)
    if not math.isfinite(number.real) or not math.isfinite(number.imag):
        raise ValueError(f"{name} must be finite")
    return number


def _bounded_probability(name: str, value: object) -> float:
    number = float(value)
    if not math.isfinite(number) or not 0.0 < number < 1.0:
        raise ValueError(f"{name} must lie between zero and one")
    return number


def _first_present_key(mapping: Mapping[str, Any], keys: tuple[str, ...]) -> Any:
    for key in keys:
        if key in mapping:
            return mapping[key]
    return None


def _first_present_attr(value: Any, names: tuple[str, ...]) -> Any:
    for name in names:
        if hasattr(value, name):
            return getattr(value, name)
    return None


__all__ = [
    "HIGGS_LFV_MODEL_V1",
    "HIGGS_LFV_INPUT_BUNDLE_V1",
    "HIGGS_LFV_YUKAWA_CONVENTION",
    "HIGGS_LFV_RS_PROXY_V1",
    "HiggsLFVInputs",
    "HiggsLFVYukawaProxyInput",
    "HiggsLFVYukawaProxy",
    "HiggsLFVBranchingResult",
    "default_higgs_lfv_inputs",
    "h_lfv_yukawa_proxy_input",
    "h_lfv_partial_width",
    "h_lfv_branching_fraction_from_yukawas",
    "h_lfv_effective_yukawa_limit",
    "h_lfv_yukawa_proxy",
    "h_lfv_branching_fraction_with_proxy",
]
