"""Adapter for coherent mu->e conversion constraints.

Constraint modules import this adapter only.  It composes the L001 lepton
dipole adapter with the reusable :mod:`quarkConstraints.mu_e_conversion`
overlap-integral core.

Tree-level light-Z vector coefficients are read from Phase-4a
``rs_ew_couplings`` up/down contacts when available.  Scalar and dipole pieces,
including their relative phase, remain deferred partial inputs.
"""

from __future__ import annotations

import math
from typing import Any, Mapping

from quarkConstraints.mu_e_conversion import (
    G_F_GEV_MINUS2,
    MU_E_CONVERSION_DIPOLE_CONVENTION,
    MU_E_CONVERSION_INPUT_BUNDLE_V1,
    MU_E_CONVERSION_MODEL_V1,
    MU_E_CONVERSION_OPERATOR_CONVENTION,
    MU_E_CONVERSION_PROXY_V1,
    MuEConversionCoefficientProxyInput,
    MuEConversionNuclearInputs,
    MuEConversionResult,
    aluminum_nuclear_inputs as _aluminum_nuclear_inputs,
    gold_nuclear_inputs as _gold_nuclear_inputs,
    mu_e_conversion_coefficients as _coefficients,
    mu_e_conversion_from_components as _from_components,
    mu_e_conversion_has_coefficient_proxy as _has_coefficient_proxy,
    mu_e_conversion_proxy_input as _proxy_input,
    nuclear_inputs_for_target as _nuclear_inputs_for_target,
    titanium_nuclear_inputs as _titanium_nuclear_inputs,
    zero_mu_e_conversion_coefficients as _zero_coefficients,
)

from .lepton import mu_to_e_gamma_from_lepton_input

__all__ = [
    "MU_E_CONVERSION_MODEL_V1",
    "MU_E_CONVERSION_INPUT_BUNDLE_V1",
    "MU_E_CONVERSION_OPERATOR_CONVENTION",
    "MU_E_CONVERSION_DIPOLE_CONVENTION",
    "MU_E_CONVERSION_PROXY_V1",
    "MU_E_CONVERSION_VECTOR_TREE_RIGOROUS_V1",
    "MU_E_CONVERSION_DEFERRED_PIECES_V1",
    "MuEConversionCoefficientProxyInput",
    "MuEConversionNuclearInputs",
    "MuEConversionResult",
    "mu_e_conversion_aluminum_nuclear_inputs",
    "mu_e_conversion_titanium_nuclear_inputs",
    "mu_e_conversion_gold_nuclear_inputs",
    "mu_e_conversion_nuclear_inputs_for_target",
    "mu_e_conversion_proxy_input",
    "mu_e_conversion_from_lepton_input",
]

MU_E_CONVERSION_VECTOR_TREE_RIGOROUS_V1 = (
    "tree-level light-Z vector mu-e conversion coefficients from "
    "rs_ew_couplings diagonal up/down llqq contacts; zero for the Phase-4a "
    "diagonal charged-lepton fit with U_e=I and universal c_L"
)
MU_E_CONVERSION_DEFERRED_PIECES_V1 = (
    "NEEDS-HUMAN-PHYSICS: loop dipole matching, dipole-contact relative phase "
    "when chiral dipoles are not supplied, scalar Higgs-lepton matching, "
    "heavy neutral exchange, and vector-scalar interference with deferred "
    "scalars remain partial"
)
_VECTOR_TREE_MISSING_ASSUMPTION = (
    "tree-level light-Z vector mu-e conversion contact not evaluated because "
    "rs_ew_couplings is absent; zero vector placeholder used while available "
    "partial scalar/dipole pieces are preserved"
)
_LEPTON_FINAL_E = 0
_LEPTON_INITIAL_MU = 1
_UP_INDEX = 0
_DOWN_INDEX = 0
_VECTOR_COEFFICIENT_KEYS = (
    "g_lv_p",
    "g_LV_p",
    "g_lv_proton",
    "g_LV_proton",
    "g_lv_n",
    "g_LV_n",
    "g_lv_neutron",
    "g_LV_neutron",
    "g_rv_p",
    "g_RV_p",
    "g_rv_proton",
    "g_RV_proton",
    "g_rv_n",
    "g_RV_n",
    "g_rv_neutron",
    "g_RV_neutron",
)
_SCALAR_COEFFICIENT_KEYS = (
    "g_ls_p",
    "g_LS_p",
    "g_ls_proton",
    "g_LS_proton",
    "g_ls_n",
    "g_LS_n",
    "g_ls_neutron",
    "g_LS_neutron",
    "g_rs_p",
    "g_RS_p",
    "g_rs_proton",
    "g_RS_proton",
    "g_rs_n",
    "g_RS_n",
    "g_rs_neutron",
    "g_RS_neutron",
)
_COEFFICIENT_DIPOLE_KEYS = (
    "dipole_amplitude_left",
    "dipole_A_L",
    "a_left",
    "A_L",
    "dipole_amplitude_right",
    "dipole_A_R",
    "a_right",
    "A_R",
)


def mu_e_conversion_aluminum_nuclear_inputs() -> MuEConversionNuclearInputs:
    """Return bundled aluminum nuclear inputs."""

    return _aluminum_nuclear_inputs()


def mu_e_conversion_titanium_nuclear_inputs() -> MuEConversionNuclearInputs:
    """Return bundled titanium nuclear inputs."""

    return _titanium_nuclear_inputs()


def mu_e_conversion_gold_nuclear_inputs() -> MuEConversionNuclearInputs:
    """Return bundled gold nuclear inputs."""

    return _gold_nuclear_inputs()


def mu_e_conversion_nuclear_inputs_for_target(
    target: str,
) -> MuEConversionNuclearInputs:
    """Return bundled nuclear inputs for a supported target."""

    return _nuclear_inputs_for_target(target)


def mu_e_conversion_proxy_input(
    *,
    g_lv_p: complex = 0.0j,
    g_lv_n: complex = 0.0j,
    g_rv_p: complex = 0.0j,
    g_rv_n: complex = 0.0j,
    g_ls_p: complex = 0.0j,
    g_ls_n: complex = 0.0j,
    g_rs_p: complex = 0.0j,
    g_rs_n: complex = 0.0j,
    dipole_amplitude_left: complex | None = None,
    dipole_amplitude_right: complex | None = None,
    m_kk_gev: float | None = None,
    source: str = "caller-supplied mu-e conversion coefficient proxy",
) -> MuEConversionCoefficientProxyInput:
    """Build a coefficient proxy accepted by this adapter."""

    return _proxy_input(
        g_lv_p=g_lv_p,
        g_lv_n=g_lv_n,
        g_rv_p=g_rv_p,
        g_rv_n=g_rv_n,
        g_ls_p=g_ls_p,
        g_ls_n=g_ls_n,
        g_rs_p=g_rs_p,
        g_rs_n=g_rs_n,
        dipole_amplitude_left=dipole_amplitude_left,
        dipole_amplitude_right=dipole_amplitude_right,
        m_kk_gev=m_kk_gev,
        source=source,
    )


def mu_e_conversion_from_lepton_input(
    lepton_input: Any,
    *,
    conversion_rate_limit: float,
    dipole_br_limit: float,
    dipole_prefactor_br: float,
    reference_scale_gev: float = 3000.0,
    m_kk_gev: float | None = None,
    nuclear_inputs: MuEConversionNuclearInputs | None = None,
) -> MuEConversionResult:
    """Evaluate ``CR(mu -> e, target)`` from dipole and coefficient proxies."""

    has_dipole = _has_dipole_input(lepton_input)
    has_rs_ew = _has_rs_ew_couplings(lepton_input)
    has_any_coefficients = _has_any_coefficient_proxy(lepton_input)
    has_partial_coefficients = _has_partial_coefficient_proxy(lepton_input)
    has_coefficients = has_rs_ew or has_partial_coefficients or (
        has_dipole and has_any_coefficients
    )
    if not has_dipole and not has_coefficients:
        raise TypeError(
            "lepton input must provide a dipole source and/or mu-e conversion "
            "scalar/dipole coefficient or rs_ew_couplings vector inputs"
        )

    dipole_parent_br = 0.0
    resolved_m_kk = m_kk_gev
    dipole_diagnostics: dict[str, Any] = {
        "dipole_input_present": bool(has_dipole),
        "dipole_component_reuses": "flavor_catalog_constraints.physics_adapters.lepton",
    }
    if has_dipole:
        dipole_source = _dipole_source(lepton_input)
        try:
            dipole = mu_to_e_gamma_from_lepton_input(
                dipole_source,
                br_limit=dipole_br_limit,
                prefactor_br=dipole_prefactor_br,
                reference_scale_gev=reference_scale_gev,
                m_kk_gev=m_kk_gev,
            )
        except (KeyError, TypeError, ValueError) as exc:
            if not has_coefficients:
                raise
            dipole_diagnostics.update(
                {
                    "dipole_evaluation_failed": True,
                    "dipole_exception_type": type(exc).__name__,
                    "dipole_exception": str(exc),
                }
            )
        else:
            dipole_parent_br = float(dipole.branching_fraction)
            resolved_m_kk = float(dipole.m_kk_gev)
            dipole_diagnostics.update(
                {
                    "dipole_evaluation_failed": False,
                    "dipole_parent_branching_fraction": float(
                        dipole.branching_fraction
                    ),
                    "dipole_parent_ratio_to_limit": float(dipole.ratio_to_limit),
                    "dipole_parent_limit": float(dipole.br_limit),
                    "dipole_lhs": float(dipole.dipole_lhs),
                    "dipole_rhs": float(dipole.dipole_rhs),
                    "dipole_off_diagonal_12": complex(dipole.off_diagonal_12),
                    "dipole_product_matrix": dipole.product_matrix,
                    "dipole_input_kind": dipole.input_kind,
                    "dipole_used_proxy": bool(dipole.used_proxy),
                }
            )

    coefficient_diagnostics: dict[str, Any]
    if has_rs_ew:
        coefficients, coefficient_diagnostics = _rigorous_vector_coefficients(
            lepton_input
        )
    elif has_any_coefficients:
        coefficients, coefficient_diagnostics = _absent_vector_coefficients(
            lepton_input
        )
    else:
        coefficients = _zero_coefficients(source="dipole-only mu-e conversion proxy")
        coefficient_diagnostics = {
            "vector_tree_rigorous": False,
            "vector_tree_source": "absent",
            "vector_tree_matching": _VECTOR_TREE_MISSING_ASSUMPTION,
            "matching_assumption": _VECTOR_TREE_MISSING_ASSUMPTION,
            "vector_tree_missing_extra": "rs_ew_couplings",
            "missing_extra": "rs_ew_couplings",
            "legacy_vector_proxy_ignored": False,
            "deferred_matching_needs_human_physics": (
                MU_E_CONVERSION_DEFERRED_PIECES_V1
            ),
        }
    if coefficients.m_kk_gev is None and resolved_m_kk is not None:
        coefficients = _proxy_input(
            g_lv_p=coefficients.g_lv_p,
            g_lv_n=coefficients.g_lv_n,
            g_rv_p=coefficients.g_rv_p,
            g_rv_n=coefficients.g_rv_n,
            g_ls_p=coefficients.g_ls_p,
            g_ls_n=coefficients.g_ls_n,
            g_rs_p=coefficients.g_rs_p,
            g_rs_n=coefficients.g_rs_n,
            dipole_amplitude_left=coefficients.dipole_amplitude_left,
            dipole_amplitude_right=coefficients.dipole_amplitude_right,
            m_kk_gev=resolved_m_kk,
            source=coefficients.source,
        )

    result = _from_components(
        dipole_parent_branching_fraction=dipole_parent_br,
        coefficients=coefficients,
        conversion_rate_limit=conversion_rate_limit,
        nuclear_inputs=nuclear_inputs,
    )
    diagnostics = dict(result.diagnostics)
    diagnostics.update(dipole_diagnostics)
    diagnostics.update(coefficient_diagnostics)
    return MuEConversionResult(
        model_label=result.model_label,
        input_bundle=result.input_bundle,
        target=result.target,
        conversion_rate=result.conversion_rate,
        conversion_rate_lower=result.conversion_rate_lower,
        conversion_rate_upper=result.conversion_rate_upper,
        sm_conversion_rate=result.sm_conversion_rate,
        np_shift_conversion_rate=result.np_shift_conversion_rate,
        conversion_width_gev=result.conversion_width_gev,
        conversion_width_lower_gev=result.conversion_width_lower_gev,
        conversion_width_upper_gev=result.conversion_width_upper_gev,
        capture_rate_gev=result.capture_rate_gev,
        dipole_component=result.dipole_component,
        scalar_component=result.scalar_component,
        vector_component=result.vector_component,
        vector_scalar_interference_component=(
            result.vector_scalar_interference_component
        ),
        contact_component=result.contact_component,
        dipole_contact_interference_component=(
            result.dipole_contact_interference_component
        ),
        dipole_contact_interference_lower=(
            result.dipole_contact_interference_lower
        ),
        dipole_contact_interference_upper=(
            result.dipole_contact_interference_upper
        ),
        dipole_contact_interference_treatment=(
            result.dipole_contact_interference_treatment
        ),
        dipole_parent_branching_fraction=result.dipole_parent_branching_fraction,
        dipole_amplitude_left=result.dipole_amplitude_left,
        dipole_amplitude_right=result.dipole_amplitude_right,
        left_nuclear_amplitude=result.left_nuclear_amplitude,
        right_nuclear_amplitude=result.right_nuclear_amplitude,
        ratio_to_limit=result.ratio_to_limit,
        ratio_to_limit_lower=result.ratio_to_limit_lower,
        ratio_to_limit_upper=result.ratio_to_limit_upper,
        conversion_rate_limit=result.conversion_rate_limit,
        passes=result.passes,
        coefficients=result.coefficients,
        nuclear_inputs=result.nuclear_inputs,
        diagnostics=diagnostics,
    )


def _has_dipole_input(value: Any) -> bool:
    if isinstance(value, Mapping):
        if "dipole" in value or "yukawa_result" in value:
            return True
        has_y = any(key in value for key in ("y_n_bar", "Y_N_bar"))
        has_u = any(key in value for key in ("pmns", "PMNS", "pmns_matrix"))
        return bool(has_y and has_u)
    return hasattr(value, "Y_N_matrix") or (
        hasattr(value, "y_n_bar") and hasattr(value, "pmns")
    )


def _dipole_source(value: Any) -> Any:
    if isinstance(value, Mapping) and "dipole" in value:
        return value["dipole"]
    if isinstance(value, Mapping) and "lepton_mass_basis_couplings" in value:
        return value["lepton_mass_basis_couplings"]
    return value


def _has_partial_coefficient_proxy(value: Any) -> bool:
    return _has_any_key_or_attr(
        value,
        (*_SCALAR_COEFFICIENT_KEYS, *_COEFFICIENT_DIPOLE_KEYS),
    )


def _has_any_coefficient_proxy(value: Any) -> bool:
    return _has_partial_coefficient_proxy(value) or _has_vector_coefficient_proxy(value)


def _has_vector_coefficient_proxy(value: Any) -> bool:
    return _has_any_key_or_attr(value, _VECTOR_COEFFICIENT_KEYS)


def _has_any_key_or_attr(value: Any, keys: tuple[str, ...]) -> bool:
    if isinstance(value, Mapping):
        if any(_is_live_complex(value[key]) for key in keys if key in value):
            return True
        nested = value.get("lepton_mass_basis_couplings")
        return bool(
            nested is not None
            and nested is not value
            and _has_any_key_or_attr(nested, keys)
        )
    if any(
        _is_live_complex(getattr(value, key))
        for key in keys
        if hasattr(value, key)
    ):
        return True
    nested = getattr(value, "lepton_mass_basis_couplings", None)
    return bool(
        nested is not None
        and nested is not value
        and _has_any_key_or_attr(nested, keys)
    )


def _is_live_complex(value: Any) -> bool:
    if value is None:
        return False
    try:
        number = complex(value)
    except (TypeError, ValueError):
        return True
    if not math.isfinite(number.real) or not math.isfinite(number.imag):
        return True
    return abs(number) > 0.0


def _has_rs_ew_couplings(value: Any) -> bool:
    return _rs_ew_couplings(value) is not None


def _rs_ew_couplings(value: Any) -> Any | None:
    if _looks_like_rs_ew_couplings(value):
        return value
    if isinstance(value, Mapping):
        candidate = value.get("rs_ew_couplings")
        return candidate if _looks_like_rs_ew_couplings(candidate) else None
    candidate = getattr(value, "rs_ew_couplings", None)
    return candidate if _looks_like_rs_ew_couplings(candidate) else None


def _looks_like_rs_ew_couplings(value: Any) -> bool:
    return hasattr(value, "contact") and hasattr(value, "contact_units")


def _rigorous_vector_coefficients(
    source: Any,
) -> tuple[MuEConversionCoefficientProxyInput, dict[str, Any]]:
    couplings = _rs_ew_couplings(source)
    if couplings is None:
        raise TypeError("rs_ew_couplings vector-contact source is absent")
    if getattr(couplings, "contact_units", None) != "GeV^-2":
        raise ValueError("rs_ew_couplings.contact_units must be 'GeV^-2'")

    base = _coefficient_base(source)
    contacts = {
        "u_LL": _contact(couplings, "u", "L", "L", _UP_INDEX),
        "u_LR": _contact(couplings, "u", "L", "R", _UP_INDEX),
        "u_RL": _contact(couplings, "u", "R", "L", _UP_INDEX),
        "u_RR": _contact(couplings, "u", "R", "R", _UP_INDEX),
        "d_LL": _contact(couplings, "d", "L", "L", _DOWN_INDEX),
        "d_LR": _contact(couplings, "d", "L", "R", _DOWN_INDEX),
        "d_RL": _contact(couplings, "d", "R", "L", _DOWN_INDEX),
        "d_RR": _contact(couplings, "d", "R", "R", _DOWN_INDEX),
    }
    norm = math.sqrt(2.0) * G_F_GEV_MINUS2
    g_lv_u = complex((contacts["u_LL"] + contacts["u_RL"]) / norm)
    g_rv_u = complex((contacts["u_LR"] + contacts["u_RR"]) / norm)
    g_lv_d = complex((contacts["d_LL"] + contacts["d_RL"]) / norm)
    g_rv_d = complex((contacts["d_LR"] + contacts["d_RR"]) / norm)
    g_lv_p = complex(2.0 * g_lv_u + g_lv_d)
    g_lv_n = complex(g_lv_u + 2.0 * g_lv_d)
    g_rv_p = complex(2.0 * g_rv_u + g_rv_d)
    g_rv_n = complex(g_rv_u + 2.0 * g_rv_d)
    coeffs = _proxy_input(
        g_lv_p=g_lv_p,
        g_lv_n=g_lv_n,
        g_rv_p=g_rv_p,
        g_rv_n=g_rv_n,
        g_ls_p=base.g_ls_p,
        g_ls_n=base.g_ls_n,
        g_rs_p=base.g_rs_p,
        g_rs_n=base.g_rs_n,
        dipole_amplitude_left=base.dipole_amplitude_left,
        dipole_amplitude_right=base.dipole_amplitude_right,
        m_kk_gev=getattr(couplings, "kk_ew_mass_gev", base.m_kk_gev),
        source="rigorous rs_ew_couplings vector contacts plus deferred scalar/dipole inputs",
    )
    vector_zero = all(abs(value) <= 1.0e-30 for value in (g_lv_p, g_lv_n, g_rv_p, g_rv_n))
    diagnostics = {
        "vector_tree_rigorous": True,
        "vector_tree_source": "rs_ew_couplings",
        "vector_tree_matching": MU_E_CONVERSION_VECTOR_TREE_RIGOROUS_V1,
        "matching_assumption": MU_E_CONVERSION_VECTOR_TREE_RIGOROUS_V1,
        "vector_tree_zero": vector_zero,
        "vector_tree_zero_for_diagonal_fit": vector_zero,
        "vector_tree_input_path": (
            "contact(u/d, q_chirality, lepton_chirality, q, q, e, mu)"
        ),
        "vector_tree_contact_units": "GeV^-2",
        "vector_tree_quark_to_kko_normalization": (
            "g_LV(q)=(C_LL+C_RL)/(sqrt(2)*G_F), "
            "g_RV(q)=(C_LR+C_RR)/(sqrt(2)*G_F)"
        ),
        "vector_tree_nucleon_map": (
            "g^p=2*g^u+g^d, g^n=g^u+2*g^d"
        ),
        "g_lv_u": g_lv_u,
        "g_lv_d": g_lv_d,
        "g_rv_u": g_rv_u,
        "g_rv_d": g_rv_d,
        "rs_ew_vector_contacts": contacts,
        "scalar_matching_status": "deferred_proxy_or_zero_NEEDS-HUMAN-PHYSICS",
        "deferred_matching_needs_human_physics": (
            MU_E_CONVERSION_DEFERRED_PIECES_V1
        ),
    }
    return coeffs, diagnostics


def _absent_vector_coefficients(
    source: Any,
) -> tuple[MuEConversionCoefficientProxyInput, dict[str, Any]]:
    base = _coefficient_base(source)
    legacy_vector_ignored = _has_vector_coefficient_proxy(source)
    coeffs = _proxy_input(
        g_lv_p=0.0j,
        g_lv_n=0.0j,
        g_rv_p=0.0j,
        g_rv_n=0.0j,
        g_ls_p=base.g_ls_p,
        g_ls_n=base.g_ls_n,
        g_rs_p=base.g_rs_p,
        g_rs_n=base.g_rs_n,
        dipole_amplitude_left=base.dipole_amplitude_left,
        dipole_amplitude_right=base.dipole_amplitude_right,
        m_kk_gev=base.m_kk_gev,
        source="absent rs_ew_couplings vector placeholder plus deferred scalar/dipole inputs",
    )
    diagnostics = {
        "vector_tree_rigorous": False,
        "vector_tree_source": "absent",
        "vector_tree_matching": _VECTOR_TREE_MISSING_ASSUMPTION,
        "matching_assumption": _VECTOR_TREE_MISSING_ASSUMPTION,
        "vector_tree_missing_extra": "rs_ew_couplings",
        "missing_extra": "rs_ew_couplings",
        "vector_tree_zero": True,
        "legacy_vector_proxy_ignored": legacy_vector_ignored,
        "scalar_matching_status": "deferred_proxy_or_zero_NEEDS-HUMAN-PHYSICS",
        "deferred_matching_needs_human_physics": (
            MU_E_CONVERSION_DEFERRED_PIECES_V1
        ),
    }
    return coeffs, diagnostics


def _coefficient_base(source: Any) -> MuEConversionCoefficientProxyInput:
    if _has_coefficient_proxy(source):
        return _coefficients(source)
    nested = None
    if isinstance(source, Mapping):
        nested = source.get("lepton_mass_basis_couplings")
    else:
        nested = getattr(source, "lepton_mass_basis_couplings", None)
    if nested is not None and nested is not source and _has_coefficient_proxy(nested):
        return _coefficients(nested)
    return _coefficients(source)


def _contact(
    couplings: Any,
    quark_sector: str,
    quark_chirality: str,
    lepton_chirality: str,
    quark_index: int,
) -> complex:
    value = complex(
        couplings.contact(
            quark_sector,
            quark_chirality,
            lepton_chirality,
            quark_index,
            quark_index,
            _LEPTON_FINAL_E,
            _LEPTON_INITIAL_MU,
        )
    )
    if not math.isfinite(value.real) or not math.isfinite(value.imag):
        raise ValueError("rs_ew_couplings mu-e vector contact must be finite")
    return value
