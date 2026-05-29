"""Coherent mu->e conversion from dipole and nucleon contact amplitudes.

This module implements the low-energy overlap-integral formula used for
coherent muon-to-electron conversion in a nucleus,

    omega_conv = 2 G_F^2 m_mu^5
                 (|A_R^* D + g_LS^p S_p + g_LS^n S_n
                   + g_LV^p V_p + g_LV^n V_n|^2
                  + |A_L^* D + g_RS^p S_p + g_RS^n S_n
                   + g_RV^p V_p + g_RV^n V_n|^2),
    CR = omega_conv / Gamma_capture.

The convention is the Kitano-Koike-Okada nucleon-coefficient convention for
dimensionless ``g`` coefficients normalized to ``G_F``.  Dipole amplitudes are
normalized by ``BR(mu -> e gamma) = 384*pi^2*(|A_L|^2 + |A_R|^2)``.

NEEDS-HUMAN-PHYSICS: full RS matching to the scalar/vector nucleon
coefficients requires charged-lepton neutral-current, Higgs/scalar, EW KK/Z/Z',
and lepton-quark matching inputs that are not on the catalog ParameterPoint.
The coefficient inputs accepted here are explicit low-energy proxies.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any, Mapping

MU_E_CONVERSION_MODEL_V1 = "mu_e_conversion_kko_overlap_proxy_v1"
MU_E_CONVERSION_INPUT_BUNDLE_V1 = "mu_e_conversion_nuclear_inputs_2026_v1"
MU_E_CONVERSION_OPERATOR_CONVENTION = (
    "omega_conv = 2*G_F^2*m_mu^5*(|A_R^*D + g_LS^p S_p + g_LS^n S_n + "
    "g_LV^p V_p + g_LV^n V_n|^2 + |A_L^*D + g_RS^p S_p + "
    "g_RS^n S_n + g_RV^p V_p + g_RV^n V_n|^2); "
    "CR = omega_conv/Gamma_capture; g coefficients are dimensionless "
    "nucleon coefficients in the Kitano-Koike-Okada convention; tabulated "
    "overlaps D,V,S are quoted in units of m_mu^(5/2)."
)
MU_E_CONVERSION_DIPOLE_CONVENTION = (
    "Dipole amplitudes use BR(mu -> e gamma) = "
    "384*pi^2*(|A_L|^2 + |A_R|^2). If explicit chiral dipole amplitudes are "
    "not available, the dipole-contact relative phase is unknown; the result "
    "reports the destructive-to-constructive conversion-rate interval and uses "
    "the lower envelope for the limit verdict."
)
MU_E_CONVERSION_PROXY_V1 = (
    "NEEDS-HUMAN-PHYSICS: scalar/vector lepton-quark RS matching is not "
    "available on ParameterPoint; caller-supplied low-energy nucleon "
    "coefficients g_LV/RV/LS/RS^(p,n) are used as documented proxies."
)

HBAR_GEV_S = 6.582119569e-25
G_F_GEV_MINUS2 = 1.1663787e-5
MUON_MASS_GEV = 0.1056583755
KKO_OVERLAP_DIMENSION_FACTOR_GEV5 = MUON_MASS_GEV**5


@dataclass(frozen=True)
class MuEConversionNuclearInputs:
    """Target-dependent overlap integrals and capture rate."""

    target: str
    isotope: str
    D: float
    V_p: float
    V_n: float
    S_p: float
    S_n: float
    capture_rate_s_inv: float
    source: str
    source_url: str
    input_bundle: str = MU_E_CONVERSION_INPUT_BUNDLE_V1

    @property
    def capture_rate_gev(self) -> float:
        return float(self.capture_rate_s_inv * HBAR_GEV_S)


@dataclass(frozen=True)
class MuEConversionCoefficientProxyInput:
    """Explicit low-energy proxy coefficients for mu->e conversion.

    The ``g_*`` fields are dimensionless proton/neutron coefficients in the
    KKO nucleon convention.  ``dipole_amplitude_left/right`` are optional
    chiral dipole amplitudes in the same normalization as the parent
    ``mu -> e gamma`` branching fraction.
    """

    g_lv_p: complex = 0.0j
    g_lv_n: complex = 0.0j
    g_rv_p: complex = 0.0j
    g_rv_n: complex = 0.0j
    g_ls_p: complex = 0.0j
    g_ls_n: complex = 0.0j
    g_rs_p: complex = 0.0j
    g_rs_n: complex = 0.0j
    dipole_amplitude_left: complex | None = None
    dipole_amplitude_right: complex | None = None
    m_kk_gev: float | None = None
    source: str = "caller-supplied mu-e conversion coefficient proxy"


@dataclass(frozen=True)
class MuEConversionResult:
    """Conversion-rate prediction normalized to the target capture rate."""

    model_label: str
    input_bundle: str
    target: str
    conversion_rate: float
    conversion_rate_lower: float
    conversion_rate_upper: float
    sm_conversion_rate: float
    np_shift_conversion_rate: float
    conversion_width_gev: float
    conversion_width_lower_gev: float
    conversion_width_upper_gev: float
    capture_rate_gev: float
    dipole_component: float
    scalar_component: float
    vector_component: float
    vector_scalar_interference_component: float
    contact_component: float
    dipole_contact_interference_component: float
    dipole_contact_interference_lower: float
    dipole_contact_interference_upper: float
    dipole_contact_interference_treatment: str
    dipole_parent_branching_fraction: float
    dipole_amplitude_left: complex | None
    dipole_amplitude_right: complex | None
    left_nuclear_amplitude: complex
    right_nuclear_amplitude: complex
    ratio_to_limit: float | None
    ratio_to_limit_lower: float | None
    ratio_to_limit_upper: float | None
    conversion_rate_limit: float | None
    passes: bool | None
    coefficients: MuEConversionCoefficientProxyInput
    nuclear_inputs: MuEConversionNuclearInputs
    diagnostics: Mapping[str, Any] = field(default_factory=dict)


def aluminum_nuclear_inputs() -> MuEConversionNuclearInputs:
    """Return Al overlap inputs reusable by L003.

    Values are the standard KKO aluminum overlap/capture inputs used in
    mu->e conversion phenomenology.
    """

    return MuEConversionNuclearInputs(
        target="Al",
        isotope="27Al",
        D=0.0362,
        V_p=0.0161,
        V_n=0.0173,
        S_p=0.0155,
        S_n=0.0167,
        capture_rate_s_inv=0.7054e6,
        source="Kitano, Koike, Okada, Phys. Rev. D66 (2002) 096002, Table I",
        source_url="https://arxiv.org/abs/hep-ph/0203110",
    )


def titanium_nuclear_inputs() -> MuEConversionNuclearInputs:
    """Return Ti overlap inputs for reuse by L005."""

    return MuEConversionNuclearInputs(
        target="Ti",
        isotope="48Ti",
        D=0.0864,
        V_p=0.0396,
        V_n=0.0468,
        S_p=0.0368,
        S_n=0.0435,
        capture_rate_s_inv=2.59e6,
        source="Kitano, Koike, Okada, Phys. Rev. D66 (2002) 096002, Table I",
        source_url="https://arxiv.org/abs/hep-ph/0203110",
    )


def gold_nuclear_inputs() -> MuEConversionNuclearInputs:
    """Return Au overlap inputs for reuse by L004."""

    return MuEConversionNuclearInputs(
        target="Au",
        isotope="197Au",
        D=0.189,
        V_p=0.0974,
        V_n=0.146,
        S_p=0.0614,
        S_n=0.0918,
        capture_rate_s_inv=13.07e6,
        source="Kitano, Koike, Okada, Phys. Rev. D66 (2002) 096002, Table I",
        source_url="https://arxiv.org/abs/hep-ph/0203110",
    )


def nuclear_inputs_for_target(target: str) -> MuEConversionNuclearInputs:
    """Return bundled nuclear inputs for a supported target label."""

    key = str(target).strip().lower()
    if key in {"al", "aluminum", "aluminium", "27al"}:
        return aluminum_nuclear_inputs()
    if key in {"ti", "titanium", "48ti"}:
        return titanium_nuclear_inputs()
    if key in {"au", "gold", "197au"}:
        return gold_nuclear_inputs()
    raise ValueError(f"unsupported mu-e conversion target {target!r}")


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
    """Build shape-checked low-energy coefficient proxy inputs."""

    return MuEConversionCoefficientProxyInput(
        g_lv_p=_finite_complex(g_lv_p, "g_lv_p"),
        g_lv_n=_finite_complex(g_lv_n, "g_lv_n"),
        g_rv_p=_finite_complex(g_rv_p, "g_rv_p"),
        g_rv_n=_finite_complex(g_rv_n, "g_rv_n"),
        g_ls_p=_finite_complex(g_ls_p, "g_ls_p"),
        g_ls_n=_finite_complex(g_ls_n, "g_ls_n"),
        g_rs_p=_finite_complex(g_rs_p, "g_rs_p"),
        g_rs_n=_finite_complex(g_rs_n, "g_rs_n"),
        dipole_amplitude_left=_optional_finite_complex(
            dipole_amplitude_left,
            "dipole_amplitude_left",
        ),
        dipole_amplitude_right=_optional_finite_complex(
            dipole_amplitude_right,
            "dipole_amplitude_right",
        ),
        m_kk_gev=None if m_kk_gev is None else _positive_float(m_kk_gev, "m_kk_gev"),
        source=str(source),
    )


def mu_e_conversion_has_coefficient_proxy(source: Any) -> bool:
    """Return true when ``source`` carries mu->e conversion coefficient data."""

    if isinstance(source, MuEConversionCoefficientProxyInput):
        return True
    if isinstance(source, Mapping):
        return any(key in source for key in _COEFFICIENT_MAPPING_KEYS)
    return any(hasattr(source, name) for name in _COEFFICIENT_ATTR_NAMES)


def mu_e_conversion_coefficients(source: Any) -> MuEConversionCoefficientProxyInput:
    """Coerce mapping/object coefficient proxy data into a typed view."""

    if isinstance(source, MuEConversionCoefficientProxyInput):
        return mu_e_conversion_proxy_input(
            g_lv_p=source.g_lv_p,
            g_lv_n=source.g_lv_n,
            g_rv_p=source.g_rv_p,
            g_rv_n=source.g_rv_n,
            g_ls_p=source.g_ls_p,
            g_ls_n=source.g_ls_n,
            g_rs_p=source.g_rs_p,
            g_rs_n=source.g_rs_n,
            dipole_amplitude_left=source.dipole_amplitude_left,
            dipole_amplitude_right=source.dipole_amplitude_right,
            m_kk_gev=source.m_kk_gev,
            source=source.source,
        )
    if isinstance(source, Mapping):
        return _coefficients_from_mapping(source)
    return _coefficients_from_object(source)


def zero_mu_e_conversion_coefficients(
    *,
    source: str = "zero mu-e conversion contact coefficients",
) -> MuEConversionCoefficientProxyInput:
    """Return explicit zero contact coefficients."""

    return mu_e_conversion_proxy_input(source=source)


def mu_e_conversion_from_components(
    *,
    dipole_parent_branching_fraction: float,
    coefficients: MuEConversionCoefficientProxyInput,
    conversion_rate_limit: float | None = None,
    nuclear_inputs: MuEConversionNuclearInputs | None = None,
    dipole_amplitude_left: complex | None = None,
    dipole_amplitude_right: complex | None = None,
) -> MuEConversionResult:
    """Evaluate normalized coherent mu->e conversion from components."""

    n = aluminum_nuclear_inputs() if nuclear_inputs is None else nuclear_inputs
    _validate_nuclear_inputs(n)
    dipole_br = _nonnegative_float(
        dipole_parent_branching_fraction,
        "dipole_parent_branching_fraction",
    )
    coeffs = coefficients
    explicit_left = _first_not_none(
        dipole_amplitude_left,
        coeffs.dipole_amplitude_left,
    )
    explicit_right = _first_not_none(
        dipole_amplitude_right,
        coeffs.dipole_amplitude_right,
    )

    scalar_left = complex(coeffs.g_ls_p * n.S_p + coeffs.g_ls_n * n.S_n)
    scalar_right = complex(coeffs.g_rs_p * n.S_p + coeffs.g_rs_n * n.S_n)
    vector_left = complex(coeffs.g_lv_p * n.V_p + coeffs.g_lv_n * n.V_n)
    vector_right = complex(coeffs.g_rv_p * n.V_p + coeffs.g_rv_n * n.V_n)
    contact_left = complex(scalar_left + vector_left)
    contact_right = complex(scalar_right + vector_right)
    contact_inner = float(abs(contact_left) ** 2 + abs(contact_right) ** 2)
    dipole_norm = math.sqrt(dipole_br / (384.0 * math.pi**2))
    dipole_inner = float(n.D**2 * dipole_norm**2)

    if explicit_left is not None or explicit_right is not None:
        a_left = _finite_complex(0.0j if explicit_left is None else explicit_left, "A_L")
        a_right = _finite_complex(0.0j if explicit_right is None else explicit_right, "A_R")
        left_nuclear = complex(a_right.conjugate() * n.D + contact_left)
        right_nuclear = complex(a_left.conjugate() * n.D + contact_right)
        total_inner = float(abs(left_nuclear) ** 2 + abs(right_nuclear) ** 2)
        dipole_inner = float(n.D**2 * (abs(a_left) ** 2 + abs(a_right) ** 2))
        interference_inner = float(total_inner - dipole_inner - contact_inner)
        lower_interference = interference_inner
        upper_interference = interference_inner
        lower_inner = total_inner
        upper_inner = total_inner
        verdict_inner = total_inner
        treatment = "explicit_chiral_dipole_amplitudes"
        phase_status = "resolved_explicit_chiral_dipole_amplitudes"
        verdict_branch = "explicit_phase"
    else:
        a_left = None
        a_right = None
        envelope_interference = float(
            2.0 * n.D * dipole_norm * math.sqrt(contact_inner)
        )
        lower_interference = -envelope_interference
        upper_interference = envelope_interference
        lower_inner = float(
            max(0.0, dipole_inner + contact_inner + lower_interference)
        )
        upper_inner = float(dipole_inner + contact_inner + upper_interference)
        verdict_inner = lower_inner
        interference_inner = lower_interference
        total_inner = verdict_inner
        left_nuclear = complex(contact_left)
        right_nuclear = complex(contact_right)
        if envelope_interference == 0.0:
            treatment = "not_applicable_zero_dipole_or_contact"
            phase_status = "not_applicable_zero_dipole_or_contact"
            verdict_branch = "single_component"
        else:
            treatment = "unknown_relative_phase_interval_NEEDS-HUMAN-PHYSICS"
            phase_status = "NEEDS-HUMAN-PHYSICS"
            verdict_branch = "lower_envelope_unknown_relative_phase"

    prefactor = 2.0 * G_F_GEV_MINUS2**2 * KKO_OVERLAP_DIMENSION_FACTOR_GEV5
    capture_gev = n.capture_rate_gev
    conversion_width = float(prefactor * total_inner)
    conversion_width_lower = float(prefactor * lower_inner)
    conversion_width_upper = float(prefactor * upper_inner)
    conversion_rate = float(conversion_width / capture_gev)
    conversion_rate_lower = float(conversion_width_lower / capture_gev)
    conversion_rate_upper = float(conversion_width_upper / capture_gev)
    dipole_component = float(prefactor * dipole_inner / capture_gev)
    contact_component = float(prefactor * contact_inner / capture_gev)
    scalar_component = float(
        prefactor * (abs(scalar_left) ** 2 + abs(scalar_right) ** 2) / capture_gev
    )
    vector_component = float(
        prefactor * (abs(vector_left) ** 2 + abs(vector_right) ** 2) / capture_gev
    )
    vector_scalar_interference = float(
        contact_component - scalar_component - vector_component
    )
    dipole_contact_interference = float(
        prefactor * interference_inner / capture_gev
    )
    dipole_contact_lower = float(prefactor * lower_interference / capture_gev)
    dipole_contact_upper = float(prefactor * upper_interference / capture_gev)

    limit = None if conversion_rate_limit is None else _bounded_probability(
        conversion_rate_limit,
        "conversion_rate_limit",
    )
    ratio = None if limit is None else float(conversion_rate / limit)
    ratio_lower = None if limit is None else float(conversion_rate_lower / limit)
    ratio_upper = None if limit is None else float(conversion_rate_upper / limit)
    passes = None if ratio is None else bool(ratio <= 1.0)
    diagnostics = {
        "model_label": MU_E_CONVERSION_MODEL_V1,
        "input_bundle": n.input_bundle,
        "operator_convention": MU_E_CONVERSION_OPERATOR_CONVENTION,
        "dipole_convention": MU_E_CONVERSION_DIPOLE_CONVENTION,
        "matching_assumption": MU_E_CONVERSION_PROXY_V1,
        "sm_conversion_rate": 0.0,
        "sm_lfv_policy": (
            "Charged-LFV coherent mu->e conversion has negligible SM rate for "
            "catalog purposes; the bound is applied to the pure NP prediction."
        ),
        "hbar_gev_s": float(HBAR_GEV_S),
        "g_fermi_gev_minus2": float(G_F_GEV_MINUS2),
        "muon_mass_gev": float(MUON_MASS_GEV),
        "kko_overlap_units": "m_mu^(5/2)",
        "kko_overlap_dimension_factor_gev5": float(
            KKO_OVERLAP_DIMENSION_FACTOR_GEV5
        ),
        "target": n.target,
        "isotope": n.isotope,
        "overlap_D": float(n.D),
        "overlap_V_p": float(n.V_p),
        "overlap_V_n": float(n.V_n),
        "overlap_S_p": float(n.S_p),
        "overlap_S_n": float(n.S_n),
        "capture_rate_s_inv": float(n.capture_rate_s_inv),
        "capture_rate_gev": float(capture_gev),
        "nuclear_input_source": n.source,
        "nuclear_input_source_url": n.source_url,
        "g_lv_p": complex(coeffs.g_lv_p),
        "g_lv_n": complex(coeffs.g_lv_n),
        "g_rv_p": complex(coeffs.g_rv_p),
        "g_rv_n": complex(coeffs.g_rv_n),
        "g_ls_p": complex(coeffs.g_ls_p),
        "g_ls_n": complex(coeffs.g_ls_n),
        "g_rs_p": complex(coeffs.g_rs_p),
        "g_rs_n": complex(coeffs.g_rs_n),
        "scalar_left_nuclear_amplitude": scalar_left,
        "scalar_right_nuclear_amplitude": scalar_right,
        "vector_left_nuclear_amplitude": vector_left,
        "vector_right_nuclear_amplitude": vector_right,
        "contact_left_nuclear_amplitude": contact_left,
        "contact_right_nuclear_amplitude": contact_right,
        "left_nuclear_amplitude": left_nuclear,
        "right_nuclear_amplitude": right_nuclear,
        "dipole_parent_branching_fraction": float(dipole_br),
        "dipole_amplitude_norm": float(dipole_norm),
        "dipole_amplitude_left": a_left,
        "dipole_amplitude_right": a_right,
        "conversion_rate_lower": float(conversion_rate_lower),
        "conversion_rate_upper": float(conversion_rate_upper),
        "conversion_rate_interval": (
            float(conversion_rate_lower),
            float(conversion_rate_upper),
        ),
        "conversion_rate_verdict_branch": verdict_branch,
        "conversion_width_lower_gev": float(conversion_width_lower),
        "conversion_width_upper_gev": float(conversion_width_upper),
        "dipole_component": float(dipole_component),
        "scalar_component": float(scalar_component),
        "vector_component": float(vector_component),
        "vector_scalar_interference_component": float(vector_scalar_interference),
        "contact_component": float(contact_component),
        "dipole_contact_interference_component": float(dipole_contact_interference),
        "dipole_contact_interference_lower": float(dipole_contact_lower),
        "dipole_contact_interference_upper": float(dipole_contact_upper),
        "dipole_contact_interference_treatment": treatment,
        "dipole_contact_relative_phase_status": phase_status,
        "dipole_contact_verdict_policy": (
            "explicit chiral A_L/A_R use the resolved phase; otherwise compare "
            "the lower conversion-rate envelope to the limit and report the "
            "upper envelope diagnostically."
        ),
        "ratio_to_limit_lower": ratio_lower,
        "ratio_to_limit_upper": ratio_upper,
        "proxy_source": coeffs.source,
        "m_kk_gev": None if coeffs.m_kk_gev is None else float(coeffs.m_kk_gev),
    }
    return MuEConversionResult(
        model_label=MU_E_CONVERSION_MODEL_V1,
        input_bundle=n.input_bundle,
        target=n.target,
        conversion_rate=conversion_rate,
        conversion_rate_lower=conversion_rate_lower,
        conversion_rate_upper=conversion_rate_upper,
        sm_conversion_rate=0.0,
        np_shift_conversion_rate=conversion_rate,
        conversion_width_gev=conversion_width,
        conversion_width_lower_gev=conversion_width_lower,
        conversion_width_upper_gev=conversion_width_upper,
        capture_rate_gev=float(capture_gev),
        dipole_component=dipole_component,
        scalar_component=scalar_component,
        vector_component=vector_component,
        vector_scalar_interference_component=vector_scalar_interference,
        contact_component=contact_component,
        dipole_contact_interference_component=dipole_contact_interference,
        dipole_contact_interference_lower=dipole_contact_lower,
        dipole_contact_interference_upper=dipole_contact_upper,
        dipole_contact_interference_treatment=treatment,
        dipole_parent_branching_fraction=float(dipole_br),
        dipole_amplitude_left=a_left,
        dipole_amplitude_right=a_right,
        left_nuclear_amplitude=left_nuclear,
        right_nuclear_amplitude=right_nuclear,
        ratio_to_limit=ratio,
        ratio_to_limit_lower=ratio_lower,
        ratio_to_limit_upper=ratio_upper,
        conversion_rate_limit=limit,
        passes=passes,
        coefficients=coeffs,
        nuclear_inputs=n,
        diagnostics=diagnostics,
    )


def _coefficients_from_mapping(
    mapping: Mapping[str, Any],
) -> MuEConversionCoefficientProxyInput:
    return mu_e_conversion_proxy_input(
        g_lv_p=_optional_complex_from_mapping(mapping, _G_LV_P_KEYS, "g_lv_p"),
        g_lv_n=_optional_complex_from_mapping(mapping, _G_LV_N_KEYS, "g_lv_n"),
        g_rv_p=_optional_complex_from_mapping(mapping, _G_RV_P_KEYS, "g_rv_p"),
        g_rv_n=_optional_complex_from_mapping(mapping, _G_RV_N_KEYS, "g_rv_n"),
        g_ls_p=_optional_complex_from_mapping(mapping, _G_LS_P_KEYS, "g_ls_p"),
        g_ls_n=_optional_complex_from_mapping(mapping, _G_LS_N_KEYS, "g_ls_n"),
        g_rs_p=_optional_complex_from_mapping(mapping, _G_RS_P_KEYS, "g_rs_p"),
        g_rs_n=_optional_complex_from_mapping(mapping, _G_RS_N_KEYS, "g_rs_n"),
        dipole_amplitude_left=_optional_complex_or_none_from_mapping(
            mapping,
            _DIPOLE_LEFT_KEYS,
            "dipole_amplitude_left",
        ),
        dipole_amplitude_right=_optional_complex_or_none_from_mapping(
            mapping,
            _DIPOLE_RIGHT_KEYS,
            "dipole_amplitude_right",
        ),
        m_kk_gev=_optional_positive_from_mapping(mapping, _M_KK_KEYS, "m_kk_gev"),
        source=str(mapping.get("source", "mapping mu-e conversion coefficient proxy")),
    )


def _coefficients_from_object(value: Any) -> MuEConversionCoefficientProxyInput:
    return mu_e_conversion_proxy_input(
        g_lv_p=_optional_complex_from_object(value, _G_LV_P_KEYS, "g_lv_p"),
        g_lv_n=_optional_complex_from_object(value, _G_LV_N_KEYS, "g_lv_n"),
        g_rv_p=_optional_complex_from_object(value, _G_RV_P_KEYS, "g_rv_p"),
        g_rv_n=_optional_complex_from_object(value, _G_RV_N_KEYS, "g_rv_n"),
        g_ls_p=_optional_complex_from_object(value, _G_LS_P_KEYS, "g_ls_p"),
        g_ls_n=_optional_complex_from_object(value, _G_LS_N_KEYS, "g_ls_n"),
        g_rs_p=_optional_complex_from_object(value, _G_RS_P_KEYS, "g_rs_p"),
        g_rs_n=_optional_complex_from_object(value, _G_RS_N_KEYS, "g_rs_n"),
        dipole_amplitude_left=_optional_complex_or_none_from_object(
            value,
            _DIPOLE_LEFT_KEYS,
            "dipole_amplitude_left",
        ),
        dipole_amplitude_right=_optional_complex_or_none_from_object(
            value,
            _DIPOLE_RIGHT_KEYS,
            "dipole_amplitude_right",
        ),
        m_kk_gev=_optional_positive_from_object(value, _M_KK_KEYS, "m_kk_gev"),
        source=str(getattr(value, "source", "object mu-e conversion coefficient proxy")),
    )


def _validate_nuclear_inputs(nuclear: MuEConversionNuclearInputs) -> None:
    for name in ("D", "V_p", "V_n", "S_p", "S_n", "capture_rate_s_inv"):
        _positive_float(getattr(nuclear, name), name)


def _first_not_none(*values: Any) -> Any:
    for value in values:
        if value is not None:
            return value
    return None


def _optional_complex_from_mapping(
    mapping: Mapping[str, Any],
    keys: tuple[str, ...],
    name: str,
) -> complex:
    value = _first_present_key(mapping, keys)
    return 0.0j if value is None else _finite_complex(value, name)


def _optional_complex_or_none_from_mapping(
    mapping: Mapping[str, Any],
    keys: tuple[str, ...],
    name: str,
) -> complex | None:
    value = _first_present_key(mapping, keys)
    return None if value is None else _finite_complex(value, name)


def _optional_positive_from_mapping(
    mapping: Mapping[str, Any],
    keys: tuple[str, ...],
    name: str,
) -> float | None:
    value = _first_present_key(mapping, keys)
    return None if value is None else _positive_float(value, name)


def _optional_complex_from_object(
    value: Any,
    names: tuple[str, ...],
    name: str,
) -> complex:
    item = _first_present_attr(value, names)
    return 0.0j if item is None else _finite_complex(item, name)


def _optional_complex_or_none_from_object(
    value: Any,
    names: tuple[str, ...],
    name: str,
) -> complex | None:
    item = _first_present_attr(value, names)
    return None if item is None else _finite_complex(item, name)


def _optional_positive_from_object(
    value: Any,
    names: tuple[str, ...],
    name: str,
) -> float | None:
    item = _first_present_attr(value, names)
    return None if item is None else _positive_float(item, name)


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


def _positive_float(value: Any, name: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return number


def _nonnegative_float(value: Any, name: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number < 0.0:
        raise ValueError(f"{name} must be nonnegative and finite")
    return number


def _bounded_probability(value: Any, name: str) -> float:
    number = float(value)
    if not math.isfinite(number) or not 0.0 < number < 1.0:
        raise ValueError(f"{name} must lie between zero and one")
    return number


def _finite_complex(value: Any, name: str) -> complex:
    number = complex(value)
    if not math.isfinite(number.real) or not math.isfinite(number.imag):
        raise ValueError(f"{name} must be finite")
    return number


def _optional_finite_complex(value: Any, name: str) -> complex | None:
    return None if value is None else _finite_complex(value, name)


_G_LV_P_KEYS = ("g_lv_p", "g_LV_p", "g_lv_proton", "g_LV_proton")
_G_LV_N_KEYS = ("g_lv_n", "g_LV_n", "g_lv_neutron", "g_LV_neutron")
_G_RV_P_KEYS = ("g_rv_p", "g_RV_p", "g_rv_proton", "g_RV_proton")
_G_RV_N_KEYS = ("g_rv_n", "g_RV_n", "g_rv_neutron", "g_RV_neutron")
_G_LS_P_KEYS = ("g_ls_p", "g_LS_p", "g_ls_proton", "g_LS_proton")
_G_LS_N_KEYS = ("g_ls_n", "g_LS_n", "g_ls_neutron", "g_LS_neutron")
_G_RS_P_KEYS = ("g_rs_p", "g_RS_p", "g_rs_proton", "g_RS_proton")
_G_RS_N_KEYS = ("g_rs_n", "g_RS_n", "g_rs_neutron", "g_RS_neutron")
_DIPOLE_LEFT_KEYS = ("dipole_amplitude_left", "dipole_A_L", "a_left", "A_L")
_DIPOLE_RIGHT_KEYS = ("dipole_amplitude_right", "dipole_A_R", "a_right", "A_R")
_M_KK_KEYS = ("m_kk_gev", "M_KK_gev", "M_KK")

_COEFFICIENT_MAPPING_KEYS = frozenset(
    {
        *_G_LV_P_KEYS,
        *_G_LV_N_KEYS,
        *_G_RV_P_KEYS,
        *_G_RV_N_KEYS,
        *_G_LS_P_KEYS,
        *_G_LS_N_KEYS,
        *_G_RS_P_KEYS,
        *_G_RS_N_KEYS,
        *_DIPOLE_LEFT_KEYS,
        *_DIPOLE_RIGHT_KEYS,
    }
)
_COEFFICIENT_ATTR_NAMES = tuple(_COEFFICIENT_MAPPING_KEYS)


__all__ = [
    "MU_E_CONVERSION_MODEL_V1",
    "MU_E_CONVERSION_INPUT_BUNDLE_V1",
    "MU_E_CONVERSION_OPERATOR_CONVENTION",
    "MU_E_CONVERSION_DIPOLE_CONVENTION",
    "MU_E_CONVERSION_PROXY_V1",
    "HBAR_GEV_S",
    "G_F_GEV_MINUS2",
    "MUON_MASS_GEV",
    "KKO_OVERLAP_DIMENSION_FACTOR_GEV5",
    "MuEConversionNuclearInputs",
    "MuEConversionCoefficientProxyInput",
    "MuEConversionResult",
    "aluminum_nuclear_inputs",
    "titanium_nuclear_inputs",
    "gold_nuclear_inputs",
    "nuclear_inputs_for_target",
    "mu_e_conversion_proxy_input",
    "mu_e_conversion_has_coefficient_proxy",
    "mu_e_conversion_coefficients",
    "zero_mu_e_conversion_coefficients",
    "mu_e_conversion_from_components",
]
