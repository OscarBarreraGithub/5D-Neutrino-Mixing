"""Adapter over :mod:`quarkConstraints.leptonic_tree`.

This is the catalog boundary for charged-current leptonic meson decays.
Constraint modules import this adapter only; the reusable tree-level formula
and the documented charged-current proxy remain isolated in ``quarkConstraints``.
"""

from __future__ import annotations

from quarkConstraints.leptonic_tree import (
    LEPTONIC_TREE_BPLUS_TAUNU_YAML_INPUT_BUNDLE_V1,
    LEPTONIC_TREE_CHARGED_CURRENT_PROXY_ASSUMPTION_V1,
    LEPTONIC_TREE_INPUT_BUNDLE_V1,
    LEPTONIC_TREE_KPLUS_EMUNU_YAML_INPUT_BUNDLE_V1,
    LEPTONIC_TREE_LFU_RATIO_OPERATOR_CONVENTION,
    LEPTONIC_TREE_LFU_RATIO_PROXY_ASSUMPTION_V1,
    LEPTONIC_TREE_MODEL_V1,
    LEPTONIC_TREE_OPERATOR_CONVENTION,
    LeptonicTreeBranchingResult,
    LeptonicTreeInputs,
    LeptonicTreeRatioInputs,
    LeptonicTreeRatioResult,
)
from quarkConstraints.leptonic_tree import (
    bplus_tau_nu_inputs_from_theory_anchors as _bplus_tau_nu_inputs_from_theory_anchors,
)
from quarkConstraints.leptonic_tree import (
    charged_current_proxy_amplitude_ratio as _charged_current_proxy_amplitude_ratio,
)
from quarkConstraints.leptonic_tree import (
    default_bplus_tau_nu_inputs as _default_bplus_tau_nu_inputs,
)
from quarkConstraints.leptonic_tree import (
    default_kplus_enu_over_munu_inputs as _default_kplus_enu_over_munu_inputs,
)
from quarkConstraints.leptonic_tree import (
    evaluate_leptonic_branching_fraction as _evaluate_leptonic_branching_fraction,
)
from quarkConstraints.leptonic_tree import (
    evaluate_leptonic_lfu_ratio as _evaluate_leptonic_lfu_ratio,
)
from quarkConstraints.leptonic_tree import (
    kplus_enu_over_munu_inputs_from_sm_ratio_anchor as _kplus_enu_over_munu_inputs_from_sm_ratio_anchor,  # noqa: E501
)
from quarkConstraints.leptonic_tree import (
    leptonic_lfu_tree_ratio_without_radiation as _leptonic_lfu_tree_ratio_without_radiation,
)
from quarkConstraints.leptonic_tree import (
    lifetime_ps_to_gev_inverse as _lifetime_ps_to_gev_inverse,
)
from quarkConstraints.leptonic_tree import (
    sm_branching_fraction as _sm_branching_fraction,
)
from quarkConstraints.leptonic_tree import (
    sm_leptonic_lfu_ratio as _sm_leptonic_lfu_ratio,
)

__all__ = [
    "LEPTONIC_TREE_MODEL_V1",
    "LEPTONIC_TREE_INPUT_BUNDLE_V1",
    "LEPTONIC_TREE_BPLUS_TAUNU_YAML_INPUT_BUNDLE_V1",
    "LEPTONIC_TREE_KPLUS_EMUNU_YAML_INPUT_BUNDLE_V1",
    "LEPTONIC_TREE_OPERATOR_CONVENTION",
    "LEPTONIC_TREE_LFU_RATIO_OPERATOR_CONVENTION",
    "LEPTONIC_TREE_CHARGED_CURRENT_PROXY_ASSUMPTION_V1",
    "LEPTONIC_TREE_LFU_RATIO_PROXY_ASSUMPTION_V1",
    "LeptonicTreeInputs",
    "LeptonicTreeBranchingResult",
    "LeptonicTreeRatioInputs",
    "LeptonicTreeRatioResult",
    "leptonic_tree_default_bplus_tau_nu_inputs",
    "leptonic_tree_bplus_tau_nu_inputs_from_theory_anchors",
    "leptonic_tree_default_kplus_enu_over_munu_inputs",
    "leptonic_tree_kplus_enu_over_munu_inputs_from_sm_ratio_anchor",
    "leptonic_tree_lifetime_ps_to_gev_inverse",
    "leptonic_tree_sm_branching_fraction",
    "leptonic_tree_lfu_tree_ratio_without_radiation",
    "leptonic_tree_sm_lfu_ratio",
    "leptonic_tree_charged_current_proxy_amplitude_ratio",
    "bplus_tau_nu_branching_fraction",
    "kplus_enu_over_munu_lfu_ratio",
]


def leptonic_tree_default_bplus_tau_nu_inputs() -> LeptonicTreeInputs:
    """Return the default ``B+ -> tau+ nu_tau`` input bundle."""
    return _default_bplus_tau_nu_inputs()


def leptonic_tree_bplus_tau_nu_inputs_from_theory_anchors(
    *,
    decay_constant_gev: float,
    ckm_abs: float,
    constants_citation: str,
) -> LeptonicTreeInputs:
    """Return ``B+ -> tau+ nu_tau`` inputs with YAML-backed theory values."""
    return _bplus_tau_nu_inputs_from_theory_anchors(
        decay_constant_gev=decay_constant_gev,
        ckm_abs=ckm_abs,
        constants_citation=constants_citation,
    )


def leptonic_tree_default_kplus_enu_over_munu_inputs() -> LeptonicTreeRatioInputs:
    """Return default ``R_K = Gamma(K -> e nu)/Gamma(K -> mu nu)`` inputs."""
    return _default_kplus_enu_over_munu_inputs()


def leptonic_tree_kplus_enu_over_munu_inputs_from_sm_ratio_anchor(
    *,
    sm_ratio: float,
    constants_citation: str,
) -> LeptonicTreeRatioInputs:
    """Return K LFU ratio inputs with radiative effects fixed by an SM anchor."""
    return _kplus_enu_over_munu_inputs_from_sm_ratio_anchor(
        sm_ratio=sm_ratio,
        constants_citation=constants_citation,
    )


def leptonic_tree_lifetime_ps_to_gev_inverse(
    tau_ps: float,
    hbar_gev_s: float,
) -> float:
    """Convert a lifetime in ps to natural units."""
    return _lifetime_ps_to_gev_inverse(tau_ps, hbar_gev_s)


def leptonic_tree_sm_branching_fraction(
    inputs: LeptonicTreeInputs | None = None,
) -> LeptonicTreeBranchingResult:
    """Evaluate the SM tree-level branching fraction."""
    return _sm_branching_fraction(inputs)


def leptonic_tree_lfu_tree_ratio_without_radiation(
    inputs: LeptonicTreeRatioInputs | None = None,
) -> float:
    """Return the tree-level LFU ratio before the radiative multiplier."""
    return _leptonic_lfu_tree_ratio_without_radiation(inputs)


def leptonic_tree_sm_lfu_ratio(
    inputs: LeptonicTreeRatioInputs | None = None,
) -> LeptonicTreeRatioResult:
    """Evaluate the SM charged-leptonic LFU ratio."""
    return _sm_leptonic_lfu_ratio(inputs)


def leptonic_tree_charged_current_proxy_amplitude_ratio(
    *,
    meson_mass_gev: float,
    m_kk_gev: float,
    proxy_strength: float = 1.0,
) -> complex:
    """Return the documented v1 charged-current proxy amplitude."""
    return _charged_current_proxy_amplitude_ratio(
        meson_mass_gev=meson_mass_gev,
        m_kk_gev=m_kk_gev,
        proxy_strength=proxy_strength,
    )


def bplus_tau_nu_branching_fraction(
    *,
    m_kk_gev: float | None = None,
    np_amplitude_ratio: complex | None = None,
    inputs: LeptonicTreeInputs | None = None,
) -> LeptonicTreeBranchingResult:
    """Evaluate ``BR(B+ -> tau+ nu_tau)`` in the SM or proxy model."""
    return _evaluate_leptonic_branching_fraction(
        m_kk_gev=m_kk_gev,
        np_amplitude_ratio=np_amplitude_ratio,
        inputs=inputs,
    )


def kplus_enu_over_munu_lfu_ratio(
    *,
    m_kk_gev: float | None = None,
    numerator_np_amplitude_ratio: complex | None = None,
    denominator_np_amplitude_ratio: complex | None = None,
    inputs: LeptonicTreeRatioInputs | None = None,
) -> LeptonicTreeRatioResult:
    """Evaluate ``R_K`` in the SM or documented lepton-nonuniversal proxy."""
    return _evaluate_leptonic_lfu_ratio(
        m_kk_gev=m_kk_gev,
        numerator_np_amplitude_ratio=numerator_np_amplitude_ratio,
        denominator_np_amplitude_ratio=denominator_np_amplitude_ratio,
        inputs=inputs,
    )
