"""Adapter over :mod:`quarkConstraints.deltaf2` (the Delta F = 2 core).

This is the import boundary the kaon / beauty / charm mixing constraints
use to reach the Delta F = 2 physics. Per the append-only convention,
constraints may *add* wrappers here but must not change the signature of
an existing one; an upstream signature change is absorbed inside the
wrapper body so callers see a stable surface.

The upstream result dataclasses (``EpsilonKResult``, ``DeltaMKResult``,
``MesonMixingResult``) are re-exported unchanged.
"""

from __future__ import annotations

from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.deltaf2 import (
    DEFAULT_DELTA_F2_INPUTS_V1,
    DELTA_M_BD_EXP,
    DELTA_M_BD_SM,
    DELTA_M_K,
    DELTA_M_BS_EXP,
    DELTA_M_BS_SM,
    B_1_BS,
    B_1_BD,
    B_1_D,
    B_4_BS,
    B_4_BD,
    B_4_D,
    B_5_BS,
    B_5_BD,
    B_5_D,
    DeltaF2WilsonCoefficients,
    DeltaMKResult,
    EpsilonKBudgetPolicy,
    EpsilonKResult,
    F_BS,
    F_BD,
    F_D,
    M_BS,
    M_BD,
    M_B_QUARK,
    M_C_QUARK,
    M_D_QUARK_BD,
    M_D0,
    M_S_QUARK_BS,
    M_U_QUARK,
    MesonMixingResult,
    compute_m12_np as _compute_meson_m12_np,
    compute_delta_f2_wilsons,
    delta_f2_epsilon_k_budget_policy as _epsilon_k_budget_policy,
    evaluate_bd_mixing_with_running as _evaluate_bd_mixing_with_running,
    evaluate_bs_mixing_with_running as _evaluate_bs_mixing_with_running,
    evaluate_delta_mk as _evaluate_delta_mk,
    evaluate_delta_mk_with_running as _evaluate_delta_mk_with_running,
    evaluate_d0_mixing_with_running as _evaluate_d0_mixing_with_running,
    evaluate_epsilon_k as _evaluate_epsilon_k,
    evaluate_epsilon_k_with_running as _evaluate_epsilon_k_with_running,
    _evolve_wilsons as _evolve_delta_f2_wilsons,
)

__all__ = [
    "QuarkMassBasisCouplings",
    "DeltaF2WilsonCoefficients",
    "EpsilonKBudgetPolicy",
    "EpsilonKResult",
    "DeltaMKResult",
    "MesonMixingResult",
    "epsilon_k_wilsons_from_couplings",
    "delta_mk_wilsons_from_couplings",
    "d0_mixing_wilsons_from_couplings",
    "bd_mixing_wilsons_from_couplings",
    "bs_mixing_wilsons_from_couplings",
    "epsilon_k_from_wilsons",
    "epsilon_k_from_wilsons_with_running",
    "epsilon_k_budget_policy",
    "delta_mk_from_wilsons_with_running",
    "delta_mk_core_inputs",
    "d0_mixing_from_wilsons_with_running",
    "d0_mixing_m12_np_from_wilsons_with_running",
    "bd_mixing_from_wilsons_with_running",
    "bd_mixing_m12_np_from_wilsons_with_running",
    "bd_mixing_core_inputs",
    "bs_mixing_from_wilsons_with_running",
    "bs_mixing_m12_np_from_wilsons_with_running",
    "bs_mixing_core_inputs",
    "epsilon_k_from_couplings",
    "delta_mk_from_couplings",
    "delta_mk_from_couplings_with_running",
]


def _kaon_wilsons(couplings: QuarkMassBasisCouplings) -> DeltaF2WilsonCoefficients:
    """Return the kaon Delta F=2 Wilson coefficients for ``couplings``.

    Centralizes the "pick the epsilon_k entry out of the default input
    bundle" step so both wrappers below share it. If upstream renames
    the bundle key, this is the one place to update.
    """
    wilsons = compute_delta_f2_wilsons(couplings, inputs=DEFAULT_DELTA_F2_INPUTS_V1)
    for w in wilsons:
        if w.input.key == "epsilon_k":
            return w
    raise RuntimeError("epsilon_k entry missing from DEFAULT_DELTA_F2_INPUTS_V1")


def _d0_wilsons(couplings: QuarkMassBasisCouplings) -> DeltaF2WilsonCoefficients:
    """Return the D0 Delta F=2 Wilson coefficients for ``couplings``."""
    wilsons = compute_delta_f2_wilsons(couplings, inputs=DEFAULT_DELTA_F2_INPUTS_V1)
    for w in wilsons:
        if w.input.key == "d":
            return w
    raise RuntimeError("d entry missing from DEFAULT_DELTA_F2_INPUTS_V1")


def _bd_wilsons(couplings: QuarkMassBasisCouplings) -> DeltaF2WilsonCoefficients:
    """Return the B_d Delta F=2 Wilson coefficients for ``couplings``."""
    wilsons = compute_delta_f2_wilsons(couplings, inputs=DEFAULT_DELTA_F2_INPUTS_V1)
    for w in wilsons:
        if w.input.key == "b_d":
            return w
    raise RuntimeError("b_d entry missing from DEFAULT_DELTA_F2_INPUTS_V1")


def _bs_wilsons(couplings: QuarkMassBasisCouplings) -> DeltaF2WilsonCoefficients:
    """Return the B_s Delta F=2 Wilson coefficients for ``couplings``."""
    wilsons = compute_delta_f2_wilsons(couplings, inputs=DEFAULT_DELTA_F2_INPUTS_V1)
    for w in wilsons:
        if w.input.key == "b_s":
            return w
    raise RuntimeError("b_s entry missing from DEFAULT_DELTA_F2_INPUTS_V1")


def epsilon_k_from_couplings(couplings: QuarkMassBasisCouplings) -> EpsilonKResult:
    """Compute the NP contribution to epsilon_K from mass-basis couplings."""
    return _evaluate_epsilon_k(_kaon_wilsons(couplings))


def delta_mk_from_couplings(couplings: QuarkMassBasisCouplings) -> DeltaMKResult:
    """Compute the NP contribution to Delta m_K from mass-basis couplings."""
    return _evaluate_delta_mk(_kaon_wilsons(couplings))


def delta_mk_from_couplings_with_running(
    couplings: QuarkMassBasisCouplings,
    *,
    mu_had: float = 2.0,
    m12_np_budget: float | None = None,
) -> DeltaMKResult:
    """Compute Delta m_K after QCD-evolving kaon Wilsons to ``mu_had``."""
    return delta_mk_from_wilsons_with_running(
        _kaon_wilsons(couplings),
        mu_had=mu_had,
        m12_np_budget=m12_np_budget,
    )


def epsilon_k_wilsons_from_couplings(
    couplings: QuarkMassBasisCouplings,
) -> DeltaF2WilsonCoefficients:
    """Return the kaon Delta F=2 Wilson coefficients used for epsilon_K."""
    return _kaon_wilsons(couplings)


def delta_mk_wilsons_from_couplings(
    couplings: QuarkMassBasisCouplings,
) -> DeltaF2WilsonCoefficients:
    """Return the kaon Delta F=2 Wilson coefficients used for Delta m_K."""
    return _kaon_wilsons(couplings)


def d0_mixing_wilsons_from_couplings(
    couplings: QuarkMassBasisCouplings,
) -> DeltaF2WilsonCoefficients:
    """Return the D0 Delta F=2 Wilson coefficients used for neutral-D mixing."""
    return _d0_wilsons(couplings)


def bd_mixing_wilsons_from_couplings(
    couplings: QuarkMassBasisCouplings,
) -> DeltaF2WilsonCoefficients:
    """Return the B_d Delta F=2 Wilson coefficients used for neutral-B_d mixing."""
    return _bd_wilsons(couplings)


def bs_mixing_wilsons_from_couplings(
    couplings: QuarkMassBasisCouplings,
) -> DeltaF2WilsonCoefficients:
    """Return the B_s Delta F=2 Wilson coefficients used for neutral-B_s mixing."""
    return _bs_wilsons(couplings)


def epsilon_k_from_wilsons(
    wilsons: DeltaF2WilsonCoefficients,
    *,
    epsilon_k_np_budget: float | None = None,
) -> EpsilonKResult:
    """Compute epsilon_K from kaon Wilsons, optionally with a catalog budget."""
    return _evaluate_epsilon_k(
        wilsons,
        epsilon_k_np_budget_override=epsilon_k_np_budget,
    )


def epsilon_k_from_wilsons_with_running(
    wilsons: DeltaF2WilsonCoefficients,
    *,
    mu_had: float = 2.0,
    epsilon_k_np_budget: float | None = None,
) -> EpsilonKResult:
    """Compute epsilon_K after QCD-evolving kaon Wilsons to ``mu_had``."""
    return _evaluate_epsilon_k_with_running(
        wilsons,
        mu_had=mu_had,
        epsilon_k_np_budget_override=epsilon_k_np_budget,
    )


def epsilon_k_budget_policy() -> EpsilonKBudgetPolicy:
    """Return the shared core epsilon_K budget policy."""

    return _epsilon_k_budget_policy()


def delta_mk_core_inputs() -> dict[str, float]:
    """Return the Delta m_K mass-splitting input hardwired in the Delta F=2 core."""
    return {
        "delta_m_k_exp_gev": float(DELTA_M_K),
        "core_m12_budget_gev": float(DELTA_M_K / 2.0),
    }


def delta_mk_from_wilsons_with_running(
    wilsons: DeltaF2WilsonCoefficients,
    *,
    mu_had: float = 2.0,
    m12_np_budget: float | None = None,
) -> DeltaMKResult:
    """Compute Delta m_K after QCD-evolving kaon Wilsons to ``mu_had``.

    ``quarkConstraints.deltaf2.evaluate_delta_mk_with_running`` supplies the
    audited running and kaon matrix-element path.  The optional budget override
    lets catalog constraints use their YAML-loaded experimental anchor while
    keeping the core prediction for ``|M12^NP|``.
    """
    result = _evaluate_delta_mk_with_running(wilsons, mu_had=mu_had)
    if m12_np_budget is None:
        return result
    if m12_np_budget <= 0.0:
        raise ValueError("m12_np_budget must be positive")
    budget = float(m12_np_budget)
    abs_m12_np = float(result.abs_m12_np)
    ratio = abs_m12_np / budget
    return DeltaMKResult(
        abs_m12_np=abs_m12_np,
        ratio_to_exp=ratio,
        passes=ratio <= 1.0,
    )


def d0_mixing_from_wilsons_with_running(
    wilsons: DeltaF2WilsonCoefficients,
    *,
    mu_had: float = 2.0,
    m12_np_budget: float | None = None,
) -> MesonMixingResult:
    """Compute D0 mixing after QCD-evolving Wilsons to ``mu_had``.

    ``quarkConstraints.deltaf2.evaluate_d0_mixing_with_running`` carries the
    audited running and matrix-element path.  The optional budget override lets
    catalog constraints use their YAML-loaded experimental anchor while keeping
    the core convention ``|M12^NP| <= Delta m_D^exp / 2``.
    """
    result = _evaluate_d0_mixing_with_running(wilsons, mu_had=mu_had)
    if m12_np_budget is None:
        return result
    if m12_np_budget <= 0.0:
        raise ValueError("m12_np_budget must be positive")
    budget = float(m12_np_budget)
    abs_m12_np = float(result.abs_m12_np)
    ratio = abs_m12_np / budget
    return MesonMixingResult(
        system=result.system,
        abs_m12_np=abs_m12_np,
        budget=budget,
        ratio_to_budget=ratio,
        passes=ratio <= 1.0,
    )


def d0_mixing_m12_np_from_wilsons_with_running(
    wilsons: DeltaF2WilsonCoefficients,
    *,
    mu_had: float = 2.0,
) -> complex:
    """Return complex ``M12^NP`` for D0 mixing after QCD running.

    This is the phase-preserving companion to
    :func:`d0_mixing_from_wilsons_with_running`.  It intentionally reuses the
    same Delta F=2 core evolution and D0 matrix-element helper that feed the
    audited magnitude evaluator, but does not collapse the result to
    ``abs(M12^NP)``.
    """
    evolved = _evolve_delta_f2_wilsons(wilsons, mu_had=mu_had)
    return complex(
        _compute_meson_m12_np(
            evolved,
            F_D,
            M_D0,
            M_C_QUARK,
            M_U_QUARK,
            B_1_D,
            B_4_D,
            B_5_D,
        )
    )


def bd_mixing_core_inputs() -> dict[str, float]:
    """Return the B_d mass-splitting inputs hardwired in the Delta F=2 core."""
    return {
        "delta_m_bd_exp_gev": float(DELTA_M_BD_EXP),
        "delta_m_bd_sm_gev": float(DELTA_M_BD_SM),
        "core_m12_budget_gev": float(
            max(DELTA_M_BD_EXP / 2.0, abs(DELTA_M_BD_EXP - DELTA_M_BD_SM) / 2.0)
        ),
    }


def bd_mixing_from_wilsons_with_running(
    wilsons: DeltaF2WilsonCoefficients,
    *,
    mu_had: float = 2.0,
    m12_np_budget: float | None = None,
) -> MesonMixingResult:
    """Compute B_d mixing after QCD-evolving Wilsons to ``mu_had``.

    ``quarkConstraints.deltaf2.evaluate_bd_mixing_with_running`` supplies the
    audited running and B_d matrix-element path.  The optional budget override
    lets catalog constraints apply an uncertainty-aware SM-vs-experiment room
    while preserving the core prediction for ``|M12^NP|``.
    """
    result = _evaluate_bd_mixing_with_running(wilsons, mu_had=mu_had)
    if m12_np_budget is None:
        return result
    if m12_np_budget <= 0.0:
        raise ValueError("m12_np_budget must be positive")
    budget = float(m12_np_budget)
    abs_m12_np = float(result.abs_m12_np)
    ratio = abs_m12_np / budget
    return MesonMixingResult(
        system=result.system,
        abs_m12_np=abs_m12_np,
        budget=budget,
        ratio_to_budget=ratio,
        passes=ratio <= 1.0,
    )


def bd_mixing_m12_np_from_wilsons_with_running(
    wilsons: DeltaF2WilsonCoefficients,
    *,
    mu_had: float = 2.0,
) -> complex:
    """Return complex ``M12^NP`` for B_d mixing after QCD running.

    This is the phase-preserving companion to
    :func:`bd_mixing_from_wilsons_with_running`.  It intentionally reuses the
    same Delta F=2 core evolution and B_d matrix-element helper that feed the
    audited magnitude evaluator, but does not collapse the result to
    ``abs(M12^NP)``.
    """
    evolved = _evolve_delta_f2_wilsons(wilsons, mu_had=mu_had)
    return complex(
        _compute_meson_m12_np(
            evolved,
            F_BD,
            M_BD,
            M_B_QUARK,
            M_D_QUARK_BD,
            B_1_BD,
            B_4_BD,
            B_5_BD,
        )
    )


def bs_mixing_core_inputs() -> dict[str, float]:
    """Return the B_s mass-splitting inputs hardwired in the Delta F=2 core."""
    return {
        "delta_m_bs_exp_gev": float(DELTA_M_BS_EXP),
        "delta_m_bs_sm_gev": float(DELTA_M_BS_SM),
        "core_m12_budget_gev": float(
            max(DELTA_M_BS_EXP / 2.0, abs(DELTA_M_BS_EXP - DELTA_M_BS_SM) / 2.0)
        ),
    }


def bs_mixing_from_wilsons_with_running(
    wilsons: DeltaF2WilsonCoefficients,
    *,
    mu_had: float = 2.0,
    m12_np_budget: float | None = None,
) -> MesonMixingResult:
    """Compute B_s mixing after QCD-evolving Wilsons to ``mu_had``.

    ``quarkConstraints.deltaf2.evaluate_bs_mixing_with_running`` supplies the
    audited running and B_s matrix-element path.  The optional budget override
    lets catalog constraints apply an uncertainty-aware SM-vs-experiment room
    while preserving the core prediction for ``|M12^NP|``.
    """
    result = _evaluate_bs_mixing_with_running(wilsons, mu_had=mu_had)
    if m12_np_budget is None:
        return result
    if m12_np_budget <= 0.0:
        raise ValueError("m12_np_budget must be positive")
    budget = float(m12_np_budget)
    abs_m12_np = float(result.abs_m12_np)
    ratio = abs_m12_np / budget
    return MesonMixingResult(
        system=result.system,
        abs_m12_np=abs_m12_np,
        budget=budget,
        ratio_to_budget=ratio,
        passes=ratio <= 1.0,
    )


def bs_mixing_m12_np_from_wilsons_with_running(
    wilsons: DeltaF2WilsonCoefficients,
    *,
    mu_had: float = 2.0,
) -> complex:
    """Return complex ``M12^NP`` for B_s mixing after QCD running.

    This is the phase-preserving companion to
    :func:`bs_mixing_from_wilsons_with_running`.  It intentionally reuses the
    same Delta F=2 core evolution and B_s matrix-element helper that feed the
    audited magnitude evaluator, but does not collapse the result to
    ``abs(M12^NP)``.
    """
    evolved = _evolve_delta_f2_wilsons(wilsons, mu_had=mu_had)
    return complex(
        _compute_meson_m12_np(
            evolved,
            F_BS,
            M_BS,
            M_B_QUARK,
            M_S_QUARK_BS,
            B_1_BS,
            B_4_BS,
            B_5_BS,
        )
    )
