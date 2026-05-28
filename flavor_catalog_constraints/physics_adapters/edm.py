"""Placeholder adapters for EDM and CP-violating-dipole constraints.

The current physics core does not implement CP-violating dipole
operators, quark/lepton EDM matching, chromo-EDM, or the Weinberg
three-gluon operator. This adapter makes that absence explicit at the
catalog-constraint boundary so EDM observables can register without
silently reusing unrelated machinery.
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = [
    "EDM_DEFERRED_REASON_LEPTON",
    "EDM_DEFERRED_REASON_HADRONIC",
    "EDM_DEFERRED_REASON_WEINBERG",
    "EdmPlaceholderResult",
    "evaluate_lepton_edm_placeholder",
    "evaluate_hadronic_edm_placeholder",
    "evaluate_chromo_edm_placeholder",
    "evaluate_weinberg_placeholder",
]


EDM_DEFERRED_REASON_LEPTON = (
    "requires CP-violating dipole Wilson matching for lepton EDMs, plus "
    "5D RS one-loop Higgs/KK-gauge contributions to imaginary parts of "
    "dipole operators, not yet in the physics core"
)

EDM_DEFERRED_REASON_HADRONIC = (
    "requires CP-violating dipole + chromo-EDM Wilson matching for "
    "quark EDMs, hadronic matrix elements (Schiff moment for atomic, "
    "QCD sum rules for neutron), and Peccei-Quinn / theta-bar treatment, "
    "not yet in the physics core"
)

EDM_DEFERRED_REASON_WEINBERG = (
    "requires Weinberg three-gluon operator coefficient C_GG~ + "
    "hadronic matrix elements for neutron EDM extraction, "
    "not yet in the physics core"
)


@dataclass(frozen=True)
class EdmPlaceholderResult:
    """Non-vetoing placeholder result for deferred EDM evaluation."""

    passes: bool = True
    predicted: float | None = None
    deferred: bool = True
    reason: str = EDM_DEFERRED_REASON_LEPTON


def evaluate_lepton_edm_placeholder(species: str = "electron") -> EdmPlaceholderResult:
    """Return an explicit deferred result for a lepton EDM (e, mu, tau)."""
    return EdmPlaceholderResult(
        reason=f"{species} EDM: {EDM_DEFERRED_REASON_LEPTON}"
    )


def evaluate_hadronic_edm_placeholder(species: str = "neutron") -> EdmPlaceholderResult:
    """Return an explicit deferred result for a hadronic/atomic EDM."""
    return EdmPlaceholderResult(
        reason=f"{species} EDM: {EDM_DEFERRED_REASON_HADRONIC}"
    )


def evaluate_chromo_edm_placeholder() -> EdmPlaceholderResult:
    """Return an explicit deferred result for quark chromo-EDM bounds."""
    return EdmPlaceholderResult(reason=EDM_DEFERRED_REASON_HADRONIC + " (chromo-EDM)")


def evaluate_weinberg_placeholder() -> EdmPlaceholderResult:
    """Return an explicit deferred result for the Weinberg 3-gluon operator."""
    return EdmPlaceholderResult(reason=EDM_DEFERRED_REASON_WEINBERG)
