"""Placeholder adapters for rare kaon decays.

The current physics core contains Delta F = 2 kaon-mixing machinery,
but no Delta S = 1 semileptonic ``s -> d`` Wilson matching or
branching-ratio evaluator. This adapter makes that absence explicit at
the catalog-constraint boundary so rare kaon modes can register without
silently reusing unrelated kaon-mixing or eps'/eps machinery.
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = [
    "KLONG_MUMU_DEFERRED_REASON",
    "KLONG_PI0EE_DEFERRED_REASON",
    "KLONG_PI0MUMU_DEFERRED_REASON",
    "KLONG_PI0_NUNU_DEFERRED_REASON",
    "KSHORT_PI0EE_DEFERRED_REASON",
    "KPLUS_PIPLUS_NUNU_DEFERRED_REASON",
    "RareKaonDecayPlaceholderResult",
    "evaluate_klong_mumu_placeholder",
    "evaluate_klong_pi0ee_placeholder",
    "evaluate_klong_pi0mumu_placeholder",
    "evaluate_klong_pi0_nunu_placeholder",
    "evaluate_kshort_pi0ee_placeholder",
    "evaluate_kplus_piplus_nunu_placeholder",
]


KPLUS_PIPLUS_NUNU_DEFERRED_REASON = (
    "requires Delta S=1 K -> pi nu nubar Wilson matching, CKM inputs, "
    "SM-NP interference conventions, and branching-ratio machinery not "
    "yet in the physics core"
)

KLONG_PI0_NUNU_DEFERRED_REASON = (
    "requires Delta S=1 K -> pi nu nubar Wilson matching, CKM inputs, "
    "SM-NP interference conventions, CP-odd K_L projection, and "
    "branching-ratio machinery not yet in the physics core"
)

KLONG_MUMU_DEFERRED_REASON = (
    "requires Delta S=1 s -> d mu+ mu- short-distance Wilson matching, "
    "CKM inputs, SM-NP interference conventions, and long-distance "
    "two-photon treatment or a documented short-distance subtraction "
    "not yet in the physics core"
)

KLONG_PI0EE_DEFERRED_REASON = (
    "requires Delta S=1 s -> d e+ e- Wilson matching, CKM inputs, "
    "SM-NP interference conventions, indirect-CP input from "
    "K_S -> pi0 e+ e-, and CP-conserving K_L -> pi0 gamma gamma "
    "long-distance treatment not yet in the physics core"
)

KSHORT_PI0EE_DEFERRED_REASON = (
    "requires Delta S=1 s -> d e+ e- Wilson matching, CKM inputs, "
    "SM-NP interference conventions, K_S -> pi0 e+ e- form-factor "
    "treatment, phase-space extrapolation, and branching-ratio machinery "
    "not yet in the physics core"
)

KLONG_PI0MUMU_DEFERRED_REASON = (
    "requires Delta S=1 s -> d mu+ mu- Wilson matching, CKM inputs, "
    "SM-NP interference conventions, indirect-CP input from "
    "K_S -> pi0 mu+ mu-, and CP-conserving photon-mediated "
    "K_L -> pi0 gamma gamma long-distance treatment not yet in the physics core"
)


@dataclass(frozen=True)
class RareKaonDecayPlaceholderResult:
    """Non-vetoing placeholder result for deferred rare-kaon evaluation."""

    passes: bool = True
    predicted: float | None = None
    deferred: bool = True
    reason: str = KPLUS_PIPLUS_NUNU_DEFERRED_REASON


def evaluate_kplus_piplus_nunu_placeholder() -> RareKaonDecayPlaceholderResult:
    """Return an explicit deferred result until rare-kaon support exists."""
    return RareKaonDecayPlaceholderResult()


def evaluate_klong_pi0_nunu_placeholder() -> RareKaonDecayPlaceholderResult:
    """Return an explicit deferred result until neutral rare-kaon support exists."""
    return RareKaonDecayPlaceholderResult(reason=KLONG_PI0_NUNU_DEFERRED_REASON)


def evaluate_klong_mumu_placeholder() -> RareKaonDecayPlaceholderResult:
    """Return an explicit deferred result until K_L -> mu mu support exists."""
    return RareKaonDecayPlaceholderResult(reason=KLONG_MUMU_DEFERRED_REASON)


def evaluate_klong_pi0ee_placeholder() -> RareKaonDecayPlaceholderResult:
    """Return an explicit deferred result until K_L -> pi0 ee support exists."""
    return RareKaonDecayPlaceholderResult(reason=KLONG_PI0EE_DEFERRED_REASON)


def evaluate_kshort_pi0ee_placeholder() -> RareKaonDecayPlaceholderResult:
    """Return an explicit deferred result until K_S -> pi0 ee support exists."""
    return RareKaonDecayPlaceholderResult(reason=KSHORT_PI0EE_DEFERRED_REASON)


def evaluate_klong_pi0mumu_placeholder() -> RareKaonDecayPlaceholderResult:
    """Return an explicit deferred result until K_L -> pi0 mumu support exists."""
    return RareKaonDecayPlaceholderResult(reason=KLONG_PI0MUMU_DEFERRED_REASON)
