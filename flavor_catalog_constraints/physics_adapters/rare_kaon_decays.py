"""Placeholder adapters for rare kaon decays.

The current physics core contains Delta F = 2 kaon-mixing machinery,
but no Delta S = 1 semileptonic ``s -> d nu nubar`` Wilson matching or
branching-ratio evaluator. This adapter makes that absence explicit at
the catalog-constraint boundary so rare kaon modes can register without
silently reusing unrelated kaon-mixing or eps'/eps machinery.
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = [
    "KPLUS_PIPLUS_NUNU_DEFERRED_REASON",
    "RareKaonDecayPlaceholderResult",
    "evaluate_kplus_piplus_nunu_placeholder",
]


KPLUS_PIPLUS_NUNU_DEFERRED_REASON = (
    "requires Delta S=1 K -> pi nu nubar Wilson matching, CKM inputs, "
    "SM-NP interference conventions, and branching-ratio machinery not "
    "yet in the physics core"
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
