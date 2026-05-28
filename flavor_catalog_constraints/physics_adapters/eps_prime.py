"""Placeholder adapter for K003, ``Re(epsilon'/epsilon)``.

The current physics core contains Delta F = 2 kaon-mixing machinery,
but no Delta S = 1 Wilson-coefficient evolution or K -> pi pi matrix
element evaluator. This adapter makes that absence explicit at the
catalog-constraint boundary so K003 can register without silently
reusing the wrong physics.
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = [
    "EPS_PRIME_DEFERRED_REASON",
    "EpsPrimePlaceholderResult",
    "evaluate_eps_prime_placeholder",
]


EPS_PRIME_DEFERRED_REASON = (
    "requires Delta S=1 Wilson-coefficient machinery not yet in the physics core"
)


@dataclass(frozen=True)
class EpsPrimePlaceholderResult:
    """Non-vetoing placeholder result for deferred eps'/eps evaluation."""

    passes: bool = True
    predicted: float | None = None
    deferred: bool = True
    reason: str = EPS_PRIME_DEFERRED_REASON


def evaluate_eps_prime_placeholder() -> EpsPrimePlaceholderResult:
    """Return an explicit deferred result until Delta S=1 support exists."""
    return EpsPrimePlaceholderResult()
