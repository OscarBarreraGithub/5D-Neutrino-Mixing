"""Placeholder adapters for charm-sector observables."""

from __future__ import annotations
from dataclasses import dataclass

__all__ = ["CharmPlaceholderResult", "evaluate_charm_placeholder"]

_REASON = "requires charm-sector Delta C=1/2 Wilson matching + 5D RS NP; not yet in core"


@dataclass(frozen=True)
class CharmPlaceholderResult:
    passes: bool = True
    predicted: float | None = None
    deferred: bool = True
    reason: str = _REASON


def evaluate_charm_placeholder(channel: str = "charm mixing") -> CharmPlaceholderResult:
    return CharmPlaceholderResult(reason=f"{channel}: {_REASON}")
