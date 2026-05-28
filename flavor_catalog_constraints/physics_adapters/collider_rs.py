"""Placeholder adapters for collider-direct RS-search bounds."""

from __future__ import annotations
from dataclasses import dataclass

__all__ = ["ColliderRsPlaceholderResult", "evaluate_collider_rs_placeholder"]

_REASON = "requires LHC direct-search recasting + 5D RS mass-coupling map; not yet in core"


@dataclass(frozen=True)
class ColliderRsPlaceholderResult:
    passes: bool = True
    predicted: float | None = None
    deferred: bool = True
    reason: str = _REASON


def evaluate_collider_rs_placeholder(channel: str = "collider RS") -> ColliderRsPlaceholderResult:
    return ColliderRsPlaceholderResult(reason=f"{channel}: {_REASON}")
