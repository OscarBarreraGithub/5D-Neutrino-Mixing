"""Placeholder adapters for B-sector observables."""

from __future__ import annotations
from dataclasses import dataclass

__all__ = ["BeautyPlaceholderResult", "evaluate_beauty_placeholder"]

_REASON = "requires B-sector ΔB=1/2 Wilson matching + 5D RS NP; not yet in core"


@dataclass(frozen=True)
class BeautyPlaceholderResult:
    passes: bool = True
    predicted: float | None = None
    deferred: bool = True
    reason: str = _REASON


def evaluate_beauty_placeholder(channel: str = "B decay") -> BeautyPlaceholderResult:
    return BeautyPlaceholderResult(reason=f"{channel}: {_REASON}")
