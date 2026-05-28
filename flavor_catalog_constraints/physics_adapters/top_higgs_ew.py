"""Placeholder adapters for top/Higgs/EW-precision observables."""

from __future__ import annotations
from dataclasses import dataclass

__all__ = ["TopHiggsEwPlaceholderResult", "evaluate_top_higgs_ew_placeholder"]

_REASON = "requires top-FCNC/EW-precision/LFV-boson machinery; not yet in core"


@dataclass(frozen=True)
class TopHiggsEwPlaceholderResult:
    passes: bool = True
    predicted: float | None = None
    deferred: bool = True
    reason: str = _REASON


def evaluate_top_higgs_ew_placeholder(channel: str = "top/EW") -> TopHiggsEwPlaceholderResult:
    return TopHiggsEwPlaceholderResult(reason=f"{channel}: {_REASON}")
