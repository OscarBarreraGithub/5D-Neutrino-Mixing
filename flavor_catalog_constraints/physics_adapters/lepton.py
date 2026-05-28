"""Placeholder adapters for charged-lepton-flavor-violating + neutrino observables."""

from __future__ import annotations
from dataclasses import dataclass

__all__ = [
    "LeptonPlaceholderResult",
    "evaluate_lfv_dipole_placeholder",
    "evaluate_lfv_4lepton_placeholder",
    "evaluate_mue_conversion_placeholder",
    "evaluate_muonium_placeholder",
    "evaluate_trident_placeholder",
]


_LFV_DIPOLE_REASON = (
    "requires LFV dipole Wilson matching (C_eγ, C_uγ); 5D RS one-loop "
    "contributions from KK gauge / Higgs to imaginary parts; not yet in core"
)
_LFV_4L_REASON = (
    "requires LFV 4-fermion Wilson matching (vector + scalar + tensor); "
    "5D RS tree- and loop-level NP contributions; not yet in core"
)
_MUE_CONV_REASON = (
    "requires LFV 4-fermion + dipole interference + nuclear-overlap "
    "form factors (Kitano-Koike-Okada); not yet in core"
)
_TRIDENT_REASON = (
    "requires Z'-mediated trident amplitude (Altmannshofer et al 1406.2332); "
    "5D RS Z' / KK contributions to ν N → ν N μμ; not yet in core"
)


@dataclass(frozen=True)
class LeptonPlaceholderResult:
    passes: bool = True
    predicted: float | None = None
    deferred: bool = True
    reason: str = _LFV_DIPOLE_REASON


def evaluate_lfv_dipole_placeholder(channel: str = "mu->e gamma") -> LeptonPlaceholderResult:
    return LeptonPlaceholderResult(reason=f"{channel}: {_LFV_DIPOLE_REASON}")

def evaluate_lfv_4lepton_placeholder(channel: str = "mu->3e") -> LeptonPlaceholderResult:
    return LeptonPlaceholderResult(reason=f"{channel}: {_LFV_4L_REASON}")

def evaluate_mue_conversion_placeholder(nucleus: str = "Al") -> LeptonPlaceholderResult:
    return LeptonPlaceholderResult(reason=f"mu->e conv on {nucleus}: {_MUE_CONV_REASON}")

def evaluate_muonium_placeholder() -> LeptonPlaceholderResult:
    return LeptonPlaceholderResult(reason=f"muonium-antimuonium: {_LFV_4L_REASON}")

def evaluate_trident_placeholder() -> LeptonPlaceholderResult:
    return LeptonPlaceholderResult(reason=_TRIDENT_REASON)
