"""Anarchic Yukawa sampling and scoring utilities.

This module provides a reproducible definition of the v1 "anarchic prior":
- i.i.d. complex entries,
- log-uniform magnitudes,
- uniform phases,
- optional score used to prioritize candidates.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Tuple

import numpy as np


@dataclass(frozen=True)
class AnarchyConfig:
    """Configuration for anarchic Yukawa sampling and scoring."""

    magnitude_min: float = 1.0 / 3.0
    magnitude_max: float = 3.0
    phase_min: float = 0.0
    phase_max: float = 2.0 * np.pi
    yN_overall_min: float = 0.01
    yN_overall_max: float = 0.2
    w_band: float = 1.0
    w_cond: float = 1.0
    w_fit: float = 0.0

    def __post_init__(self):
        if self.magnitude_min <= 0:
            raise ValueError("magnitude_min must be positive")
        if self.magnitude_max <= self.magnitude_min:
            raise ValueError("magnitude_max must be greater than magnitude_min")
        if self.phase_max <= self.phase_min:
            raise ValueError("phase_max must be greater than phase_min")
        if self.yN_overall_min <= 0:
            raise ValueError("yN_overall_min must be positive")
        if self.yN_overall_max <= self.yN_overall_min:
            raise ValueError("yN_overall_max must be greater than yN_overall_min")


def _log_uniform(rng: np.random.Generator, low: float, high: float, size: Tuple[int, ...]) -> np.ndarray:
    """Sample a log-uniform random array over [low, high]."""
    log_low = np.log(low)
    log_high = np.log(high)
    return np.exp(rng.uniform(log_low, log_high, size=size))


def sample_complex_matrix(
    rng: np.random.Generator,
    shape: Tuple[int, int],
    magnitude_min: float,
    magnitude_max: float,
    phase_min: float,
    phase_max: float,
) -> np.ndarray:
    """Sample a complex matrix with log-uniform magnitudes and uniform phases."""
    magnitudes = _log_uniform(rng, magnitude_min, magnitude_max, size=shape)
    phases = rng.uniform(phase_min, phase_max, size=shape)
    return magnitudes * np.exp(1j * phases)


def band_penalty(matrix: np.ndarray, magnitude_min: float, magnitude_max: float) -> float:
    """Penalty for matrix entries outside [magnitude_min, magnitude_max]."""
    abs_entries = np.abs(np.asarray(matrix, dtype=complex))
    upper = np.maximum(0.0, np.log(abs_entries) - np.log(magnitude_max))
    lower = np.maximum(0.0, np.log(magnitude_min) - np.log(abs_entries))
    return float(np.sum(upper * upper + lower * lower))


def compute_anarchy_score(
    ytilde_n: np.ndarray,
    p_band: float,
    w_band: float,
    w_cond: float,
    w_fit: float,
    chi2_total: float = 0.0,
) -> Tuple[float, float]:
    """Compute the v1 anarchic score and condition penalty.

    Returns
    -------
    (score, condition_penalty)
    """
    cond_val = float(np.linalg.cond(ytilde_n))
    cond_penalty = float(np.log(cond_val) ** 2)
    score = -w_band * p_band - w_cond * cond_penalty - w_fit * float(chi2_total)
    return float(score), cond_penalty


def sample_anarchy_state(
    rng: np.random.Generator,
    config: AnarchyConfig,
    chi2_total: float = 0.0,
) -> Dict[str, object]:
    """Sample anarchic Yukawa structures and compute score metadata."""
    ytilde_e = sample_complex_matrix(
        rng,
        shape=(3, 3),
        magnitude_min=config.magnitude_min,
        magnitude_max=config.magnitude_max,
        phase_min=config.phase_min,
        phase_max=config.phase_max,
    )
    ytilde_n = sample_complex_matrix(
        rng,
        shape=(3, 3),
        magnitude_min=config.magnitude_min,
        magnitude_max=config.magnitude_max,
        phase_min=config.phase_min,
        phase_max=config.phase_max,
    )
    yN_overall = float(_log_uniform(rng, config.yN_overall_min, config.yN_overall_max, size=(1,))[0])

    p_band = band_penalty(ytilde_n, config.magnitude_min, config.magnitude_max)
    score, cond_penalty = compute_anarchy_score(
        ytilde_n,
        p_band=p_band,
        w_band=config.w_band,
        w_cond=config.w_cond,
        w_fit=config.w_fit,
        chi2_total=chi2_total,
    )

    return {
        "Ytilde_E": ytilde_e,
        "Ytilde_N": ytilde_n,
        "yN_overall": yN_overall,
        "band_penalty": p_band,
        "condition_penalty": cond_penalty,
        "score": score,
        "w_band": float(config.w_band),
        "w_cond": float(config.w_cond),
        "w_fit": float(config.w_fit),
    }
