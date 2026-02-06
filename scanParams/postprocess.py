"""Post-processing utilities for scan result reclassification."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional, Tuple

import numpy as np

from neutrinos.neutrinoValues import get_pmns

from .anarchy import AnarchyConfig, score_anarchy_from_matrix


def _to_float(value: Any, name: str) -> float:
    """Parse a scalar into float with a stable error message."""
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Could not parse {name}='{value}' as float") from exc


@dataclass(frozen=True)
class ReclassifyConfig:
    """Configuration for post-hoc row classification."""

    max_Y_bar: float = 4.0
    naturalness_range: Tuple[float, float] = (0.1, 4.0)
    require_lfv: bool = True
    anarchy: Optional[AnarchyConfig] = None
    anarchy_min_score: Optional[float] = None

    def __post_init__(self):
        if self.max_Y_bar <= 0:
            raise ValueError("max_Y_bar must be positive")
        lo, hi = self.naturalness_range
        if lo <= 0 or hi <= 0 or hi < lo:
            raise ValueError("naturalness_range must satisfy 0 < lo <= hi")
        if self.anarchy_min_score is not None and not np.isfinite(self.anarchy_min_score):
            raise ValueError("anarchy_min_score must be finite")


def reconstruct_y_n_bar_matrix(row: Dict[str, Any]) -> np.ndarray:
    """Reconstruct Ybar_N matrix from row-level eigenvalues and PMNS metadata."""
    y_n_bar = np.array(
        [
            _to_float(row["Y_N_bar_1"], "Y_N_bar_1"),
            _to_float(row["Y_N_bar_2"], "Y_N_bar_2"),
            _to_float(row["Y_N_bar_3"], "Y_N_bar_3"),
        ],
        dtype=float,
    )
    ordering = str(row.get("ordering", "normal"))
    alpha = _to_float(row.get("majorana_alpha", 0.0), "majorana_alpha")
    beta = _to_float(row.get("majorana_beta", 0.0), "majorana_beta")
    pmns = get_pmns(ordering=ordering, alpha=alpha, beta=beta)
    return pmns @ np.diag(y_n_bar.astype(complex))


def classify_row(row: Dict[str, Any], config: ReclassifyConfig) -> Dict[str, Any]:
    """Classify one row and return reclassification metadata."""
    y_vals = np.array(
        [
            _to_float(row["Y_E_bar_1"], "Y_E_bar_1"),
            _to_float(row["Y_E_bar_2"], "Y_E_bar_2"),
            _to_float(row["Y_E_bar_3"], "Y_E_bar_3"),
            _to_float(row["Y_N_bar_1"], "Y_N_bar_1"),
            _to_float(row["Y_N_bar_2"], "Y_N_bar_2"),
            _to_float(row["Y_N_bar_3"], "Y_N_bar_3"),
        ],
        dtype=float,
    )
    abs_y = np.abs(y_vals)
    max_y = float(np.max(abs_y))
    min_y = float(np.min(abs_y))

    reclass_perturbative = bool(max_y < config.max_Y_bar)
    lo, hi = config.naturalness_range
    reclass_natural = bool((min_y >= lo) and (max_y <= hi))

    reclass_lfv_passes = True
    if config.require_lfv:
        if "lfv_ratio" in row and row["lfv_ratio"] not in ("", None):
            ratio = _to_float(row["lfv_ratio"], "lfv_ratio")
            reclass_lfv_passes = bool(ratio <= 1.0)
        elif (
            "lfv_lhs" in row
            and row["lfv_lhs"] not in ("", None)
            and "lfv_rhs" in row
            and row["lfv_rhs"] not in ("", None)
        ):
            lhs = _to_float(row["lfv_lhs"], "lfv_lhs")
            rhs = _to_float(row["lfv_rhs"], "lfv_rhs")
            reclass_lfv_passes = bool(lhs <= rhs)
        else:
            raise ValueError("Row is missing LFV fields (need lfv_ratio or lfv_lhs/lfv_rhs)")

    reasons = []
    if not reclass_perturbative:
        reasons.append("perturbativity")
    if not reclass_natural:
        reasons.append("naturalness")
    if not reclass_lfv_passes:
        reasons.append("mu_to_e_gamma")

    anarchy_score = np.nan
    anarchy_band_penalty = np.nan
    anarchy_condition_penalty = np.nan
    anarchy_yN_overall = np.nan
    if config.anarchy is not None:
        y_n_bar_matrix = reconstruct_y_n_bar_matrix(row)
        anarchy_state = score_anarchy_from_matrix(
            y_n_bar_matrix=y_n_bar_matrix,
            config=config.anarchy,
        )
        anarchy_score = float(anarchy_state["score"])
        anarchy_band_penalty = float(anarchy_state["band_penalty"])
        anarchy_condition_penalty = float(anarchy_state["condition_penalty"])
        anarchy_yN_overall = float(anarchy_state["yN_overall"])
        if config.anarchy_min_score is not None and anarchy_score < config.anarchy_min_score:
            reasons.append("anarchy_score")

    return {
        "reclass_max_Y_bar_observed": max_y,
        "reclass_min_Y_bar_observed": min_y,
        "reclass_perturbative": reclass_perturbative,
        "reclass_natural": reclass_natural,
        "reclass_lfv_passes": reclass_lfv_passes,
        "reclass_anarchy_score": anarchy_score,
        "reclass_anarchy_band_penalty": anarchy_band_penalty,
        "reclass_anarchy_condition_penalty": anarchy_condition_penalty,
        "reclass_anarchy_yN_overall": anarchy_yN_overall,
        "reclass_passes_all": len(reasons) == 0,
        "reclass_reject_reason": ";".join(reasons),
    }
