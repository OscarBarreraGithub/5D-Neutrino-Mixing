#!/usr/bin/env python3
"""Numerical audit for repo-v1 Delta F=2 Wilson running.

HONESTY BANNER: this script is a SELF-CONSISTENCY / linear-algebra check, not
physics validation against external paper targets. Its "independent LO
textbook" path shares the codebase alpha_s routine and the same LO anomalous
dimension block, so it verifies propagation algebra and basis bookkeeping only.
A true validation requires external Buras/UTfit/literature reference numbers
with a pass/fail tolerance, which are not embedded here. Treat this as a
diagnostic for the quarantined Lane C / Delta-F=2 stack, not production
validation; see `docs/audits/full_repo_audit_2026-07/FIX_LEDGER.md`.

Related reporting caveat: the Perez-Randall "O(3-4)x tension" conclusion is
convention-dependent; it shrinks to about 1.8x under plausible (seesaw
prefactor, v) choices and should not be treated as a robust physics-validation
failure.

The script compares the public in-code evolution with an independent LO
textbook calculation in the scalar LR basis used by ``deltaf2.py``:

    Q1_LR^BMU = -2 O5_LR,   Q2_LR^BMU = O4_LR.

Run from the repo root:

    python scripts/audit_wilson_rg.py
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from quarkConstraints.qcd_running import evolve_deltaf2_wilsons, run_alpha_s


MU_HIGH_GEV = 3000.0
MU_LOW_GEV = 2.0
M_T_GEV = 163.5
M_B_GEV = 4.18
M_C_GEV = 1.27

UNIT_WILSONS: tuple[tuple[str, tuple[complex, complex, complex, complex]], ...] = (
    ("C1_VLL", (1.0 + 0.0j, 0.0 + 0.0j, 0.0 + 0.0j, 0.0 + 0.0j)),
    ("C1_VRR", (0.0 + 0.0j, 1.0 + 0.0j, 0.0 + 0.0j, 0.0 + 0.0j)),
    ("C4_LR", (0.0 + 0.0j, 0.0 + 0.0j, 1.0 + 0.0j, 0.0 + 0.0j)),
    ("C5_LR", (0.0 + 0.0j, 0.0 + 0.0j, 0.0 + 0.0j, 1.0 + 0.0j)),
)


@dataclass(frozen=True)
class Segment:
    mu_upper: float
    mu_lower: float
    n_f: int
    alpha_upper: float
    alpha_lower: float


def beta0(n_f: int) -> float:
    return (33.0 - 2.0 * n_f) / 3.0


def n_f_for_upper_scale(mu: float) -> int:
    if mu > M_T_GEV:
        return 6
    if mu > M_B_GEV:
        return 5
    if mu > M_C_GEV:
        return 4
    return 3


def build_segments(mu_high: float, mu_low: float) -> tuple[Segment, ...]:
    thresholds = sorted(
        [mass for mass in (M_T_GEV, M_B_GEV, M_C_GEV) if mu_low < mass < mu_high],
        reverse=True,
    )
    boundaries = [mu_high, *thresholds, mu_low]
    segments: list[Segment] = []
    for upper, lower in zip(boundaries[:-1], boundaries[1:], strict=True):
        segments.append(
            Segment(
                mu_upper=upper,
                mu_lower=lower,
                n_f=n_f_for_upper_scale(upper),
                alpha_upper=run_alpha_s(upper),
                alpha_lower=run_alpha_s(lower),
            )
        )
    return tuple(segments)


def scalar_lr_segment_matrix(segment: Segment) -> np.ndarray:
    """Closed-form LO scalar LR matrix for one fixed-n_f segment.

    This is the BMU LR LO coefficient block ``[[2, 0], [12, -16]]``
    conjugated to the code's conventional scalar ``[C4_LR, C5_LR]`` order,
    giving ``[[-16, -6], [0, 2]]``.

    Defensive guard (R04-I1, C03 cleanup, 2026-05-25): the closed-form
    propagator below assumes the scalar LR coefficient ADM is strictly upper
    triangular (``gamma_54 == 0``), which is what conjugating the BMU LR
    block into the conventional ``[C4_LR, C5_LR]`` order produces.  If the
    pinned anomalous-dimension entries are ever edited in a way that breaks
    that property the shortcut becomes silently wrong, so we cross-check at
    runtime and assert the lower-left entry is exactly zero.
    """
    gamma44 = -16.0
    gamma45 = -6.0
    gamma54 = 0.0  # required by the upper-tri closed-form below
    gamma55 = 2.0
    gamma_lr_local = np.array([[gamma44, gamma45], [gamma54, gamma55]], dtype=float)
    if gamma_lr_local[1, 0] != 0.0 or not np.allclose(
        gamma_lr_local - np.triu(gamma_lr_local), 0.0
    ):
        raise ValueError(
            "audit_wilson_rg.scalar_lr_segment_matrix requires an upper-triangular "
            "scalar LR ADM (gamma_54 == 0); update the closed-form propagator if "
            "the anomalous-dimension block is changed."
        )
    eta = segment.alpha_upper / segment.alpha_lower
    f44 = eta ** (gamma44 / (2.0 * beta0(segment.n_f)))
    f55 = eta ** (gamma55 / (2.0 * beta0(segment.n_f)))
    f45 = gamma45 * (f44 - f55) / (gamma44 - gamma55)
    return np.array([[f44, f45], [0.0, f55]], dtype=float)


def reference_evolution() -> tuple[complex, complex, complex, complex]:
    segments = build_segments(MU_HIGH_GEV, MU_LOW_GEV)

    vll_factor = 1.0
    lr_matrix = np.eye(2, dtype=float)
    for segment in segments:
        eta = segment.alpha_upper / segment.alpha_lower
        vll_factor *= eta ** (4.0 / (2.0 * beta0(segment.n_f)))
        lr_matrix = scalar_lr_segment_matrix(segment) @ lr_matrix

    return (
        complex(vll_factor),
        complex(vll_factor),
        complex(lr_matrix[0, 0]),
        complex(lr_matrix[1, 0]),
    )


def main() -> None:
    segments = build_segments(MU_HIGH_GEV, MU_LOW_GEV)
    print(f"Wilson RG audit: mu_high={MU_HIGH_GEV:g} GeV, mu_low={MU_LOW_GEV:g} GeV")
    print("Segments:")
    for segment in segments:
        print(
            "  "
            f"{segment.mu_upper:g}->{segment.mu_lower:g} GeV, "
            f"n_f={segment.n_f}, "
            f"alpha_s={segment.alpha_upper:.12g}->{segment.alpha_lower:.12g}"
        )

    print("\nIn-code unit-vector evolution:")
    observed_vectors: list[np.ndarray] = []
    for name, vector in UNIT_WILSONS:
        evolved = evolve_deltaf2_wilsons(*vector, MU_HIGH_GEV, MU_LOW_GEV)
        observed = np.array(evolved, dtype=np.complex128)
        observed_vectors.append(observed)
        formatted = ", ".join(f"{value.real:.12g}{value.imag:+.12g}j" for value in observed)
        print(f"  {name}: ({formatted})")

    expected_c1, _, expected_c4_from_c4, expected_c5_from_c4 = reference_evolution()
    expected_c5_from_c5 = complex(
        np.linalg.multi_dot([scalar_lr_segment_matrix(segment) for segment in segments[::-1]])[
            1, 1
        ]
    )
    expected_c4_from_c5 = complex(
        np.linalg.multi_dot([scalar_lr_segment_matrix(segment) for segment in segments[::-1]])[
            0, 1
        ]
    )
    expected_vectors = [
        np.array([expected_c1, 0.0, 0.0, 0.0], dtype=np.complex128),
        np.array([0.0, expected_c1, 0.0, 0.0], dtype=np.complex128),
        np.array([0.0, 0.0, expected_c4_from_c4, expected_c5_from_c4], dtype=np.complex128),
        np.array([0.0, 0.0, expected_c4_from_c5, expected_c5_from_c5], dtype=np.complex128),
    ]

    max_relative = 0.0
    for observed, expected in zip(observed_vectors, expected_vectors, strict=True):
        scale = np.maximum(np.abs(expected), 1.0e-15)
        max_relative = max(max_relative, float(np.max(np.abs(observed - expected) / scale)))
    print(f"\nMax relative discrepancy vs closed-form LO reference: {max_relative:.3e}")


if __name__ == "__main__":
    main()
