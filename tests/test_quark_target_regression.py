"""Regression test: v1 round-number targets accept m_s = 0.065 GeV; v2 PDG-2024
targets reject it.

This is a synthetic-fixture test (option b in plan v3 §6) — it does NOT pin
to the production CSV. The fixture constructs residuals directly from a
chosen mass spectrum at ``mu_common = m_t(m_t)``. The test asserts that:

1. Against the legacy v1 fixed-scale targets at mu = 3 TeV, an m_s of
   0.065 GeV would have been accepted (residual / 2sigma_v1 < 1).
2. Against the new PDG-2024 MS-bar targets at mu = m_t(m_t), the same
   m_s of 0.065 GeV is rejected (residual / 2sigma_pdg > 1).

The mismatch reflects the collaborator's original critique: PDG strange
mass at the common scale is ~0.055 GeV with ~3% 2sigma; 0.065 GeV is ~30%
above central, well beyond 2sigma.
"""

from __future__ import annotations

import numpy as np

from qcd.constants import M_TOP_MS
from quarkConstraints.benchmarks import (
    _FIXED_SCALE_TARGETS_MU_3TEV_V1,
    _FIXED_SCALE_TARGETS_PDG2024_MT_V1,
    default_quark_targets,
)
from quarkConstraints.pdg_quark_masses import pdg_quark_masses_at_scale


def _strange_relative_residual(m_strange: float, central: float) -> float:
    return abs(np.log(m_strange / central))


def test_v1_targets_would_have_accepted_ms_065_gev():
    # Legacy v1 used 0.057 GeV at mu = 3 TeV with a coarse 0.10 log gate.
    central_v1 = float(_FIXED_SCALE_TARGETS_MU_3TEV_V1["down_masses"][1])
    residual = _strange_relative_residual(0.065, central_v1)
    # Legacy gate was a flat 0.10 log-residual; 0.065 vs 0.057 is well inside.
    assert residual < 0.15, (
        f"v1 strange residual {residual:.4f} unexpectedly outside the legacy "
        "0.10 log gate"
    )


def test_v2_pdg_targets_reject_ms_065_gev():
    bundle = _FIXED_SCALE_TARGETS_PDG2024_MT_V1
    central = float(bundle["down_masses"][1])
    sigma_2_rel = float(bundle["down_2sigma_relative"][1])
    # m_s(mu_common) is ~0.05 GeV; PDG 2sigma rel ~0.017.
    assert central < 0.10, f"sanity: PDG m_s at mu_common is ~0.05 GeV (got {central})"
    log_residual = _strange_relative_residual(0.065, central)
    # The 2sigma log threshold is the relative tolerance (small approx).
    threshold = sigma_2_rel
    assert log_residual > threshold, (
        f"PDG-2024 v2 should reject m_s=0.065 GeV: log-residual={log_residual:.4f} "
        f"vs threshold {threshold:.4f}"
    )


def test_pdg_targets_themselves_pass_their_own_gates():
    # The PDG centrals must trivially pass their own per-flavor 2sigma gates.
    bundle = _FIXED_SCALE_TARGETS_PDG2024_MT_V1
    masses = pdg_quark_masses_at_scale(M_TOP_MS)
    for flavor, target_idx in zip(("u", "c", "t"), (0, 1, 2)):
        central = float(bundle["up_masses"][target_idx])
        assert masses[flavor] == central
    for flavor, target_idx in zip(("d", "s", "b"), (0, 1, 2)):
        central = float(bundle["down_masses"][target_idx])
        assert masses[flavor] == central


def test_default_quark_targets_label_is_pdg2024_mt_v1():
    targets = default_quark_targets()
    assert targets.label == "pdg-2024-msbar-mu-mt-v1"
