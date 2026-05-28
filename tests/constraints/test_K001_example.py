"""Per-constraint test for the K001 throwaway example.

Demonstrates the per-constraint test flow. If the orchestrator deletes
the K001 example constraint, these tests skip cleanly rather than fail.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder

_PID = "K001"

pytestmark = pytest.mark.skipif(
    _PID not in fcc.all_constraints(),
    reason="K001 example constraint not registered",
)


def test_anchor_matches_yaml():
    c = fcc.get(_PID)
    # K001.yaml canonical_experimental_value.value = 0.002228
    assert c.anchor.value == pytest.approx(0.002228)
    assert c.anchor.uncertainty == pytest.approx(0.000011)
    # SM reference pulled separately.
    assert c.sm_value == pytest.approx(0.002161)
    assert c.anchor.source_url  # provenance carried through, not defaulted


def test_evaluate_without_input_degrades_gracefully():
    r = fcc.get(_PID).evaluate(point_builder.empty_point())
    assert r.process_id == _PID
    assert r.passes is True  # no input -> no veto
    assert r.experimental == pytest.approx(0.002228)
    assert "absent" in r.notes
    assert r.diagnostics["missing_extra"] == "quark_mass_basis_couplings"


def _sd_couplings(left: complex, right: complex, M_KK: float = 3000.0):
    """Minimal QuarkMassBasisCouplings with only the (s,d) slot populated.

    Mirrors the helper in tests/test_epsilon_k_physics.py so the K001
    adapter->physics path is exercised on a real, valid input object.
    """
    from quarkConstraints.couplings import QuarkMassBasisCouplings

    zeros = np.zeros((3, 3), dtype=np.complex128)
    left_down = zeros.copy()
    right_down = zeros.copy()
    left_down[0, 1] = left
    left_down[1, 0] = np.conj(left)
    right_down[0, 1] = right
    right_down[1, 0] = np.conj(right)
    return QuarkMassBasisCouplings(
        M_KK=M_KK,
        xi_KK=1.0,
        alpha_s=0.09,
        g_s=1.0,
        left_overlap=zeros,
        right_up_overlap=zeros,
        right_down_overlap=zeros,
        left_up=zeros,
        left_down=left_down,
        right_up=zeros,
        right_down=right_down,
    )


def test_evaluate_with_input():
    """End-to-end: adapter -> evaluate_epsilon_k -> ConstraintResult.

    Guards the physics path against regression and confirms the numeric
    fields are real, finite floats (no complex leakage out of the adapter).
    """
    # Tiny complex couplings so epsilon_K^NP is finite and well below 1.
    couplings = _sd_couplings(left=0.001 + 0.0005j, right=0.001 + 0.0002j)
    point = point_builder.build_from_quark_couplings(couplings)
    r = fcc.get(_PID).evaluate(point)

    assert r.process_id == _PID
    for value in (r.predicted, r.ratio, r.budget, r.sm_prediction, r.experimental):
        assert isinstance(value, float)
        assert math.isfinite(value)
