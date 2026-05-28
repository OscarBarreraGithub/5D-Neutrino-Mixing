"""TEMPLATE for a per-constraint test. Copy to test_<ID>.py and fill in.

The global contract test already checks protocol/metadata/sidecar/
isolation for ALL constraints. A per-constraint test adds the things
only the author knows: the anchor value matches the yaml, and evaluate()
produces sensible results on a hand-built point.

This template module itself is inert (the example body is commented and
skipped) so it neither fails nor double-counts before you copy it.
"""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.skip(reason="template; copy to test_<ID>.py and fill in")

# import flavor_catalog_constraints as fcc
# from flavor_catalog_constraints import point_builder
#
# _PID = "XX000"
#
#
# def test_anchor_matches_yaml():
#     c = fcc.get(_PID)
#     # Pin the typed anchor back to the raw yaml value the author expects.
#     assert c.anchor.value == pytest.approx(EXPECTED_VALUE)
#
#
# def test_evaluate_without_input_degrades_gracefully():
#     r = fcc.get(_PID).evaluate(point_builder.empty_point())
#     assert r.process_id == _PID
#     assert "absent" in r.notes
#
#
# def test_evaluate_with_input():
#     couplings = ...  # build a QuarkMassBasisCouplings (see quarkConstraints)
#     point = point_builder.build_from_quark_couplings(couplings)
#     r = fcc.get(_PID).evaluate(point)
#     assert isinstance(r.ratio, float)
