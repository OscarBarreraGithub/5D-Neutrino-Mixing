"""Tests for deterministic quark-sector benchmarks."""

import numpy as np

from quarkConstraints import default_quark_benchmark, rough_sm_targets


def test_rough_sm_targets_have_expected_shapes():
    targets = rough_sm_targets()

    assert targets.up_masses.shape == (3,)
    assert targets.down_masses.shape == (3,)
    assert targets.ckm.shape == (3, 3)


def test_default_benchmark_has_expected_metadata():
    bench = default_quark_benchmark()

    assert bench.name == "repo-local-mfv-benchmark"
    assert bench.point.metadata["preferred_r_window"] == (0.1, 0.4)
    assert np.isclose(bench.point.metadata["overall_y"], 2.8)
