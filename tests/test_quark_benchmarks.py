"""Tests for deterministic quark-sector benchmarks."""

import warnings

import numpy as np
import pytest

from qcd.constants import M_TOP_MS
from quarkConstraints.benchmarks import (
    _FIXED_SCALE_TARGETS_PDG2024_MT_V1,
    default_quark_benchmark,
    default_quark_targets,
    rough_sm_targets,
)
from quarkConstraints.pdg_quark_masses import (
    PDG_QUARK_MASSES_EDITION,
    pdg_up_down_arrays_at_scale,
)
from quarkConstraints.scales import (
    DEFAULT_QUARK_BENCHMARK_H_RS_MAX,
    DEFAULT_QUARK_BENCHMARK_XI_KK,
    DEFAULT_QUARK_FIT_SCALE_GEV,
    DEFAULT_QUARK_TARGET_SCALE_GEV,
    DEFAULT_QUARK_XI_KK,
    GAUGE_KK_ROOT_NN,
    default_quark_m_kk_from_lambda_ir,
)


def test_rough_sm_targets_have_expected_shapes_and_emit_deprecation_warning():
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        targets = rough_sm_targets()
    assert any(issubclass(w.category, DeprecationWarning) for w in captured)
    assert targets.up_masses.shape == (3,)
    assert targets.down_masses.shape == (3,)
    assert targets.ckm.shape == (3, 3)
    assert targets.label == "rough-sm-like-compatibility"


def test_default_quark_targets_use_pdg2024_mt_label():
    targets = default_quark_targets()

    assert targets.label == "pdg-2024-msbar-mu-mt-v1"
    # Mass-target scale moves to mu = m_t(m_t).
    assert np.isclose(DEFAULT_QUARK_FIT_SCALE_GEV, M_TOP_MS)
    assert np.isclose(DEFAULT_QUARK_FIT_SCALE_GEV, 163.5)
    # WC matching scale stays at 3 TeV (orthogonal).
    assert np.isclose(DEFAULT_QUARK_TARGET_SCALE_GEV, 3000.0)


def test_default_quark_targets_match_pdg2024_evolved_arrays():
    targets = default_quark_targets()
    expected_up, expected_down = pdg_up_down_arrays_at_scale(DEFAULT_QUARK_FIT_SCALE_GEV)
    np.testing.assert_allclose(targets.up_masses, expected_up, rtol=1e-12)
    np.testing.assert_allclose(targets.down_masses, expected_down, rtol=1e-12)


def test_pdg2024_target_bundle_provenance_assertions():
    bundle = _FIXED_SCALE_TARGETS_PDG2024_MT_V1
    assert bundle["label"] == "pdg-2024-msbar-mu-mt-v1"
    assert bundle["edition"] == PDG_QUARK_MASSES_EDITION
    assert np.isclose(bundle["scale_GeV"], 163.5)
    # Per-flavor 2sigma relative arrays are present and positive.
    assert bundle["up_2sigma_relative"].shape == (3,)
    assert bundle["down_2sigma_relative"].shape == (3,)
    assert (bundle["up_2sigma_relative"] > 0).all()
    assert (bundle["down_2sigma_relative"] > 0).all()


def test_default_quark_m_kk_from_lambda_ir_respects_explicit_convention():
    assert np.isclose(default_quark_m_kk_from_lambda_ir(3000.0), 3000.0)
    assert np.isclose(default_quark_m_kk_from_lambda_ir(3000.0, xi_KK=2.5), 7500.0)
    assert np.isclose(GAUGE_KK_ROOT_NN, 2.448687135269161)
    assert np.isclose(DEFAULT_QUARK_BENCHMARK_XI_KK, DEFAULT_QUARK_XI_KK)
    assert np.isclose(DEFAULT_QUARK_BENCHMARK_H_RS_MAX, 1.0)


def test_default_benchmark_has_expected_metadata():
    bench = default_quark_benchmark()

    assert bench.name == "repo-local-mfv-benchmark"
    assert bench.point.metadata["preferred_r_window"] == (0.1, 0.4)
    assert np.isclose(bench.point.metadata["overall_y"], 2.8)
    assert np.isclose(bench.point.metadata["default_target_scale_GeV"], 163.5)
    assert bench.point.metadata["default_target_label"] == "pdg-2024-msbar-mu-mt-v1"
