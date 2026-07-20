import numpy as np
import pytest

from warpConfig.wavefuncs import f_IR, f_UV

EPSILON_RS = 1.0e-15


def test_overlap_c_half_limit_matches_logarithmic_result():
    expected = np.sqrt(1.0 / (-2.0 * np.log(EPSILON_RS)))

    assert float(f_IR(0.5, EPSILON_RS)) == pytest.approx(expected)
    assert float(f_UV(0.5, EPSILON_RS)) == pytest.approx(expected)


def test_overlap_c_zero_exact_values():
    ir_sq = float(f_IR(0.0, EPSILON_RS)) ** 2
    uv_sq = float(f_UV(0.0, EPSILON_RS)) ** 2

    assert ir_sq == pytest.approx(0.5 / (1.0 - EPSILON_RS), rel=1e-14)
    assert uv_sq == pytest.approx(
        EPSILON_RS * 0.5 / (1.0 - EPSILON_RS), rel=1e-14
    )


def test_overlap_localization_limits_are_finite_without_overflow():
    c_uv_localized = 10.0
    c_ir_localized = -9.5

    with np.errstate(over="raise", divide="raise", invalid="raise"):
        uv_localized_ir = float(f_IR(c_uv_localized, EPSILON_RS))
        uv_localized_uv = float(f_UV(c_uv_localized, EPSILON_RS))
        ir_localized_ir = float(f_IR(c_ir_localized, EPSILON_RS))
        ir_localized_uv = float(f_UV(c_ir_localized, EPSILON_RS))

    values = np.array(
        [uv_localized_ir, uv_localized_uv, ir_localized_ir, ir_localized_uv]
    )
    assert np.all(np.isfinite(values))

    assert uv_localized_ir < 1e-120
    assert uv_localized_uv == pytest.approx(np.sqrt(c_uv_localized - 0.5))

    assert ir_localized_ir == pytest.approx(np.sqrt(0.5 - c_ir_localized))
    assert ir_localized_uv < 1e-140


def test_ir_uv_overlap_identity():
    c_values = np.array([-0.75, -0.25, 0.0, 0.25, 0.5, 0.75, 1.25])

    ir = f_IR(c_values, EPSILON_RS)
    uv = f_UV(c_values, EPSILON_RS)

    assert np.allclose(
        uv,
        (EPSILON_RS ** (0.5 - c_values)) * ir,
        rtol=1e-12,
        atol=1e-300,
    )
