"""Production tests for B002 (S_psiK_S / sin(2 beta))."""

from __future__ import annotations

import cmath
from dataclasses import replace
import math
from pathlib import Path

import numpy as np
import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.deltaf2 import (
    B_1_BD,
    B_4_BD,
    B_5_BD,
    DEFAULT_DELTA_F2_INPUTS_V1,
    DELTA_M_BD_SM,
    F_BD,
    M_BD,
    M_B_QUARK,
    M_D_QUARK_BD,
    compute_delta_f2_wilsons,
    compute_m12_np,
    evaluate_bd_mixing_with_running,
    _evolve_wilsons,
)
from quarkConstraints.ckm_extraction import (
    repo_default_ckm_matrix,
    repo_default_ckm_phases,
)

_PID = "B002"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "beauty" / "B002.yaml"


def _yaml():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)


def _yaml_pdg_block():
    return _yaml()["pdg_or_equivalent"]


def _bd_couplings(
    left: complex,
    right: complex,
    M_KK: float = 3000.0,
    ckm_matrix: np.ndarray | None = None,
) -> QuarkMassBasisCouplings:
    """Minimal valid mass-basis couplings with only the d-b slot populated."""
    zeros = np.zeros((3, 3), dtype=np.complex128)
    left_down = zeros.copy()
    right_down = zeros.copy()
    left_down[0, 2] = left
    left_down[2, 0] = np.conj(left)
    right_down[0, 2] = right
    right_down[2, 0] = np.conj(right)
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
        ckm_matrix=ckm_matrix,
        ckm_source=("B002 test CKM" if ckm_matrix is not None else None),
    )


def _rephase_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    up_phases: np.ndarray,
    down_phases: np.ndarray,
) -> QuarkMassBasisCouplings:
    """Apply a simultaneous Dirac mass-eigenstate rephasing."""

    up = np.diag(up_phases)
    down = np.diag(down_phases)
    def transform_up(matrix: np.ndarray) -> np.ndarray:
        return up.conjugate() @ matrix @ up

    def transform_down(matrix: np.ndarray) -> np.ndarray:
        return down.conjugate() @ matrix @ down

    return replace(
        couplings,
        left_overlap=transform_down(couplings.left_overlap),
        right_up_overlap=transform_up(couplings.right_up_overlap),
        right_down_overlap=transform_down(couplings.right_down_overlap),
        left_up=transform_up(couplings.left_up),
        left_down=transform_down(couplings.left_down),
        right_up=transform_up(couplings.right_up),
        right_down=transform_down(couplings.right_down),
        ckm_matrix=(
            up.conjugate()
            @ np.asarray(couplings.ckm_matrix, dtype=np.complex128)
            @ down
        ),
        ckm_source="B002 random rephasing",
    )


def _bd_wilsons(couplings: QuarkMassBasisCouplings):
    wilsons = compute_delta_f2_wilsons(
        couplings,
        inputs=DEFAULT_DELTA_F2_INPUTS_V1,
    )
    return next(w for w in wilsons if w.input.key == "b_d")


def _budget_from_yaml_and_core():
    pdg = _yaml_pdg_block()
    exp = pdg["canonical_experimental_average"]
    beta = pdg["beta_physical_solution"]
    penguin = pdg["penguin_phase_bound"]
    ckm_phase = repo_default_ckm_phases()

    two_beta_rad = ckm_phase.two_beta
    sigma_two_beta_rad = math.radians(
        2.0
        * max(
            abs(float(beta["upper_uncertainty_degrees"])),
            abs(float(beta["lower_uncertainty_degrees"])),
        )
    )
    cos_two_beta_abs = abs(math.cos(two_beta_rad))
    sm_sigma = cos_two_beta_abs * sigma_two_beta_rad
    penguin_sigma = cos_two_beta_abs * math.radians(
        abs(float(penguin["bound_abs_delta_phi_d_degrees"]))
    )
    budget = math.sqrt(
        float(exp["uncertainty"]) * float(exp["uncertainty"])
        + sm_sigma * sm_sigma
        + penguin_sigma * penguin_sigma
    )
    sm_sin2beta = ckm_phase.sin_2beta
    return {
        "two_beta_rad": two_beta_rad,
        "two_beta_deg": ckm_phase.two_beta_degrees,
        "beta_rad": ckm_phase.beta,
        "beta_deg": ckm_phase.beta_degrees,
        "ckm_phase_source": ckm_phase.source,
        "sigma_two_beta_rad": sigma_two_beta_rad,
        "sm_sin2beta": sm_sin2beta,
        "central_residual": abs(float(exp["value"]) - sm_sin2beta),
        "sm_sigma": sm_sigma,
        "penguin_sigma": penguin_sigma,
        "budget": budget,
        "m12_sm": DELTA_M_BD_SM / 2.0,
    }


def _core_spsi_ks_from_running_m12(
    couplings: QuarkMassBasisCouplings,
    constraint,
):
    wilsons = _bd_wilsons(couplings)
    evolved = _evolve_wilsons(wilsons, mu_had=4.18)
    m12_np = compute_m12_np(
        evolved,
        F_BD,
        M_BD,
        M_B_QUARK,
        M_D_QUARK_BD,
        B_1_BD,
        B_4_BD,
        B_5_BD,
    )
    magnitude = evaluate_bd_mixing_with_running(wilsons, mu_had=4.18)
    ckm = (
        repo_default_ckm_matrix()
        if couplings.ckm_matrix is None
        else np.asarray(couplings.ckm_matrix, dtype=np.complex128)
    )
    ckm_factor = (np.conjugate(ckm[2, 0]) * ckm[2, 2]) ** 2
    m12_sm = (
        constraint.anchor.budget_band.m12_sm_gev
        * ckm_factor
        / abs(ckm_factor)
    )
    m12_ratio = m12_np / m12_sm
    phi_d_np = cmath.phase(1.0 + m12_ratio)
    predicted = math.sin(constraint.anchor.budget_band.two_beta_radians + phi_d_np)
    residual = abs(predicted - constraint.anchor.value)
    ratio = residual / constraint.anchor.budget
    return m12_np, magnitude, phi_d_np, predicted, residual, ratio


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "beauty"
    assert constraint.observable == "S_psiK_S"


def test_anchor_matches_yaml():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    exp = pdg["canonical_experimental_average"]
    mode = pdg["mode_specific_jpsi_ks"]
    beta = pdg["beta_physical_solution"]
    penguin = pdg["penguin_phase_bound"]
    budget = _budget_from_yaml_and_core()

    assert constraint.anchor.experimental.value == pytest.approx(exp["value"])
    assert constraint.anchor.experimental.uncertainty == pytest.approx(
        exp["uncertainty"]
    )
    assert constraint.anchor.experimental.source_url == exp["source_url"]
    assert constraint.anchor.mode_specific_jpsi_ks.value == pytest.approx(
        mode["value"]
    )
    assert constraint.anchor.mode_specific_jpsi_ks.uncertainty == pytest.approx(
        mode["uncertainty"]
    )
    assert constraint.anchor.beta_solution.beta_degrees == pytest.approx(
        beta["value_degrees"]
    )
    assert constraint.anchor.beta_solution.upper_uncertainty_degrees == pytest.approx(
        beta["upper_uncertainty_degrees"]
    )
    assert constraint.anchor.beta_solution.lower_uncertainty_degrees == pytest.approx(
        beta["lower_uncertainty_degrees"]
    )
    assert constraint.anchor.penguin_phase_bound.bound_abs_delta_phi_d_degrees == (
        pytest.approx(penguin["bound_abs_delta_phi_d_degrees"])
    )

    assert constraint.anchor.sm_value == pytest.approx(budget["sm_sin2beta"])
    assert constraint.anchor.budget_band.two_beta_radians == pytest.approx(
        budget["two_beta_rad"]
    )
    assert constraint.anchor.budget_band.two_beta_degrees == pytest.approx(
        budget["two_beta_deg"]
    )
    assert constraint.anchor.budget_band.beta_radians == pytest.approx(
        budget["beta_rad"]
    )
    assert constraint.anchor.budget_band.beta_degrees == pytest.approx(
        budget["beta_deg"]
    )
    assert constraint.anchor.budget_band.ckm_phase_source == budget["ckm_phase_source"]
    assert constraint.anchor.budget_band.sigma_two_beta_radians == pytest.approx(
        budget["sigma_two_beta_rad"]
    )
    assert constraint.anchor.budget_band.central_residual == pytest.approx(
        budget["central_residual"]
    )
    assert constraint.anchor.budget_band.sm_sin2beta_sigma == pytest.approx(
        budget["sm_sigma"]
    )
    assert constraint.anchor.budget_band.penguin_phase_sigma == pytest.approx(
        budget["penguin_sigma"]
    )
    assert constraint.anchor.budget == pytest.approx(budget["budget"])
    assert constraint.anchor.budget == pytest.approx(0.017723477116003937)
    assert constraint.anchor.budget_band.m12_sm_gev == pytest.approx(
        budget["m12_sm"]
    )


def test_evaluate_without_input_degrades_gracefully():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(constraint.anchor.value)
    assert result.sm_prediction == pytest.approx(constraint.anchor.sm_value)
    assert result.budget == pytest.approx(constraint.anchor.budget)
    assert result.diagnostics["missing_extra"] == "quark_mass_basis_couplings"
    assert "needs_human_physics" not in result.diagnostics
    assert "NEEDS-HUMAN-PHYSICS" not in result.notes
    assert result.diagnostics["ckm_phase_source"] == repo_default_ckm_phases().source


def test_sm_limit_uses_in_core_ckm_phase():
    constraint = fcc.get(_PID)
    phase = repo_default_ckm_phases()
    couplings = _bd_couplings(left=0.0j, right=0.0j)
    result = constraint.evaluate(point_builder.build_from_quark_couplings(couplings))

    assert result.predicted == pytest.approx(phase.sin_2beta)
    assert result.predicted == pytest.approx(0.7083388934693238)
    assert result.predicted == pytest.approx(result.sm_prediction)
    assert result.passes is True
    assert result.diagnostics["phi_d_np_rad"] == pytest.approx(0.0)
    assert result.diagnostics["two_beta_deg"] == pytest.approx(phase.two_beta_degrees)
    assert result.diagnostics["ckm_phase_source"] == phase.source
    assert "needs_human_physics" not in result.diagnostics
    assert "NEEDS-HUMAN-PHYSICS" not in result.notes


def test_evaluate_runs_end_to_end_with_real_couplings_and_real_finite_fields():
    couplings = _bd_couplings(left=(1.0e-3 + 1.0e-3j), right=0.0j)
    point = point_builder.build_from_quark_couplings(couplings)
    result = fcc.get(_PID).evaluate(point)

    assert result.process_id == _PID
    for value in (
        result.predicted,
        result.ratio,
        result.budget,
        result.sm_prediction,
        result.experimental,
    ):
        assert isinstance(value, float)
        assert math.isfinite(value)
    assert isinstance(result.diagnostics["m12_np_gev"], complex)
    assert isinstance(result.diagnostics["m12_np_over_m12_sm"], complex)
    assert isinstance(result.diagnostics["m12_sm_box_gev"], complex)
    assert isinstance(result.diagnostics["m12_sm_box_ckm_factor"], complex)
    assert isinstance(result.diagnostics["left_db_coupling"], complex)
    assert isinstance(result.diagnostics["wilson_coefficients"]["C1_VLL"], complex)
    for key in (
        "abs_m12_np_gev",
        "core_abs_m12_np_gev",
        "m12_sm_gev",
        "m12_sm_box_phase_rad",
        "re_m12_np_over_m12_sm",
        "im_m12_np_over_m12_sm",
        "phi_d_np_rad",
        "phi_d_np_deg",
        "s_psi_ks_total",
        "s_psi_ks_residual",
        "hadronic_scale_gev",
        "matching_scale_gev",
        "m_kk_gev",
        "central_sm_exp_residual",
        "hard_veto_budget",
        "budget_experimental_sigma",
        "budget_sm_sin2beta_sigma",
        "budget_penguin_phase_sigma",
        "beta_rad",
        "beta_deg",
        "two_beta_rad",
        "two_beta_deg",
        "sigma_two_beta_rad",
        "penguin_bound_rad",
        "delta_m_bd_sm_gev",
        "core_delta_m_bd_exp_gev",
        "core_legacy_m12_budget",
        "core_legacy_ratio_to_budget",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert result.diagnostics["qcd_running_applied"] is True
    assert result.diagnostics["hadronic_scale_gev"] == pytest.approx(4.18)
    assert result.diagnostics["core_input_key"] == "b_d"
    assert result.diagnostics["down_sector_indices"] == (0, 2)
    assert result.diagnostics["phase_uses_complex_m12_not_abs"] is True
    assert abs(result.diagnostics["m12_sm_box_gev"].imag) > 1.0e-18
    assert result.diagnostics["ckm_phase_source"] == repo_default_ckm_phases().source
    assert "needs_human_physics" not in result.diagnostics
    assert "NEEDS-HUMAN-PHYSICS" not in result.notes


def test_numbers_match_direct_running_complex_m12_phase_evaluator():
    couplings = _bd_couplings(left=(1.0e-3 + 1.0e-3j), right=0.0j)
    constraint = fcc.get(_PID)
    point = point_builder.build_from_quark_couplings(couplings)
    result = constraint.evaluate(point)
    m12_np, magnitude, phi_d_np, predicted, residual, ratio = (
        _core_spsi_ks_from_running_m12(couplings, constraint)
    )

    assert result.predicted == pytest.approx(predicted)
    assert result.ratio == pytest.approx(ratio)
    assert result.budget == pytest.approx(constraint.anchor.budget)
    assert result.diagnostics["s_psi_ks_residual"] == pytest.approx(residual)
    assert result.diagnostics["phi_d_np_rad"] == pytest.approx(phi_d_np)
    assert result.diagnostics["m12_np_gev"].real == pytest.approx(m12_np.real)
    assert result.diagnostics["m12_np_gev"].imag == pytest.approx(m12_np.imag)
    assert result.diagnostics["abs_m12_np_gev"] == pytest.approx(abs(m12_np))
    assert result.diagnostics["core_abs_m12_np_gev"] == pytest.approx(
        magnitude.abs_m12_np
    )
    assert abs(m12_np) == pytest.approx(magnitude.abs_m12_np)
    # Re-pinned after restoring the complex SM-box phase.
    assert result.predicted == pytest.approx(0.7126711557121855)
    assert result.ratio == pytest.approx(0.1507128479757265)
    assert result.diagnostics["phi_d_np_deg"] == pytest.approx(
        0.3527422342297369
    )


@pytest.mark.parametrize(
    ("couplings", "expected_pass"),
    [
        (_bd_couplings(left=(1.0e-3 + 1.0e-3j), right=0.0j), True),
        (_bd_couplings(left=(5.0e-3 + 5.0e-3j), right=0.0j), False),
    ],
)
def test_safe_point_passes_and_large_np_point_fails(
    couplings: QuarkMassBasisCouplings,
    expected_pass: bool,
):
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point_builder.build_from_quark_couplings(couplings))

    assert result.passes is expected_pass
    if expected_pass:
        assert result.ratio <= 1.0
    else:
        assert result.ratio > 1.0
        # Re-pinned after restoring the complex SM-box phase.
        assert result.diagnostics["phi_d_np_deg"] == pytest.approx(
            7.636816740620537
        )


def test_evaluate_is_pure_and_deterministic():
    couplings = _bd_couplings(left=(1.0e-3 + 1.0e-3j), right=0.0j)
    before_left_down = couplings.left_down.copy()
    before_right_down = couplings.right_down.copy()
    point = point_builder.build_from_quark_couplings(couplings)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    np.testing.assert_array_equal(couplings.left_down, before_left_down)
    np.testing.assert_array_equal(couplings.right_down, before_right_down)


def test_prediction_is_invariant_under_random_quark_field_rephasings():
    rng = np.random.default_rng(20022026)
    couplings = _bd_couplings(
        left=1.0e-3 + 2.0e-3j,
        right=-0.7e-3 + 0.4e-3j,
        ckm_matrix=repo_default_ckm_matrix(),
    )
    constraint = fcc.get(_PID)
    baseline = constraint.evaluate(point_builder.build_from_quark_couplings(couplings))

    for _ in range(16):
        up_phases = np.exp(1j * rng.uniform(-math.pi, math.pi, size=3))
        down_phases = np.exp(1j * rng.uniform(-math.pi, math.pi, size=3))
        rephased = _rephase_couplings(
            couplings,
            up_phases=up_phases,
            down_phases=down_phases,
        )
        result = constraint.evaluate(point_builder.build_from_quark_couplings(rephased))

        assert result.predicted == pytest.approx(baseline.predicted, abs=2e-14)
        assert result.ratio == pytest.approx(baseline.ratio, abs=2e-12)
        assert result.diagnostics["phi_d_np_rad"] == pytest.approx(
            baseline.diagnostics["phi_d_np_rad"],
            abs=2e-14,
        )
        assert result.diagnostics["m12_np_over_m12_sm"] == pytest.approx(
            baseline.diagnostics["m12_np_over_m12_sm"],
            abs=2e-14,
        )
