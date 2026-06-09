"""Production tests for L001 (mu -> e gamma)."""

from __future__ import annotations

import importlib.util
import math
import sys
from dataclasses import fields
from pathlib import Path

import numpy as np
import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import anchors, point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.physics_adapters.lepton import (
    LMFVLeptonParameters,
    lmfv_lepton_parameters_from_yukawa_result,
    mu_to_e_gamma_from_lepton_input,
    mu_to_e_gamma_proxy_input,
)
from yukawa import compute_all_yukawas

_PID = "L001"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "charged_lepton" / "L001.yaml"
_HARNESS_PATH = _REPO_ROOT / "scripts" / "run_full_catalog_scan.py"

_ORACLE_BR = 1.5508276601368708e-10
_ORACLE_BR_RATIO = 1033.8851067579139
_ORACLE_OFFDIAG = complex(-0.034773700046830024, 0.05165132075170267)
_ORACLE_LHS = 0.062266115587389717
_ORACLE_DIP_RATIO_PAPER_C = 3.1133057793694858


def _yaml_pdg_block():
    with _SIDECAR.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _load_harness():
    spec = importlib.util.spec_from_file_location("run_full_catalog_scan_for_l001", _HARNESS_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _benchmark_yukawa():
    return compute_all_yukawas(
        Lambda_IR=3000.0,
        c_L=0.58,
        c_E=[0.75, 0.60, 0.50],
        c_N=0.27,
        M_N=1.22e18,
        lightest_nu_mass=0.002,
        ordering="normal",
    )


def _benchmark_carrier(*, m_kk_gev: float = 3000.0) -> LMFVLeptonParameters:
    return lmfv_lepton_parameters_from_yukawa_result(
        _benchmark_yukawa(),
        m_kk_gev=m_kk_gev,
    )


def _carrier_point(*, m_kk_gev: float = 3000.0, include_mass_extra: bool = False):
    carrier = _benchmark_carrier(m_kk_gev=m_kk_gev)
    extras = {"lepton_lmfv_parameters": carrier}
    if include_mass_extra:
        extras["kk_ew_mass_gev"] = m_kk_gev
    return point_builder.make_point(**extras)


def _rotation_pmns() -> np.ndarray:
    theta = 0.5
    c = math.cos(theta)
    s = math.sin(theta)
    return np.asarray(
        [
            [c, s, 0.0],
            [-s, c, 0.0],
            [0.0, 0.0, 1.0],
        ],
        dtype=complex,
    )


def _forged_carrier(base: LMFVLeptonParameters, **updates) -> LMFVLeptonParameters:
    forged = object.__new__(LMFVLeptonParameters)
    for field in fields(LMFVLeptonParameters):
        object.__setattr__(forged, field.name, updates.get(field.name, getattr(base, field.name)))
    return forged


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "charged_lepton"
    assert constraint.observable == "BR(mu -> e gamma)"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    exp = pdg["primary_current_limit"]
    repo = pdg["repo_default"]

    assert constraint.anchor.experimental.value == pytest.approx(exp["value"])
    assert constraint.anchor.experimental.source_url == exp["source_url"]
    assert constraint.anchor.budget == pytest.approx(1.5e-13)
    assert constraint.anchor.repo_default_br_limit.value == pytest.approx(
        repo["br_limit"]["value"]
    )
    assert constraint.anchor.prefactor_br.value == pytest.approx(
        repo["prefac_br"]["value"]
    )
    assert constraint.anchor.lfv_c.value == pytest.approx(repo["lfv_C"]["value"])
    assert constraint.anchor.lfv_c.value == pytest.approx(
        math.sqrt(exp["value"] / repo["prefac_br"]["value"])
    )
    assert constraint.anchor.c_paper.value == pytest.approx(
        repo["c_paper"]["value"]
    )

    with pytest.raises(anchors.AnchorError):
        anchors.load_anchor(
            _PID,
            family="charged_lepton",
            candidates=("no_such_block",),
        )


def test_evaluate_without_input_degrades_gracefully():
    result = fcc.get(_PID).evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.sm_prediction is None
    assert result.experimental == pytest.approx(fcc.get(_PID).anchor.value)
    assert result.budget == pytest.approx(fcc.get(_PID).anchor.budget)
    assert result.notes == (
        "NOT EVALUATED - no lepton dipole prediction available "
        "(LMFV lepton carrier not on ParameterPoint)"
    ).replace(" - ", " \u2014 ")
    assert "pass" not in result.notes.lower()
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["unevaluated_reason"] == (
        "no lepton dipole prediction available "
        "(LMFV lepton carrier not on ParameterPoint)"
    )
    assert "non-vetoing only" in result.diagnostics["passes_semantics"]
    assert result.diagnostics["missing_extra"] == "lepton_lmfv_parameters"
    assert "needs_human_physics" not in result.diagnostics


def test_invalid_lepton_input_is_unevaluated_not_real_pass():
    result = fcc.get(_PID).evaluate(
        point_builder.make_point(lepton_lmfv_parameters={"y_n_bar": [1.0, 2.0]})
    )

    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.sm_prediction is None
    assert result.notes.startswith("NOT EVALUATED")
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["invalid_extra"] == "lepton_lmfv_parameters"
    assert result.diagnostics["exception_type"] == "TypeError"
    assert "needs_human_physics" not in result.diagnostics


def test_generated_carrier_oracle_is_rigorous_and_br_vetoed():
    result = fcc.get(_PID).evaluate(_carrier_point(include_mass_extra=True))
    tag, matching_status, needs_human, proxy_flags = _load_harness().tag_result(result)

    assert result.predicted == pytest.approx(_ORACLE_BR)
    assert result.ratio == pytest.approx(_ORACLE_BR_RATIO)
    assert result.passes is False
    assert result.sm_prediction == pytest.approx(0.0)
    assert result.diagnostics["dipole_lhs"] == pytest.approx(_ORACLE_LHS)
    assert result.diagnostics["dipole_rhs"] == pytest.approx(0.02)
    assert result.diagnostics["dipole_ratio_to_bound"] == pytest.approx(
        _ORACLE_DIP_RATIO_PAPER_C
    )
    assert result.diagnostics["off_diagonal_12"] == pytest.approx(_ORACLE_OFFDIAG)
    assert result.diagnostics["lfv_coefficient"] == pytest.approx(0.02)
    assert result.diagnostics["c_lfv_role"] == "dipole_rhs_diagnostic_only"
    assert result.diagnostics["br_limit"] == pytest.approx(1.5e-13)
    assert result.diagnostics["br_ratio_to_limit"] == pytest.approx(_ORACLE_BR_RATIO)
    assert result.diagnostics["core_passes"] is False
    assert result.diagnostics["input_kind"] == "LMFVLeptonParameters"
    assert result.diagnostics["extra_used"] == "lepton_lmfv_parameters"
    assert result.diagnostics["used_proxy"] is False
    assert result.diagnostics["lmfv_model"] == "Perez-Randall LMFV NDA"
    assert result.diagnostics["kk_ew_mass_extra_present"] is True
    assert result.diagnostics["kk_ew_mass_extra_used"] is False
    assert result.diagnostics["kk_ew_mass_extra_consistent"] is True
    assert "needs_human_physics" not in result.diagnostics
    assert tag == "rigorous"
    assert matching_status == "locked_lmfv_nda_carrier"
    assert needs_human is None
    assert proxy_flags == {}


def test_l001_mass_extra_reports_consistency_without_overriding_carrier_mass():
    carrier = _benchmark_carrier(m_kk_gev=3000.0)
    point = point_builder.make_point(
        lepton_lmfv_parameters=carrier,
        kk_ew_mass_gev=6000.0,
    )

    result = fcc.get(_PID).evaluate(point)

    assert result.predicted == pytest.approx(_ORACLE_BR)
    assert result.diagnostics["m_kk_gev"] == pytest.approx(3000.0)
    assert result.diagnostics["kk_ew_mass_extra_present"] is True
    assert result.diagnostics["kk_ew_mass_extra_used"] is False
    assert result.diagnostics["kk_ew_mass_extra_value"] == pytest.approx(6000.0)
    assert result.diagnostics["kk_ew_mass_extra_consistent"] is False


def test_pass_fail_is_pinned_to_br_and_independent_of_c_lfv():
    carrier = _benchmark_carrier()
    constraint = fcc.get(_PID)

    paper = mu_to_e_gamma_from_lepton_input(
        carrier,
        br_limit=constraint.anchor.budget,
        prefactor_br=constraint.anchor.prefactor_br.value,
        c_lfv=0.02,
    )
    loose = mu_to_e_gamma_from_lepton_input(
        carrier,
        br_limit=constraint.anchor.budget,
        prefactor_br=constraint.anchor.prefactor_br.value,
        c_lfv=100.0,
    )

    assert paper.dipole_rhs != pytest.approx(loose.dipole_rhs)
    assert paper.diagnostics["core_passes"] is False
    assert loose.diagnostics["core_passes"] is True
    assert paper.branching_fraction == pytest.approx(loose.branching_fraction)
    assert paper.ratio_to_limit == pytest.approx(loose.ratio_to_limit)
    assert paper.passes is loose.passes is False


def test_generated_carrier_emits_no_truthy_proxy_or_recast_diagnostics():
    result = fcc.get(_PID).evaluate(_carrier_point())
    tag, _, _, proxy_flags = _load_harness().tag_result(result)

    for key, value in result.diagnostics.items():
        key_l = str(key).lower()
        if "proxy" in key_l or "recast" in key_l:
            assert value is False or value is None
    assert tag == "rigorous"
    assert proxy_flags == {}


def test_legacy_adapter_proxy_and_mapping_inputs_remain_marked_proxy():
    constraint = fcc.get(_PID)
    proxy = mu_to_e_gamma_proxy_input(
        (0.10, 0.20, 0.30),
        _rotation_pmns(),
        3000.0,
        source="L001 adapter-level proxy test",
    )
    mapping = {
        "Y_N_bar": (0.10, 0.20, 0.30),
        "pmns": _rotation_pmns(),
        "M_KK": 3000.0,
        "source": "L001 adapter-level mapping proxy test",
    }

    for lepton_input in (proxy, mapping):
        result = mu_to_e_gamma_from_lepton_input(
            lepton_input,
            br_limit=constraint.anchor.budget,
            prefactor_br=constraint.anchor.prefactor_br.value,
        )
        assert result.input_kind == "MuToEGammaProxyInput"
        assert result.used_proxy is True
        assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


@pytest.mark.parametrize(
    "updates",
    [
        lambda c: {"Y_N_bar": np.array([np.nan, c.Y_N_bar[1], c.Y_N_bar[2]])},
        lambda c: {"pmns": np.eye(3, dtype=np.complex128) * 2.0},
        lambda c: {"M_KK_gev": 0.0},
        lambda c: {"Y_N_bar_matrix": c.Y_N_bar_matrix + np.eye(3) * 1.0e-3},
    ],
)
def test_malformed_carrier_degrades_to_invalid_extra(updates):
    base = _benchmark_carrier()
    malformed = _forged_carrier(base, **updates(base))

    result = fcc.get(_PID).evaluate(
        point_builder.make_point(lepton_lmfv_parameters=malformed)
    )

    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["invalid_extra"] == "lepton_lmfv_parameters"
