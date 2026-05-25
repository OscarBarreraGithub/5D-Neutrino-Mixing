"""Tests for the ``--epsilon-k-budget`` CLI flag on ``scripts/run_rs_anarchy.py``.

Coverage (cleanup unit C02a-code, R03-I3 task 2):

1. The flag round-trips: parsed args, ``EnsembleConfig``, and the resolved
   numerical override on the ``deltaf2`` side all agree.
2. The default (``central``) reproduces the pre-flag code path bit-for-bit
   by passing ``None`` for ``epsilon_k_np_budget_override`` — i.e. the
   ``deltaf2`` default ~6.7e-5 is used.
3. The 'low' and 'high' edges map to the documented numerical values
   (1e-5 and 3e-4 respectively) and propagate into ``evaluate_epsilon_k``.
4. Invalid edges are rejected by argparse.

No SLURM scan is executed here; the actual three-edge RUNA reruns are
deferred to cleanup unit C02c. This test exercises only the CLI plumbing.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO / "scripts" / "run_rs_anarchy.py"


@pytest.fixture(scope="module")
def rs_module():
    name = "_test_run_rs_anarchy_module"
    spec = importlib.util.spec_from_file_location(name, SCRIPT_PATH)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_budget_edge_table_has_three_entries(rs_module):
    """The {central, low, high} edge contract is fixed by the C02a plan."""
    assert set(rs_module.EPSILON_K_BUDGET_EDGES.keys()) == {
        "central",
        "low",
        "high",
    }
    assert rs_module.EPSILON_K_BUDGET_EDGES["central"] is None
    assert rs_module.EPSILON_K_BUDGET_EDGES["low"] == pytest.approx(1.0e-5)
    assert rs_module.EPSILON_K_BUDGET_EDGES["high"] == pytest.approx(3.0e-4)
    assert rs_module.DEFAULT_EPSILON_K_BUDGET == "central"


def test_argparser_default_is_central(rs_module):
    """Omitting the flag must reproduce the pre-flag behaviour."""
    parser = rs_module._build_argparser()
    args = parser.parse_args(["--output-dir", "/tmp/_rs_anarchy_dummy"])
    assert args.epsilon_k_budget == "central"


@pytest.mark.parametrize(
    "edge, expected_override",
    [
        ("central", None),
        ("low", 1.0e-5),
        ("high", 3.0e-4),
    ],
)
def test_argparser_accepts_three_edges(rs_module, edge, expected_override):
    """Each documented edge label parses and maps to the expected value."""
    parser = rs_module._build_argparser()
    args = parser.parse_args(
        ["--output-dir", "/tmp/_rs_anarchy_dummy", "--epsilon-k-budget", edge]
    )
    assert args.epsilon_k_budget == edge
    assert rs_module.EPSILON_K_BUDGET_EDGES[args.epsilon_k_budget] == (
        pytest.approx(expected_override)
        if expected_override is not None
        else None
    )


def test_argparser_rejects_unknown_edge(rs_module):
    """Anything outside the contract must be rejected, not silently ignored."""
    parser = rs_module._build_argparser()
    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                "--output-dir",
                "/tmp/_rs_anarchy_dummy",
                "--epsilon-k-budget",
                "fifth-edge",
            ]
        )


def test_ensemble_config_round_trip(rs_module):
    """EnsembleConfig stores both the edge label and the resolved override."""
    cfg = rs_module.EnsembleConfig(
        mkk_values_GeV=(5000.0,),
        n_draws_per_tile=1,
        epsilon_k_budget_edge="low",
        epsilon_k_np_budget_override=rs_module.EPSILON_K_BUDGET_EDGES["low"],
    )
    assert cfg.epsilon_k_budget_edge == "low"
    assert cfg.epsilon_k_np_budget_override == pytest.approx(1.0e-5)

    cfg_default = rs_module.EnsembleConfig(
        mkk_values_GeV=(5000.0,), n_draws_per_tile=1
    )
    assert cfg_default.epsilon_k_budget_edge == "central"
    assert cfg_default.epsilon_k_np_budget_override is None


# ---------------------------------------------------------------------------
# Cross-check: the override propagates correctly into the deltaf2 evaluator.
# ---------------------------------------------------------------------------

def _make_kaon_wilsons(c_imag: float = 1.0):
    """Build a minimal DeltaF2WilsonCoefficients for the epsilon_K key.

    Pure-imaginary VLL coefficient so Im(M_12) is non-trivial; the LR
    coefficients are zero so the budget rescaling is the only handle on
    the ratio.
    """
    from quarkConstraints.deltaf2 import (
        DEFAULT_DELTA_F2_INPUTS_V1,
        DeltaF2WilsonCoefficients,
    )

    eps_input = next(
        item for item in DEFAULT_DELTA_F2_INPUTS_V1 if item.key == "epsilon_k"
    )
    return DeltaF2WilsonCoefficients(
        input=eps_input,
        M_KK=5000.0,
        matching_scale=5000.0,
        left_coupling=0j,
        right_coupling=0j,
        c1_vll=complex(0.0, c_imag),
        c1_vrr=0j,
        c4_lr=0j,
        c5_lr=0j,
    )


def test_evaluate_epsilon_k_default_budget_matches_implicit():
    """Passing ``None`` reproduces the historical default exactly."""
    from quarkConstraints.deltaf2 import (
        EPSILON_K_EXP,
        EPSILON_K_SM,
        evaluate_epsilon_k,
    )

    wilsons = _make_kaon_wilsons()
    implicit = evaluate_epsilon_k(wilsons)
    explicit = evaluate_epsilon_k(
        wilsons,
        epsilon_k_np_budget_override=None,
    )
    assert implicit.epsilon_k_np_budget == pytest.approx(
        abs(EPSILON_K_EXP - EPSILON_K_SM)
    )
    assert explicit.epsilon_k_np_budget == pytest.approx(
        implicit.epsilon_k_np_budget
    )
    assert explicit.ratio_to_budget == pytest.approx(implicit.ratio_to_budget)


@pytest.mark.parametrize("override", [1.0e-5, 3.0e-4])
def test_evaluate_epsilon_k_override_rescales_ratio(override):
    """Replacing the budget rescales ratio_to_budget by 1/override exactly."""
    from quarkConstraints.deltaf2 import evaluate_epsilon_k

    wilsons = _make_kaon_wilsons()
    baseline = evaluate_epsilon_k(wilsons)
    overridden = evaluate_epsilon_k(
        wilsons, epsilon_k_np_budget_override=override
    )
    assert overridden.epsilon_k_np_budget == pytest.approx(override)
    # epsilon_k_np is the numerator and must NOT change.
    assert overridden.epsilon_k_np == pytest.approx(baseline.epsilon_k_np)
    expected_ratio = baseline.epsilon_k_np / override
    assert overridden.ratio_to_budget == pytest.approx(expected_ratio)


def test_evaluate_epsilon_k_rejects_nonpositive_override():
    from quarkConstraints.deltaf2 import evaluate_epsilon_k

    wilsons = _make_kaon_wilsons()
    with pytest.raises(ValueError, match="must be positive"):
        evaluate_epsilon_k(wilsons, epsilon_k_np_budget_override=0.0)
    with pytest.raises(ValueError, match="must be positive"):
        evaluate_epsilon_k(wilsons, epsilon_k_np_budget_override=-1.0e-5)


# ---------------------------------------------------------------------------
# End-to-end: the multiprocessing-worker config rehydration in
# ``_worker_init`` must round-trip the epsilon_K override. This guards against
# the (very real) regression class where new EnsembleConfig fields are
# threaded through main() but silently dropped when rebuilt inside workers,
# producing a scan that silently ignores ``--epsilon-k-budget``.
# ---------------------------------------------------------------------------

def test_worker_init_round_trips_budget_override(rs_module, monkeypatch):
    """`_worker_init` must rebuild EnsembleConfig with the override intact."""
    monkeypatch.setattr(rs_module, "_load_pdg_targets", lambda: {})
    cfg_dict = {
        "mkk_values_GeV": [5000.0],
        "n_draws_per_tile": 1,
        "c_Q": list(rs_module.DEFAULT_C_Q),
        "c_u": list(rs_module.DEFAULT_C_U),
        "c_d": list(rs_module.DEFAULT_C_D),
        "xi_KK": rs_module.DEFAULT_XI_KK,
        "k_GeV": rs_module.DEFAULT_K_GEV,
        "v_GeV": rs_module.DEFAULT_V_GEV,
        "y_half_range": rs_module.DEFAULT_Y_HALF_RANGE,
        "y_floor": rs_module.DEFAULT_Y_FLOOR,
        "y_prior": rs_module.DEFAULT_Y_PRIOR,
        "y_sigma": rs_module.DEFAULT_Y_SIGMA,
        "y_trunc_sigma": rs_module.DEFAULT_Y_TRUNC_SIGMA,
        "pdg_mass_factor": rs_module.PDG_MASS_FACTOR_TOL,
        "pdg_ckm_factor": rs_module.PDG_CKM_FACTOR_TOL,
        "pdg_j_factor": rs_module.PDG_J_FACTOR_TOL,
        "base_seed": 20260506,
        "epsilon_k_budget_edge": "high",
        "epsilon_k_np_budget_override": 3.0e-4,
    }
    rs_module._worker_init(cfg_dict)
    cfg = rs_module._GLOBAL_CFG
    assert cfg is not None
    assert cfg.epsilon_k_budget_edge == "high"
    assert cfg.epsilon_k_np_budget_override == pytest.approx(3.0e-4)

    # And the default-edge path stays None (the legacy bit-for-bit case).
    cfg_dict["epsilon_k_budget_edge"] = "central"
    cfg_dict["epsilon_k_np_budget_override"] = None
    rs_module._worker_init(cfg_dict)
    cfg = rs_module._GLOBAL_CFG
    assert cfg is not None
    assert cfg.epsilon_k_budget_edge == "central"
    assert cfg.epsilon_k_np_budget_override is None
