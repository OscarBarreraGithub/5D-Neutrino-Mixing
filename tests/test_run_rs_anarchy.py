"""Tests for the ``--epsilon-k-budget`` CLI flag on ``scripts/run_rs_anarchy.py``.

Coverage (cleanup unit C02a-code, R03-I3 task 2):

1. The flag round-trips: parsed args, ``EnsembleConfig``, and the resolved
   numerical override on the ``deltaf2`` side all agree.
2. The default (``central``) delegates to the shared deltaf2 sign-aware
   one-sigma policy by passing ``None`` for ``epsilon_k_np_budget_override``.
3. The 'low' and 'high' edges map to the shared policy's lower/upper signed
   budgets and propagate into ``evaluate_epsilon_k``.
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
    """The {central, low, high} edge contract follows the shared policy."""
    from quarkConstraints.deltaf2 import delta_f2_epsilon_k_budget_policy

    policy = delta_f2_epsilon_k_budget_policy()
    assert set(rs_module.EPSILON_K_BUDGET_EDGES.keys()) == {
        "central",
        "low",
        "high",
    }
    assert rs_module.EPSILON_K_BUDGET_EDGES["central"] is None
    assert rs_module.EPSILON_K_BUDGET_EDGES["low"] == pytest.approx(
        policy.budget_lowers_epsilon_k
    )
    assert rs_module.EPSILON_K_BUDGET_EDGES["high"] == pytest.approx(
        policy.budget_raises_epsilon_k
    )
    assert rs_module.DEFAULT_EPSILON_K_BUDGET == "central"


def test_argparser_default_is_central(rs_module):
    """Omitting the flag must select the shared-policy default."""
    parser = rs_module._build_argparser()
    args = parser.parse_args(["--output-dir", "/tmp/_rs_anarchy_dummy"])
    assert args.epsilon_k_budget == "central"


@pytest.mark.parametrize(
    "edge",
    [
        "central",
        "low",
        "high",
    ],
)
def test_argparser_accepts_three_edges(rs_module, edge):
    """Each documented edge label parses and maps to the expected value."""
    parser = rs_module._build_argparser()
    args = parser.parse_args(
        ["--output-dir", "/tmp/_rs_anarchy_dummy", "--epsilon-k-budget", edge]
    )
    assert args.epsilon_k_budget == edge
    if edge == "central":
        assert rs_module.EPSILON_K_BUDGET_EDGES[edge] is None
    else:
        from quarkConstraints.deltaf2 import delta_f2_epsilon_k_budget_policy

        policy = delta_f2_epsilon_k_budget_policy()
        expected = (
            policy.budget_lowers_epsilon_k
            if edge == "low"
            else policy.budget_raises_epsilon_k
        )
        assert rs_module.EPSILON_K_BUDGET_EDGES[edge] == pytest.approx(expected)


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
    from quarkConstraints.deltaf2 import delta_f2_epsilon_k_budget_policy

    policy = delta_f2_epsilon_k_budget_policy()
    cfg = rs_module.EnsembleConfig(
        mkk_values_GeV=(5000.0,),
        n_draws_per_tile=1,
        epsilon_k_budget_edge="low",
        epsilon_k_np_budget_override=rs_module.EPSILON_K_BUDGET_EDGES["low"],
    )
    assert cfg.epsilon_k_budget_edge == "low"
    assert cfg.epsilon_k_np_budget_override == pytest.approx(
        policy.budget_lowers_epsilon_k
    )

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
    """Passing ``None`` uses the shared sign-aware default exactly."""
    from quarkConstraints.deltaf2 import (
        delta_f2_epsilon_k_budget_policy,
        evaluate_epsilon_k,
    )

    wilsons = _make_kaon_wilsons()
    implicit = evaluate_epsilon_k(wilsons)
    explicit = evaluate_epsilon_k(
        wilsons,
        epsilon_k_np_budget_override=None,
    )
    expected_budget, expected_direction = (
        delta_f2_epsilon_k_budget_policy().selected_signed_budget(
            implicit.epsilon_k_np_signed
        )
    )
    assert implicit.epsilon_k_np_budget == pytest.approx(
        expected_budget
    )
    assert implicit.selected_budget_direction == expected_direction
    assert explicit.epsilon_k_np_budget == pytest.approx(
        implicit.epsilon_k_np_budget
    )
    assert explicit.ratio_to_budget == pytest.approx(implicit.ratio_to_budget)


@pytest.mark.parametrize("edge", ["low", "high"])
def test_evaluate_epsilon_k_override_rescales_ratio(rs_module, edge):
    """Replacing the budget rescales ratio_to_budget by 1/override exactly."""
    from quarkConstraints.deltaf2 import evaluate_epsilon_k

    override = rs_module.EPSILON_K_BUDGET_EDGES[edge]
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
        "epsilon_k_np_budget_override": rs_module.EPSILON_K_BUDGET_EDGES["high"],
    }
    rs_module._worker_init(cfg_dict)
    cfg = rs_module._GLOBAL_CFG
    assert cfg is not None
    assert cfg.epsilon_k_budget_edge == "high"
    assert cfg.epsilon_k_np_budget_override == pytest.approx(
        rs_module.EPSILON_K_BUDGET_EDGES["high"]
    )

    # And the default-edge path stays None so the shared core policy is used.
    cfg_dict["epsilon_k_budget_edge"] = "central"
    cfg_dict["epsilon_k_np_budget_override"] = None
    rs_module._worker_init(cfg_dict)
    cfg = rs_module._GLOBAL_CFG
    assert cfg is not None
    assert cfg.epsilon_k_budget_edge == "central"
    assert cfg.epsilon_k_np_budget_override is None


# ---------------------------------------------------------------------------
# C06 / R07-I3: git_sha embedded in tile_summary.json so manifest provenance
# can be verified programmatically rather than inferred from directory mtime.
# ---------------------------------------------------------------------------


def test_resolve_git_sha_inside_repo_returns_hex(rs_module):
    """Inside this checkout, the helper returns the repo HEAD SHA."""
    import re
    sha = rs_module._resolve_git_sha()
    # In-repo execution always lands in a git tree; either we get the HEAD
    # SHA (40 hex chars) or — if git is somehow unavailable on the runner —
    # the documented "unknown" sentinel. Anything else is a regression.
    assert sha == "unknown" or re.fullmatch(r"[0-9a-f]{40}", sha) is not None


def test_resolve_git_sha_outside_git_tree_returns_unknown(rs_module, tmp_path, monkeypatch):
    """Outside a git tree the helper returns the ``"unknown"`` sentinel."""
    # Point the helper's repo root at a non-git directory and confirm the
    # fallback fires without raising. The helper hard-codes _REPO_ROOT, so
    # monkeypatch the module attribute.
    monkeypatch.setattr(rs_module, "_REPO_ROOT", tmp_path)
    assert rs_module._resolve_git_sha() == "unknown"


def test_tile_summary_embeds_git_sha(rs_module, tmp_path):
    """A freshly-written tile_summary.json contains the ``git_sha`` field.

    Drives ``main()`` with the minimum payload (one tile, one draw, single
    worker) and inspects the on-disk file. The intent is provenance, not
    physics — we only assert the field is present and shaped correctly.
    """
    import json
    import re

    output_dir = tmp_path / "rs_anarchy_smoke_git_sha"
    rc = rs_module.main(
        [
            "--output-dir",
            str(output_dir),
            "--n-draws",
            "1",
            "--m-kk-tev",
            "10",
            "--n-workers",
            "1",
            "--base-seed",
            "20260525",
        ]
    )
    assert rc == 0
    summary_path = output_dir / "tile_summary.json"
    assert summary_path.exists()
    payload = json.loads(summary_path.read_text())
    assert "git_sha" in payload, (
        "tile_summary.json must embed `git_sha` (C06 / R07-I3)"
    )
    sha = payload["git_sha"]
    assert isinstance(sha, str) and sha, "git_sha must be a non-empty string"
    # In-repo execution: either a 40-char hex SHA, or the "unknown" sentinel
    # if `git` is unavailable on the runner.
    assert sha == "unknown" or re.fullmatch(r"[0-9a-f]{40}", sha) is not None


def test_legacy_tile_summary_without_git_sha_still_loads(tmp_path):
    """Backward-compat: loaders must tolerate tile_summary.json without git_sha.

    Mirrors the (minimal) shape that older writers produced and confirms a
    plain ``json.load`` round-trip still works, with ``.get('git_sha',
    'unknown')`` defaulting cleanly for re-runners that compare against a
    manifest's ``code_sha_at_run_time`` field.
    """
    import json

    legacy_path = tmp_path / "tile_summary.json"
    legacy_payload = {
        "config": {"mkk_values_GeV": [10000.0], "n_draws_per_tile": 1},
        "elapsed_seconds": 0.0,
        "tiles": [
            {
                "M_KK_GeV": 10000.0,
                "n_ok": 1,
                "n_pdg_pass": 0,
                "pdg_pass_fraction": 0.0,
            }
        ],
        "schema": "rs_anarchy_acps_v1",
    }
    legacy_path.write_text(json.dumps(legacy_payload))
    loaded = json.loads(legacy_path.read_text())
    assert "git_sha" not in loaded  # legacy shape, by construction
    # Loader contract: missing field defaults to "unknown" via .get().
    assert loaded.get("git_sha", "unknown") == "unknown"
    # All previously-required top-level keys still present.
    for key in ("config", "tiles", "schema"):
        assert key in loaded
