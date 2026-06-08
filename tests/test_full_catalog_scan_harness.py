from __future__ import annotations

import importlib.util
import ast
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from flavor_catalog_constraints.base import ConstraintResult, Severity


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "scripts" / "run_full_catalog_scan.py"


def _load_harness():
    spec = importlib.util.spec_from_file_location("run_full_catalog_scan", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_hard_veto_policy_separates_rigorous_proxy_partial_stub_and_advisory():
    harness = _load_harness()
    payload = harness._classify_results(
        {
            "R001": ConstraintResult(
                process_id="R001",
                severity=Severity.HARD,
                passes=False,
                ratio=2.0,
                diagnostics={"evaluated": True, "matching_status": "tree_level_exact"},
            ),
            "P001": ConstraintResult(
                process_id="P001",
                severity=Severity.HARD,
                passes=False,
                ratio=3.0,
                diagnostics={
                    "evaluated": True,
                    "used_proxy": True,
                    "needs_human_physics": "documented mass-limit proxy",
                },
            ),
            "H001": ConstraintResult(
                process_id="H001",
                severity=Severity.HARD,
                passes=False,
                ratio=4.0,
                diagnostics={
                    "evaluated": True,
                    "needs_human_physics": "NEEDS-HUMAN-PHYSICS: partial EFT matching",
                },
            ),
            "S001": ConstraintResult(
                process_id="S001",
                severity=Severity.HARD,
                passes=False,
                diagnostics={"evaluated": False, "missing_extra": "rs_loop_input"},
            ),
            "SOFT1": ConstraintResult(
                process_id="SOFT1",
                severity=Severity.SOFT,
                passes=False,
                ratio=9.0,
            ),
            "INFO1": ConstraintResult(
                process_id="INFO1",
                severity=Severity.INFO,
                passes=False,
                ratio=9.0,
            ),
        }
    )

    assert payload["excluded_by_rigorous"] == ["R001"]
    assert payload["excluded_by_proxy"] == ["P001"]
    assert payload["hard_not_evaluated"] == ["H001", "S001"]
    assert payload["survives_all_HARD_strict"] is False
    assert payload["survives_all_HARD_inclusive"] is False
    assert "SOFT:SOFT1" in payload["advisory_flags"]
    assert "INFO:INFO1" in payload["advisory_flags"]
    assert payload["constraints"]["H001"]["tag"] == "partial"
    assert payload["constraints"]["S001"]["tag"] == "stub"


def test_evaluated_proxy_with_unevaluated_subobservable_note_stays_proxy():
    harness = _load_harness()
    result = ConstraintResult(
        process_id="B013",
        severity=Severity.HARD,
        passes=False,
        ratio=2.0,
        notes="A_Delta and S_phi gamma are not evaluated by this branching proxy.",
        diagnostics={
            "evaluated": True,
            "mass_proxy": "exclusive C7 normalization proxy",
            "needs_human_physics": "documented C7 proxy",
        },
    )

    tag, _, _, proxy_flags = harness.tag_result(result)

    assert tag == "proxy"
    assert proxy_flags["mass_proxy"] == "exclusive C7 normalization proxy"


def test_resolved_proxy_wording_in_matching_status_remains_rigorous():
    harness = _load_harness()
    result = ConstraintResult(
        process_id="B022",
        severity=Severity.HARD,
        passes=False,
        ratio=1.1,
        diagnostics={
            "evaluated": True,
            "rs_semileptonic_nunu_matching_status": (
                "Phase-4d rigorous active-neutrino Wilson; old one-Z-like "
                "proxy resolved and proxy is not used."
            ),
        },
    )

    tag, matching_status, needs_human, proxy_flags = harness.tag_result(result)

    assert tag == "rigorous"
    assert "proxy resolved" in matching_status
    assert needs_human is None
    assert proxy_flags == {}


def test_deferred_custodial_refinement_without_needs_human_stays_rigorous():
    harness = _load_harness()
    result = ConstraintResult(
        process_id="T010",
        severity=Severity.HARD,
        passes=False,
        ratio=2.0,
        diagnostics={
            "evaluated": True,
            "minimal_rs_tree_complete": True,
            "minimal_rs_tree_veto_ready": True,
            "custodial_variant_deferred": True,
            "custodial_variant_deferred_note": "Deferred refinement: custodial variant.",
        },
    )

    tag, matching_status, needs_human, proxy_flags = harness.tag_result(result)

    assert tag == "rigorous"
    assert matching_status is None
    assert needs_human is None
    assert proxy_flags == {}


def test_tile_seed_stride_and_config_hash_are_deterministic():
    harness = _load_harness()
    cfg = harness.ScanConfig(
        mkk_values_gev=(1000.0, 3000.0),
        n_draws_per_tile=5,
        xi_kk=2.0,
        base_seed=11,
        tile_seed_stride=17,
    )

    tiles = harness._build_tiles(cfg)
    same_hash = harness._config_hash(cfg)
    changed_hash = harness._config_hash(
        harness.ScanConfig(mkk_values_gev=(1000.0, 3000.0), n_draws_per_tile=6)
    )

    assert [tile.seed for tile in tiles] == [11, 28]
    assert [tile.lambda_ir_gev for tile in tiles] == [500.0, 1500.0]
    assert same_hash == harness._config_hash(cfg)
    assert same_hash != changed_hash


def test_default_tile_seed_stride_covers_large_smoke_tiles():
    harness = _load_harness()
    cfg = harness.ScanConfig(
        mkk_values_gev=(1000.0, 3000.0, 6000.0, 10000.0),
        n_draws_per_tile=10_000,
    )
    tiles = harness._build_tiles(cfg)

    assert harness.DEFAULT_TILE_SEED_STRIDE > cfg.n_draws_per_tile
    ranges = [
        set(range(tile.seed, tile.seed + tile.n_draws))
        for tile in tiles
    ]
    assert sum(len(seed_range) for seed_range in ranges) == len(set().union(*ranges))


def test_completed_summary_resume_requires_complete_matching_hash_and_draws(tmp_path):
    harness = _load_harness()
    path = tmp_path / "tile-00000.summary.json"
    path.write_text(
        json.dumps({"complete": True, "config_hash": "abc", "n_requested": 10, "n_rows": 10}),
        encoding="utf-8",
    )

    assert harness._completed_summary(path, expected_hash="abc", expected_draws=10) is not None
    assert harness._completed_summary(path, expected_hash="def", expected_draws=10) is None
    assert harness._completed_summary(path, expected_hash="abc", expected_draws=9) is None

    path.write_text(
        json.dumps({"complete": True, "config_hash": "abc", "n_requested": 10, "n_rows": 9}),
        encoding="utf-8",
    )
    assert harness._completed_summary(path, expected_hash="abc", expected_draws=10) is None

    path.write_text(
        json.dumps({"complete": False, "config_hash": "abc", "n_requested": 10, "n_rows": 10}),
        encoding="utf-8",
    )
    assert harness._completed_summary(path, expected_hash="abc", expected_draws=10) is None


def test_quark_only_config_payload_preserves_full_mode_hash_surface():
    harness = _load_harness()
    full_cfg = harness.ScanConfig(
        mkk_values_gev=(1000.0, 3000.0),
        n_draws_per_tile=5,
        xi_kk=2.0,
        base_seed=11,
        tile_seed_stride=17,
    )
    quark_cfg = harness.ScanConfig(
        mkk_values_gev=(1000.0, 3000.0),
        n_draws_per_tile=5,
        quark_only=True,
        xi_kk=2.0,
        base_seed=11,
        tile_seed_stride=17,
    )

    assert "quark_only" not in harness._config_payload(full_cfg)
    assert harness._config_payload(quark_cfg)["quark_only"] is True
    assert harness._config_hash(full_cfg) != harness._config_hash(quark_cfg)


def test_quark_only_allowlist_matches_candidate_verification_and_drops_deferred_leptons():
    harness = _load_harness()
    expected_collider_in = (
        "CR001",
        "CR002",
        "CR003",
        "CR004",
        "CR007",
        "CR008",
        "CR010",
        "CR012",
        "CR013",
    )
    expected_collider_deferred = ("CR005", "CR006", "CR009", "CR011", "CR014")

    assert set(harness.QUARK_ONLY_CANDIDATE_IDS) == (
        set(harness.QUARK_ONLY_ALLOWLIST_IDS)
        | set(harness.QUARK_ONLY_DEFERRED_LEPTON_FOLLOWUP)
    )
    assert len(harness.QUARK_ONLY_ALLOWLIST_IDS) == 46
    assert set(expected_collider_in) <= harness.QUARK_ONLY_ALLOWLIST_SET
    assert set(expected_collider_deferred).isdisjoint(harness.QUARK_ONLY_ALLOWLIST_SET)
    assert harness.QUARK_ONLY_DEFERRED_LEPTON_FOLLOWUP == (
        "EW002",
        *expected_collider_deferred,
    )
    assert "EW002" not in harness.QUARK_ONLY_ALLOWLIST_SET

    actual = _candidate_get_extra_usage(harness)
    for pid in harness.QUARK_ONLY_ALLOWLIST_IDS:
        assert tuple(actual[pid]) == harness.QUARK_ONLY_ALLOWLIST_EXTRAS[pid]
        forbidden = set(actual[pid]) & harness.QUARK_ONLY_FORBIDDEN_EXTRAS
        if pid not in harness.QUARK_ONLY_OPTIONAL_LEPTON_DIAGNOSTIC_IDS:
            assert forbidden == set()
    for pid in harness.QUARK_ONLY_DEFERRED_LEPTON_FOLLOWUP:
        assert actual[pid] == harness.QUARK_ONLY_DEFERRED_EXTRAS[pid]
    for pid in expected_collider_in:
        assert not set(actual[pid]) & harness.QUARK_ONLY_FORBIDDEN_EXTRAS


def test_quark_only_collider_allowlist_is_active_proxy_and_monotonic():
    from flavor_catalog_constraints import point_builder

    harness = _load_harness()
    collider_ids = (
        "CR001",
        "CR002",
        "CR003",
        "CR004",
        "CR007",
        "CR008",
        "CR010",
        "CR012",
        "CR013",
    )

    def point(m_kk_gev: float):
        return point_builder.make_point(
            kk_ew_mass_gev=m_kk_gev,
            kk_gluon_mass_gev=m_kk_gev,
            quark_mass_basis_couplings=SimpleNamespace(M_KK=m_kk_gev, xi_KK=1.0),
        )

    low_results = harness._evaluate_constraint_ids(point(1000.0), collider_ids)
    high_results = harness._evaluate_constraint_ids(point(50000.0), collider_ids)
    low_payload = harness._classify_results(low_results)
    high_payload = harness._classify_results(high_results)

    assert low_payload["excluded_by_proxy"]
    assert low_payload["survives_all_HARD_inclusive"] is False
    assert high_payload["excluded_by_proxy"] == []
    assert high_payload["survives_all_HARD_inclusive"] is True

    for pid in collider_ids:
        result = low_results[pid]
        high_result = high_results[pid]
        diagnostics = dict(result.diagnostics)
        assert result.severity is Severity.HARD
        assert low_payload["constraints"][pid]["evaluated"] is True
        assert low_payload["constraints"][pid]["active"] is True
        assert low_payload["constraints"][pid]["tag"] == "proxy"
        assert "missing_extra" not in diagnostics
        assert "missing_extras" not in diagnostics
        assert result.ratio is not None
        assert high_result.ratio is not None
        assert result.ratio > high_result.ratio
        assert high_result.passes is True
        if pid in {"CR001", "CR007", "CR012", "CR013"}:
            assert diagnostics["mass_source"] in {
                "kk_gluon_mass_gev",
                "kk_ew_mass_gev",
            }
        else:
            assert diagnostics["mass_source"] == "quark_mass_basis_couplings.M_KK"


def test_quark_only_evaluate_draw_skips_leptons_filters_allowlist_and_is_deterministic(monkeypatch):
    harness = _load_harness()
    cfg = harness.ScanConfig(
        mkk_values_gev=(3000.0,),
        n_draws_per_tile=1,
        quark_only=True,
        quark_fit_max_nfev=1,
    )
    tile = harness.TileSpec(
        tile_id=0,
        mkk_gev=3000.0,
        lambda_ir_gev=3000.0 / cfg.xi_kk,
        n_draws=1,
        seed=123,
    )
    spectrum = object()

    @dataclass
    class _Cache:
        max_a_rel_err: float = 0.0

    @dataclass
    class _BulkState:
        c_Q: np.ndarray
        c_u: np.ndarray
        c_d: np.ndarray

    @dataclass
    class _Point:
        Y_u: np.ndarray
        Y_d: np.ndarray

    @dataclass
    class _FitResult:
        bulk_state: _BulkState
        point: _Point
        score: float
        residual_norm: float
        mass_residuals_up: np.ndarray
        mass_residuals_down: np.ndarray
        ckm_residuals: np.ndarray

    @dataclass
    class _Solution:
        result: _FitResult
        success: bool = True
        message: str = "ok"
        nfev: int = 1
        initial_score: float = 1.0

    def fake_fit(*args, **kwargs):
        assert kwargs["fit_orientation"] is True
        assert kwargs["r"] == cfg.quark_fit_r
        fit_result = _FitResult(
            bulk_state=_BulkState(
                c_Q=np.array([0.61, 0.62, 0.63]),
                c_u=np.array([0.64, 0.65, 0.66]),
                c_d=np.array([0.67, 0.68, 0.69]),
            ),
            point=_Point(
                Y_u=np.ones((3, 3), dtype=np.complex128),
                Y_d=np.ones((3, 3), dtype=np.complex128),
            ),
            score=0.0,
            residual_norm=0.0,
            mass_residuals_up=np.zeros(3),
            mass_residuals_down=np.zeros(3),
            ckm_residuals=np.zeros(4),
        )
        return _Solution(result=fit_result)

    def forbidden(*_args, **_kwargs):
        raise AssertionError("lepton stage must not run in quark-only mode")

    def fake_build(fit_result, **kwargs):
        assert kwargs["spectrum"] is spectrum
        assert kwargs["rs_ew_cache"] is cache
        assert kwargs["lepton_yukawa_result"] is None
        assert kwargs["include_charged_current"] is False
        assert kwargs["include_fermion_kk_mixing"] is True
        assert kwargs["include_higgs_yukawas"] is False
        return harness.point_builder.make_point(
            raw=kwargs["raw"],
            kk_ew_mass_gev=3000.0,
            rs_ew_spectrum=spectrum,
            rs_ew_couplings=object(),
        )

    @dataclass
    class _Couplings:
        M_KK: float

    def fake_couplings(_fit_result, **kwargs):
        return _Couplings(M_KK=float(kwargs["M_KK"]))

    captured_ids: list[tuple[str, ...]] = []

    def fake_evaluate_ids(_point, process_ids):
        ids = tuple(process_ids)
        captured_ids.append(ids)
        return {
            pid: ConstraintResult(
                process_id=pid,
                severity=Severity.HARD,
                passes=pid != "K001",
                ratio=2.0 if pid == "K001" else 0.1,
                diagnostics={"evaluated": True},
            )
            for pid in ids
        }

    cache = _Cache()
    monkeypatch.setattr(harness, "fit_quark_sector", fake_fit)
    monkeypatch.setattr(harness, "_draw_lepton_inputs", forbidden)
    monkeypatch.setattr(harness, "compute_all_yukawas", forbidden)
    monkeypatch.setattr(harness, "_require_perturbative_leptons", forbidden)
    monkeypatch.setattr(harness.point_builder, "build_from_rs_ew_inputs", fake_build)
    monkeypatch.setattr(harness, "compute_quark_kk_gluon_couplings", fake_couplings)
    monkeypatch.setattr(harness, "_evaluate_constraint_ids", fake_evaluate_ids)

    kwargs = dict(
        tile=tile,
        draw_idx=0,
        draw_seed=123,
        cfg=cfg,
        spectrum=spectrum,
        overlap_cache=cache,
        provenance={"mode": "quark_only", "quark_fit_r": cfg.quark_fit_r},
        registry_count=harness.EXPECTED_REGISTRY_COUNT,
        config_hash="abc",
    )
    row1 = harness._evaluate_draw(np.random.default_rng(123), **kwargs)
    row2 = harness._evaluate_draw(np.random.default_rng(123), **kwargs)

    assert row1["skipped"] is False
    assert row1["mode"] == "quark_only"
    assert row1["quark_fit_r"] == cfg.quark_fit_r
    assert row1["params"]["quark_fit_r"] == cfg.quark_fit_r
    assert row1["provenance"]["quark_fit_r"] == cfg.quark_fit_r
    assert row1["provenance"]["config_hash"] == "abc"
    assert "lepton_inputs" not in row1["params"]
    assert "fitted_up_yukawa_singular_values" in row1["fit_diagnostics"]
    assert "fitted_down_yukawa_singular_values" in row1["fit_diagnostics"]
    np.testing.assert_allclose(
        row1["fit_diagnostics"]["fitted_up_yukawa_singular_values"],
        [0.0, 0.0, 3.0],
        atol=1.0e-12,
    )
    assert set(row1["constraints"]) == harness.QUARK_ONLY_ALLOWLIST_SET
    assert "EW002" not in row1["constraints"]
    assert captured_ids[-1] == harness.QUARK_ONLY_ALLOWLIST_IDS
    assert row1["params"] == row2["params"]
    assert row1["constraints"] == row2["constraints"]


def test_quark_only_constraint_tallies_count_evaluated_active_failed_and_vetoed():
    harness = _load_harness()
    tallies = harness._empty_constraint_tallies(("K001", "EW003"))
    counters = {}
    from collections import Counter

    row = {
        "skipped": False,
        "survives_all_HARD_strict": False,
        "survives_all_HARD_inclusive": False,
        "excluded_by_rigorous": ["K001"],
        "excluded_by_proxy": [],
        "hard_not_evaluated": [],
        "constraints": {
            "K001": {
                "passes": False,
                "severity": "HARD",
                "ratio": 2.0,
                "active": True,
                "evaluated": True,
                "tag": "rigorous",
            },
            "EW003": {
                "passes": False,
                "severity": "SOFT",
                "ratio": 1.2,
                "active": True,
                "evaluated": True,
                "tag": "partial",
            },
        },
    }
    harness._accumulate_row(
        row,
        counters=Counter(counters),
        hard_veto_rigorous=Counter(),
        hard_veto_proxy=Counter(),
        hard_not_evaluated=Counter(),
        tag_counts=Counter(),
        exception_ids=Counter(),
        constraint_tallies=tallies,
    )
    final = harness._finalize_constraint_tallies(tallies)

    assert final["K001"]["evaluated"] == 1
    assert final["K001"]["active"] == 1
    assert final["K001"]["failed"] == 1
    assert final["K001"]["vetoed"] == 1
    assert final["EW003"]["failed"] == 1
    assert final["EW003"]["vetoed"] == 0


def test_deltaf2_bucket1_ratios_relax_with_high_mkk():
    from flavor_catalog_constraints import point_builder, registry
    from quarkConstraints.benchmarks import default_quark_targets, default_spurion_seed
    from quarkConstraints.couplings import compute_quark_kk_gluon_couplings
    from quarkConstraints.fit import fit_quark_sector

    fit = fit_quark_sector(
        default_quark_targets(),
        seed=default_spurion_seed(),
        overall_scale=3.0,
        max_nfev=120,
    ).result
    low = compute_quark_kk_gluon_couplings(fit, M_KK=1500.0, xi_KK=1.0, g_s_star=None)
    high = compute_quark_kk_gluon_couplings(fit, M_KK=15000.0, xi_KK=1.0, g_s_star=None)
    low_point = point_builder.make_point(
        quark_mass_basis_couplings=low,
        kk_gluon_mass_gev=1500.0,
    )
    high_point = point_builder.make_point(
        quark_mass_basis_couplings=high,
        kk_gluon_mass_gev=15000.0,
    )

    k001_low = registry.get("K001").evaluate(low_point)
    k001_high = registry.get("K001").evaluate(high_point)
    k002_low = registry.get("K002").evaluate(low_point)
    k002_high = registry.get("K002").evaluate(high_point)

    assert k001_low.ratio > k001_high.ratio
    assert k002_low.ratio > k002_high.ratio
    assert k001_low.passes is False
    assert k001_high.passes is True


def _candidate_get_extra_usage(harness):
    constants_by_id = {}
    usage = {}
    for pid in harness.QUARK_ONLY_CANDIDATE_IDS:
        paths = list((REPO_ROOT / "flavor_catalog_constraints").glob(f"**/{pid}.py"))
        assert len(paths) == 1
        tree = ast.parse(paths[0].read_text(encoding="utf-8"))
        constants = {}
        for node in tree.body:
            if isinstance(node, ast.Assign):
                for target in node.targets:
                    if isinstance(target, ast.Name):
                        constants[target.id] = _literal_or_none(node.value)
        constants_by_id[pid] = constants
        extras: list[str] = []

        class _Visitor(ast.NodeVisitor):
            def __init__(self):
                self.env = {}

            def visit_For(self, node):
                if isinstance(node.target, ast.Name) and isinstance(node.iter, ast.Name):
                    values = constants.get(node.iter.id)
                    if isinstance(values, tuple):
                        old = self.env.get(node.target.id, None)
                        had_old = node.target.id in self.env
                        for value in values:
                            self.env[node.target.id] = (str(value),)
                            for child in node.body:
                                self.visit(child)
                        if had_old:
                            self.env[node.target.id] = old
                        else:
                            self.env.pop(node.target.id, None)
                        for child in node.orelse:
                            self.visit(child)
                        return
                self.generic_visit(node)

            def visit_Call(self, node):
                if isinstance(node.func, ast.Attribute) and node.func.attr == "get_extra":
                    if node.args:
                        extras.extend(_resolve_extra_arg(node.args[0], constants, self.env))
                self.generic_visit(node)

        _Visitor().visit(tree)
        usage[pid] = tuple(dict.fromkeys(extras))
    return usage


def _literal_or_none(node):
    try:
        return ast.literal_eval(node)
    except (ValueError, SyntaxError):
        return None


def _resolve_extra_arg(node, constants, env):
    if isinstance(node, ast.Constant):
        return [str(node.value)]
    if isinstance(node, ast.Name):
        if node.id in env:
            return list(env[node.id])
        value = constants.get(node.id)
        if isinstance(value, str):
            return [value]
        if isinstance(value, tuple):
            return [str(item) for item in value]
        return [node.id]
    return [ast.unparse(node)]
