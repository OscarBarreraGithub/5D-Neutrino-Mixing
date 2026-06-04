"""Quark-sector flavor constraints for warped Randall-Sundrum models.

The top-level package uses lazy attribute loading so importing
`quarkConstraints.paper_0710_1869` does not eagerly pull in the repo-v1
observable stack as a side effect.
"""

from __future__ import annotations

from importlib import import_module

_EXPORTS = {
    "AlignmentDiagnostics": (".proxies", "AlignmentDiagnostics"),
    "BulkMassMap": (".model", "BulkMassMap"),
    "DELTA_F2_INPUT_BUNDLE": (".deltaf2", "DELTA_F2_INPUT_BUNDLE"),
    "DELTA_F2_OPERATOR_CONVENTION": (".deltaf2", "DELTA_F2_OPERATOR_CONVENTION"),
    "DEFAULT_QUARK_BENCHMARK_H_RS_MAX": (".scales", "DEFAULT_QUARK_BENCHMARK_H_RS_MAX"),
    "DEFAULT_QUARK_BENCHMARK_XI_KK": (".scales", "DEFAULT_QUARK_BENCHMARK_XI_KK"),
    "DEFAULT_QUARK_MAX_MISALIGNMENT_STRESS": (
        ".scales",
        "DEFAULT_QUARK_MAX_MISALIGNMENT_STRESS",
    ),
    "DEFAULT_QUARK_PAPER_H_RS_MAX": (".scales", "DEFAULT_QUARK_PAPER_H_RS_MAX"),
    "DEFAULT_QUARK_TARGET_SCALE_GEV": (".scales", "DEFAULT_QUARK_TARGET_SCALE_GEV"),
    "DEFAULT_QUARK_XI_KK": (".scales", "DEFAULT_QUARK_XI_KK"),
    "DeltaF2ConstraintSummary": (".deltaf2", "DeltaF2ConstraintSummary"),
    "DeltaF2Input": (".deltaf2", "DeltaF2Input"),
    "DeltaF2ObservableSummary": (".deltaf2", "DeltaF2ObservableSummary"),
    "DeltaF2WilsonSet": (".deltaf2", "DeltaF2WilsonSet"),
    "GAUGE_KK_ROOT_NN": (".scales", "GAUGE_KK_ROOT_NN"),
    "SPIN2_GRAVITON_KK_ROOT": (".scales", "SPIN2_GRAVITON_KK_ROOT"),
    "ProxySummary": (".proxies", "ProxySummary"),
    "QuarkBenchmark": (".benchmarks", "QuarkBenchmark"),
    "QuarkBulkState": (".model", "QuarkBulkState"),
    "QuarkFitResult": (".fit", "QuarkFitResult"),
    "QuarkFitSeed": (".fit", "QuarkFitSeed"),
    "QuarkFitSolution": (".fit", "QuarkFitSolution"),
    "QuarkMassBasisCouplings": (".couplings", "QuarkMassBasisCouplings"),
    "QuarkScanConfig": (".scan", "QuarkScanConfig"),
    "QuarkSpurionPoint": (".model", "QuarkSpurionPoint"),
    "QuarkTargets": (".benchmarks", "QuarkTargets"),
    "RSHiggsYukawaCouplings": (".rs_higgs_yukawas", "RSHiggsYukawaCouplings"),
    "RSEWSpectrum": (".rs_ew_spectrum", "RSEWSpectrum"),
    "RotationParameters": (".model", "RotationParameters"),
    "benchmark_fit_summary": (".validation", "benchmark_fit_summary"),
    "benchmark_plot_data": (".validation", "benchmark_plot_data"),
    "benchmark_solution": (".validation", "benchmark_solution"),
    "benchmark_spurion_input": (".benchmarks", "benchmark_spurion_input"),
    "bulk_mass_map_comparison_data": (".validation", "bulk_mass_map_comparison_data"),
    "build_mass_matrices": (".fit", "build_mass_matrices"),
    "build_mfv_point_from_singular_values": (".model", "build_mfv_point_from_singular_values"),
    "build_rs_higgs_yukawas": (".rs_higgs_yukawas", "build_rs_higgs_yukawas"),
    "build_spurion_matrix": (".model", "build_spurion_matrix"),
    "ckm_like_unitary": (".model", "ckm_like_unitary"),
    "ckm_observables": (".fit", "ckm_observables"),
    "compute_alignment_diagnostics": (".proxies", "compute_alignment_diagnostics"),
    "compute_delta_f2_wilsons": (".deltaf2", "compute_delta_f2_wilsons"),
    "compute_proxy_summary": (".proxies", "compute_proxy_summary"),
    "compute_quark_kk_gluon_couplings": (".couplings", "compute_quark_kk_gluon_couplings"),
    "default_delta_f2_inputs": (".deltaf2", "default_delta_f2_inputs"),
    "default_quark_benchmark": (".benchmarks", "default_quark_benchmark"),
    "default_quark_m_kk_from_lambda_ir": (".scales", "default_quark_m_kk_from_lambda_ir"),
    "default_quark_targets": (".benchmarks", "default_quark_targets"),
    "default_spurion_seed": (".benchmarks", "default_spurion_seed"),
    "derive_bulk_state": (".model", "derive_bulk_state"),
    "evaluate_delta_f2_constraints": (".deltaf2", "evaluate_delta_f2_constraints"),
    "evaluate_quark_fit": (".fit", "evaluate_quark_fit"),
    "fit_quark_sector": (".fit", "fit_quark_sector"),
    "spin2_graviton_mass_from_lambda_ir": (
        ".scales",
        "spin2_graviton_mass_from_lambda_ir",
    ),
    "jarlskog_invariant": (".fit", "jarlskog_invariant"),
    "r_sweep_plot_data": (".validation", "r_sweep_plot_data"),
    "rough_sm_targets": (".benchmarks", "rough_sm_targets"),
    "run_quark_scan": (".scan", "run_quark_scan"),
    "spurion_svd_summary": (".model", "spurion_svd_summary"),
    "summarize_flavor_diagnostics": (".proxies", "summarize_flavor_diagnostics"),
}

__all__ = list(_EXPORTS)


def __getattr__(name: str):
    try:
        module_name, attr_name = _EXPORTS[name]
    except KeyError as exc:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from exc
    module = import_module(module_name, __name__)
    value = getattr(module, attr_name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
