"""Export a self-contained current-code collaborator CSV near 5 TeV.

This exporter intentionally reads the small current-code rerun in
``scan_outputs/collaborator_5tev_current_20260518T232900`` rather than the
older May 6 dense accepted table.  The rerun was made after the PDG-mass,
benchmark-seed, and Delta-F=2/Wilson updates, so the exported ratios are
recomputed at the physical first-KK-gluon mass with the current backend and
QCD evolution from the physical matching scale.

The output includes the full complex CKM matrix, fitted masses, target values,
bulk masses, zero-mode overlaps, bulk-basis Yukawas, mass matrices, fit
diagnostics, physical-mass constraint ratios, and compact convention metadata.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace
from typing import Iterable

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from quarkConstraints.benchmarks import default_quark_targets, default_spurion_seed  # noqa: E402
from quarkConstraints.fit import fit_quark_sector, jarlskog_invariant  # noqa: E402
from quarkConstraints.modern.bridge_artifacts import (
    build_modern_point_bridge_artifact,  # noqa: E402
)
from quarkConstraints.modern.couplings import build_modern_point_couplings  # noqa: E402
from quarkConstraints.modern.evaluation import MODERN_POINT_EVALUATION_SCHEMA_ID  # noqa: E402
from quarkConstraints.modern.matching import build_modern_point_matching  # noqa: E402
from quarkConstraints.modern.phenomenology import (
    build_modern_point_phenomenology_artifact,  # noqa: E402
)
from quarkConstraints.scales import (  # noqa: E402
    DEFAULT_QUARK_FIT_SCALE_GEV,
    DEFAULT_QUARK_TARGET_SCALE_GEV,
    GAUGE_KK_ROOT_NN,
)

RUN_DIR = REPO_ROOT / "scan_outputs/collaborator_5tev_current_20260518T232900"
MERGED_JSONL = RUN_DIR / "merged/results.jsonl"
CONFIG_JSON = RUN_DIR / "config.json"
MERGE_MANIFEST = RUN_DIR / "merged/manifest.json"
VERIFICATION_JSON = RUN_DIR / "merged/verification.json"
DEFAULT_OUTPUT = REPO_ROOT / "artifacts/collaborator_5tev_points.csv"

DEFAULT_COORDINATES = (
    (0.022995139907954716, 1.5),
    (0.021949975309861126, 1.736842105263158),
    (0.02409007080517565, 6.0),
    (0.025237137661320411, 6.0),
    (0.026438822969320579, 1.736842105263158),
)

UP_LABELS = ("u", "c", "t")
DOWN_LABELS = ("d", "s", "b")
CKM_ROWS = ("u", "c", "t")
CKM_COLS = ("d", "s", "b")
SYSTEMS = ("epsilon_K", "K", "B_d", "B_s", "D0")
SYSTEM_COLUMN_LABELS = {
    "epsilon_K": "epsilon_K",
    "K": "DeltaM_K",
    "B_d": "DeltaM_B_d",
    "B_s": "DeltaM_B_s",
    "D0": "DeltaM_D0",
}
OPERATORS = ("C1_VLL", "C1_VRR", "C4_LR", "C5_LR")
MATRIX_INDICES = ("1", "2", "3")


def _vector_column_names(prefix: str, labels: tuple[str, ...]) -> list[str]:
    return [f"{prefix}_{label}" for label in labels]


def _complex_matrix_column_names(
    prefix: str,
    *,
    row_labels: tuple[str, ...] = MATRIX_INDICES,
    col_labels: tuple[str, ...] = MATRIX_INDICES,
    include_abs: bool = False,
) -> list[str]:
    columns: list[str] = []
    for row_label in row_labels:
        for col_label in col_labels:
            stem = f"{prefix}_{row_label}{col_label}"
            columns.extend([f"{stem}_re", f"{stem}_im"])
            if include_abs:
                columns.append(f"{stem}_abs")
    return columns


def _constraint_column_names() -> list[str]:
    columns = [
        "constraint_M_KK_GeV",
        "constraint_xi_KK",
        "constraint_matching_scale_GeV",
        "constraint_rg_low_scale_GeV",
        "constraint_evaluation",
        "constraint_operator_basis",
        "constraint_weight_policy",
        "binding_system_physical_mgkk",
        "max_ratio_to_bound_physical_mgkk",
    ]
    for system in SYSTEMS:
        label = SYSTEM_COLUMN_LABELS[system]
        columns.extend(
            [
                f"ratio_{label}_physical_mgkk",
                f"passes_{label}_physical_mgkk",
                f"bound_{label}",
                f"dominant_operator_{label}",
                f"dominant_operator_size_{label}",
            ]
        )
        columns.extend(f"weighted_size_{label}_{operator}" for operator in OPERATORS)
    return columns


COLLABORATOR_COLUMNS = [
    "point",
    "claim_level",
    "sample_type",
    "representative_sample",
    "fpr_literal_reproduction",
    "r",
    "overall_scale",
    "Lambda_IR_scan_GeV",
    "M_KK_scan_GeV",
    "xi_KK_scan",
    "m_gkk_publication_TeV",
    "xi_kk_publication",
    "k_GeV",
    "g_s_star",
    "v_GeV",
    "target_label",
    "mass_scheme",
    "mass_target_source",
    "ckm_target_source",
    "mass_target_scale_GeV",
    "wilson_reference_scale_GeV",
    "fit_success",
    "fit_score",
    "residual_norm",
    "max_abs_log_mass_residual",
    "max_abs_ckm_observable_residual",
    *_constraint_column_names(),
    "Y_basis",
    "matrix_index_order",
    "mass_matrix_formula",
    "ckm_phase_convention",
    "max_abs_Y_u_bulk",
    "max_abs_Y_d_bulk",
    "max_abs_Y_bulk",
    *_vector_column_names("m_fit_GeV", (*UP_LABELS, *DOWN_LABELS)),
    *_vector_column_names("target_m_GeV", (*UP_LABELS, *DOWN_LABELS)),
    *_vector_column_names("log_mass_residual", (*UP_LABELS, *DOWN_LABELS)),
    *_vector_column_names("ckm_observable", ("Vus", "Vcb", "Vub", "J")),
    *_vector_column_names("target_ckm_observable", ("Vus", "Vcb", "Vub", "J")),
    *_vector_column_names(
        "ckm_observable_relative_residual",
        ("Vus", "Vcb", "Vub", "J"),
    ),
    "Jarlskog_refit",
    *_complex_matrix_column_names(
        "CKM",
        row_labels=CKM_ROWS,
        col_labels=CKM_COLS,
        include_abs=True,
    ),
    *_complex_matrix_column_names(
        "target_CKM",
        row_labels=CKM_ROWS,
        col_labels=CKM_COLS,
        include_abs=True,
    ),
    *_vector_column_names("c_Q", MATRIX_INDICES),
    *_vector_column_names("c_u", MATRIX_INDICES),
    *_vector_column_names("c_d", MATRIX_INDICES),
    *_vector_column_names("F_Q", MATRIX_INDICES),
    *_vector_column_names("F_u", MATRIX_INDICES),
    *_vector_column_names("F_d", MATRIX_INDICES),
    *_complex_matrix_column_names("Y_u_bulk"),
    *_complex_matrix_column_names("Y_d_bulk"),
    *_complex_matrix_column_names("M_u_GeV"),
    *_complex_matrix_column_names("M_d_GeV"),
]


def _format_float(value: float) -> str:
    numeric = float(value)
    if not math.isfinite(numeric):
        raise ValueError(f"non-finite float {value!r}")
    return format(numeric, ".17g")


def _format_bool(value: object) -> str:
    return "true" if bool(value) else "false"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _load_scan_rows(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def _coordinate_key(r: float, overall_scale: float) -> tuple[str, str]:
    return (_format_float(r), _format_float(overall_scale))


def _select_rows(rows: Iterable[dict], coordinates: Iterable[tuple[float, float]]) -> list[dict]:
    by_coord = {
        _coordinate_key(row["r"], row["overall_scale"]): row
        for row in rows
    }
    selected: list[dict] = []
    missing: list[tuple[float, float]] = []
    for r_value, overall_scale in coordinates:
        key = _coordinate_key(r_value, overall_scale)
        row = by_coord.get(key)
        if row is None:
            missing.append((r_value, overall_scale))
        else:
            selected.append(row)
    if missing:
        raise ValueError(f"{MERGED_JSONL} is missing coordinate pairs: {missing}")
    return selected


def _add_vector(row: dict[str, str], prefix: str, labels: tuple[str, ...], values: np.ndarray) -> None:
    for label, value in zip(labels, np.asarray(values, dtype=float), strict=True):
        row[f"{prefix}_{label}"] = _format_float(value)


def _add_complex_matrix(
    row: dict[str, str],
    prefix: str,
    matrix: np.ndarray,
    *,
    row_labels: tuple[str, ...] = MATRIX_INDICES,
    col_labels: tuple[str, ...] = MATRIX_INDICES,
    include_abs: bool = False,
) -> None:
    arr = np.asarray(matrix, dtype=np.complex128)
    if arr.shape != (3, 3):
        raise ValueError(f"{prefix} has shape {arr.shape}, expected (3, 3)")
    for i, row_label in enumerate(row_labels):
        for j, col_label in enumerate(col_labels):
            stem = f"{prefix}_{row_label}{col_label}"
            row[f"{stem}_re"] = _format_float(arr[i, j].real)
            row[f"{stem}_im"] = _format_float(arr[i, j].imag)
            if include_abs:
                row[f"{stem}_abs"] = _format_float(abs(arr[i, j]))


def _physical_mgkk_phenomenology(
    result,
    *,
    point_id: str,
    point_label: str,
    m_gkk_gev: float,
    g_s_star: float,
):
    couplings = build_modern_point_couplings(
        result,
        point_id=point_id,
        point_label=point_label,
        M_KK=m_gkk_gev,
        g_s_star=g_s_star,
    )
    matching = build_modern_point_matching(couplings)
    source = SimpleNamespace(
        schema_id=MODERN_POINT_EVALUATION_SCHEMA_ID,
        point_id=point_id,
        point_label=point_label,
        input_bundle_schema_id=couplings.input_bundle_schema_id,
        input_bundle_id=couplings.input_bundle_id,
        input_provenance_id=couplings.input_provenance_id,
        input_resolution_policy_id=couplings.input_resolution_policy_id,
        qcd_metadata_id=couplings.qcd_metadata_id,
        alpha_s_policy_id=couplings.alpha_s_policy_id,
        coupling_schema_id=couplings.schema_id,
        matching_schema_id=matching.schema_id,
        couplings=couplings.as_dict(),
        matching=matching.as_dict(),
    )
    bridge = build_modern_point_bridge_artifact(source)
    return build_modern_point_phenomenology_artifact(bridge), matching


def _make_row(scan_row: dict, config: dict) -> dict[str, str]:
    targets = default_quark_targets()
    lambda_ir = float(scan_row["Lambda_IR"])
    m_gkk_gev = lambda_ir * GAUGE_KK_ROOT_NN
    ratio_rescale = 1.0 / (GAUGE_KK_ROOT_NN**2)
    raw_ratios = dict(scan_row["ratio_to_bound_by_system"])
    publication_ratios = {
        system: float(raw_ratios[system]) * ratio_rescale
        for system in SYSTEMS
    }
    binding_publication = max(publication_ratios, key=publication_ratios.get)

    solution = fit_quark_sector(
        targets,
        r=float(scan_row["r"]),
        overall_scale=float(scan_row["overall_scale"]),
        seed=default_spurion_seed(),
        k=float(scan_row["k"]),
        Lambda_IR=lambda_ir,
        max_nfev=int(config.get("max_nfev", 120)),
        fit_orientation=bool(config.get("fit_orientation", True)),
    )
    result = solution.result
    state = result.state
    g_s_star = float(config.get("g_s_star", 3.0))
    physical_phenomenology, physical_matching = _physical_mgkk_phenomenology(
        result,
        point_id=str(scan_row["point_id"]),
        point_label=str(scan_row["point_label"]),
        m_gkk_gev=m_gkk_gev,
        g_s_star=g_s_star,
    )
    physical_results = {
        system_result.system_id: system_result
        for system_result in physical_phenomenology.system_results
    }
    physical_ratios = {
        system: float(physical_results[system].ratio_to_bound)
        for system in SYSTEMS
    }
    binding_physical = max(physical_ratios, key=physical_ratios.get)
    y_u_abs = np.abs(state.Y_u_bulk_basis)
    y_d_abs = np.abs(state.Y_d_bulk_basis)

    row: dict[str, str] = {
        "point_id": str(scan_row["point_id"]),
        "point_index": str(scan_row["point_index"]),
        "point_label": str(scan_row["point_label"]),
        "source_run_dir": str(RUN_DIR),
        "source_results_jsonl": str(MERGED_JSONL),
        "source_results_sha256": _sha256(MERGED_JSONL),
        "source_config_json": str(CONFIG_JSON),
        "source_config_hash": str(scan_row["config_hash"]),
        "source_claim_level": str(scan_row["claim_level"]),
        "source_git_commit": str(scan_row["git_commit"]),
        "source_dirty_tree": _format_bool(scan_row.get("dirty_tree", False)),
        "source_input_bundle_id": str(scan_row["input_bundle_id"]),
        "source_input_provenance_id": str(scan_row["input_provenance_id"]),
        "source_lane_id": str(scan_row["lane_id"]),
        "artifact_path": str(RUN_DIR / scan_row["artifact_path"]),
        "bridge_artifact_path": str(RUN_DIR / scan_row["bridge_artifact_path"]),
        "phenomenology_artifact_path": str(RUN_DIR / scan_row["phenomenology_artifact_path"]),
        "verifier_ok": _format_bool(scan_row["verifier_ok"]),
        "bridge_verifier_ok": _format_bool(scan_row["bridge_verifier_ok"]),
        "phenomenology_verifier_ok": _format_bool(scan_row["phenomenology_verifier_ok"]),
        "accepted_raw_scan_convention": _format_bool(scan_row["accepted"]),
        "accepted_publication_convention": _format_bool(max(publication_ratios.values()) <= 1.0),
        "accepted_physical_mgkk_convention": _format_bool(physical_phenomenology.non_cp_passes),
        "r": _format_float(scan_row["r"]),
        "overall_scale": _format_float(scan_row["overall_scale"]),
        "Lambda_IR_scan_GeV": _format_float(lambda_ir),
        "M_KK_scan_GeV": _format_float(scan_row["M_KK"]),
        "xi_KK_scan": _format_float(scan_row["xi_KK"]),
        "xi_kk_publication": _format_float(GAUGE_KK_ROOT_NN),
        "m_gkk_GeV_publication": _format_float(m_gkk_gev),
        "m_gkk_TeV_publication": _format_float(m_gkk_gev / 1000.0),
        "ratio_rescale_factor_publication_over_raw": _format_float(ratio_rescale),
        "g_s_star_scan": _format_float(g_s_star),
        "k_GeV": _format_float(scan_row["k"]),
        "v_GeV": _format_float(state.point.v),
        "target_label": targets.label,
        "mass_scheme": "PDG-2024 MSbar masses evolved to mu=m_t(m_t)",
        "mass_target_source": "quarkConstraints.benchmarks.default_quark_targets",
        "ckm_target_source": "quarkConstraints.benchmarks.ckm_like_unitary PDG-like central target",
        "mass_target_scale_GeV": _format_float(DEFAULT_QUARK_FIT_SCALE_GEV),
        "wilson_reference_scale_GeV": _format_float(DEFAULT_QUARK_TARGET_SCALE_GEV),
        "fit_success_source": _format_bool(scan_row["fit_success"]),
        "fit_success_refit": _format_bool(solution.success),
        "fit_message_source": str(scan_row["fit_message"]),
        "fit_message_refit": solution.message,
        "fit_nfev_source": str(scan_row["fit_nfev"]),
        "fit_nfev_refit": str(solution.nfev),
        "fit_score_source": _format_float(scan_row["fit_score"]),
        "fit_score_refit": _format_float(result.score),
        "fit_score_abs_delta": _format_float(abs(result.score - float(scan_row["fit_score"]))),
        "residual_norm_source": _format_float(scan_row["residual_norm"]),
        "residual_norm_refit": _format_float(result.residual_norm),
        "max_abs_log_mass_residual": _format_float(
            np.max(np.abs(np.concatenate([result.mass_residuals_up, result.mass_residuals_down])))
        ),
        "max_abs_ckm_observable_residual": _format_float(np.max(np.abs(result.ckm_residuals))),
        "binding_system_publication": binding_publication,
        "max_ratio_to_bound_publication": _format_float(max(publication_ratios.values())),
        "max_ratio_to_bound_raw_scan": _format_float(scan_row["max_ratio_to_bound"]),
        "max_non_cp_ratio_to_bound_raw_scan": _format_float(scan_row["max_non_cp_ratio_to_bound"]),
        "constraint_M_KK_GeV": _format_float(physical_phenomenology.M_KK),
        "constraint_xi_KK": _format_float(physical_phenomenology.xi_KK),
        "constraint_matching_scale_GeV": _format_float(physical_matching.M_KK),
        "constraint_rg_low_scale_GeV": _format_float(2.0),
        "constraint_evaluation": "full_current_backend_at_physical_mgkk_with_qcd_rg",
        "constraint_operator_basis": physical_matching.operator_basis_id,
        "constraint_weight_policy": physical_matching.weight_policy_id,
        "binding_system_physical_mgkk": binding_physical,
        "max_ratio_to_bound_physical_mgkk": _format_float(max(physical_ratios.values())),
        "Y_basis": "bulk_mass_eigenbasis",
        "matrix_index_order": "1,2,3 light_to_heavy",
        "mass_matrix_formula": "M=2*v*diag(F_Q)*Y_bulk*diag(F_u_or_F_d)",
        "ckm_phase_convention": "raw SVD phase convention; compare abs entries and J unless rephased",
        "max_abs_Y_u_bulk": _format_float(np.max(y_u_abs)),
        "max_abs_Y_d_bulk": _format_float(np.max(y_d_abs)),
        "max_abs_Y_bulk": _format_float(max(np.max(y_u_abs), np.max(y_d_abs))),
    }

    for system in SYSTEMS:
        publication_ratio = publication_ratios[system]
        raw_ratio = float(raw_ratios[system])
        physical_result = physical_results[system]
        label = SYSTEM_COLUMN_LABELS[system]
        row[f"ratio_{system}_publication"] = _format_float(publication_ratio)
        row[f"passes_{system}_publication"] = _format_bool(publication_ratio <= 1.0)
        row[f"ratio_{system}_raw_scan"] = _format_float(raw_ratio)
        row[f"passes_{system}_raw_scan"] = _format_bool(raw_ratio <= 1.0)
        row[f"ratio_{label}_physical_mgkk"] = _format_float(physical_result.ratio_to_bound)
        row[f"passes_{label}_physical_mgkk"] = _format_bool(physical_result.passes)
        row[f"bound_{label}"] = _format_float(physical_result.bound)
        row[f"dominant_operator_{label}"] = str(physical_result.dominant_operator)
        row[f"dominant_operator_size_{label}"] = _format_float(
            physical_result.dominant_operator_size
        )
        for operator in OPERATORS:
            row[f"weighted_size_{label}_{operator}"] = _format_float(
                physical_result.weighted_operator_sizes[operator]
            )

    _add_vector(row, "m_fit_GeV", UP_LABELS, result.masses_up)
    _add_vector(row, "m_fit_GeV", DOWN_LABELS, result.masses_down)
    _add_vector(row, "target_m_GeV", UP_LABELS, targets.up_masses)
    _add_vector(row, "target_m_GeV", DOWN_LABELS, targets.down_masses)
    _add_vector(row, "log_mass_residual", UP_LABELS, result.mass_residuals_up)
    _add_vector(row, "log_mass_residual", DOWN_LABELS, result.mass_residuals_down)

    ckm_obs_labels = ("Vus", "Vcb", "Vub", "J")
    _add_vector(row, "ckm_observable", ckm_obs_labels, result.ckm_observables)
    _add_vector(row, "target_ckm_observable", ckm_obs_labels, targets.ckm_observables)
    _add_vector(row, "ckm_observable_relative_residual", ckm_obs_labels, result.ckm_residuals)
    row["Jarlskog_refit"] = _format_float(jarlskog_invariant(result.ckm))

    _add_complex_matrix(row, "CKM", result.ckm, row_labels=CKM_ROWS, col_labels=CKM_COLS, include_abs=True)
    _add_complex_matrix(row, "target_CKM", targets.ckm, row_labels=CKM_ROWS, col_labels=CKM_COLS, include_abs=True)
    _add_complex_matrix(row, "M_u_GeV", result.M_u)
    _add_complex_matrix(row, "M_d_GeV", result.M_d)
    _add_complex_matrix(row, "Y_u_bulk", state.Y_u_bulk_basis)
    _add_complex_matrix(row, "Y_d_bulk", state.Y_d_bulk_basis)

    _add_vector(row, "c_Q", MATRIX_INDICES, state.c_Q)
    _add_vector(row, "c_u", MATRIX_INDICES, state.c_u)
    _add_vector(row, "c_d", MATRIX_INDICES, state.c_d)
    _add_vector(row, "F_Q", MATRIX_INDICES, state.F_Q)
    _add_vector(row, "F_u", MATRIX_INDICES, state.F_u)
    _add_vector(row, "F_d", MATRIX_INDICES, state.F_d)
    return row


def _collaborator_row(full_row: dict[str, str], point_number: int) -> dict[str, str]:
    row = {
        "point": str(point_number),
        "claim_level": full_row["source_claim_level"],
        "sample_type": "focused_fitted_envelope_existence_point",
        "representative_sample": "false",
        "fpr_literal_reproduction": "false",
        "r": full_row["r"],
        "overall_scale": full_row["overall_scale"],
        "Lambda_IR_scan_GeV": full_row["Lambda_IR_scan_GeV"],
        "M_KK_scan_GeV": full_row["M_KK_scan_GeV"],
        "xi_KK_scan": full_row["xi_KK_scan"],
        "m_gkk_publication_TeV": full_row["m_gkk_TeV_publication"],
        "xi_kk_publication": full_row["xi_kk_publication"],
        "k_GeV": full_row["k_GeV"],
        "g_s_star": full_row["g_s_star_scan"],
        "v_GeV": full_row["v_GeV"],
        "target_label": full_row["target_label"],
        "mass_scheme": full_row["mass_scheme"],
        "mass_target_source": full_row["mass_target_source"],
        "ckm_target_source": full_row["ckm_target_source"],
        "mass_target_scale_GeV": full_row["mass_target_scale_GeV"],
        "wilson_reference_scale_GeV": full_row["wilson_reference_scale_GeV"],
        "fit_success": full_row["fit_success_refit"],
        "fit_score": full_row["fit_score_refit"],
        "residual_norm": full_row["residual_norm_refit"],
        "max_abs_log_mass_residual": full_row["max_abs_log_mass_residual"],
        "max_abs_ckm_observable_residual": full_row["max_abs_ckm_observable_residual"],
    }

    for column in COLLABORATOR_COLUMNS:
        if column in row:
            continue
        row[column] = full_row[column]
    return row


def _parse_coordinate(value: str) -> tuple[float, float]:
    try:
        r_text, scale_text = value.split(",", 1)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("coordinate must be R,OVERALL_SCALE") from exc
    return float(r_text), float(scale_text)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--coordinate",
        action="append",
        type=_parse_coordinate,
        help="Coordinate pair R,OVERALL_SCALE to export. Repeat to override defaults.",
    )
    args = parser.parse_args()

    config_payload = _load_json(CONFIG_JSON)
    config = config_payload.get("config", config_payload)
    merge_manifest = _load_json(MERGE_MANIFEST)
    verification = _load_json(VERIFICATION_JSON)
    selected_rows = _select_rows(
        _load_scan_rows(MERGED_JSONL),
        args.coordinate or DEFAULT_COORDINATES,
    )
    full_rows = [_make_row(scan_row, config) for scan_row in selected_rows]
    rows = [
        _collaborator_row(full_row, point_number=index + 1)
        for index, full_row in enumerate(full_rows)
    ]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLLABORATOR_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    provenance_path = args.output.with_suffix(".provenance.json")
    provenance = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "generator_script": "scripts/export_collaborator_5tev_points.py",
        "generator_script_sha256": _sha256(Path(__file__).resolve()),
        "output_csv": str(args.output),
        "output_csv_sha256": _sha256(args.output),
        "row_count": len(rows),
        "point_ids": [row["point_id"] for row in full_rows],
        "coordinates": [
            {"r": row["r"], "overall_scale": row["overall_scale"]}
            for row in rows
        ],
        "csv_columns": COLLABORATOR_COLUMNS,
        "source_run_dir": str(RUN_DIR),
        "source_results_jsonl": str(MERGED_JSONL),
        "source_results_sha256": _sha256(MERGED_JSONL),
        "source_config_json": str(CONFIG_JSON),
        "source_config_sha256": _sha256(CONFIG_JSON),
        "source_merge_manifest": str(MERGE_MANIFEST),
        "source_merge_manifest_sha256": _sha256(MERGE_MANIFEST),
        "source_verification_json": str(VERIFICATION_JSON),
        "source_verification_ok": bool(verification.get("ok", False)),
        "source_git_commit": full_rows[0]["source_git_commit"],
        "source_dirty_tree": full_rows[0]["source_dirty_tree"],
        "source_dirty_tree_note": (
            "The focused rerun recorded dirty_tree=true because collaborator "
            "export artifacts/scripts were untracked in the working tree; "
            "git status showed no tracked physics-source modifications."
        ),
        "source_config_hash": merge_manifest.get("config_hash"),
        "source_total_point_count": merge_manifest.get("total_point_count"),
        "source_accepted_point_count": merge_manifest.get("accepted_point_count"),
        "selection_note": (
            "Five accepted points from a current-code 15-point focused rerun on "
            "the nearest physical first-KK-gluon shell to 5 TeV "
            f"(m_gkk_TeV={rows[0]['m_gkk_publication_TeV']}). The May 6 dense "
            "accepted CSV is not used as a data source for these ratios."
        ),
        "conventions": {
            "scan_M_KK_GeV": "Lambda_IR_scan_GeV because xi_KK_scan=1 in the scan",
            "publication_m_gkk_GeV": "Lambda_IR_scan_GeV * xi_kk_publication",
            "xi_kk_publication": GAUGE_KK_ROOT_NN,
            "constraint_ratios_in_csv": (
                "full current-backend recomputation at physical m_gkk with "
                "Wilson evolution from mu_high=m_gkk to mu_low=2 GeV"
            ),
            "mass_target_scale_GeV": DEFAULT_QUARK_FIT_SCALE_GEV,
            "wilson_reference_scale_GeV": DEFAULT_QUARK_TARGET_SCALE_GEV,
            "wilson_reference_scale_note": (
                "3000 GeV is the repository input/reference convention; the "
                "exported constraint matching scale is constraint_M_KK_GeV."
            ),
            "ckm_phase_convention": (
                "Complex CKM entries are in the raw SVD phase convention; compare "
                "magnitudes and J unless a rephasing convention is imposed."
            ),
            "Y_basis": "bulk_mass_eigenbasis",
            "mass_matrix_formula": "M=2*v*diag(F_Q)*Y_bulk*diag(F_u_or_F_d)",
            "claim_level": config.get("claim_level", "operational_scan_only"),
            "warning": (
                "This is a focused fitted-envelope rerun for collaborator "
                "inspection, not an anarchic typicality claim or a full dense "
                "allowed-region replacement."
            ),
        },
    }
    provenance_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(f"Wrote {args.output} ({len(rows)} rows)")
    print(f"Wrote {provenance_path}")


if __name__ == "__main__":
    main()
