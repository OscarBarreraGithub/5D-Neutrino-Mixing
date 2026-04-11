"""Contract tests for the paper-mode artifact and verifier boundary."""

from __future__ import annotations

import ast
import hashlib
import json
import subprocess
import sys
from dataclasses import asdict, replace
from pathlib import Path

import pytest

import quarkConstraints.paper_0710_1869.verifier as verifier
from quarkConstraints.paper_0710_1869.artifacts import (
    DEFAULT_KAON_SUPPORTED_OPERATORS,
    DEFAULT_KAON_UNSUPPORTED_OPERATORS,
    HADRONIC_BUNDLE_SCHEMA,
    OBSERVABLE_BUNDLE_SCHEMA,
    PAPER_MODE,
    PROVENANCE_BUNDLE_SCHEMA,
    WILSON_BUNDLE_SCHEMA,
    ArtifactMetadata,
    ArtifactScale,
    ArtifactSchemaError,
    ArtifactSourceRecord,
    ComplexValue,
    HadronicArtifactBundleV1,
    ObservableArtifactBundleV1,
    ObservableRecord,
    ProvenanceBundleV1,
    ProvenanceRecord,
    WilsonArtifactBundleV1,
    WilsonCoefficientRecord,
    artifact_from_dict,
    build_default_paper_0710_1869_kaon_artifact_export_set,
    read_artifact,
    write_default_paper_0710_1869_kaon_artifact_exports,
)
from quarkConstraints.paper_0710_1869.verifier import (
    TolerancePolicy,
    read_verifier_inputs,
    verify_artifact_paths,
    verify_inputs,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
BENCHMARK_SCRIPT = REPO_ROOT / "scripts" / "benchmark_quark_0710_1869.py"
GOLDEN_ARTIFACT_DIR = (
    REPO_ROOT / "tests" / "golden" / "paper_0710_1869" / "default_kaon_np_only"
)
EXPECTED_ARTIFACT_FILENAMES = {
    "wilsons": "wilsons.json",
    "hadronic": "hadronic.json",
    "observables": "observables.json",
    "provenance": "provenance.json",
}


def _export_default_artifacts(tmp_path: Path) -> dict[str, object]:
    completed = subprocess.run(
        [
            sys.executable,
            str(BENCHMARK_SCRIPT),
            "--emit-json",
            "--export-artifacts-dir",
            str(tmp_path),
        ],
        check=True,
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )
    return json.loads(completed.stdout)


def _verify_exported_artifacts(tmp_path: Path, *, check: bool = True) -> dict[str, object]:
    completed = subprocess.run(
        [
            sys.executable,
            str(BENCHMARK_SCRIPT),
            "--verify-artifacts-dir",
            str(tmp_path),
        ],
        check=check,
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )
    return json.loads(completed.stdout)


def _canonical_tolerance_policy() -> dict[str, object]:
    return asdict(TolerancePolicy())


def _sample_scales() -> tuple[ArtifactScale, ...]:
    return (
        ArtifactScale(name="mu_match", role="matching scale", value_gev=3000.0),
        ArtifactScale(name="mu_had", role="hadronic evaluation scale", value_gev=2.0),
        ArtifactScale(name="m_g1", role="KK-gluon propagator mass", value_gev=3000.0),
    )


def _sample_q1_matrix_element(*, m_K0_GeV: float, f_K_GeV: float, B_K_mu_had: float) -> float:
    return (8.0 / 3.0) * (f_K_GeV**2) * (m_K0_GeV**2) * B_K_mu_had


def _sample_quartet() -> tuple[
    WilsonArtifactBundleV1,
    HadronicArtifactBundleV1,
    ObservableArtifactBundleV1,
    ProvenanceBundleV1,
]:
    benchmark_id = "benchmark-point-001"
    point_id = benchmark_id
    operator_basis = "kk_gluon_tree_np_only.v1"
    operator_normalization = "paper_0710_1869.deltaf2.kk_gluon_tree_color_normalization.v1"
    scheme = "bmu.hep-ph-0005183.ndr-ms.lo.v1"
    scales = _sample_scales()

    provenance = ProvenanceBundleV1(
        metadata=ArtifactMetadata.create(
            schema_name=PROVENANCE_BUNDLE_SCHEMA,
            bundle_id="provenance-bundle-001",
            point_id=point_id,
        ),
        records=(
            ProvenanceRecord(
                record_id="inputs.deltaf2",
                category="paper",
                label="paper conventions",
                version="arxiv:0710.1869",
                source="doi:10.48550/arXiv.0710.1869",
            ),
            ProvenanceRecord(
                record_id="hadronic.bundle.v1",
                category="hadronic-input-bundle",
                label="kaon hadronic bundle",
                version="paper-defaults:v1",
                source="paper default hadronic bundle",
            ),
            ProvenanceRecord(
                record_id="pdg.k0.mass.v1",
                category="particle-property",
                label="neutral kaon mass",
                version="PDG 2024",
                source="PDG neutral kaon listing",
            ),
            ProvenanceRecord(
                record_id="pdg.fk.v1",
                category="review-average",
                label="kaon decay constant",
                version="PDG 2024",
                source="PDG leptonic decays review",
            ),
            ProvenanceRecord(
                record_id="flag.bk.v1",
                category="derived-from-rgi",
                label="kaon bag parameter",
                version="FLAG/PDG 2024",
                source="PDG lattice review",
            ),
        ),
    )

    wilson = WilsonArtifactBundleV1(
        metadata=ArtifactMetadata.create(
            schema_name=WILSON_BUNDLE_SCHEMA,
            bundle_id="wilson-bundle-001",
            point_id=point_id,
        ),
        operator_basis=operator_basis,
        operator_normalization=operator_normalization,
        renormalization_scheme=scheme,
        scales=scales,
        coefficients=(
            WilsonCoefficientRecord(
                sector="down",
                system="kaon",
                operator="Q1_VLL",
                value=ComplexValue(real=1.25e-12, imag=3.5e-13),
            ),
            WilsonCoefficientRecord(
                sector="down",
                system="kaon",
                operator="Q1_VRR",
                value=ComplexValue(real=-4.0e-13, imag=1.5e-13),
            ),
            WilsonCoefficientRecord(
                sector="down",
                system="kaon",
                operator="Q4_LR",
                value=ComplexValue(real=0.0, imag=0.0),
            ),
            WilsonCoefficientRecord(
                sector="down",
                system="kaon",
                operator="Q5_LR",
                value=ComplexValue(real=0.0, imag=0.0),
            ),
        ),
        coefficient_scale_name="mu_had",
        matching_scale_name="mu_match",
        provenance_ids=("inputs.deltaf2",),
    )

    hadronic = HadronicArtifactBundleV1(
        metadata=ArtifactMetadata.create(
            schema_name=HADRONIC_BUNDLE_SCHEMA,
            bundle_id="hadronic-bundle-001",
            point_id=point_id,
        ),
        benchmark_id=benchmark_id,
        system="kaon",
        operator_basis=operator_basis,
        operator_normalization=operator_normalization,
        renormalization_scheme=scheme,
        scales=(ArtifactScale(name="mu_had", role="hadronic evaluation scale", value_gev=2.0),),
        matrix_element_formula_id="kaon.q1_vll_vrr.8over3_fk2_mk2_bk_mu.v1",
        hamiltonian_convention_id="heff.sum_ci_qi.no_hc_factor.v1",
        parity_relation_id="kaon.q1_vll_equals_q1_vrr.by_parity.v1",
        supported_operator_names=DEFAULT_KAON_SUPPORTED_OPERATORS,
        unsupported_operator_names=DEFAULT_KAON_UNSUPPORTED_OPERATORS,
        mu_had_GeV=2.0,
        m_K0_GeV=0.497611,
        f_K_GeV=0.1557,
        B_K_mu_had=0.7625,
        q1_matrix_element_GeV4=_sample_q1_matrix_element(
            m_K0_GeV=0.497611,
            f_K_GeV=0.1557,
            B_K_mu_had=0.7625,
        ),
        source_id="hadronic.bundle.v1",
        hadronic_source_id="hadronic.bundle.v1",
        mass_source_id="pdg.k0.mass.v1",
        decay_constant_source_id="pdg.fk.v1",
        bag_parameter_source_id="flag.bk.v1",
        source_wilson_bundle_id=wilson.metadata.bundle_id,
        mass_source=ArtifactSourceRecord(
            source_id="pdg.k0.mass.v1",
            source_kind="particle-property",
            citation="PDG neutral kaon listing",
            locator_label="K0 mass listing",
            year=2024,
        ),
        decay_constant_source=ArtifactSourceRecord(
            source_id="pdg.fk.v1",
            source_kind="review-average",
            citation="PDG leptonic decays review",
            locator_label="Eq. (72.14)",
            year=2024,
        ),
        bag_parameter_source=ArtifactSourceRecord(
            source_id="flag.bk.v1",
            source_kind="derived-from-rgi",
            citation="PDG lattice review",
            locator_label="Sec. 17.3.2",
            year=2024,
            renormalization_scheme_id="rgi.flag21.average.v1",
            transformation_id="hat_bk_to_bk_mu_had.lo_inverse_q1_running.v1",
        ),
        provenance_ids=(
            "hadronic.bundle.v1",
            "pdg.k0.mass.v1",
            "pdg.fk.v1",
            "flag.bk.v1",
        ),
    )

    m12_value = (
        (
            wilson.coefficients[0].value.to_complex()
            + wilson.coefficients[1].value.to_complex()
        )
        * hadronic.q1_matrix_element_GeV4
        / (2.0 * hadronic.m_K0_GeV)
    )
    observable = ObservableArtifactBundleV1(
        metadata=ArtifactMetadata.create(
            schema_name=OBSERVABLE_BUNDLE_SCHEMA,
            bundle_id="observable-bundle-001",
            point_id=point_id,
        ),
        source_wilson_bundle_id=wilson.metadata.bundle_id,
        source_hadronic_bundle_id=hadronic.metadata.bundle_id,
        interpretation="kaon.np_only.v1",
        operator_basis=operator_basis,
        operator_normalization=operator_normalization,
        renormalization_scheme=scheme,
        scales=scales,
        observables=(
            ObservableRecord(
                name="M12_K_NP.re",
                system="kaon",
                value=float(m12_value.real),
                units="GeV",
            ),
            ObservableRecord(
                name="M12_K_NP.im",
                system="kaon",
                value=float(m12_value.imag),
                units="GeV",
            ),
            ObservableRecord(
                name="Delta_m_K_NP",
                system="kaon",
                value=float(2.0 * m12_value.real),
                units="GeV",
            ),
        ),
        supported_operator_names=DEFAULT_KAON_SUPPORTED_OPERATORS,
        unsupported_operator_names=DEFAULT_KAON_UNSUPPORTED_OPERATORS,
        provenance_ids=(
            "inputs.deltaf2",
            "hadronic.bundle.v1",
            "pdg.k0.mass.v1",
            "pdg.fk.v1",
            "flag.bk.v1",
        ),
        observable_scope_id="kaon.np_only.m12.v1",
        m12_observable_id="M12_K_NP",
        delta_m_observable_id="Delta_m_K_NP",
        delta_m_relation_id="delta_m_k_np.equals.2_re_m12_k_np.v1",
        m12_value=ComplexValue.from_complex(m12_value),
    )

    return wilson, hadronic, observable, provenance


def _replace_coefficient_value(
    bundle: WilsonArtifactBundleV1,
    *,
    operator: str,
    value: complex,
) -> WilsonArtifactBundleV1:
    return replace(
        bundle,
        coefficients=tuple(
            replace(record, value=ComplexValue.from_complex(value))
            if record.operator == operator
            else record
            for record in bundle.coefficients
        ),
    )


def _sha256_digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.mark.parametrize(
    ("artifact_name", "filename"),
    (
        ("wilsons", "wilsons.json"),
        ("hadronic", "hadronic.json"),
        ("observables", "observables.json"),
        ("provenance", "provenance.json"),
    ),
)
def test_artifacts_round_trip_through_disk(
    tmp_path: Path,
    artifact_name: str,
    filename: str,
) -> None:
    wilson, hadronic, observable, provenance = _sample_quartet()
    artifacts = {
        "wilsons": wilson,
        "hadronic": hadronic,
        "observables": observable,
        "provenance": provenance,
    }
    artifact = artifacts[artifact_name]
    path = tmp_path / filename

    artifact.write_json(path)
    loaded = read_artifact(path)

    assert loaded == artifact
    assert path.read_text(encoding="utf-8").endswith("\n")


def test_schema_helpers_cover_round_trip_and_invalid_payload_shapes(tmp_path: Path) -> None:
    wilson, hadronic, observable, provenance = _sample_quartet()

    assert ComplexValue.from_complex(1.5 - 0.25j).to_complex() == complex(1.5, -0.25)

    source_record = ArtifactSourceRecord(
        source_id="flag.hat_bk",
        source_kind="derived-from-rgi",
        citation="FLAG hat_B_K source",
        locator_label="Sec. 17.3.2",
        year=2024,
        renormalization_scheme_id="NDR",
        scale_GeV=2.0,
        transformation_id="rgi_to_mu_had.v1",
        notes="Converted to the hadronic scheme at 2 GeV.",
    )
    assert ArtifactSourceRecord.from_dict(source_record.to_dict()) == source_record

    for artifact in (wilson, hadronic, observable, provenance):
        assert artifact_from_dict(artifact.to_dict()) == artifact

    wilson_path = wilson.write_json(tmp_path / "wilsons.json")
    hadronic_path = hadronic.write_json(tmp_path / "hadronic.json")
    observable_path = observable.write_json(tmp_path / "observables.json")
    provenance_path = provenance.write_json(tmp_path / "provenance.json")

    assert WilsonArtifactBundleV1.read_json(wilson_path) == wilson
    assert HadronicArtifactBundleV1.read_json(hadronic_path) == hadronic
    assert ObservableArtifactBundleV1.read_json(observable_path) == observable
    assert ProvenanceBundleV1.read_json(provenance_path) == provenance

    with pytest.raises(ArtifactSchemaError):
        artifact_from_dict("not-a-mapping")  # type: ignore[arg-type]
    with pytest.raises(ArtifactSchemaError):
        ArtifactMetadata.from_dict(
            {
                "schema_name": WILSON_BUNDLE_SCHEMA,
                "schema_version": "1",
                "bundle_id": "bundle",
                "point_id": "point",
                "mode": PAPER_MODE,
            }
        )
    with pytest.raises(ArtifactSchemaError):
        ArtifactScale.from_dict({"name": "mu_match", "role": "matching", "value_gev": "bad"})
    with pytest.raises(ArtifactSchemaError):
        ObservableRecord.from_dict(
            {
                "name": "Delta_m_K_NP",
                "system": "kaon",
                "value": float("nan"),
                "units": "GeV",
            }
        )


def test_hadronic_bundle_validators_reject_invalid_schema() -> None:
    wilson, hadronic, observable, provenance = _sample_quartet()

    with pytest.raises(ArtifactSchemaError):
        replace(hadronic, metadata=replace(hadronic.metadata, schema_name="wrong"))
    with pytest.raises(ArtifactSchemaError):
        replace(hadronic, supported_operator_names=())
    with pytest.raises(ArtifactSchemaError):
        HadronicArtifactBundleV1.from_dict(
            {
                **hadronic.to_dict(),
                "unsupported_operator_names": ["Q4_LR", "Q4_LR"],
            }
        )
    with pytest.raises(ArtifactSchemaError):
        replace(hadronic, mass_source="bad")  # type: ignore[arg-type]

    assert observable.source_hadronic_bundle_id == hadronic.metadata.bundle_id
    assert provenance.metadata.schema_name == PROVENANCE_BUNDLE_SCHEMA
    assert wilson.metadata.schema_name == WILSON_BUNDLE_SCHEMA


def test_read_verifier_inputs_requires_four_typed_paths(tmp_path: Path) -> None:
    wilson, hadronic, observable, provenance = _sample_quartet()
    wilson_path = wilson.write_json(tmp_path / "wilsons.json")
    hadronic_path = hadronic.write_json(tmp_path / "hadronic.json")
    observable_path = observable.write_json(tmp_path / "observables.json")
    provenance_path = provenance.write_json(tmp_path / "provenance.json")

    typed_inputs = read_verifier_inputs(
        wilson_path,
        hadronic_path,
        observable_path,
        provenance_path,
    )

    assert typed_inputs.wilson_bundle == wilson
    assert typed_inputs.hadronic_bundle == hadronic
    assert typed_inputs.observable_bundle == observable
    assert typed_inputs.provenance_bundle == provenance

    with pytest.raises(ArtifactSchemaError, match="expected HadronicArtifactBundleV1"):
        read_verifier_inputs(
            wilson_path,
            provenance_path,
            observable_path,
            provenance_path,
        )


def test_verifier_accepts_consistent_artifact_quartet(tmp_path: Path) -> None:
    wilson, hadronic, observable, provenance = _sample_quartet()
    wilson_path = wilson.write_json(tmp_path / "wilsons.json")
    hadronic_path = hadronic.write_json(tmp_path / "hadronic.json")
    observable_path = observable.write_json(tmp_path / "observables.json")
    provenance_path = provenance.write_json(tmp_path / "provenance.json")

    report = verify_artifact_paths(
        wilson_path,
        hadronic_path,
        observable_path,
        provenance_path,
    )

    assert set(report.issue_codes) <= {
        "import_isolation_failed",
        "import_isolation_runtime_violation",
    }
    assert report.point_id == wilson.metadata.point_id
    assert report.wilson_bundle_id == wilson.metadata.bundle_id
    assert report.hadronic_bundle_id == hadronic.metadata.bundle_id
    assert report.observable_bundle_id == observable.metadata.bundle_id
    assert report.provenance_bundle_id == provenance.metadata.bundle_id
    assert report.coefficient_count == len(wilson.coefficients)
    assert report.observable_count == len(observable.observables)


def test_verifier_reports_cross_bundle_contract_failures() -> None:
    wilson, hadronic, observable, provenance = _sample_quartet()
    broken_hadronic = replace(
        hadronic,
        source_wilson_bundle_id="other-wilson-bundle",
        operator_basis="JMS",
        renormalization_scheme="HV",
        provenance_ids=("hadronic.bundle.v1", "missing.record"),
    )
    object.__setattr__(
        broken_hadronic,
        "scales",
        (ArtifactScale(name="mu_had", role="hadronic", value_gev=3.0),),
    )
    object.__setattr__(broken_hadronic, "mu_had_GeV", 3.0)
    broken_observable = replace(
        observable,
        source_wilson_bundle_id="other-wilson-bundle",
        source_hadronic_bundle_id="other-hadronic-bundle",
        operator_basis="JMS",
        renormalization_scheme="HV",
        scales=observable.scales[:-1],
        provenance_ids=("inputs.deltaf2", "missing.record"),
    )

    report = verify_inputs(
        verifier.VerifierInputSet(
            wilson_bundle=wilson,
            hadronic_bundle=broken_hadronic,
            observable_bundle=broken_observable,
            provenance_bundle=provenance,
        )
    )

    issue_codes = {issue.code for issue in report.issues}
    assert not report.ok
    assert "source_wilson_bundle_id_mismatch" in issue_codes
    assert "hadronic_source_wilson_bundle_id_mismatch" in issue_codes
    assert "source_hadronic_bundle_id_mismatch" in issue_codes
    assert "operator_basis_mismatch" in issue_codes
    assert "hadronic_operator_basis_mismatch" in issue_codes
    assert "renormalization_scheme_mismatch" in issue_codes
    assert "hadronic_renormalization_scheme_mismatch" in issue_codes
    assert "scale_name_set_mismatch" in issue_codes
    assert "hadronic_scale_value_mismatch" in issue_codes
    assert "missing_provenance_record" in issue_codes


def test_verifier_enforces_q1_matrix_element_identity_and_numeric_reconstruction() -> None:
    wilson, hadronic, observable, provenance = _sample_quartet()

    broken_hadronic = replace(hadronic)
    object.__setattr__(
        broken_hadronic,
        "q1_matrix_element_GeV4",
        hadronic.q1_matrix_element_GeV4 * 1.01,
    )
    identity_report = verify_inputs(
        verifier.VerifierInputSet(
            wilson_bundle=wilson,
            hadronic_bundle=broken_hadronic,
            observable_bundle=observable,
            provenance_bundle=provenance,
        )
    )
    assert "q1_matrix_element_mismatch" in {issue.code for issue in identity_report.issues}

    broken_wilson = _replace_coefficient_value(
        wilson,
        operator="Q1_VLL",
        value=wilson.coefficients[0].value.to_complex() * 1.05,
    )
    numeric_report = verify_inputs(
        verifier.VerifierInputSet(
            wilson_bundle=broken_wilson,
            hadronic_bundle=hadronic,
            observable_bundle=observable,
            provenance_bundle=provenance,
        )
    )
    numeric_codes = {issue.code for issue in numeric_report.issues}
    assert "m12_reconstruction_mismatch" in numeric_codes
    assert "delta_m_reconstruction_mismatch" in numeric_codes


def test_verifier_rejects_nonzero_lr_coefficients() -> None:
    wilson, hadronic, observable, provenance = _sample_quartet()
    broken_wilson = _replace_coefficient_value(
        wilson,
        operator="Q4_LR",
        value=1.0e-12 + 0.0j,
    )

    report = verify_inputs(
        verifier.VerifierInputSet(
            wilson_bundle=broken_wilson,
            hadronic_bundle=hadronic,
            observable_bundle=observable,
            provenance_bundle=provenance,
        )
    )

    assert {issue.code for issue in report.issues} >= {"lr_coefficients_nonzero"}


def test_verifier_imports_only_artifact_schema_helpers() -> None:
    source = Path(verifier.__file__).read_text(encoding="utf-8")
    tree = ast.parse(source)
    allowed_relative_imports = {"artifacts", "conventions"}
    allowed_absolute_imports = {
        "quarkConstraints.paper_0710_1869.artifacts",
        "quarkConstraints.paper_0710_1869.conventions",
    }

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                assert not alias.name.startswith("quarkConstraints")
        if isinstance(node, ast.ImportFrom):
            if node.level == 1:
                assert node.module in allowed_relative_imports
            elif node.module is not None and node.module.startswith("quarkConstraints"):
                assert node.module in allowed_absolute_imports


def test_verifier_runtime_import_keeps_repo_v1_modules_unloaded() -> None:
    script = """
import importlib
import json
import sys

importlib.import_module("quarkConstraints.paper_0710_1869.verifier")

forbidden = sorted(
    name
    for name in sys.modules
    if name in {
        "quarkConstraints.deltaf2",
        "quarkConstraints.couplings",
        "quarkConstraints.model",
    }
)
print(json.dumps(forbidden))
"""
    completed = subprocess.run(
        [sys.executable, "-c", script],
        check=True,
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )
    assert json.loads(completed.stdout) == []


def test_default_kaon_artifact_writer_outputs_are_deterministic(tmp_path: Path) -> None:
    writer_a = tmp_path / "writer-a"
    writer_b = tmp_path / "writer-b"
    manual = tmp_path / "manual"

    paths_a = write_default_paper_0710_1869_kaon_artifact_exports(writer_a)
    paths_b = write_default_paper_0710_1869_kaon_artifact_exports(writer_b)

    export_set = build_default_paper_0710_1869_kaon_artifact_export_set()
    manual.mkdir(parents=True, exist_ok=True)
    export_set.wilson_bundle.write_json(manual / paths_a.wilson_path.name)
    export_set.hadronic_bundle.write_json(manual / paths_a.hadronic_path.name)
    export_set.observable_bundle.write_json(manual / paths_a.observable_path.name)
    export_set.provenance_bundle.write_json(manual / paths_a.provenance_path.name)

    assert paths_a.wilson_path.name == "wilsons.json"
    assert paths_a.hadronic_path.name == "hadronic.json"
    assert paths_a.observable_path.name == "observables.json"
    assert paths_a.provenance_path.name == "provenance.json"
    assert _sha256_digest(paths_a.wilson_path) == _sha256_digest(paths_b.wilson_path)
    assert _sha256_digest(paths_a.hadronic_path) == _sha256_digest(paths_b.hadronic_path)
    assert _sha256_digest(paths_a.observable_path) == _sha256_digest(paths_b.observable_path)
    assert _sha256_digest(paths_a.provenance_path) == _sha256_digest(paths_b.provenance_path)
    assert _sha256_digest(paths_a.wilson_path) == _sha256_digest(manual / paths_a.wilson_path.name)
    assert _sha256_digest(paths_a.hadronic_path) == _sha256_digest(
        manual / paths_a.hadronic_path.name
    )
    assert _sha256_digest(paths_a.observable_path) == _sha256_digest(
        manual / paths_a.observable_path.name
    )
    assert _sha256_digest(paths_a.provenance_path) == _sha256_digest(
        manual / paths_a.provenance_path.name
    )
    assert isinstance(read_artifact(paths_a.hadronic_path), HadronicArtifactBundleV1)


def test_benchmark_exports_default_artifact_set_and_standalone_verifier_passes(
    tmp_path: Path,
) -> None:
    summary = _export_default_artifacts(tmp_path)
    verifier_payload = _verify_exported_artifacts(tmp_path)

    assert summary["artifacts"]["status"] == "ok"
    for filename in EXPECTED_ARTIFACT_FILENAMES.values():
        assert (tmp_path / filename).exists(), filename

    checks = summary["artifacts"]["checks"]
    assert checks["observable_source_wilson_bundle_linked"] is True
    assert checks["observable_source_hadronic_bundle_linked"] is True
    assert checks["hadronic_mu_had_matches_wilson_scale"] is True
    assert checks["hadronic_operator_normalization_matches_wilson"] is True
    assert checks["observable_operator_normalization_matches_wilson"] is True
    assert checks["writer_outputs_are_deterministic"] is True

    assert verifier_payload["ok"] is True
    assert verifier_payload["import_isolation_ok"] is True
    assert verifier_payload["schema_ok"] is True
    assert verifier_payload["numeric_match_ok"] is True
    assert verifier_payload["scope_ok"] is True
    assert verifier_payload["loaded_forbidden_modules"] == []
    assert verifier_payload["tolerance_policy"] == _canonical_tolerance_policy()
    assert verifier_payload["tolerances"] == _canonical_tolerance_policy()


def test_checked_in_golden_hadronic_artifact_freezes_hat_bk_conversion_fields(
    tmp_path: Path,
) -> None:
    _export_default_artifacts(tmp_path)

    golden_hadronic_path = GOLDEN_ARTIFACT_DIR / EXPECTED_ARTIFACT_FILENAMES["hadronic"]
    golden_provenance_path = GOLDEN_ARTIFACT_DIR / EXPECTED_ARTIFACT_FILENAMES["provenance"]
    exported_hadronic_path = tmp_path / EXPECTED_ARTIFACT_FILENAMES["hadronic"]

    golden_hadronic = json.loads(golden_hadronic_path.read_text(encoding="utf-8"))
    golden_provenance = json.loads(golden_provenance_path.read_text(encoding="utf-8"))
    exported_hadronic = json.loads(exported_hadronic_path.read_text(encoding="utf-8"))

    assert golden_hadronic["hat_B_K_rgi_source_value"] == pytest.approx(0.7625)
    assert golden_hadronic["B_K_mu_had"] < golden_hadronic["hat_B_K_rgi_source_value"]
    assert golden_hadronic["bag_parameter_source_scheme_id"] == "rgi.flag21.average.v1"
    assert (
        golden_hadronic["bag_parameter_source"]["renormalization_scheme_id"]
        == golden_hadronic["bag_parameter_source_scheme_id"]
    )
    assert (
        golden_hadronic["bag_parameter_transformation_id"]
        == "hat_bk_to_bk_mu_had.lo_inverse_q1_running.v1"
    )
    assert (
        golden_hadronic["bag_parameter_source"]["transformation_id"]
        == golden_hadronic["bag_parameter_transformation_id"]
    )
    assert golden_hadronic["bag_parameter_conversion_alpha_s_mu_had"] > 0.0
    assert golden_hadronic["bag_parameter_conversion_beta0"] == pytest.approx(25.0 / 3.0)
    assert golden_hadronic["bag_parameter_conversion_exponent"] == pytest.approx(0.24)
    assert (
        golden_hadronic["bag_parameter_conversion_formula_id"]
        == "b_k_mu_had.equals.hat_b_k_rgi_times_alpha_s_mu_had_pow_2_over_beta0_lo.v1"
    )
    assert golden_hadronic["bag_parameter_conversion_n_f"] == 4
    assert golden_hadronic["alpha_s_policy_id"] == "qcd.alpha_s.lo.beta0.thresholds_continuous.v1"
    assert golden_hadronic["bag_parameter_source_id"] in {
        record["record_id"] for record in golden_provenance["records"]
    }

    for field_name in (
        "B_K_mu_had",
        "hat_B_K_rgi_source_value",
        "bag_parameter_source_id",
        "bag_parameter_source_scheme_id",
        "bag_parameter_transformation_id",
        "bag_parameter_conversion_alpha_s_mu_had",
        "bag_parameter_conversion_beta0",
        "bag_parameter_conversion_exponent",
        "bag_parameter_conversion_formula_id",
        "bag_parameter_conversion_n_f",
        "alpha_s_policy_id",
    ):
        assert exported_hadronic[field_name] == golden_hadronic[field_name]
    assert exported_hadronic["bag_parameter_source"] == golden_hadronic["bag_parameter_source"]


def test_standalone_verifier_rejects_structural_tampering(tmp_path: Path) -> None:
    _export_default_artifacts(tmp_path)
    observable_path = tmp_path / EXPECTED_ARTIFACT_FILENAMES["observables"]
    payload = json.loads(observable_path.read_text(encoding="utf-8"))
    payload["source_hadronic_bundle_id"] = "tampered-hadronic-bundle"
    observable_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    verifier_payload = _verify_exported_artifacts(tmp_path, check=False)

    assert verifier_payload["ok"] is False
    assert verifier_payload["schema_ok"] is False
    assert "source_hadronic_bundle_id_mismatch" in verifier_payload["issue_codes"]


def test_standalone_verifier_rejects_numeric_tampering(tmp_path: Path) -> None:
    _export_default_artifacts(tmp_path)
    wilson_path = tmp_path / EXPECTED_ARTIFACT_FILENAMES["wilsons"]
    payload = json.loads(wilson_path.read_text(encoding="utf-8"))
    for record in payload["coefficients"]:
        if record["operator"] == "Q1_VLL":
            record["value"]["real"] *= 1.05
            break
    wilson_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    verifier_payload = _verify_exported_artifacts(tmp_path, check=False)

    assert verifier_payload["ok"] is False
    assert verifier_payload["schema_ok"] is True
    assert verifier_payload["numeric_match_ok"] is False
    assert "m12_reconstruction_mismatch" in verifier_payload["issue_codes"]


def test_standalone_verifier_rejects_scope_creep_from_epsilon_k_and_lr_support(
    tmp_path: Path,
) -> None:
    _export_default_artifacts(tmp_path)

    observable_path = tmp_path / EXPECTED_ARTIFACT_FILENAMES["observables"]
    observable_payload = json.loads(observable_path.read_text(encoding="utf-8"))
    observable_payload["observables"].append(
        {
            "name": "epsilon_K_NP",
            "system": "kaon",
            "units": "dimensionless",
            "value": 1.0e-6,
        }
    )
    observable_path.write_text(
        json.dumps(observable_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    wilson_path = tmp_path / EXPECTED_ARTIFACT_FILENAMES["wilsons"]
    wilson_payload = json.loads(wilson_path.read_text(encoding="utf-8"))
    wilson_payload.pop("supported_operator_names", None)
    wilson_payload.pop("unsupported_operator_names", None)
    wilson_payload["coefficients"].append(
        {
            "operator": "Q4_LR",
            "sector": "down",
            "system": "kaon",
            "value": {"real": 1.0e-9, "imag": 0.0},
        }
    )
    wilson_path.write_text(
        json.dumps(wilson_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    verifier_payload = _verify_exported_artifacts(tmp_path, check=False)

    assert verifier_payload["ok"] is False
    assert verifier_payload["schema_ok"] is True
    assert verifier_payload["scope_ok"] is False
    assert "lr_coefficients_nonzero" in verifier_payload["issue_codes"]
    assert "unexpected_observable_row" in verifier_payload["issue_codes"]
