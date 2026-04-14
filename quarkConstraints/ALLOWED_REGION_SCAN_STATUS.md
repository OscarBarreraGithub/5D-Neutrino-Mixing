# Allowed Region Scan Status

Updated: 2026-04-14

## Goal

Build the missing end-to-end machinery for an honest quark-sector allowed-parameter-space scan of the Yukawas in the `modern` lane, with a production-quality, shardable, resumable SLURM pipeline from quark point evaluation through artifact export, verification, merge/reduction, and runbooked pilot/production execution.

## Current Repo Truth

- `quarkConstraints/paper_0710_1869/` is the frozen paper-facing audit lane. Its published claim boundary is complete on that boundary and must not be silently widened.
- The broader program contract in [`QUARK_COMPLETION_PROGRAM.md`](./QUARK_COMPLETION_PROGRAM.md) says QS1 exists in the working tree, but QS2 through QS7 still block any honest full allowed-quark-Yukawa claim.
- The current `modern` lane already contains lane conventions/export surface, explicit versioned modern default input payloads and provenance, schemaed pointwise coupling/matching bridge objects consumed by modern evaluation, deterministic per-point evaluation contracts, per-point artifact export, artifact-only verification, and an operational scan engine with deterministic shards/manifests, warm-start seed caching, merge/reduce, and quark-specific SLURM wrappers.
- The current `modern` lane now contains explicit bridge and phenomenology sidecars with artifact-only verifier/export paths, but it still does not contain a fully modern `Delta F = 2` backend or a distinct non-CP kaon observable. The repo is therefore still below a final allowed-region claim.
- The repo-level [`scan.py`](./scan.py) is exploratory `repo_v1` machinery and is not sufficient by itself for the final allowed-region claim.
- Existing scripts under `scripts/` are useful operational templates, but the current `run_scan_sharded.sbatch`, `reclassify_scan_shards.sbatch`, and `submit_scan_pipeline.sh` are not yet the final modern quark allowed-region pipeline.
- The working tree is dirty and already contains unrelated user changes, especially in the paper lane. Those changes must not be reverted.

## Agreed Milestone Order

Reconciled from the required author/reviewer planning loop on 2026-04-13:

1. Keep lane separation strict and preserve the frozen `paper_0710_1869` boundary exactly as stated in [`paper_0710_1869/CURRENT_STATUS.md`](./paper_0710_1869/CURRENT_STATUS.md).
2. Close the practical QS4 layer first in the claim-bearing `modern` lane: real modern input bundles with numeric payloads, provenance ids, and interpretation policy ids.
3. Close the practical QS2 bridge in `modern`: fitted quark point -> KK-gluon couplings -> EFT Wilsons/matching metadata with explicit scale and bundle tags.
4. Close the practical QS5 phenomenology layer in `modern`: explicit `K`, `B_d`, `B_s`, conservative `D0`, and honest `epsilon_K` handling whenever kaon viability is claimed.
5. Close the practical QS7 production scan layer in `modern`: deterministic point ids/configs, shard partitioning, resumability, warm-start/cache plumbing, per-point artifacts, manifests, merge/reduce, and verifier-backed acceptance.
6. Only after the modern scan engine exists, replace the current exploratory/lepton-flavored SLURM scripts with a quark-specific shard -> merge -> verify pipeline and run local smoke plus pilot validation.

## Completed Work

- Read the mandated starting files listed in the takeover request.
- Created this durable handoff/status file for the allowed-region scan program.
- Captured the dirty working-tree warning so later work does not accidentally trample user-owned edits.
- Completed the required two-subagent planning loop:
  - author draft: practical implementation slices centered on modern inputs, coupling/matching, phenomenology, and the production scan engine
  - independent review: flagged the current exploratory `quarkConstraints/scan.py` and existing `scripts/*.sbatch` pipeline as non-claim-bearing and unsuitable for honest modern allowed-region claims
- Reconciled one agreed milestone order: modern numeric input bundles and modern phenomenology must land before scan scale-out and before any allowed-region claim.
- Landed the first major `modern` scan-engine slice in [`modern/scan.py`](./modern/scan.py):
  - deterministic canonical point enumeration over `r`, `overall_scale`, and `Lambda_IR`
  - stable `config_hash` and per-point `point_id`
  - deterministic shard partitioning
  - resumable shard outputs with per-point result rows
  - per-point artifact export through the existing modern artifact contract
  - subprocess verifier execution so the artifact-only verifier boundary stays intact
  - shard manifests, merged reduction manifests, and merged verification payloads
  - CLI entrypoint via `python -m quarkConstraints.modern.scan`
- Preserved and exposed the scan surface through the lazy [`modern/__init__.py`](./modern/__init__.py) export layer without forcing `modern.scan` to import during unrelated modern-lane imports.
- Added focused scan tests in [`tests/test_modern_scan.py`](../tests/test_modern_scan.py) for:
  - deterministic point/shard enumeration
  - resume idempotence
  - shard-plus-merge smoke flow with real artifact writing and verifier subprocesses
- Added quark-specific SLURM wrappers around the landed modern scan CLI without touching the old exploratory/lepton-flavored scripts:
  - [`scripts/run_modern_quark_scan_shard.sbatch`](../scripts/run_modern_quark_scan_shard.sbatch)
  - [`scripts/merge_modern_quark_scan.sbatch`](../scripts/merge_modern_quark_scan.sbatch)
  - [`scripts/verify_modern_quark_scan.sbatch`](../scripts/verify_modern_quark_scan.sbatch)
  - [`scripts/submit_modern_quark_scan_pipeline.sh`](../scripts/submit_modern_quark_scan_pipeline.sh)
- Added wrapper-level config materialization for `smoke` and `pilot` presets, plus env-driven config freezing via `MODERN_SCAN_CONFIG_PATH`.
- Validated the new wrapper path locally outside SLURM:
  - one-point smoke shard run completed under `/tmp/modern_quark_scan_smoke2/`
  - merge completed and wrote `/tmp/modern_quark_scan_smoke2/merged/manifest.json`
  - standalone verify completed and wrote `/tmp/modern_quark_scan_smoke2/merged/verification.json` with `ok = true`
- Landed the practical QS4 slice in [`modern/inputs.py`](./modern/inputs.py) and [`modern/evaluation.py`](./modern/evaluation.py):
  - replaced the registry-only `ModernDefaultInputs` placeholder with explicit versioned modern payload blocks for neutral-meson targets, operator-weight policy, CKM targets, quark-mass targets, QCD metadata, and provenance records
  - kept `strict_paper_inputs` frozen and lane-separated while preserving `paper_inputs = ()` for the modern default family
  - rewired `evaluate_modern_point()` to pass explicit modern-derived `DeltaF2Input` payloads into the backend instead of silently falling back to the repo-owned default bundle
  - stamped modern per-point evaluation records with `input_bundle_schema_id`, `input_bundle_id`, `input_provenance_id`, and `input_resolution_policy_id`
- Added QS4-focused tests in [`tests/test_modern_input_registry.py`](../tests/test_modern_input_registry.py) and [`tests/test_modern_point_evaluation.py`](../tests/test_modern_point_evaluation.py):
  - deterministic numeric bundle serialization and mutation rejection
  - explicit backend-input wiring from the modern bundle
  - `xi_KK` consistency against the bundle QCD metadata
- Fixed a runbook regression in [`modern/scan.py`](./modern/scan.py) and [`tests/test_modern_scan.py`](../tests/test_modern_scan.py):
  - `write-preset` configs written directly to `RUN_DIR/config.json` are now upgraded in place to run-config snapshots during `run-shard`
  - added a regression test for `write-preset -> run-shard` with `config.json` at the run root
- Landed the first practical QS2 subset in [`modern/couplings.py`](./modern/couplings.py), [`modern/matching.py`](./modern/matching.py), and [`modern/evaluation.py`](./modern/evaluation.py):
  - added schemaed modern pointwise coupling objects with bundle/provenance, QCD-metadata, target-id, and mass-basis matrix payloads
  - added schemaed modern tree-level matching objects with explicit `system_id`, `observable_id`, `backend_key`, `backend_system_id`, `generations`, raw left/right couplings, and raw `C1_VLL`, `C1_VRR`, `C4_LR`, `C5_LR`
  - rewired modern evaluation to consume the explicit bridge objects and to derive the current exclusion surrogate from them instead of hiding the coupling/matching hop in one backend call
- Added focused QS2 tests:
  - [`tests/test_modern_couplings.py`](../tests/test_modern_couplings.py) validates exact agreement with the current repo coupling formula at a benchmark point
  - [`tests/test_modern_matching.py`](../tests/test_modern_matching.py) validates exact agreement with the current repo Wilson formulas at a benchmark point
  - [`tests/test_modern_point_evaluation.py`](../tests/test_modern_point_evaluation.py) now asserts the bridge objects are present and remain numerically aligned with the current benchmark backend
- Closed the practical QS2 export slice end to end in the `modern` lane:
  - introduced a dedicated artifact-only bridge sidecar module in [`modern/bridge_artifacts.py`](./modern/bridge_artifacts.py) and shifted the scan/verifier path onto it
  - restored the live point-verifier import-isolation boundary by keeping the point-artifact import path backend-free
  - extended [`modern/verifier.py`](./modern/verifier.py) with explicit bridge-sidecar verification while keeping the verifier on artifact-only imports
  - wired [`modern/scan.py`](./modern/scan.py) to write one point artifact plus one bridge sidecar per point, to run both verifier subprocesses, and to record both artifact paths and verifier statuses in scan rows/manifests/merged verification
  - tightened scan resumability so existing rows cannot be skipped unless both the point artifact and the bridge sidecar still exist on disk
  - propagated canonical scan `point_id` and `point_label` into the modern evaluation path so bridge sidecars produced during scans match the deterministic shard point identity
- Added focused bridge-export tests:
  - [`tests/test_modern_bridge_artifacts.py`](../tests/test_modern_bridge_artifacts.py) validates deterministic bridge export, round-trip stability, bridge-verifier acceptance, tampering rejection, and artifact-only import isolation
  - [`tests/test_modern_scan.py`](../tests/test_modern_scan.py) now asserts bridge artifact paths/verifier bits are present and that resume fails if the bridge sidecar is missing
- Landed the first practical QS5 slice end to end in the `modern` lane for the honest first non-CP release scope:
  - [`modern/phenomenology.py`](./modern/phenomenology.py) now exports an explicit phenomenology sidecar over the bridge artifact
  - the sidecar keeps `epsilon_K` diagnostic-only, keeps generic `K` explicitly blocked, and drives release-scope acceptance from `B_d`, `B_s`, and conservative `D0`
  - [`modern/verifier.py`](./modern/verifier.py) now independently verifies the phenomenology sidecar under the same artifact-only import boundary
  - [`modern/scan.py`](./modern/scan.py) now writes, resumes, merges, and verifies point, bridge, and phenomenology artifacts together and records release-scope acceptance semantics in scan rows
  - [`tests/test_modern_phenomenology.py`](../tests/test_modern_phenomenology.py) and [`tests/test_modern_scan.py`](../tests/test_modern_scan.py) now cover the additive QS5 path
- Fixed the broken bridge-sidecar artifact path in [`modern/artifacts.py`](./modern/artifacts.py) by restoring the missing canonical complex-matrix helper used by the artifact-side bridge records; without that repair, the landed bridge and phenomenology scan path did not actually run.
- Tightened the scan-row contract after the first QS5 slice so the row schema now mirrors the phenomenology sidecar explicitly:
  - `accepted` and `phenomenology_passes` now remain tied to the non-CP release scope only
  - rows now carry `phenomenology_release_scope_id`, explicit non-CP acceptance-system ids, explicit diagnostic-only-system ids, explicit blocked-system ids, and separate ratio/failure payloads for acceptance-bearing vs diagnostic systems
  - legacy `ratio_to_bound_by_system` and `failing_system_ids` are retained only as diagnostic mirrors of the evaluated systems, not as the claim-bearing acceptance contract
  - [`tests/test_modern_scan.py`](../tests/test_modern_scan.py) now includes regressions where a diagnostic-only `epsilon_K` failure leaves `accepted = true` and where a non-CP `B_d` failure flips `accepted = false`

## Work In Progress

- Scan-engine slice is landed and locally validated.
- The SLURM-wrapper slice is also landed and locally smoke-tested as operational infrastructure over the current scan CLI.
- The practical QS4 slice is landed and locally validated.
- The practical QS2 bridge is now landed and locally validated end to end, including sidecar export, artifact-only verification, scan row/manifest integration, and resume enforcement on both artifacts.
- The additive QS5 phenomenology slice is now also landed and locally validated end to end for the honest first non-CP release scope:
  - keep the frozen point artifact unchanged and treat it as legacy diagnostic export only
  - drive release-scope acceptance from `B_d`, `B_s`, and conservative `D0`
  - export `epsilon_K` as explicit diagnostic-only metadata
  - keep generic `K` blocked until a distinct non-CP kaon observable exists in the modern bundle/backend
- Claim-bearing work still remaining:
  - distinct non-CP kaon observable/backend so generic `K` can move from blocked to evaluated
  - broader backend/provenance work before any final modern allowed-region claim
- Honesty note remains unchanged: this operational machinery does not by itself upgrade the repo to a final allowed-region claim while those later milestones remain open.

## Next Concrete Steps

1. Add a distinct non-CP kaon observable/backend so the current blocked `K` row can become evaluated instead of explicitly unavailable.
2. Run a small multi-shard SLURM pilot with the full point + bridge + phenomenology artifact path and record the exact output root, manifests, merged verification payload, and scheduler job ids here.
3. Reassess whether the current transitional repo-v1 coupling/matching backend is an acceptable wrapper boundary for a first non-CP release or whether more backend provenance must be surfaced before any claim upgrade.
4. If the pilot stays clean, scale the current first-release non-CP scan scope out across SLURM while keeping generic kaon claims blocked in docs and outputs.

## Open Blockers/Risks

- The current `modern` package is policy- and contract-heavy; the claim-bearing numerical backend still depends on top-level exploratory modules and needs careful adaptation instead of naive reuse.
- The working tree already has substantial paper-lane modifications and untracked modern-lane files; integrations must avoid mutating the frozen paper boundary or overwriting user work.
- The practical QS4 bundle is now explicit and versioned, the practical QS2 bridge objects now exist, and the additive QS5 phenomenology sidecar is now wired through scan/export/verifier, but the modern lane still projects through the compact top-level `deltaf2.py` semantics for the evaluated surrogate rows. That remains acceptable only as a transitional backend while the broader backend milestones remain open.
- The current modern policy surface explicitly lists both `epsilon_K` and `K`, and the landed phenomenology sidecar now makes that difference explicit: `epsilon_K` is diagnostic-only and generic `K` is blocked. A distinct modern non-CP kaon observable is still missing, so the broader kaon-capable QS5 scope is not yet closed.
- The first honest claim-bearing scope under the current numeric bundle is narrower than a generic kaon-capable `Delta F = 2` release: `B_d`, `B_s`, and conservative `D0` can drive acceptance, `epsilon_K` can only be diagnostic, and generic `K` must stay blocked until a distinct modern non-CP kaon observable lands.
- Pre-QS5 shard outputs are not resume-compatible with the current scan-row schema, because `phenomenology_artifact_path` is now required for deterministic restart/integrity checks.
- The current `scripts/run_scan_sharded.sbatch`, `scripts/reclassify_scan_shards.sbatch`, and `scripts/submit_scan_pipeline.sh` reference lepton/neutrino scan axes and a reclassification flow, not modern quark artifact verification.
- The scan-engine slice can be operationally complete while still being scientifically non-claim-bearing, because it still sits on a transitional repo-v1 coupling/matching backend until the later QS2/QS5 slices land.
- The landed scan engine is deterministic and resumable, but the smoke validation runs still returned `accepted_point_count = 0` with `max_nfev = 1`; that was an operational smoke choice, not a physics statement about viability.
- The merged-run verifier now exists in `python -m quarkConstraints.modern.scan verify`, but it only checks scan/output integrity and recorded per-point verifier outcomes. It is not a substitute for the still-open claim-bearing physics milestones.

## Validation Status

- Modern-lane validation suite passed after the bridge-sidecar export slice:
  - `pytest -q tests/test_modern_scan.py tests/test_modern_point_evaluation.py tests/test_modern_point_artifacts.py tests/test_modern_bridge_artifacts.py tests/test_modern_input_registry.py tests/test_modern_phenomenology.py tests/test_modern_couplings.py tests/test_modern_matching.py`
  - result after the bridge-sidecar export slice: `72 passed in 12.43s`
- Modern-lane validation suite passed again after the QS5 scan/verifier integration repair:
  - `pytest -q tests/test_modern_scan.py tests/test_modern_point_evaluation.py tests/test_modern_point_artifacts.py tests/test_modern_bridge_artifacts.py tests/test_modern_input_registry.py tests/test_modern_phenomenology.py tests/test_modern_couplings.py tests/test_modern_matching.py`
  - current result after the scan-row/manifests honesty cleanup: `77 passed in 16.01s`
- Focused QS5 validation passed after the bridge-artifact helper repair:
  - `pytest -q tests/test_modern_phenomenology.py tests/test_modern_scan.py`
  - current result after the scan-row semantics cleanup: `35 passed in 14.22s`
- Import/compile checks passed:
  - `python -m py_compile quarkConstraints/modern/__init__.py quarkConstraints/modern/scan.py tests/test_modern_scan.py`
  - `python -m py_compile quarkConstraints/modern/inputs.py quarkConstraints/modern/evaluation.py tests/test_modern_input_registry.py tests/test_modern_point_evaluation.py`
- Import/compile checks passed for the bridge-sidecar export slice:
  - `python -m py_compile quarkConstraints/modern/artifacts.py quarkConstraints/modern/bridge_artifacts.py quarkConstraints/modern/verifier.py quarkConstraints/modern/scan.py quarkConstraints/modern/__init__.py tests/test_modern_point_artifacts.py tests/test_modern_bridge_artifacts.py tests/test_modern_scan.py`
- Import/compile checks passed after the QS5 scan/verifier integration repair:
  - `python -m py_compile quarkConstraints/modern/artifacts.py quarkConstraints/modern/phenomenology.py quarkConstraints/modern/verifier.py quarkConstraints/modern/scan.py tests/test_modern_scan.py tests/test_modern_phenomenology.py`
- Direct CLI smoke path passed end to end on a fresh temporary run directory:
  - `python -m quarkConstraints.modern.scan write-preset smoke /tmp/modern-scan-smoke2/config.json`
  - `python -m quarkConstraints.modern.scan run-shard --config /tmp/modern-scan-smoke2/config.json --output-dir /tmp/modern-scan-smoke2 --shard-index 0 --shard-count 1`
  - `python -m quarkConstraints.modern.scan merge --run-dir /tmp/modern-scan-smoke2`
  - `python -m quarkConstraints.modern.scan verify --run-dir /tmp/modern-scan-smoke2 --config /tmp/modern-scan-smoke2/config.json --total-shards 1`
  - produced complete shard, merged, and verification manifests under `/tmp/modern-scan-smoke2/`
- Direct CLI smoke path was revalidated after the QS4 slice and the config-snapshot compatibility fix:
  - `python -m quarkConstraints.modern.scan write-preset smoke /tmp/modern-scan-smoke3/config.json`
  - `python -m quarkConstraints.modern.scan run-shard --config /tmp/modern-scan-smoke3/config.json --output-dir /tmp/modern-scan-smoke3 --shard-index 0 --shard-count 1`
  - `python -m quarkConstraints.modern.scan merge --run-dir /tmp/modern-scan-smoke3`
  - `python -m quarkConstraints.modern.scan verify --run-dir /tmp/modern-scan-smoke3 --config /tmp/modern-scan-smoke3/config.json --total-shards 1`
  - produced complete shard, merged, and verification manifests under `/tmp/modern-scan-smoke3/`
- Direct CLI smoke path was revalidated again after the explicit QS2 bridge subset:
  - `python -m quarkConstraints.modern.scan write-preset smoke /tmp/modern-scan-smoke4/config.json`
  - `python -m quarkConstraints.modern.scan run-shard --config /tmp/modern-scan-smoke4/config.json --output-dir /tmp/modern-scan-smoke4 --shard-index 0 --shard-count 1`
  - `python -m quarkConstraints.modern.scan merge --run-dir /tmp/modern-scan-smoke4`
  - `python -m quarkConstraints.modern.scan verify --run-dir /tmp/modern-scan-smoke4 --config /tmp/modern-scan-smoke4/config.json --total-shards 1`
  - produced complete shard, merged, and verification manifests under `/tmp/modern-scan-smoke4/`
- Direct CLI smoke path was revalidated again after the QS5 scan/verifier integration repair:
  - `python -m quarkConstraints.modern.scan write-preset smoke /tmp/modern-scan-smoke7/config.json`
  - `python -m quarkConstraints.modern.scan run-shard --config /tmp/modern-scan-smoke7/config.json --output-dir /tmp/modern-scan-smoke7 --shard-index 0 --shard-count 1`
  - `python -m quarkConstraints.modern.scan merge --run-dir /tmp/modern-scan-smoke7`
  - `python -m quarkConstraints.modern.scan verify --run-dir /tmp/modern-scan-smoke7 --config /tmp/modern-scan-smoke7/config.json --total-shards 1`
  - produced complete shard, merged, and verification manifests under `/tmp/modern-scan-smoke7/`
  - the result row now records `phenomenology_release_scope_id = quarkConstraints.modern.phenomenology.release.non_cp_deltaf2.v1`, `non_cp_acceptance_system_ids = [B_d, B_s, D0]`, `diagnostic_only_system_ids = [epsilon_K]`, and `blocked_system_ids = [K]`
  - merged verification returned `ok = true`, `bridge_verifier_failed_point_count = 0`, and `phenomenology_verifier_failed_point_count = 0`
- Initial inherited truth: the frozen paper-facing lane is complete on its published boundary, while the broader allowed-region program remains blocked by QS2-QS7.
- Planning-loop validation:
  - author/reviewer agreement: the exploratory `quarkConstraints/scan.py` path is not the production engine
  - author/reviewer agreement: modern bundles and modern phenomenology must precede production scan claims
- Current scan-engine honesty status:
  - operational scan machinery: implemented
  - practical QS4 modern input bundle: implemented
  - practical QS2 bridge export/verifier/scan integration: implemented
  - practical first-release QS5 phenomenology/export/verifier/scan semantics: implemented
  - broader kaon-capable modern allowed-region claim: still blocked by the distinct non-CP `K` backend and later milestones
- SLURM-wrapper validation passed on 2026-04-13:
  - `bash -n scripts/run_modern_quark_scan_shard.sbatch`
  - `bash -n scripts/merge_modern_quark_scan.sbatch`
  - `bash -n scripts/verify_modern_quark_scan.sbatch`
  - `bash -n scripts/submit_modern_quark_scan_pipeline.sh`
  - `MODERN_SCAN_PRESET=smoke MODERN_SCAN_OUTPUT_ROOT=/tmp/modern_quark_scan_smoke2 MODERN_SCAN_TOTAL_SHARDS=1 bash scripts/run_modern_quark_scan_shard.sbatch`
  - `MODERN_SCAN_OUTPUT_ROOT=/tmp/modern_quark_scan_smoke2 bash scripts/merge_modern_quark_scan.sbatch`
  - `MODERN_SCAN_OUTPUT_ROOT=/tmp/modern_quark_scan_smoke2 MODERN_SCAN_TOTAL_SHARDS=1 bash scripts/verify_modern_quark_scan.sbatch`
  - `MODERN_SCAN_SUBMIT_DRY_RUN=1 MODERN_SCAN_PRESET=pilot MODERN_SCAN_TOTAL_SHARDS=4 MODERN_SCAN_OUTPUT_ROOT=/tmp/modern_quark_scan_pilot2 bash scripts/submit_modern_quark_scan_pipeline.sh`
- Wrapper validation outcome:
  - direct smoke shard -> merge -> verify path completed successfully
  - verification report returned `ok: true`
  - smoke run remained operational only and accepted `0/1` points, which is fine for a wrapper sanity test

## Runbook / Commands

Initial commands executed for takeover:

```bash
sed -n '1,260p' quarkConstraints/QUARK_COMPLETION_PROGRAM.md
sed -n '1,260p' quarkConstraints/paper_0710_1869/CURRENT_STATUS.md
sed -n '1,260p' quarkConstraints/README.md
sed -n '1,260p' quarkConstraints/scan.py
sed -n '1,260p' quarkConstraints/modern/__init__.py
sed -n '1,260p' quarkConstraints/modern/inputs.py
sed -n '1,260p' quarkConstraints/modern/phenomenology.py
sed -n '1,260p' quarkConstraints/modern/evaluation.py
sed -n '1,260p' quarkConstraints/modern/artifacts.py
sed -n '1,260p' quarkConstraints/modern/verifier.py
sed -n '1,260p' scripts/run_scan_sharded.sbatch
sed -n '1,260p' scripts/reclassify_scan_shards.sbatch
sed -n '1,260p' scripts/submit_scan_pipeline.sh
git status --short --untracked-files=all
sed -n '261,520p' quarkConstraints/scan.py
sed -n '261,520p' quarkConstraints/modern/evaluation.py
sed -n '1,240p' quarkConstraints/deltaf2.py
sed -n '241,420p' quarkConstraints/deltaf2.py
sed -n '1,220p' quarkConstraints/couplings.py
sed -n '1,260p' quarkConstraints/fit.py
sed -n '1,240p' quarkConstraints/paper_0710_1869/scan.py
rg -n "def compute_quark_kk_gluon_couplings|class QuarkMassBasisCouplings" quarkConstraints/couplings.py
rg -n "def evaluate_delta_f2_constraints|class DeltaF2ConstraintSummary|class DeltaF2ObservableSummary" quarkConstraints/deltaf2.py
rg -n "class QuarkFitResult|class QuarkFitSolution|def fit_quark_sector" quarkConstraints/fit.py
rg -n "modern|scan|artifact|verifier|manifest" tests -g 'test_modern*' -g 'test_quark_fit.py'
```

Validation commands executed after landing the scan and wrapper slices:

```bash
pytest -q tests/test_modern_scan.py \
  tests/test_modern_point_evaluation.py \
  tests/test_modern_point_artifacts.py \
  tests/test_modern_input_registry.py \
  tests/test_modern_phenomenology.py
python -m py_compile \
  quarkConstraints/modern/__init__.py \
  quarkConstraints/modern/scan.py \
  tests/test_modern_scan.py
python -m quarkConstraints.modern.scan write-preset smoke /tmp/modern-scan-smoke2/config.json
python -m quarkConstraints.modern.scan run-shard \
  --config /tmp/modern-scan-smoke2/config.json \
  --output-dir /tmp/modern-scan-smoke2 \
  --shard-index 0 \
  --shard-count 1
python -m quarkConstraints.modern.scan merge --run-dir /tmp/modern-scan-smoke2
python -m quarkConstraints.modern.scan verify \
  --run-dir /tmp/modern-scan-smoke2 \
  --config /tmp/modern-scan-smoke2/config.json \
  --total-shards 1
bash -n scripts/run_modern_quark_scan_shard.sbatch
bash -n scripts/merge_modern_quark_scan.sbatch
bash -n scripts/verify_modern_quark_scan.sbatch
bash -n scripts/submit_modern_quark_scan_pipeline.sh
MODERN_SCAN_PRESET=smoke \
MODERN_SCAN_OUTPUT_ROOT=/tmp/modern_quark_scan_smoke2 \
MODERN_SCAN_TOTAL_SHARDS=1 \
bash scripts/run_modern_quark_scan_shard.sbatch
MODERN_SCAN_OUTPUT_ROOT=/tmp/modern_quark_scan_smoke2 \
bash scripts/merge_modern_quark_scan.sbatch
MODERN_SCAN_OUTPUT_ROOT=/tmp/modern_quark_scan_smoke2 \
MODERN_SCAN_TOTAL_SHARDS=1 \
bash scripts/verify_modern_quark_scan.sbatch
MODERN_SCAN_SUBMIT_DRY_RUN=1 \
MODERN_SCAN_PRESET=pilot \
MODERN_SCAN_TOTAL_SHARDS=4 \
MODERN_SCAN_OUTPUT_ROOT=/tmp/modern_quark_scan_pilot2 \
bash scripts/submit_modern_quark_scan_pipeline.sh
```

Validation commands executed after landing the QS4 slice and the config-snapshot fix:

```bash
python -m py_compile \
  quarkConstraints/modern/inputs.py \
  quarkConstraints/modern/evaluation.py \
  tests/test_modern_input_registry.py \
  tests/test_modern_point_evaluation.py
pytest -q tests/test_modern_input_registry.py \
  tests/test_modern_point_evaluation.py \
  tests/test_modern_point_artifacts.py
pytest -q tests/test_modern_scan.py
pytest -q tests/test_modern_scan.py \
  tests/test_modern_point_evaluation.py \
  tests/test_modern_point_artifacts.py \
  tests/test_modern_input_registry.py \
  tests/test_modern_phenomenology.py
rm -rf /tmp/modern-scan-smoke3
python -m quarkConstraints.modern.scan write-preset smoke /tmp/modern-scan-smoke3/config.json
python -m quarkConstraints.modern.scan run-shard \
  --config /tmp/modern-scan-smoke3/config.json \
  --output-dir /tmp/modern-scan-smoke3 \
  --shard-index 0 \
  --shard-count 1
python -m quarkConstraints.modern.scan merge --run-dir /tmp/modern-scan-smoke3
python -m quarkConstraints.modern.scan verify \
  --run-dir /tmp/modern-scan-smoke3 \
  --config /tmp/modern-scan-smoke3/config.json \
  --total-shards 1
```

Validation commands executed after landing the explicit QS2 bridge subset:

```bash
python -m py_compile \
  quarkConstraints/modern/couplings.py \
  quarkConstraints/modern/matching.py \
  quarkConstraints/modern/evaluation.py \
  tests/test_modern_couplings.py \
  tests/test_modern_matching.py \
  tests/test_modern_point_evaluation.py
pytest -q tests/test_modern_couplings.py \
  tests/test_modern_matching.py \
  tests/test_modern_point_evaluation.py
pytest -q tests/test_modern_point_artifacts.py tests/test_modern_scan.py
pytest -q tests/test_modern_scan.py \
  tests/test_modern_point_evaluation.py \
  tests/test_modern_point_artifacts.py \
  tests/test_modern_input_registry.py \
  tests/test_modern_phenomenology.py \
  tests/test_modern_couplings.py \
  tests/test_modern_matching.py
rm -rf /tmp/modern-scan-smoke4
python -m quarkConstraints.modern.scan write-preset smoke /tmp/modern-scan-smoke4/config.json
python -m quarkConstraints.modern.scan run-shard \
  --config /tmp/modern-scan-smoke4/config.json \
  --output-dir /tmp/modern-scan-smoke4 \
  --shard-index 0 \
  --shard-count 1
python -m quarkConstraints.modern.scan merge --run-dir /tmp/modern-scan-smoke4
python -m quarkConstraints.modern.scan verify \
  --run-dir /tmp/modern-scan-smoke4 \
  --config /tmp/modern-scan-smoke4/config.json \
  --total-shards 1
```

Validation commands executed after landing the bridge-sidecar export slice:

```bash
python -m py_compile \
  quarkConstraints/modern/artifacts.py \
  quarkConstraints/modern/bridge_artifacts.py \
  quarkConstraints/modern/verifier.py \
  quarkConstraints/modern/scan.py \
  quarkConstraints/modern/__init__.py \
  tests/test_modern_point_artifacts.py \
  tests/test_modern_bridge_artifacts.py \
  tests/test_modern_scan.py
pytest -q \
  tests/test_modern_point_artifacts.py::test_modern_verifier_has_artifact_only_import_isolation \
  tests/test_modern_bridge_artifacts.py \
  tests/test_modern_scan.py
pytest -q tests/test_modern_scan.py \
  tests/test_modern_point_evaluation.py \
  tests/test_modern_point_artifacts.py \
  tests/test_modern_bridge_artifacts.py \
  tests/test_modern_input_registry.py \
  tests/test_modern_phenomenology.py \
  tests/test_modern_couplings.py \
  tests/test_modern_matching.py
rm -rf /tmp/modern-scan-smoke5
python -m quarkConstraints.modern.scan write-preset smoke /tmp/modern-scan-smoke5/config.json
python -m quarkConstraints.modern.scan run-shard \
  --config /tmp/modern-scan-smoke5/config.json \
  --output-dir /tmp/modern-scan-smoke5 \
  --shard-index 0 \
  --shard-count 1
python -m quarkConstraints.modern.scan merge --run-dir /tmp/modern-scan-smoke5
python -m quarkConstraints.modern.scan verify \
  --run-dir /tmp/modern-scan-smoke5 \
  --config /tmp/modern-scan-smoke5/config.json \
  --total-shards 1
```

Validation commands executed after the scan-note honesty cleanup:

```bash
pytest -q tests/test_modern_scan.py \
  tests/test_modern_point_evaluation.py \
  tests/test_modern_point_artifacts.py \
  tests/test_modern_bridge_artifacts.py \
  tests/test_modern_input_registry.py \
  tests/test_modern_phenomenology.py \
  tests/test_modern_couplings.py \
  tests/test_modern_matching.py
rm -rf /tmp/modern-scan-smoke7
python -m quarkConstraints.modern.scan write-preset smoke /tmp/modern-scan-smoke7/config.json
python -m quarkConstraints.modern.scan run-shard \
  --config /tmp/modern-scan-smoke7/config.json \
  --output-dir /tmp/modern-scan-smoke7 \
  --shard-index 0 \
  --shard-count 1
python -m quarkConstraints.modern.scan merge --run-dir /tmp/modern-scan-smoke7
python -m quarkConstraints.modern.scan verify \
  --run-dir /tmp/modern-scan-smoke7 \
  --config /tmp/modern-scan-smoke7/config.json \
  --total-shards 1
```

Pilot and production commands ready now:

```bash
# Direct local smoke run without SLURM
MODERN_SCAN_PRESET=smoke \
MODERN_SCAN_OUTPUT_ROOT=/tmp/modern_quark_scan_smoke3 \
MODERN_SCAN_TOTAL_SHARDS=1 \
bash scripts/run_modern_quark_scan_shard.sbatch

MODERN_SCAN_OUTPUT_ROOT=/tmp/modern_quark_scan_smoke3 \
bash scripts/merge_modern_quark_scan.sbatch

MODERN_SCAN_OUTPUT_ROOT=/tmp/modern_quark_scan_smoke3 \
MODERN_SCAN_TOTAL_SHARDS=1 \
bash scripts/verify_modern_quark_scan.sbatch

# Dry-run the SLURM submission pipeline
MODERN_SCAN_SUBMIT_DRY_RUN=1 \
MODERN_SCAN_PRESET=pilot \
MODERN_SCAN_TOTAL_SHARDS=4 \
MODERN_SCAN_OUTPUT_ROOT=/tmp/modern_quark_scan_pilot3 \
bash scripts/submit_modern_quark_scan_pipeline.sh

# Real SLURM pilot submission
MODERN_SCAN_PRESET=pilot \
MODERN_SCAN_TOTAL_SHARDS=4 \
MODERN_SCAN_OUTPUT_ROOT=scan_outputs/modern_quark_pilot \
bash scripts/submit_modern_quark_scan_pipeline.sh
```

## Produced Artifacts / Manifest Locations

- Durable handoff/status file: [`ALLOWED_REGION_SCAN_STATUS.md`](./ALLOWED_REGION_SCAN_STATUS.md)
- Quark-specific wrapper entry points:
  - [`scripts/run_modern_quark_scan_shard.sbatch`](../scripts/run_modern_quark_scan_shard.sbatch)
  - [`scripts/merge_modern_quark_scan.sbatch`](../scripts/merge_modern_quark_scan.sbatch)
  - [`scripts/verify_modern_quark_scan.sbatch`](../scripts/verify_modern_quark_scan.sbatch)
  - [`scripts/submit_modern_quark_scan_pipeline.sh`](../scripts/submit_modern_quark_scan_pipeline.sh)
- Temporary smoke run produced during this slice:
  - run config: `/tmp/modern-scan-smoke2/config.json`
  - shard manifest: `/tmp/modern-scan-smoke2/shards/shard-00000-of-00001/manifest.json`
  - shard results: `/tmp/modern-scan-smoke2/shards/shard-00000-of-00001/results.jsonl`
  - merged manifest: `/tmp/modern-scan-smoke2/merged/manifest.json`
  - merged results: `/tmp/modern-scan-smoke2/merged/results.jsonl`
  - merged verification: `/tmp/modern-scan-smoke2/merged/verification.json`
- Revalidated smoke run after the QS4 slice and config-snapshot fix:
  - preset config: `/tmp/modern-scan-smoke3/config.json`
  - shard manifest: `/tmp/modern-scan-smoke3/shards/shard-00000-of-00001/manifest.json`
  - merged manifest: `/tmp/modern-scan-smoke3/merged/manifest.json`
  - merged verification: `/tmp/modern-scan-smoke3/merged/verification.json`
- Revalidated smoke run after the explicit QS2 bridge subset:
  - preset config: `/tmp/modern-scan-smoke4/config.json`
  - shard manifest: `/tmp/modern-scan-smoke4/shards/shard-00000-of-00001/manifest.json`
  - merged manifest: `/tmp/modern-scan-smoke4/merged/manifest.json`
  - merged verification: `/tmp/modern-scan-smoke4/merged/verification.json`
- Revalidated smoke run after the QS5 scan/verifier integration repair:
  - preset config: `/tmp/modern-scan-smoke7/config.json`
  - shard manifest: `/tmp/modern-scan-smoke7/shards/shard-00000-of-00001/manifest.json`
  - point artifacts: `/tmp/modern-scan-smoke7/shards/shard-00000-of-00001/artifacts/*.json`
  - bridge artifacts: `/tmp/modern-scan-smoke7/shards/shard-00000-of-00001/bridge_artifacts/*.bridge.json`
  - phenomenology artifacts: `/tmp/modern-scan-smoke7/shards/shard-00000-of-00001/phenomenology_artifacts/*.phenomenology.json`
  - merged manifest: `/tmp/modern-scan-smoke7/merged/manifest.json`
  - merged verification: `/tmp/modern-scan-smoke7/merged/verification.json`
- Revalidated smoke run after the bridge-sidecar export slice:
  - preset config: `/tmp/modern-scan-smoke5/config.json`
  - shard manifest: `/tmp/modern-scan-smoke5/shards/shard-00000-of-00001/manifest.json`
  - point artifacts: `/tmp/modern-scan-smoke5/shards/shard-00000-of-00001/artifacts/*.json`
  - bridge artifacts: `/tmp/modern-scan-smoke5/shards/shard-00000-of-00001/bridge_artifacts/*.bridge.json`
  - merged manifest: `/tmp/modern-scan-smoke5/merged/manifest.json`
  - merged verification: `/tmp/modern-scan-smoke5/merged/verification.json`
- Wrapper validation smoke run produced during this slice:
  - frozen config: `/tmp/modern_quark_scan_smoke2/config.json`
  - preset config source: `/tmp/modern_quark_scan_smoke2/configs/smoke.json`
  - shard manifest: `/tmp/modern_quark_scan_smoke2/shards/shard-00000-of-00001/manifest.json`
  - merged manifest: `/tmp/modern_quark_scan_smoke2/merged/manifest.json`
  - merged verification: `/tmp/modern_quark_scan_smoke2/merged/verification.json`
- Expected reusable run layout for future scans:
  - `RUN_DIR/config.json`
  - `RUN_DIR/shards/shard-XXXXX-of-YYYYY/manifest.json`
  - `RUN_DIR/shards/shard-XXXXX-of-YYYYY/results.jsonl`
  - `RUN_DIR/shards/shard-XXXXX-of-YYYYY/artifacts/*.json`
  - `RUN_DIR/shards/shard-XXXXX-of-YYYYY/bridge_artifacts/*.bridge.json`
  - `RUN_DIR/shards/shard-XXXXX-of-YYYYY/phenomenology_artifacts/*.phenomenology.json`
  - `RUN_DIR/shards/shard-XXXXX-of-YYYYY/cache/*.seed.json`
  - `RUN_DIR/merged/manifest.json`
  - `RUN_DIR/merged/results.jsonl`
  - `RUN_DIR/merged/verification.json`
