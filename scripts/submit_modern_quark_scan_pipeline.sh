#!/bin/bash
# Submit the modern quark scan shard -> merge -> verify pipeline.
#
# Smoke submission:
#   MODERN_SCAN_PRESET=smoke \
#   MODERN_SCAN_OUTPUT_ROOT=scan_outputs/modern_quark_smoke \
#   scripts/submit_modern_quark_scan_pipeline.sh
#
# Pilot submission:
#   MODERN_SCAN_PRESET=pilot \
#   MODERN_SCAN_TOTAL_SHARDS=4 \
#   MODERN_SCAN_OUTPUT_ROOT=scan_outputs/modern_quark_pilot \
#   scripts/submit_modern_quark_scan_pipeline.sh

set -euo pipefail

ROOT_DIR="${PWD}"
cd "$ROOT_DIR"

SCAN_SCRIPT="${MODERN_SCAN_RUN_SCRIPT:-scripts/run_modern_quark_scan_shard.sbatch}"
MERGE_SCRIPT="${MODERN_SCAN_MERGE_SCRIPT:-scripts/merge_modern_quark_scan.sbatch}"
VERIFY_SCRIPT="${MODERN_SCAN_VERIFY_SCRIPT:-scripts/verify_modern_quark_scan.sbatch}"

for script in "$SCAN_SCRIPT" "$MERGE_SCRIPT" "$VERIFY_SCRIPT"; do
  if [[ ! -f "$script" ]]; then
    echo "ERROR: Missing pipeline script: $script"
    exit 1
  fi
done

MODERN_SCAN_PRESET="${MODERN_SCAN_PRESET:-${PRESET:-}}"
MODERN_SCAN_OUTPUT_ROOT="${MODERN_SCAN_OUTPUT_ROOT:-${OUTPUT_ROOT:-scan_outputs/modern_quark_scan}}"
MODERN_SCAN_CONFIG_PATH="${MODERN_SCAN_CONFIG_PATH:-${CONFIG_PATH:-}}"
MODERN_SCAN_REQUIRE_COMPLETE="${MODERN_SCAN_REQUIRE_COMPLETE:-1}"
MODERN_SCAN_SUBMIT_DRY_RUN="${MODERN_SCAN_SUBMIT_DRY_RUN:-${DRY_RUN:-0}}"

if [[ -z "${MODERN_SCAN_TOTAL_SHARDS:-}" ]]; then
  case "${MODERN_SCAN_PRESET}" in
    smoke)
      MODERN_SCAN_TOTAL_SHARDS=1
      ;;
    pilot)
      MODERN_SCAN_TOTAL_SHARDS=4
      ;;
    *)
      MODERN_SCAN_TOTAL_SHARDS=16
      ;;
  esac
fi

if ! [[ "$MODERN_SCAN_TOTAL_SHARDS" =~ ^[0-9]+$ ]] || [[ "$MODERN_SCAN_TOTAL_SHARDS" -lt 1 ]]; then
  echo "ERROR: MODERN_SCAN_TOTAL_SHARDS must be a positive integer, got '$MODERN_SCAN_TOTAL_SHARDS'."
  exit 1
fi

mkdir -p "${MODERN_SCAN_OUTPUT_ROOT}/logs"

if [[ -z "$MODERN_SCAN_CONFIG_PATH" ]]; then
  if [[ -n "$MODERN_SCAN_PRESET" ]]; then
    MODERN_SCAN_CONFIG_PATH="${MODERN_SCAN_OUTPUT_ROOT}/configs/${MODERN_SCAN_PRESET}.json"
  else
    MODERN_SCAN_CONFIG_PATH="${MODERN_SCAN_OUTPUT_ROOT}/configs/env.json"
  fi
fi

EXPORT_VARS=(
  "MODERN_SCAN_OUTPUT_ROOT=${MODERN_SCAN_OUTPUT_ROOT}"
  "MODERN_SCAN_TOTAL_SHARDS=${MODERN_SCAN_TOTAL_SHARDS}"
  "MODERN_SCAN_CONFIG_PATH=${MODERN_SCAN_CONFIG_PATH}"
  "MODERN_SCAN_PRESET=${MODERN_SCAN_PRESET}"
  "MODERN_SCAN_REQUIRE_COMPLETE=${MODERN_SCAN_REQUIRE_COMPLETE}"
  "MODERN_SCAN_MAX_POINTS=${MODERN_SCAN_MAX_POINTS:-}"
  "MODERN_SCAN_PROGRESS_EVERY=${MODERN_SCAN_PROGRESS_EVERY:-0}"
  "MODERN_SCAN_R_VALUES=${MODERN_SCAN_R_VALUES:-}"
  "MODERN_SCAN_OVERALL_SCALE_VALUES=${MODERN_SCAN_OVERALL_SCALE_VALUES:-}"
  "MODERN_SCAN_LAMBDA_IR_VALUES=${MODERN_SCAN_LAMBDA_IR_VALUES:-}"
  "MODERN_SCAN_XI_KK=${MODERN_SCAN_XI_KK:-}"
  "MODERN_SCAN_K=${MODERN_SCAN_K:-}"
  "MODERN_SCAN_MAX_NFEV=${MODERN_SCAN_MAX_NFEV:-}"
  "MODERN_SCAN_FIT_ORIENTATION=${MODERN_SCAN_FIT_ORIENTATION:-}"
  "MODERN_SCAN_RECORD_GIT_METADATA=${MODERN_SCAN_RECORD_GIT_METADATA:-}"
  "MODERN_SCAN_MERGED_OUTPUT_DIR=${MODERN_SCAN_MERGED_OUTPUT_DIR:-}"
  "MODERN_SCAN_MERGED_DIR=${MODERN_SCAN_MERGED_DIR:-}"
  "MODERN_SCAN_MERGED_MANIFEST_PATH=${MODERN_SCAN_MERGED_MANIFEST_PATH:-}"
  "MODERN_SCAN_MERGED_RESULTS_PATH=${MODERN_SCAN_MERGED_RESULTS_PATH:-}"
  "MODERN_SCAN_VERIFICATION_DIR=${MODERN_SCAN_VERIFICATION_DIR:-}"
  "MODERN_SCAN_VERIFICATION_REPORT_PATH=${MODERN_SCAN_VERIFICATION_REPORT_PATH:-}"
)

EXPORT_ARG="ALL,$(IFS=,; echo "${EXPORT_VARS[*]}")"
ARRAY_SPEC="0-$((MODERN_SCAN_TOTAL_SHARDS - 1))"

if [[ "$MODERN_SCAN_SUBMIT_DRY_RUN" == "1" ]]; then
  echo "Dry run only. Equivalent commands:"
  echo "sbatch --parsable --array=${ARRAY_SPEC} --export=${EXPORT_ARG} ${SCAN_SCRIPT}"
  echo "sbatch --parsable --dependency=afterok:<scan_job_id> --export=${EXPORT_ARG} ${MERGE_SCRIPT}"
  echo "sbatch --parsable --dependency=afterok:<merge_job_id> --export=${EXPORT_ARG} ${VERIFY_SCRIPT}"
  exit 0
fi

if ! command -v sbatch >/dev/null 2>&1; then
  echo "ERROR: sbatch not found in PATH."
  exit 1
fi

scan_job_id=$(
  sbatch \
    --parsable \
    --array="${ARRAY_SPEC}" \
    --export="${EXPORT_ARG}" \
    "${SCAN_SCRIPT}"
)

merge_job_id=$(
  sbatch \
    --parsable \
    --dependency="afterok:${scan_job_id}" \
    --export="${EXPORT_ARG}" \
    "${MERGE_SCRIPT}"
)

verify_job_id=$(
  sbatch \
    --parsable \
    --dependency="afterok:${merge_job_id}" \
    --export="${EXPORT_ARG}" \
    "${VERIFY_SCRIPT}"
)

echo "Submitted modern scan array job: ${scan_job_id} (array ${ARRAY_SPEC})"
echo "Submitted modern merge job:      ${merge_job_id} (afterok:${scan_job_id})"
echo "Submitted modern verify job:     ${verify_job_id} (afterok:${merge_job_id})"
echo "Output root: ${MODERN_SCAN_OUTPUT_ROOT}"
echo "Config path: ${MODERN_SCAN_CONFIG_PATH}"
echo "Track with: squeue -j ${scan_job_id},${merge_job_id},${verify_job_id}"
