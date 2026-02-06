#!/bin/bash
# Submit scan shards, then submit dependent reclassification job.
#
# Example:
#   TOTAL_SHARDS=32 N_C_L=121 N_C_N=121 N_LAMBDA=9 N_MLIGHTEST=41 \
#   RECLASS_ANARCHY_MIN_SCORE=-2.5 \
#   scripts/submit_scan_pipeline.sh

set -euo pipefail

ROOT_DIR="${PWD}"
cd "$ROOT_DIR"

if ! command -v sbatch >/dev/null 2>&1; then
  echo "ERROR: sbatch not found in PATH."
  exit 1
fi

TOTAL_SHARDS="${TOTAL_SHARDS:-16}"
if ! [[ "$TOTAL_SHARDS" =~ ^[0-9]+$ ]] || [[ "$TOTAL_SHARDS" -lt 1 ]]; then
  echo "ERROR: TOTAL_SHARDS must be a positive integer, got '$TOTAL_SHARDS'."
  exit 1
fi

SCAN_SCRIPT="${SCAN_SCRIPT:-scripts/run_scan_sharded.sbatch}"
RECLASS_SCRIPT="${RECLASS_SCRIPT:-scripts/reclassify_scan_shards.sbatch}"

if [[ ! -f "$SCAN_SCRIPT" ]]; then
  echo "ERROR: Missing scan script: $SCAN_SCRIPT"
  exit 1
fi
if [[ ! -f "$RECLASS_SCRIPT" ]]; then
  echo "ERROR: Missing reclass script: $RECLASS_SCRIPT"
  exit 1
fi

array_spec="0-$((TOTAL_SHARDS - 1))"
scan_job_id=$(
  sbatch \
    --parsable \
    --array="$array_spec" \
    --export="ALL,TOTAL_SHARDS=${TOTAL_SHARDS}" \
    "$SCAN_SCRIPT"
)

reclass_job_id=$(
  sbatch \
    --parsable \
    --dependency="afterok:${scan_job_id}" \
    --export="ALL,TOTAL_SHARDS=${TOTAL_SHARDS}" \
    "$RECLASS_SCRIPT"
)

echo "Submitted scan array job: ${scan_job_id} (array ${array_spec})"
echo "Submitted reclass job:   ${reclass_job_id} (afterok:${scan_job_id})"
echo "Track with: squeue -j ${scan_job_id},${reclass_job_id}"
