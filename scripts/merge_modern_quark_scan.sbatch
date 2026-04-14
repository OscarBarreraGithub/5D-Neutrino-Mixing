#!/bin/bash
# Merge wrapper for the operational modern-quark scan pipeline.
#
# Required exports:
#   one of RUN_DIR or MODERN_SCAN_OUTPUT_ROOT

#SBATCH --job-name=modern-qmerge
#SBATCH --account=randall_lab
#SBATCH --partition=sapphire
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=slurm-modern-qmerge-%j.out
#SBATCH --error=slurm-modern-qmerge-%j.err

set -euo pipefail

ROOT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$ROOT_DIR"

RUN_DIR="${RUN_DIR:-${MODERN_SCAN_OUTPUT_ROOT:-}}"
if [[ -z "$RUN_DIR" ]]; then
  echo "ERROR: RUN_DIR or MODERN_SCAN_OUTPUT_ROOT must be set."
  exit 1
fi
REQUIRE_COMPLETE="${MODERN_SCAN_REQUIRE_COMPLETE:-1}"
MERGED_OUTPUT_DIR="${MODERN_SCAN_MERGED_OUTPUT_DIR:-}"
LOG_DIR="${RUN_DIR}/logs"

mkdir -p "$LOG_DIR"

exec > >(tee -a "${LOG_DIR}/merge_${SLURM_JOB_ID:-manual}.out")
exec 2> >(tee -a "${LOG_DIR}/merge_${SLURM_JOB_ID:-manual}.err" >&2)

echo "RUN_DIR=${RUN_DIR}"

merge_args=(
  merge
  --run-dir "${RUN_DIR}"
)
if [[ -n "$MERGED_OUTPUT_DIR" ]]; then
  merge_args+=(--output-dir "${MERGED_OUTPUT_DIR}")
fi
if [[ "$REQUIRE_COMPLETE" == "0" ]]; then
  merge_args+=(--no-require-complete)
else
  merge_args+=(--require-complete)
fi

python -m quarkConstraints.modern.scan "${merge_args[@]}"
