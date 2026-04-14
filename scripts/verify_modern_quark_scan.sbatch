#!/bin/bash
# Verification wrapper for the operational modern-quark scan pipeline.
#
# Required exports:
#   one of RUN_DIR or MODERN_SCAN_OUTPUT_ROOT
#   one of TOTAL_SHARDS or MODERN_SCAN_TOTAL_SHARDS
#
# Optional exports:
#   one of SCAN_CONFIG or MODERN_SCAN_CONFIG_PATH (defaults to RUN_DIR/config.json)

#SBATCH --job-name=modern-qverify
#SBATCH --account=randall_lab
#SBATCH --partition=sapphire
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=slurm-modern-qverify-%j.out
#SBATCH --error=slurm-modern-qverify-%j.err

set -euo pipefail

ROOT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$ROOT_DIR"

RUN_DIR="${RUN_DIR:-${MODERN_SCAN_OUTPUT_ROOT:-}}"
TOTAL_SHARDS="${TOTAL_SHARDS:-${MODERN_SCAN_TOTAL_SHARDS:-}}"
if [[ -z "$RUN_DIR" ]]; then
  echo "ERROR: RUN_DIR or MODERN_SCAN_OUTPUT_ROOT must be set."
  exit 1
fi
if [[ -z "$TOTAL_SHARDS" ]]; then
  echo "ERROR: TOTAL_SHARDS or MODERN_SCAN_TOTAL_SHARDS must be set."
  exit 1
fi
SCAN_CONFIG="${SCAN_CONFIG:-${MODERN_SCAN_CONFIG_PATH:-${RUN_DIR}/config.json}}"
LOG_DIR="${RUN_DIR}/logs"

mkdir -p "$LOG_DIR"

exec > >(tee -a "${LOG_DIR}/verify_${SLURM_JOB_ID:-manual}.out")
exec 2> >(tee -a "${LOG_DIR}/verify_${SLURM_JOB_ID:-manual}.err" >&2)

echo "RUN_DIR=${RUN_DIR}"
echo "SCAN_CONFIG=${SCAN_CONFIG}"
echo "TOTAL_SHARDS=${TOTAL_SHARDS}"

python -m quarkConstraints.modern.scan \
  verify \
  --run-dir "${RUN_DIR}" \
  --config "${SCAN_CONFIG}" \
  --total-shards "${TOTAL_SHARDS}"
