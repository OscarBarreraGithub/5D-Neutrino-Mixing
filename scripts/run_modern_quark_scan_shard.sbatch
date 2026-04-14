#!/bin/bash
# Operational modern-quark scan shard wrapper.
#
# Required exports:
#   one of RUN_DIR or MODERN_SCAN_OUTPUT_ROOT
#   one of TOTAL_SHARDS or MODERN_SCAN_TOTAL_SHARDS
#
# Optional exports:
#   one of SCAN_CONFIG or MODERN_SCAN_CONFIG_PATH
#   SHARD_INDEX
#   MODERN_SCAN_PRESET=smoke|pilot

#SBATCH --job-name=modern-qscan
#SBATCH --account=randall_lab
#SBATCH --partition=sapphire
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=slurm-modern-qscan-%A_%a.out
#SBATCH --error=slurm-modern-qscan-%A_%a.err

set -euo pipefail

ROOT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$ROOT_DIR"

RUN_DIR="${RUN_DIR:-${MODERN_SCAN_OUTPUT_ROOT:-}}"
TOTAL_SHARDS="${TOTAL_SHARDS:-${MODERN_SCAN_TOTAL_SHARDS:-}}"
SCAN_CONFIG="${SCAN_CONFIG:-${MODERN_SCAN_CONFIG_PATH:-}}"
PRESET="${MODERN_SCAN_PRESET:-}"
if [[ -z "$RUN_DIR" ]]; then
  echo "ERROR: RUN_DIR or MODERN_SCAN_OUTPUT_ROOT must be set."
  exit 1
fi
if [[ -z "$TOTAL_SHARDS" ]]; then
  echo "ERROR: TOTAL_SHARDS or MODERN_SCAN_TOTAL_SHARDS must be set."
  exit 1
fi
SHARD_INDEX="${SLURM_ARRAY_TASK_ID:-${SHARD_INDEX:-0}}"
LOG_DIR="${RUN_DIR}/logs"

mkdir -p "$LOG_DIR"

if [[ -z "$SCAN_CONFIG" ]]; then
  if [[ -n "$PRESET" ]]; then
    SCAN_CONFIG="${RUN_DIR}/configs/${PRESET}.json"
  else
    echo "ERROR: SCAN_CONFIG/MODERN_SCAN_CONFIG_PATH not set and no MODERN_SCAN_PRESET given."
    exit 1
  fi
fi

if [[ ! -f "$SCAN_CONFIG" ]]; then
  if [[ -n "$PRESET" ]]; then
    mkdir -p "$(dirname "$SCAN_CONFIG")"
    python -m quarkConstraints.modern.scan write-preset "$PRESET" "$SCAN_CONFIG"
  else
    echo "ERROR: missing scan config: $SCAN_CONFIG"
    exit 1
  fi
fi

exec > >(tee -a "${LOG_DIR}/shard_${SLURM_ARRAY_JOB_ID:-manual}_${SHARD_INDEX}.out")
exec 2> >(tee -a "${LOG_DIR}/shard_${SLURM_ARRAY_JOB_ID:-manual}_${SHARD_INDEX}.err" >&2)

echo "RUN_DIR=${RUN_DIR}"
echo "SCAN_CONFIG=${SCAN_CONFIG}"
echo "TOTAL_SHARDS=${TOTAL_SHARDS}"
echo "SHARD_INDEX=${SHARD_INDEX}"

python -m quarkConstraints.modern.scan \
  run-shard \
  --config "${SCAN_CONFIG}" \
  --output-dir "${RUN_DIR}" \
  --shard-index "${SHARD_INDEX}" \
  --shard-count "${TOTAL_SHARDS}"
