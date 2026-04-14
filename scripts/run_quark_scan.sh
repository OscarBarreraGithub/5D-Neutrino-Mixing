#!/bin/bash
# Quark-sector allowed-region scan: write config, submit shard array, then merge+verify.
#
# Usage:
#   bash scripts/run_quark_scan.sh exploratory   # 64 points, test partition
#   bash scripts/run_quark_scan.sh production     # 1200 points, serial_requeue
#
set -euo pipefail

PRESET="${1:-exploratory}"
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TIMESTAMP="$(date +%Y%m%dT%H%M%S)"

ACCOUNT="randall_lab"

if [[ "$PRESET" == "production" ]]; then
    TOTAL_SHARDS=120
    PARTITION="serial_requeue"
    TIME="01:00:00"
    MEM="4G"
    RUN_DIR="${REPO_ROOT}/scan_outputs/production_${TIMESTAMP}"
elif [[ "$PRESET" == "dense" || "$PRESET" == "dense_wide_yukawa" ]]; then
    TOTAL_SHARDS=200
    PARTITION="serial_requeue"
    TIME="02:00:00"
    MEM="4G"
    RUN_DIR="${REPO_ROOT}/scan_outputs/${PRESET}_${TIMESTAMP}"
elif [[ "$PRESET" == "exploratory" ]]; then
    TOTAL_SHARDS=8
    PARTITION="serial_requeue"
    TIME="00:30:00"
    MEM="4G"
    RUN_DIR="${REPO_ROOT}/scan_outputs/exploratory_${TIMESTAMP}"
else
    echo "Unknown preset: $PRESET (use 'exploratory' or 'production')"
    exit 1
fi

echo "=== Quark scan: preset=$PRESET, shards=$TOTAL_SHARDS, partition=$PARTITION ==="
echo "=== Run directory: $RUN_DIR ==="

# Write the frozen config
mkdir -p "$RUN_DIR"
python -m quarkConstraints.modern.scan write-preset "$PRESET" "$RUN_DIR/config.json"
echo "Config written to $RUN_DIR/config.json"

# Submit shard array job
SHARD_JOB=$(sbatch --parsable \
    --partition="$PARTITION" \
    --account="$ACCOUNT" \
    --job-name="qscan_shard_${PRESET}" \
    --array="0-$((TOTAL_SHARDS - 1))" \
    --time="$TIME" \
    --mem="$MEM" \
    --cpus-per-task=1 \
    --output="$RUN_DIR/slurm_shard_%a.out" \
    --error="$RUN_DIR/slurm_shard_%a.err" \
    --wrap="python -m quarkConstraints.modern.scan run-shard \
        --config $RUN_DIR/config.json \
        --output-dir $RUN_DIR \
        --shard-index \$SLURM_ARRAY_TASK_ID \
        --shard-count $TOTAL_SHARDS")
echo "Shard array job submitted: $SHARD_JOB (tasks 0-$((TOTAL_SHARDS - 1)))"

# Submit merge job (depends on all shards)
MERGE_JOB=$(sbatch --parsable \
    --partition="$PARTITION" \
    --account="$ACCOUNT" \
    --job-name="qscan_merge_${PRESET}" \
    --dependency="afterok:${SHARD_JOB}" \
    --time="00:10:00" \
    --mem="4G" \
    --cpus-per-task=1 \
    --output="$RUN_DIR/slurm_merge.out" \
    --error="$RUN_DIR/slurm_merge.err" \
    --wrap="python -m quarkConstraints.modern.scan merge --run-dir $RUN_DIR")
echo "Merge job submitted: $MERGE_JOB (depends on $SHARD_JOB)"

# Submit verify job (depends on merge)
VERIFY_JOB=$(sbatch --parsable \
    --partition="$PARTITION" \
    --account="$ACCOUNT" \
    --job-name="qscan_verify_${PRESET}" \
    --dependency="afterok:${MERGE_JOB}" \
    --time="00:10:00" \
    --mem="4G" \
    --cpus-per-task=1 \
    --output="$RUN_DIR/slurm_verify.out" \
    --error="$RUN_DIR/slurm_verify.err" \
    --wrap="python -m quarkConstraints.modern.scan verify \
        --run-dir $RUN_DIR \
        --config $RUN_DIR/config.json \
        --total-shards $TOTAL_SHARDS")
echo "Verify job submitted: $VERIFY_JOB (depends on $MERGE_JOB)"

echo ""
echo "=== Pipeline submitted ==="
echo "  Shards: $SHARD_JOB (array 0-$((TOTAL_SHARDS - 1)))"
echo "  Merge:  $MERGE_JOB"
echo "  Verify: $VERIFY_JOB"
echo "  Run dir: $RUN_DIR"
echo ""
echo "Monitor with: squeue -u \$USER -j ${SHARD_JOB},${MERGE_JOB},${VERIFY_JOB}"
echo "Results will be in: $RUN_DIR/merged/results.jsonl"
