#!/bin/bash
# Reclassify scan shard outputs after the scan array completes.
#
# Typical usage (usually submitted by submit_scan_pipeline.sh):
#   sbatch --export=ALL,TOTAL_SHARDS=16 scripts/reclassify_scan_shards.sbatch

#SBATCH --job-name=rs-reclass
#SBATCH --account=iaifi_lab
#SBATCH --partition=shared
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/reclass_%j.out
#SBATCH --error=logs/reclass_%j.err

set -euo pipefail

ROOT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$ROOT_DIR"
mkdir -p logs scan_outputs

TOTAL_SHARDS="${TOTAL_SHARDS:?TOTAL_SHARDS must be set (e.g. 16)}"
INPUT_DIR="${INPUT_DIR:-scan_outputs}"
OUTPUT_DIR="${OUTPUT_DIR:-scan_outputs}"
OUTPUT_TAG="${OUTPUT_TAG:-reclass}"

RECLASS_MAX_Y_BAR="${RECLASS_MAX_Y_BAR:-4.0}"
RECLASS_NATURAL_MIN="${RECLASS_NATURAL_MIN:-0.1}"
RECLASS_NATURAL_MAX="${RECLASS_NATURAL_MAX:-4.0}"
RECLASS_REQUIRE_LFV="${RECLASS_REQUIRE_LFV:-1}"        # 1 => require LFV, 0 => skip LFV
RECLASS_COMPUTE_ANARCHY="${RECLASS_COMPUTE_ANARCHY:-1}"  # 1 => include anarchy scoring
RECLASS_ANARCHY_MIN_SCORE="${RECLASS_ANARCHY_MIN_SCORE:-}"
RECLASS_ANARCHY_MAG_MIN="${RECLASS_ANARCHY_MAG_MIN:-0.3333333333333333}"
RECLASS_ANARCHY_MAG_MAX="${RECLASS_ANARCHY_MAG_MAX:-3.0}"
RECLASS_ANARCHY_PHASE_MIN="${RECLASS_ANARCHY_PHASE_MIN:-0.0}"
RECLASS_ANARCHY_PHASE_MAX="${RECLASS_ANARCHY_PHASE_MAX:-6.283185307179586}"
RECLASS_ANARCHY_YN_MIN="${RECLASS_ANARCHY_YN_MIN:-0.01}"
RECLASS_ANARCHY_YN_MAX="${RECLASS_ANARCHY_YN_MAX:-0.2}"
RECLASS_ANARCHY_W_BAND="${RECLASS_ANARCHY_W_BAND:-1.0}"
RECLASS_ANARCHY_W_COND="${RECLASS_ANARCHY_W_COND:-1.0}"
RECLASS_ANARCHY_W_FIT="${RECLASS_ANARCHY_W_FIT:-0.0}"

echo "Reclassifying shard outputs"
echo "TOTAL_SHARDS=${TOTAL_SHARDS}"
echo "INPUT_DIR=${INPUT_DIR}"
echo "OUTPUT_DIR=${OUTPUT_DIR}"
echo "OUTPUT_TAG=${OUTPUT_TAG}"

reclass_args=(
  --max-y-bar "$RECLASS_MAX_Y_BAR"
  --natural-min "$RECLASS_NATURAL_MIN"
  --natural-max "$RECLASS_NATURAL_MAX"
)

if [[ "$RECLASS_REQUIRE_LFV" == "0" ]]; then
  reclass_args+=(--skip-lfv)
fi

if [[ "$RECLASS_COMPUTE_ANARCHY" == "1" ]]; then
  reclass_args+=(
    --compute-anarchy
    --anarchy-magnitude-min "$RECLASS_ANARCHY_MAG_MIN"
    --anarchy-magnitude-max "$RECLASS_ANARCHY_MAG_MAX"
    --anarchy-phase-min "$RECLASS_ANARCHY_PHASE_MIN"
    --anarchy-phase-max "$RECLASS_ANARCHY_PHASE_MAX"
    --anarchy-yn-overall-min "$RECLASS_ANARCHY_YN_MIN"
    --anarchy-yn-overall-max "$RECLASS_ANARCHY_YN_MAX"
    --anarchy-w-band "$RECLASS_ANARCHY_W_BAND"
    --anarchy-w-cond "$RECLASS_ANARCHY_W_COND"
    --anarchy-w-fit "$RECLASS_ANARCHY_W_FIT"
  )
  if [[ -n "$RECLASS_ANARCHY_MIN_SCORE" ]]; then
    reclass_args+=(--anarchy-min-score "$RECLASS_ANARCHY_MIN_SCORE")
  fi
elif [[ -n "$RECLASS_ANARCHY_MIN_SCORE" ]]; then
  echo "ERROR: RECLASS_ANARCHY_MIN_SCORE is set but RECLASS_COMPUTE_ANARCHY=0"
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

for shard_id in $(seq 0 $((TOTAL_SHARDS - 1))); do
  in_csv=$(printf "%s/scan_shard_%03d_of_%03d.csv" "$INPUT_DIR" "$shard_id" "$TOTAL_SHARDS")
  out_csv=$(printf "%s/scan_shard_%03d_of_%03d_%s.csv" "$OUTPUT_DIR" "$shard_id" "$TOTAL_SHARDS" "$OUTPUT_TAG")

  if [[ ! -f "$in_csv" ]]; then
    echo "ERROR: Missing shard file: $in_csv"
    exit 1
  fi

  echo "Reclassifying shard ${shard_id}/${TOTAL_SHARDS}"
  python scripts/reclassify_scan.py "$in_csv" "$out_csv" "${reclass_args[@]}"
done

echo "Reclassification complete for ${TOTAL_SHARDS} shard files."
