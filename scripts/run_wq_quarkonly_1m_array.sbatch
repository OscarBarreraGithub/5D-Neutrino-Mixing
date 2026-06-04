#!/bin/bash
#SBATCH --job-name=wq_qonly_1m
#SBATCH --partition=serial_requeue
#SBATCH --account=randall_lab
#SBATCH --cpus-per-task=48
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --array=0-49%50
#SBATCH --output=logs/wq_quarkonly_1M_%A_%a.out
#SBATCH --error=logs/wq_quarkonly_1M_%A_%a.err

# WQ-1M quark-only r x M_KK array.
#
# Physics plan:
#   r = 0.05,0.10,0.25,0.50,1.00
#   M_KK = 1,2,3,5,7,10,15,20,30,50 TeV
#   20k draws per (r,M_KK): 10 M_KK x 5 r x 20k = 1,000,000 draws
#
# Sharding:
#   10 draw shards per r; each array task runs all 10 M_KK tiles for one
#   (r, draw-shard) pair with 2k draws/M_KK/tile.
#
# Seed proof:
#   per-shard base = 202606040000 + 20000000 * task_id
#   tile seed      = base + 1000003 * tile_id
#   draw seed      = tile seed + draw_idx
#   chosen block 20000000 > stride*n_tiles + n_draws = 1000003*10 + 2000.
#   Therefore no draw seed interval overlaps within a shard or across shards.
#
# Merge/analysis:
#   No physical merge is needed. After all array tasks finish, run:
#     python scripts/analyze_wq_quarkonly.py --input scan_outputs/wq_quarkonly_1M_<array_job_id>
#   The analysis reads tile-*.jsonl recursively under every r*/shard-* directory
#   and writes plots/ plus analysis_report.md at the output root.
#
# Do not submit this array until the independent codex and Opus reviewers both approve.

set -euo pipefail

REPO_ROOT="/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing"
cd "${REPO_ROOT}"

source "${HOME}/.bashrc"
conda activate ising_bootstrap

mkdir -p "${REPO_ROOT}/logs"

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
eval "$("${REPO_ROOT}/scripts/wq_quarkonly_1m_plan.py" task-env --task-id "${TASK_ID}")"

WQ_OUTPUT_RUN_ID="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"
OUTPUT_ROOT="${WQ_OUTPUT_ROOT:-${REPO_ROOT}/scan_outputs/wq_quarkonly_1M_${WQ_OUTPUT_RUN_ID}}"
OUTPUT_DIR="${OUTPUT_ROOT}/${WQ_R_LABEL}/shard-${WQ_DRAW_SHARD_PADDED}"
N_WORKERS="${N_WORKERS:-${SLURM_CPUS_PER_TASK:-48}}"

mkdir -p "${OUTPUT_ROOT}" "${OUTPUT_DIR}"

if [[ "${TASK_ID}" == "0" ]]; then
    "${REPO_ROOT}/scripts/wq_quarkonly_1m_plan.py" summary > "${OUTPUT_ROOT}/scan_plan.json"
    "${REPO_ROOT}/scripts/wq_quarkonly_1m_plan.py" seed-proof > "${OUTPUT_ROOT}/seed_disjointness.md"
fi

echo "[wq_quarkonly_1M] OUTPUT_ROOT = ${OUTPUT_ROOT}"
echo "[wq_quarkonly_1M] OUTPUT_DIR  = ${OUTPUT_DIR}"
echo "[wq_quarkonly_1M] task        = ${WQ_TASK_ID}/${WQ_TOTAL_TASKS}"
echo "[wq_quarkonly_1M] r           = ${WQ_R_VALUE}"
echo "[wq_quarkonly_1M] draw shard  = ${WQ_DRAW_SHARD} (${WQ_DRAW_START}:${WQ_DRAW_STOP})"
echo "[wq_quarkonly_1M] base seed   = ${WQ_BASE_SEED}"
echo "[wq_quarkonly_1M] M_KK (TeV)  = ${WQ_M_KK_TEV}"
echo "[wq_quarkonly_1M] workers     = ${N_WORKERS}"
echo "[wq_quarkonly_1M] git HEAD    = $(git rev-parse HEAD)"

python "${REPO_ROOT}/scripts/run_full_catalog_scan.py" \
    --quark-only \
    --output-dir "${OUTPUT_DIR}" \
    --n-draws "${WQ_DRAWS_PER_MKK_PER_SHARD}" \
    --m-kk-tev "${WQ_M_KK_TEV}" \
    --base-seed "${WQ_BASE_SEED}" \
    --tile-seed-stride "${WQ_TILE_SEED_STRIDE}" \
    --quark-fit-r "${WQ_R_VALUE}" \
    --n-workers "${N_WORKERS}"

echo "[wq_quarkonly_1M] shard complete: ${OUTPUT_DIR}"
