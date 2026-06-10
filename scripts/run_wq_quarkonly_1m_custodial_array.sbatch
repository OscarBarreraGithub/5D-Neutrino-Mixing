#!/bin/bash
#SBATCH --job-name=wq_qonly_1m_cust
#SBATCH --partition=serial_requeue
#SBATCH --account=randall_lab
#SBATCH --cpus-per-task=48
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --array=0-49%50
#SBATCH --output=logs/wq_quarkonly_1M_custodial_%A_%a.out
#SBATCH --error=logs/wq_quarkonly_1M_custodial_%A_%a.err

# WQ-1M quark-only custodial r x M_KK array.
#
# Physics plan:
#   r = 0.05,0.10,0.25,0.50,1.00
#   M_KK = 1,2,3,5,7,10,15,20,30,50 TeV
#   20k draws per (r,M_KK): 10 M_KK x 5 r x 20k = 1,000,000 draws
#   EW model = custodial_rs_plr, tree-level custodial protection included.
#   W8 top-partner loops are deferred for this W9 scan.
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
#     python scripts/analyze_wq_quarkonly.py --input scan_outputs/wq_quarkonly_1M_custodial_<array_job_id> --output-dir scan_outputs/wq_quarkonly_1M_custodial_<array_job_id>/analysis
#     python scripts/build_wq_quarkonly_comparison.py --minimal scan_outputs/wq_quarkonly_1M_20128400 --custodial scan_outputs/wq_quarkonly_1M_custodial_<array_job_id> --output scan_outputs/wq_quarkonly_1M_custodial_<array_job_id>/comparison
#
# Do not submit this array until the W9 implementation is reviewed and committed.

set -euo pipefail

REPO_ROOT="/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing"
cd "${REPO_ROOT}"

source "${HOME}/.bashrc"
conda activate ising_bootstrap

mkdir -p "${REPO_ROOT}/logs"

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
eval "$("${REPO_ROOT}/scripts/wq_quarkonly_1m_plan.py" task-env --task-id "${TASK_ID}")"

WQ_OUTPUT_RUN_ID="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"
OUTPUT_ROOT="${WQ_OUTPUT_ROOT:-${REPO_ROOT}/scan_outputs/wq_quarkonly_1M_custodial_${WQ_OUTPUT_RUN_ID}}"
OUTPUT_DIR="${OUTPUT_ROOT}/${WQ_R_LABEL}/shard-${WQ_DRAW_SHARD_PADDED}"
N_WORKERS="${N_WORKERS:-${SLURM_CPUS_PER_TASK:-48}}"

mkdir -p "${OUTPUT_ROOT}" "${OUTPUT_DIR}"

if [[ "${TASK_ID}" == "0" ]]; then
    "${REPO_ROOT}/scripts/wq_quarkonly_1m_plan.py" summary > "${OUTPUT_ROOT}/scan_plan.json"
    "${REPO_ROOT}/scripts/wq_quarkonly_1m_plan.py" seed-proof > "${OUTPUT_ROOT}/seed_disjointness.md"
    MANIFEST_TMP="${OUTPUT_ROOT}/custodial_scan_manifest.json.tmp.$$"
    cat > "${MANIFEST_TMP}" <<EOF
{
  "schema": "wq_quarkonly_custodial_scan_manifest_v1",
  "ew_model": "custodial_rs_plr",
  "custodial_tree_status": "included",
  "custodial_top_partner_loop_status": "deferred",
  "custodial_fcnc_mode": "pr1_minimal_offdiag",
  "output_root": "${OUTPUT_ROOT}",
  "baseline_minimal_root": "${REPO_ROOT}/scan_outputs/wq_quarkonly_1M_20128400",
  "grid_source": "scripts/wq_quarkonly_1m_plan.py",
  "git_sha": "$(git rev-parse HEAD)"
}
EOF
    mv "${MANIFEST_TMP}" "${OUTPUT_ROOT}/custodial_scan_manifest.json"
fi

echo "[wq_quarkonly_1M_custodial] OUTPUT_ROOT = ${OUTPUT_ROOT}"
echo "[wq_quarkonly_1M_custodial] OUTPUT_DIR  = ${OUTPUT_DIR}"
echo "[wq_quarkonly_1M_custodial] task        = ${WQ_TASK_ID}/${WQ_TOTAL_TASKS}"
echo "[wq_quarkonly_1M_custodial] r           = ${WQ_R_VALUE}"
echo "[wq_quarkonly_1M_custodial] draw shard  = ${WQ_DRAW_SHARD} (${WQ_DRAW_START}:${WQ_DRAW_STOP})"
echo "[wq_quarkonly_1M_custodial] base seed   = ${WQ_BASE_SEED}"
echo "[wq_quarkonly_1M_custodial] M_KK (TeV)  = ${WQ_M_KK_TEV}"
echo "[wq_quarkonly_1M_custodial] ew model    = custodial_rs_plr"
echo "[wq_quarkonly_1M_custodial] workers     = ${N_WORKERS}"
echo "[wq_quarkonly_1M_custodial] git HEAD    = $(git rev-parse HEAD)"

python "${REPO_ROOT}/scripts/run_full_catalog_scan.py" \
    --quark-only \
    --ew-model custodial_rs_plr \
    --output-dir "${OUTPUT_DIR}" \
    --n-draws "${WQ_DRAWS_PER_MKK_PER_SHARD}" \
    --m-kk-tev "${WQ_M_KK_TEV}" \
    --base-seed "${WQ_BASE_SEED}" \
    --tile-seed-stride "${WQ_TILE_SEED_STRIDE}" \
    --quark-fit-r "${WQ_R_VALUE}" \
    --n-workers "${N_WORKERS}"

echo "[wq_quarkonly_1M_custodial] shard complete: ${OUTPUT_DIR}"
