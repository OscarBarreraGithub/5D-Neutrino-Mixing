#!/bin/bash
#SBATCH --job-name=rs_anarchy_c02c
#SBATCH --partition=serial_requeue
#SBATCH --account=randall_lab
#SBATCH --requeue
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=logs/rs_anarchy_c02c_%x_%j.out
#SBATCH --error=logs/rs_anarchy_c02c_%x_%j.err

# C02c: epsilon_K NP-budget band reruns.
#
# Submit one job per budget edge:
#   for e in low central high; do
#     sbatch --export=ALL,BUDGET_EDGE=$e scripts/run_rs_anarchy_c02c.sbatch
#   done
#
# Required env:  BUDGET_EDGE  (one of: low, central, high)
# Optional envs: N_DRAWS      (default 1600000)
#                M_KK_TEV     (default "3,5,7,10,15,20,30,50")
#                N_WORKERS    (default SLURM_CPUS_PER_TASK or 32)

set -euo pipefail

REPO_ROOT="/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing"
cd "${REPO_ROOT}"

# ── validate required env ────────────────────────────────────────────
if [[ -z "${BUDGET_EDGE:-}" ]]; then
    echo "[c02c] ERROR: BUDGET_EDGE env var is required (low|central|high)" >&2
    exit 1
fi
case "${BUDGET_EDGE}" in
    low|central|high) ;;
    *)
        echo "[c02c] ERROR: BUDGET_EDGE='${BUDGET_EDGE}' invalid; must be low|central|high" >&2
        exit 1
        ;;
esac

# ── defaults ─────────────────────────────────────────────────────────
N_DRAWS="${N_DRAWS:-1600000}"
M_KK_TEV="${M_KK_TEV:-3,5,7,10,15,20,30,50}"
N_WORKERS="${N_WORKERS:-${SLURM_CPUS_PER_TASK:-32}}"

# ── thread pinning (one thread per worker process) ───────────────────
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# ── conda activation ────────────────────────────────────────────────
source /n/home09/obarrera/.conda/etc/profile.d/conda.sh 2>/dev/null \
    || source "${HOME}/.bashrc"
conda activate ising_bootstrap

# ── output directory ─────────────────────────────────────────────────
mkdir -p "${REPO_ROOT}/logs"
TIMESTAMP="$(date -u +%Y%m%dT%H%M%S)"
OUTPUT_DIR="${REPO_ROOT}/scan_outputs/rs_anarchy_c02c_budget_${BUDGET_EDGE}_${TIMESTAMP}"
mkdir -p "${OUTPUT_DIR}"

# ── diagnostics ──────────────────────────────────────────────────────
echo "[c02c] BUDGET_EDGE = ${BUDGET_EDGE}"
echo "[c02c] OUTPUT_DIR  = ${OUTPUT_DIR}"
echo "[c02c] N_DRAWS     = ${N_DRAWS}"
echo "[c02c] M_KK_TEV    = ${M_KK_TEV}"
echo "[c02c] N_WORKERS   = ${N_WORKERS}"
echo "[c02c] SLURM_JOB_ID= ${SLURM_JOB_ID:-local}"
echo "[c02c] hostname     = $(hostname)"
echo "[c02c] date         = $(date -u +%Y-%m-%dT%H:%M:%SZ)"

# ── run ──────────────────────────────────────────────────────────────
python "${REPO_ROOT}/scripts/run_rs_anarchy.py" \
    --output-dir "${OUTPUT_DIR}" \
    --n-draws "${N_DRAWS}" \
    --m-kk-tev "${M_KK_TEV}" \
    --epsilon-k-budget "${BUDGET_EDGE}" \
    --n-workers "${N_WORKERS}"

echo "[c02c] DONE: ${OUTPUT_DIR}"
