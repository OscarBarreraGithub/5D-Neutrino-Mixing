#!/bin/bash
#SBATCH --job-name=rs_anarchy_runA
#SBATCH --partition=serial_requeue
#SBATCH --account=randall_lab
#SBATCH --requeue
#SBATCH --cpus-per-task=32
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --output=logs/rs_anarchy_runA_%j.out
#SBATCH --error=logs/rs_anarchy_runA_%j.err

# Run A: high-statistics baseline (1M draws/tile, 8 tiles).
#
# Halt trigger: at the M_KK=3 TeV tile we expect n_pdg_pass = 163,600 +/- 2000
# (binomial 5 sigma around the baseline pdg_pass_fraction = 0.1636). The
# orchestrator should grep tile_summary.json after completion to verify.

set -euo pipefail

REPO_ROOT="/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing"
cd "${REPO_ROOT}"

N_DRAWS="${N_DRAWS:-1000000}"
M_KK_TEV="${M_KK_TEV:-3,5,7,10,15,20,30,50}"
N_WORKERS="${N_WORKERS:-${SLURM_CPUS_PER_TASK:-32}}"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

source /n/home09/obarrera/.conda/etc/profile.d/conda.sh 2>/dev/null \
    || source "${HOME}/.bashrc"
conda activate ising_bootstrap

mkdir -p "${REPO_ROOT}/logs"
TIMESTAMP="$(date -u +%Y%m%dT%H%M%S)"
OUTPUT_DIR="${REPO_ROOT}/scan_outputs/rs_anarchy_runA_${TIMESTAMP}"
mkdir -p "${OUTPUT_DIR}"

echo "[runA] OUTPUT_DIR = ${OUTPUT_DIR}"
echo "[runA] N_DRAWS    = ${N_DRAWS}"
echo "[runA] M_KK_TEV   = ${M_KK_TEV}"
echo "[runA] N_WORKERS  = ${N_WORKERS}"

python "${REPO_ROOT}/scripts/run_rs_anarchy.py" \
    --output-dir "${OUTPUT_DIR}" \
    --n-draws "${N_DRAWS}" \
    --m-kk-tev "${M_KK_TEV}" \
    --n-workers "${N_WORKERS}"

echo "[runA] DONE: ${OUTPUT_DIR}"
