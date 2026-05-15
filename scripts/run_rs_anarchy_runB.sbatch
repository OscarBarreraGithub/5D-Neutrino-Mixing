#!/bin/bash
#SBATCH --job-name=rs_anarchy_runB
#SBATCH --partition=serial_requeue
#SBATCH --account=randall_lab
#SBATCH --requeue
#SBATCH --cpus-per-task=32
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=logs/rs_anarchy_runB_%x_%j.out
#SBATCH --error=logs/rs_anarchy_runB_%x_%j.err

# Run B: Y-prior sensitivity. Env-var dispatcher.
#
# Required:
#   PRIOR_TAG     one of {narrow_uniform, wide_uniform, gaussian_3sigma}
#
# Submission examples (orchestrator):
#   sbatch --export=ALL,PRIOR_TAG=narrow_uniform   scripts/run_rs_anarchy_runB.sbatch
#   sbatch --export=ALL,PRIOR_TAG=wide_uniform     scripts/run_rs_anarchy_runB.sbatch
#   sbatch --export=ALL,PRIOR_TAG=gaussian_3sigma  scripts/run_rs_anarchy_runB.sbatch

set -euo pipefail

REPO_ROOT="/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing"
cd "${REPO_ROOT}"

: "${PRIOR_TAG:?PRIOR_TAG is required (narrow_uniform / wide_uniform / gaussian_3sigma)}"

case "${PRIOR_TAG}" in
    narrow_uniform)
        PRIOR_FLAGS="--y-half-range 1.0 --y-floor 0.1"
        ;;
    wide_uniform)
        PRIOR_FLAGS="--y-half-range 3.0 --y-floor 0.1"
        ;;
    gaussian_3sigma)
        PRIOR_FLAGS="--y-prior gaussian --y-sigma 1.0 --y-trunc-sigma 3.0 --y-floor 0.1"
        ;;
    *)
        echo "[runB] ERROR: unknown PRIOR_TAG=${PRIOR_TAG}" >&2
        exit 2
        ;;
esac

N_DRAWS="${N_DRAWS:-200000}"
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
OUTPUT_DIR="${REPO_ROOT}/scan_outputs/rs_anarchy_runB_${PRIOR_TAG}_${TIMESTAMP}"
mkdir -p "${OUTPUT_DIR}"

echo "[runB] PRIOR_TAG  = ${PRIOR_TAG}"
echo "[runB] PRIOR_FLAGS= ${PRIOR_FLAGS}"
echo "[runB] OUTPUT_DIR = ${OUTPUT_DIR}"
echo "[runB] N_DRAWS    = ${N_DRAWS}"

python "${REPO_ROOT}/scripts/run_rs_anarchy.py" \
    --output-dir "${OUTPUT_DIR}" \
    --n-draws "${N_DRAWS}" \
    --m-kk-tev "${M_KK_TEV}" \
    ${PRIOR_FLAGS} \
    --n-workers "${N_WORKERS}"

echo "[runB] DONE: ${OUTPUT_DIR}"
