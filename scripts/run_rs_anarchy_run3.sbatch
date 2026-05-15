#!/bin/bash
#SBATCH --job-name=rs_anarchy_run3
#SBATCH --partition=serial_requeue
#SBATCH --account=randall_lab
#SBATCH --requeue
#SBATCH --cpus-per-task=32
#SBATCH --mem=24G
#SBATCH --time=01:30:00
#SBATCH --output=logs/rs_anarchy_run3_%x_%j.out
#SBATCH --error=logs/rs_anarchy_run3_%x_%j.err

# Run 3 dispatcher: c-value sensitivity sweep.
#
# Reads pattern from env vars:
#   PATTERN_TAG    e.g. "baseline", "qtop_shifted", "moreUV", "moreIR"
#   C_Q            comma-separated triple, e.g. "0.63,0.57,0.20"
#   C_U            comma-separated triple
#   C_D            comma-separated triple
#
# Optional:
#   N_DRAWS        default 200000
#   M_KK_TEV       default "3,5,7,10,15,20,30,50"
#
# Submission examples (orchestrator):
#   sbatch --export=ALL,PATTERN_TAG=baseline,C_Q="0.63,0.57,0.20",\
#       C_U="0.66,0.50,-0.50",C_D="0.66,0.61,0.55" \
#       scripts/run_rs_anarchy_run3.sbatch

set -euo pipefail

REPO_ROOT="/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing"
cd "${REPO_ROOT}"

# Required env vars
: "${PATTERN_TAG:?PATTERN_TAG is required (e.g. baseline / qtop_shifted / moreUV / moreIR)}"
: "${C_Q:?C_Q is required, e.g. \"0.63_0.57_0.20\"}"
: "${C_U:?C_U is required, e.g. \"0.66_0.50_-0.50\"}"
: "${C_D:?C_D is required, e.g. \"0.66_0.61_0.55\"}"

# sbatch --export parses commas as env-var separators, so we accept c-values
# with underscores and convert them to commas here.
C_Q="${C_Q//_/,}"
C_U="${C_U//_/,}"
C_D="${C_D//_/,}"

# Optional
N_DRAWS="${N_DRAWS:-200000}"
M_KK_TEV="${M_KK_TEV:-3,5,7,10,15,20,30,50}"
N_WORKERS="${N_WORKERS:-${SLURM_CPUS_PER_TASK:-32}}"

# Performance hygiene: numpy/scipy threadpools off (we parallelize at the
# Python level via mp.Pool).
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Conda env
source /n/home09/obarrera/.conda/etc/profile.d/conda.sh 2>/dev/null \
    || source "${HOME}/.bashrc"
conda activate ising_bootstrap

mkdir -p "${REPO_ROOT}/logs"
TIMESTAMP="$(date -u +%Y%m%dT%H%M%S)"
OUTPUT_DIR="${REPO_ROOT}/scan_outputs/rs_anarchy_run3_${PATTERN_TAG}_${TIMESTAMP}"
mkdir -p "${OUTPUT_DIR}"

echo "[run3] PATTERN_TAG = ${PATTERN_TAG}"
echo "[run3] C_Q         = ${C_Q}"
echo "[run3] C_U         = ${C_U}"
echo "[run3] C_D         = ${C_D}"
echo "[run3] N_DRAWS     = ${N_DRAWS}"
echo "[run3] M_KK_TEV    = ${M_KK_TEV}"
echo "[run3] OUTPUT_DIR  = ${OUTPUT_DIR}"

python "${REPO_ROOT}/scripts/run_rs_anarchy.py" \
    --output-dir "${OUTPUT_DIR}" \
    --n-draws "${N_DRAWS}" \
    --m-kk-tev "${M_KK_TEV}" \
    --c-Q "${C_Q}" \
    --c-u "${C_U}" \
    --c-d "${C_D}" \
    --n-workers "${N_WORKERS}"

echo "[run3] DONE: ${OUTPUT_DIR}"
