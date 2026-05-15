#!/bin/bash
#SBATCH --job-name=rs_anarchy_runC
#SBATCH --partition=serial_requeue
#SBATCH --account=randall_lab
#SBATCH --requeue
#SBATCH --cpus-per-task=32
#SBATCH --mem=24G
#SBATCH --time=01:30:00
#SBATCH --output=logs/rs_anarchy_runC_%j.out
#SBATCH --error=logs/rs_anarchy_runC_%j.err

# Run C: tighter PDG gate + tighter Y-floor (CFW-like setup).
#   - baseline c-pattern
#   - --y-half-range 3.0 --y-floor 0.5  (mimics |m~| >= 0.5 in CFW)
#   - --pdg-mass-factor 1.5 --pdg-ckm-factor 1.5 --pdg-j-factor 2.5
#   - 500k draws/tile, 8 tiles
#
# We CANNOT fully reproduce CFW's exact c+m~ scan (the original PDF is
# inaccessible). The plotter overlays the CFW abstract benchmarks at
# 21 TeV (generic anarchic) and 33 TeV (PGB Higgs) for context.

set -euo pipefail

REPO_ROOT="/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing"
cd "${REPO_ROOT}"

N_DRAWS="${N_DRAWS:-500000}"
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
OUTPUT_DIR="${REPO_ROOT}/scan_outputs/rs_anarchy_runC_${TIMESTAMP}"
mkdir -p "${OUTPUT_DIR}"

echo "[runC] OUTPUT_DIR = ${OUTPUT_DIR}"
echo "[runC] N_DRAWS    = ${N_DRAWS}"
echo "[runC] M_KK_TEV   = ${M_KK_TEV}"

python "${REPO_ROOT}/scripts/run_rs_anarchy.py" \
    --output-dir "${OUTPUT_DIR}" \
    --n-draws "${N_DRAWS}" \
    --m-kk-tev "${M_KK_TEV}" \
    --y-half-range 3.0 \
    --y-floor 0.5 \
    --pdg-mass-factor 1.5 \
    --pdg-ckm-factor 1.5 \
    --pdg-j-factor 2.5 \
    --n-workers "${N_WORKERS}"

echo "[runC] DONE: ${OUTPUT_DIR}"
