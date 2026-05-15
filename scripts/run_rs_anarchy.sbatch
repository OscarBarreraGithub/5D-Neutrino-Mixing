#!/bin/bash
#SBATCH --job-name=rs_anarchy
#SBATCH --partition=serial_requeue
#SBATCH --account=randall_lab
#SBATCH --cpus-per-task=32
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --output=logs/rs_anarchy_%j.out
#SBATCH --error=logs/rs_anarchy_%j.err

# RS flavor anarchy (ACPS) ensemble.
#
# Forward evaluation only (NO optimizer).  Fixed c-value pattern; anarchic
# O(1) Yukawas.  Per draw, all five Delta-F=2 systems are evaluated through
# the repo's tree-level KK-gluon -> Wilson-coefficient pipeline.
#
# DO NOT execute by hand without the orchestrator's go-ahead.  Each tile is
# 100k draws x 8 M_KK values = 800k total evaluations.  Empirically ~2 ms
# per draw on this cluster, so ~25 min total wall time at 32 cores.

set -euo pipefail

REPO_ROOT="/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing"
cd "${REPO_ROOT}"

# ------------------------------------------------------------------
# Conda/Python environment
# ------------------------------------------------------------------
# Match the rest of the repo's SLURM scripts: rely on the user's profile to
# expose conda, then activate the project env.
source "${HOME}/.bashrc"
conda activate ising_bootstrap

# ------------------------------------------------------------------
# Output directory (timestamped)
# ------------------------------------------------------------------
TIMESTAMP="$(date -u +%Y%m%dT%H%M%S)"
OUTPUT_DIR="${REPO_ROOT}/scan_outputs/rs_anarchy_${TIMESTAMP}"
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${REPO_ROOT}/logs"

# ------------------------------------------------------------------
# Run
# ------------------------------------------------------------------
N_DRAWS="${N_DRAWS:-100000}"
M_KK_TEV="${M_KK_TEV:-3,5,7,10,15,20,30,50}"
N_WORKERS="${N_WORKERS:-${SLURM_CPUS_PER_TASK:-32}}"

echo "[rs_anarchy] OUTPUT_DIR  = ${OUTPUT_DIR}"
echo "[rs_anarchy] M_KK (TeV)  = ${M_KK_TEV}"
echo "[rs_anarchy] N_DRAWS     = ${N_DRAWS}"
echo "[rs_anarchy] N_WORKERS   = ${N_WORKERS}"

# Forward-only smoke check (defensive grep — abort if optimizer ever leaks in).
if grep -qE "fit_quark_sector\(|scipy\.optimize|least_squares\(|leastsq\(|minimize\(|evaluate_quark_fit\(" \
    "${REPO_ROOT}/scripts/run_rs_anarchy.py"; then
    echo "[rs_anarchy] ERROR: optimizer reference detected in run_rs_anarchy.py" >&2
    exit 1
fi

python "${REPO_ROOT}/scripts/run_rs_anarchy.py" \
    --output-dir "${OUTPUT_DIR}" \
    --n-draws "${N_DRAWS}" \
    --m-kk-tev "${M_KK_TEV}" \
    --n-workers "${N_WORKERS}"

echo "[rs_anarchy] SLURM job complete: ${OUTPUT_DIR}"
