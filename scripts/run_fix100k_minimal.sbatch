#!/bin/bash
#SBATCH --job-name=fix100k_min
#SBATCH --partition=serial_requeue
#SBATCH --account=randall_lab
#SBATCH --cpus-per-task=48
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --output=logs/fix100k_minimal_%j.out
#SBATCH --error=logs/fix100k_minimal_%j.err

# QUARK-FIX-2026-06-21 verification scan.
#
# 100k MINIMAL (no custodial) FULL-CATALOG scan with all six audit fixes
# (B1 Zbb, B2 epsilon_K rephasing, B3 dF=2 O4/O5, M1 T010 gate, M2 EW001, M5 tag)
# in the working tree. Full-catalog mode (NO --quark-only) so the lepton sector
# is drawn and mu->e gamma (L001) is evaluated alongside the quark/EW IDs.
#
# Constraint set the user asked for (the review_local set): SM masses + CKM
# (per-draw fit gate), epsilon_K (K001), D0 (C001/C002), dm_d (B001), dm_s
# (B003), mu->e gamma (L001), Z->bb (T010/T011), S,T,U (EW001), plus a simple
# collider cutoff m1 >= 5.5 TeV applied post-hoc on M_KK.
#
# 10 M_KK tiles x 10000 draws = 100000 total, single r = 0.25.
# ~0.3-0.45 core-s/evaluated point => ~9-13 core-hours => ~11-16 min wall at 48 cores.

set -euo pipefail

REPO_ROOT="/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing"
cd "${REPO_ROOT}"

source "${HOME}/.bashrc"
conda activate ising_bootstrap

TIMESTAMP="$(date -u +%Y%m%dT%H%M%S)"
OUTPUT_DIR="${REPO_ROOT}/scan_outputs/fix100k_minimal_${TIMESTAMP}"
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${REPO_ROOT}/logs"

N_DRAWS="${N_DRAWS:-10000}"
M_KK_TEV="${M_KK_TEV:-1,2,3,5,7,10,15,20,30,50}"
N_WORKERS="${N_WORKERS:-${SLURM_CPUS_PER_TASK:-48}}"

echo "[fix100k] OUTPUT_DIR = ${OUTPUT_DIR}"
echo "[fix100k] M_KK (TeV) = ${M_KK_TEV}"
echo "[fix100k] N_DRAWS/tile = ${N_DRAWS}"
echo "[fix100k] N_WORKERS  = ${N_WORKERS}"
echo "[fix100k] ew-model   = minimal_rs (NO custodial), full-catalog (lepton on)"
echo "[fix100k] git HEAD   = $(git rev-parse HEAD) (+ uncommitted QUARK-FIX working tree)"

# NOTE: deliberately NO --quark-only, so L001 (mu->e gamma) is evaluated.
python "${REPO_ROOT}/scripts/run_full_catalog_scan.py" \
    --ew-model minimal_rs \
    --output-dir "${OUTPUT_DIR}" \
    --n-draws "${N_DRAWS}" \
    --m-kk-tev "${M_KK_TEV}" \
    --quark-fit-r 0.25 \
    --n-workers "${N_WORKERS}" \
    --smoke-report "${OUTPUT_DIR}/run_report.json" \
    --smoke-report-md "${OUTPUT_DIR}/run_report.md"

echo "[fix100k] SLURM job complete: ${OUTPUT_DIR}"
echo "${OUTPUT_DIR}" > "${REPO_ROOT}/.orchestration/runs/QUARK-FIX-2026-06-21/.last_scan_dir"
