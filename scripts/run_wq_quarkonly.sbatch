#!/bin/bash
#SBATCH --job-name=wq_quarkonly
#SBATCH --partition=serial_requeue
#SBATCH --account=randall_lab
#SBATCH --cpus-per-task=48
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --output=logs/wq_quarkonly_%j.out
#SBATCH --error=logs/wq_quarkonly_%j.err

# WQ quark-only Bucket-1 mini-scan (dual-approved mode, commit 070cdb6).
#
# Fits the quark sector per draw (anarchic O(1) Yukawas) and evaluates ONLY the
# 46 Bucket-1 constraints that need ZERO swept-lepton data. The lepton sector is
# DROPPED (not rigorous) and tagged as such. Purpose: first look at the RS quark
# sector — M_KK reach vs anarchic Yukawa structure + per-constraint constraining
# power.
#
# NOTE: unlike scripts/run_rs_anarchy.sbatch (forward-only), this harness USES
# fit_quark_sector by design — do NOT add the forward-only optimizer guard here.
#
# 10 M_KK tiles x 10000 draws = 100000 total. ~0.28 s/evaluated point post-cache
# => ~7.7 core-hours => ~10-15 min wall at 48 cores.

set -euo pipefail

REPO_ROOT="/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing"
cd "${REPO_ROOT}"

source "${HOME}/.bashrc"
conda activate ising_bootstrap

TIMESTAMP="$(date -u +%Y%m%dT%H%M%S)"
OUTPUT_DIR="${REPO_ROOT}/scan_outputs/wq_quarkonly_${TIMESTAMP}"
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${REPO_ROOT}/logs"

N_DRAWS="${N_DRAWS:-10000}"
M_KK_TEV="${M_KK_TEV:-1,2,3,5,7,10,15,20,30,50}"
N_WORKERS="${N_WORKERS:-${SLURM_CPUS_PER_TASK:-48}}"

echo "[wq_quarkonly] OUTPUT_DIR = ${OUTPUT_DIR}"
echo "[wq_quarkonly] M_KK (TeV) = ${M_KK_TEV}"
echo "[wq_quarkonly] N_DRAWS/tile = ${N_DRAWS}"
echo "[wq_quarkonly] N_WORKERS  = ${N_WORKERS}"
echo "[wq_quarkonly] git HEAD   = $(git rev-parse HEAD)"

python "${REPO_ROOT}/scripts/run_full_catalog_scan.py" \
    --quark-only \
    --output-dir "${OUTPUT_DIR}" \
    --n-draws "${N_DRAWS}" \
    --m-kk-tev "${M_KK_TEV}" \
    --n-workers "${N_WORKERS}" \
    --smoke-report "${OUTPUT_DIR}/run_report.json" \
    --smoke-report-md "${OUTPUT_DIR}/run_report.md"

echo "[wq_quarkonly] SLURM job complete: ${OUTPUT_DIR}"
