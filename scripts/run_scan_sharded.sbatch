#!/bin/bash
# Example:
#   sbatch --array=0-15 scripts/run_scan_sharded.sbatch

#SBATCH --job-name=rs-scan
#SBATCH --account=iaifi_lab
#SBATCH --partition=shared
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --output=logs/scan_%A_%a.out
#SBATCH --error=logs/scan_%A_%a.err

set -euo pipefail

ROOT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$ROOT_DIR"

mkdir -p logs scan_outputs

# Use array count when available; fallback to TOTAL_SHARDS env or 1.
# Export so the Python heredoc subprocess can see them via os.environ.
export TOTAL_SHARDS="${SLURM_ARRAY_TASK_COUNT:-${TOTAL_SHARDS:-1}}"
export SHARD_ID="${SLURM_ARRAY_TASK_ID:-0}"

# Optional knobs via environment variables.
export BASE_SEED="${BASE_SEED:-20260205}"
export N_C_L="${N_C_L:-61}"
export N_C_N="${N_C_N:-61}"
export N_MLIGHTEST="${N_MLIGHTEST:-25}"

echo "Running shard ${SHARD_ID}/${TOTAL_SHARDS}"
echo "ROOT_DIR=${ROOT_DIR}"

python - <<'PY'
import os
import numpy as np

from scanParams import ScanConfig, run_scan

total_shards = int(os.environ["TOTAL_SHARDS"])
shard_id = int(os.environ["SHARD_ID"])

base_seed = int(os.environ["BASE_SEED"])
n_c_l = int(os.environ["N_C_L"])
n_c_n = int(os.environ["N_C_N"])
n_mlightest = int(os.environ["N_MLIGHTEST"])

# Global grids (from scanner theory priors).
all_c_L = np.linspace(0.52, 0.70, n_c_l)
all_c_N = np.linspace(0.15, 0.45, n_c_n)
all_lambda_ir = np.array([3000.0, 4000.0, 5000.0, 6000.0, 7000.0], dtype=float)
all_m_lightest = np.logspace(-5, np.log10(3e-2), n_mlightest)

# Shard over c_L so each task gets a deterministic subset.
c_l_shards = np.array_split(all_c_L, total_shards)
c_l_local = c_l_shards[shard_id]

if c_l_local.size == 0:
    print(f"[shard {shard_id}] No c_L points assigned; exiting cleanly.")
    raise SystemExit(0)

config = ScanConfig(
    # Geometry and priors
    k=1.2209e19,
    Lambda_IR_values=all_lambda_ir,
    xi_KK=1.0,
    max_fL_ratio=1.1,
    # Bulk masses
    c_L_values=c_l_local,
    c_N_values=all_c_N,
    c_E_fixed=[0.75, 0.60, 0.50],
    # UV Majorana
    MN_mode="fixed_ratio",
    MN_over_k=0.1,
    # Neutrino mass prior
    lightest_nu_mass_values=all_m_lightest,
    ordering="normal",
    majorana_alpha=0.0,
    majorana_beta=0.0,
    # LFV (MEG II 2025)
    br_limit=1.5e-13,
    prefac_br=4.0e-8,
    lfv_reference_scale=3000.0,
    # Reproducibility
    rng_seed_global=base_seed + shard_id,
    record_git_metadata=True,
)

out_csv = f"scan_outputs/scan_shard_{shard_id:03d}_of_{total_shards:03d}.csv"
rows = run_scan(config, output_csv=out_csv, progress_every=500)

n_total = len(rows)
n_pass = sum(1 for r in rows if r["passes_all"])
print(f"[shard {shard_id}] wrote {out_csv}")
print(f"[shard {shard_id}] accepted {n_pass}/{n_total}")
PY

