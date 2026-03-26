#!/bin/bash
# Production scan — Stage 1 (c_E grid, M_N/k scan, Λ_IR 3–15 TeV)
#
# Launch:
#   sbatch --array=0-1052 scripts/run_scan_sharded.sbatch
#
# Sharding: flattened (c_L, c_N) grid = 81 × 13 = 1053 pairs.
# Each shard iterates over the remaining axes (c_N_inner, c_E, M_N/k, m_lightest).
# Expected: ~91M total points, ~87k per shard, ~8 min per shard at 180 pts/s.

#SBATCH --job-name=rs-scan
#SBATCH --account=randall_lab
#SBATCH --partition=sapphire
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
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
BASE_SEED="${BASE_SEED:-20260325}"
K_GEV="${K_GEV:-1.2209e19}"
N_C_L="${N_C_L:-81}"
N_C_N="${N_C_N:-51}"
N_LAMBDA="${N_LAMBDA:-13}"
N_MLIGHTEST="${N_MLIGHTEST:-21}"
C_L_MIN="${C_L_MIN:-0.50}"
C_L_MAX="${C_L_MAX:-0.60}"
C_N_MIN="${C_N_MIN:-0.15}"
C_N_MAX="${C_N_MAX:-0.40}"
LAMBDA_MIN="${LAMBDA_MIN:-3000.0}"
LAMBDA_MAX="${LAMBDA_MAX:-15000.0}"
MLIGHTEST_MIN="${MLIGHTEST_MIN:-1e-4}"
MLIGHTEST_MAX="${MLIGHTEST_MAX:-2e-2}"
# M_N/k scan points (space-separated).
MN_OVER_K_VALUES="${MN_OVER_K_VALUES:-0.01 0.05 0.1}"
# c_E grid (space-separated per generation; 3x3x3 = 27 combinations).
C_E1_VALUES="${C_E1_VALUES:-0.72 0.74 0.76}"
C_E2_VALUES="${C_E2_VALUES:-0.56 0.58 0.60}"
C_E3_VALUES="${C_E3_VALUES:-0.40 0.45 0.50}"
export TOTAL_SHARDS SHARD_ID BASE_SEED K_GEV N_C_L N_C_N N_LAMBDA N_MLIGHTEST
export C_L_MIN C_L_MAX C_N_MIN C_N_MAX LAMBDA_MIN LAMBDA_MAX MLIGHTEST_MIN MLIGHTEST_MAX
export MN_OVER_K_VALUES C_E1_VALUES C_E2_VALUES C_E3_VALUES

echo "Running shard ${SHARD_ID}/${TOTAL_SHARDS}"
echo "ROOT_DIR=${ROOT_DIR}"

python - <<'PY'
import os
import numpy as np

from scanParams import ScanConfig, run_scan

total_shards = int(os.environ["TOTAL_SHARDS"])
shard_id = int(os.environ["SHARD_ID"])

base_seed = int(os.environ["BASE_SEED"])
k_gev = float(os.environ["K_GEV"])
n_c_l = int(os.environ["N_C_L"])
n_c_n = int(os.environ["N_C_N"])
n_lambda = int(os.environ["N_LAMBDA"])
n_mlightest = int(os.environ["N_MLIGHTEST"])
c_l_min = float(os.environ["C_L_MIN"])
c_l_max = float(os.environ["C_L_MAX"])
c_n_min = float(os.environ["C_N_MIN"])
c_n_max = float(os.environ["C_N_MAX"])
lambda_min = float(os.environ["LAMBDA_MIN"])
lambda_max = float(os.environ["LAMBDA_MAX"])
mlightest_min = float(os.environ["MLIGHTEST_MIN"])
mlightest_max = float(os.environ["MLIGHTEST_MAX"])
mn_over_k_values = [float(x) for x in os.environ["MN_OVER_K_VALUES"].split()]
c_e1_values = [float(x) for x in os.environ["C_E1_VALUES"].split()]
c_e2_values = [float(x) for x in os.environ["C_E2_VALUES"].split()]
c_e3_values = [float(x) for x in os.environ["C_E3_VALUES"].split()]

if any(n < 1 for n in (n_c_l, n_c_n, n_lambda, n_mlightest)):
    raise ValueError("N_C_L, N_C_N, N_LAMBDA, N_MLIGHTEST must all be >= 1")
if c_l_max < c_l_min:
    raise ValueError("C_L_MAX must be >= C_L_MIN")
if c_n_max < c_n_min:
    raise ValueError("C_N_MAX must be >= C_N_MIN")
if lambda_max < lambda_min:
    raise ValueError("LAMBDA_MAX must be >= LAMBDA_MIN")
if mlightest_min <= 0 or mlightest_max <= 0 or mlightest_max < mlightest_min:
    raise ValueError("Require 0 < MLIGHTEST_MIN <= MLIGHTEST_MAX")

# Global grids.
all_c_L = np.linspace(c_l_min, c_l_max, n_c_l)
all_c_N = np.linspace(c_n_min, c_n_max, n_c_n)
all_lambda_ir = np.linspace(lambda_min, lambda_max, n_lambda)
all_m_lightest = np.logspace(np.log10(mlightest_min), np.log10(mlightest_max), n_mlightest)

# Shard over flattened (c_L, Lambda_IR) pairs for fine-grained parallelism.
# 81 c_L × 13 Lambda_IR = 1053 pairs — one pair per shard at default settings.
pairs = [(i, j) for i in range(len(all_c_L)) for j in range(len(all_lambda_ir))]
shard_pairs = np.array_split(pairs, total_shards)
my_pairs = shard_pairs[shard_id]

if len(my_pairs) == 0:
    print(f"[shard {shard_id}] No (c_L, Lambda_IR) pairs assigned; exiting cleanly.")
    raise SystemExit(0)

# Group by unique c_L and Lambda_IR indices in this shard.
c_l_indices = sorted(set(int(p[0]) for p in my_pairs))
lambda_indices = sorted(set(int(p[1]) for p in my_pairs))
c_l_local = all_c_L[c_l_indices]
lambda_local = all_lambda_ir[lambda_indices]

n_cE = len(c_e1_values) * len(c_e2_values) * len(c_e3_values)
n_pts = (
    len(lambda_local) * len(c_l_local) * len(all_c_N)
    * n_cE * len(mn_over_k_values) * len(all_m_lightest)
)
print(f"[shard {shard_id}] c_L points: {len(c_l_local)}, "
      f"Lambda_IR points: {len(lambda_local)}, "
      f"c_E combos: {n_cE}, M_N/k points: {len(mn_over_k_values)}, "
      f"total: {n_pts}")

config = ScanConfig(
    # Geometry and priors
    k=k_gev,
    Lambda_IR_values=lambda_local,
    xi_KK=1.0,
    max_fL_ratio=1.1,
    # Bulk masses
    c_L_values=c_l_local,
    c_N_values=all_c_N,
    c_E_fixed=None,
    c_E_grid=[
        np.array(c_e1_values),
        np.array(c_e2_values),
        np.array(c_e3_values),
    ],
    # UV Majorana
    MN_mode="scan_ratio",
    MN_over_k_values=np.array(mn_over_k_values),
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

out_csv = f"scan_outputs/scan_shard_{shard_id:04d}_of_{total_shards:04d}.csv"
rows = run_scan(config, output_csv=out_csv, progress_every=5000)

n_total = len(rows)
n_pass = sum(1 for r in rows if r["passes_all"])
print(f"[shard {shard_id}] wrote {out_csv}")
print(f"[shard {shard_id}] accepted {n_pass}/{n_total}")
PY
