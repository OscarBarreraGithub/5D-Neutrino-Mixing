#!/bin/bash
#SBATCH --job-name=bauer_deband
#SBATCH --partition=serial_requeue
#SBATCH --account=randall_lab
#SBATCH --array=0-2
#SBATCH --cpus-per-task=32
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --output=logs/bauer_deband_%a_%j.out
#SBATCH --error=logs/bauer_deband_%a_%j.err

# De-band the S1 anarchic eps_K cloud (and its m12 / Z->bb companions) by
# regenerating on a FINE 45-tile M_KK grid (1->20 TeV) with high per-tile counts,
# so the mass+CKM-gated (~12% pass) subset becomes dense and continuous instead of
# the discrete 10-tile banding in the prior figures. Same scripts + schemas as the
# files the report renderer already consumes; only the grid + draw counts change.

source ~/.bashrc
conda activate ising_bootstrap
export LD_LIBRARY_PATH=$HOME/.conda/envs/ising_bootstrap/lib:$LD_LIBRARY_PATH
cd /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing

GRID=1.0,1.07,1.15,1.23,1.31,1.41,1.5,1.61,1.72,1.85,1.98,2.11,2.26,2.42,2.59,2.78,2.97,3.18,3.41,3.65,3.9,4.18,4.47,4.79,5.12,5.49,5.87,6.29,6.73,7.2,7.71,8.25,8.83,9.46,10.12,10.84,11.6,12.42,13.29,14.23,15.23,16.31,17.45,18.68,20.0

case "$SLURM_ARRAY_TASK_ID" in
  0)
    # Main grey eps_K cloud + 5/50/95 quantiles (backdrop in solo_epsK_cloud).
    python scripts/anarchic_bauer_s1.py \
      --scenario S1 --per-tile 22000 --n-workers "${SLURM_CPUS_PER_TASK}" \
      --m-kk-tev "$GRID" \
      --out scan_outputs/anarchic_reproduction/anarchic_bauer_S1.parquet
    ;;
  1)
    # Complex-M12 companion (C_Bq-phi, D-funnel, Re/Im(M12) figures).
    python scripts/anarchic_bauer_s1.py \
      --scenario S1 --per-tile 22000 --n-workers "${SLURM_CPUS_PER_TASK}" --emit-m12 \
      --m-kk-tev "$GRID" \
      --out scan_outputs/anarchic_reproduction/anarchic_bauer_s1_m12.parquet
    ;;
  2)
    # Z->bb companion (blue "with Z->bb" / orange "with |eps_K|" overlays).
    python scripts/anarchic_bauer_s1_zbb.py \
      --scenario S1 --per-tile 5000 --workers "${SLURM_CPUS_PER_TASK}" \
      --m-kk-tev "$GRID" \
      --out scan_outputs/anarchic_reproduction/anarchic_bauer_s1_zbb.parquet
    ;;
esac
