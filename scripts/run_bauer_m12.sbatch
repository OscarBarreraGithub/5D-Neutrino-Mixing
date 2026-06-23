#!/bin/bash
#SBATCH --job-name=bauer_m12
#SBATCH --partition=serial_requeue
#SBATCH --account=randall_lab
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/bauer_m12_%j.out
#SBATCH --error=logs/bauer_m12_%j.err

source ~/.bashrc
conda activate ising_bootstrap
export LD_LIBRARY_PATH=$HOME/.conda/envs/ising_bootstrap/lib:$LD_LIBRARY_PATH
cd /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing

# S1-matched ensemble WITH complex M12 emission, for the C_Bq-phi / D-funnel /
# Re-Im(M12) figures.  Tiles mirror the existing anarchic_complex_m12 set.
python scripts/anarchic_bauer_s1.py \
  --scenario S1 --per-tile 25000 --n-workers "${SLURM_CPUS_PER_TASK}" --emit-m12 \
  --m-kk-tev 1,1.5,2,3,5,7,10 \
  --out scan_outputs/anarchic_reproduction/anarchic_bauer_s1_m12.parquet
