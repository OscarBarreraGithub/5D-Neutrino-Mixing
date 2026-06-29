#!/bin/bash
#SBATCH --job-name=zbb_dense
#SBATCH --partition=serial_requeue
#SBATCH --account=randall_lab
#SBATCH --cpus-per-task=48
#SBATCH --mem=24G
#SBATCH --time=03:30:00
#SBATCH --output=logs/zbb_deband_%j.out
#SBATCH --error=logs/zbb_deband_%j.err

# Fast single-wave de-band of the Z->bb overlay: 45 tiles on 48 cores (1 tile/worker,
# no second wave), 1500 draws/tile (~165 gated/tile at 11% -> ~7400 overlay points,
# plenty). Heavy RS-EW spectrum per draw (~1.4 s) => ~35 min single wave.
source ~/.bashrc
conda activate ising_bootstrap
export LD_LIBRARY_PATH=$HOME/.conda/envs/ising_bootstrap/lib:$LD_LIBRARY_PATH
cd /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing

GRID=1.0,1.07,1.15,1.23,1.31,1.41,1.5,1.61,1.72,1.85,1.98,2.11,2.26,2.42,2.59,2.78,2.97,3.18,3.41,3.65,3.9,4.18,4.47,4.79,5.12,5.49,5.87,6.29,6.73,7.2,7.71,8.25,8.83,9.46,10.12,10.84,11.6,12.42,13.29,14.23,15.23,16.31,17.45,18.68,20.0

python scripts/anarchic_bauer_s1_zbb.py \
  --scenario S1 --per-tile 4500 --workers "${SLURM_CPUS_PER_TASK}" \
  --m-kk-tev "$GRID" \
  --out scan_outputs/anarchic_reproduction/anarchic_bauer_s1_zbb.parquet
