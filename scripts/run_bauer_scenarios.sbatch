#!/bin/bash
#SBATCH --job-name=bauer_anarch
#SBATCH --partition=serial_requeue
#SBATCH --account=randall_lab
#SBATCH --array=0-3
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/bauer_anarch_%a_%j.out
#SBATCH --error=logs/bauer_anarch_%a_%j.err

source ~/.bashrc
conda activate ising_bootstrap
export LD_LIBRARY_PATH=$HOME/.conda/envs/ising_bootstrap/lib:$LD_LIBRARY_PATH
cd /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing

SCEN=(S1 S2 S3 S4)
SC=${SCEN[$SLURM_ARRAY_TASK_ID]}
PT=20000
if [ "$SC" != "S1" ]; then PT=15000; fi

python scripts/anarchic_bauer_s1.py \
  --scenario "$SC" --per-tile "$PT" --n-workers "${SLURM_CPUS_PER_TASK}" \
  --m-kk-tev 1,1.5,2,2.5,3,3.5,4,5,6,7,8,9,10 \
  --out scan_outputs/anarchic_reproduction/anarchic_bauer_${SC}.parquet
