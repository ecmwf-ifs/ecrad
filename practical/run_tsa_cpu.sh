#!/bin/bash -l
#SBATCH --job-name="ecrad_cpu"
#SBATCH --account="s83"
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=debug
#SBATCH --constraint=gpu
#SBATCH --exclusive

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export CRAY_CUDA_MPS=1

srun ../bin/ecrad config.nam era5slice.nc tt.nc
