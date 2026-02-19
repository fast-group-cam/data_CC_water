#!/bin/bash
#SBATCH --job-name=aims
#SBATCH -p gen
#SBATCH -C icelake
#SBATCH --nodes=1
#SBATCH --tasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=3:00:00
    module purge
    module load modules/2.4-20250724
    module load gcc/13.3.0 openmpi/4.1.8 intel-oneapi-mkl/2024.2.2 hdf5/1.12.3 boost/1.87.0 python/3.11.11 cmake/3.31.6 cuda/12.5.1 nccl/2.23.4-1
    source /mnt/home/bshi1/Programs/venvs/mace-torch-lammps/bin/activate

EXEC="/mnt/home/bshi1/ceph/Programs/cpu/lammps-symmetrix/lammps/build/lmp"
srun --overlap --oversubscribe --ntasks=64 --unbuffered ${EXEC} -in lammps-full.inp -v runID 0
