#!/usr/bin/bash
#SBATCH -p gpu
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1
#SBATCH --mem=200G
#SBATCH --constraint=a100  # if you want a particular type of GPU
#SBATCH --time=24:00:00

module load modules/2.4-20250724
module load gcc/13.3.0 
module load openmpi/4.1.8
module load intel-oneapi-mkl/2024.2.2
module load hdf5/1.12.3
module load boost/1.87.0
module load python/3.11.11
module load cmake/3.31.6
module load cuda/12.5.1
module load nccl/2.23.4-1
source /mnt/home/bshi1/Programs/venvs/mace-torch-lammps/bin/activate


# Load the LAMMPSMBX module
#module load gcc openmpi fftw
#export MBX_HOME="/mnt/home/bshi1/Programs/cpu/MBX_LAMMPS"
#export LD_LIBRARY_PATH=$MBX_HOME/lib/:$LD_LIBRARY_PATH
#export LAMMPS_HOME="/mnt/home/bshi1/Programs/cpu/LAMMPS_MBX"
export EXEC="/mnt/home/bshi1/ceph/Programs/gpu/lammps-symmetrix/lammps-symmetrix/lammps/build-gpu/lmp"

# Avoid any unintentional OpenMP threading by setting OMP_NUM_THREADS
export OMP_NUM_THREADS=1 #$SLURM_CPUS_PER_TASK

# Ensure the cpus-per-task option is propagated to srun commands
export SRUN_CPUS_PER_TASK=1 #$SLURM_CPUS_PER_TASK
$EXEC -k on g 1 -sf kk -in lammps-full.inp -v runID 0


#mpirun --bind-to none $EXEC -in lammps-full.inp -v runID 0
