#!/bin/bash

#SBATCH --job-name=
#SBATCH --account=campphyre
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=256G
#SBATCH --time=2:00:00 # Max runtime in DD-HH:MM:SS format.

#SBATCH --export=all
#SBATCH --output=%x_%j.out # where STDOUT goes
#SBATCH --error=%x_%j.err  # where STDERR goes

# Modules to use (optional).
source /mmfs1/home/zzh0903/.bash_profile

# init_flow input
srun ./nga.dp.gnu.opt.mpi.exe -i input

# We are Done
exit

