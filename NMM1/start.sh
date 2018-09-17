#!/bin/sh
#SBATCH -p short
#SBATCH -n12
#SBATCH -C alpha
mpirun.openmpi python3 NMM1_mpi.py