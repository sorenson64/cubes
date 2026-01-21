#!/bin/bash

#PBS -q normal
#PBS -P qk9
#PBS -l ncpus=480
#PBS -l mem=400GB
#PBS -l jobfs=10GB
#PBS -l walltime=20:00:00
#PBS -l storage=scratch/qk9
#PBS -l wd

module load openmpi
mpirun code/parmain
