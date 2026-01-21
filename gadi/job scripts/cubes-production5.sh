#!/bin/bash

#PBS -q normal
#PBS -P qk9
#PBS -l ncpus=1680
#PBS -l mem=500GB
#PBS -l jobfs=10GB
#PBS -l walltime=10:00:00
#PBS -l storage=scratch/qk9
#PBS -l wd

module load openmpi
mpirun code/parmain
