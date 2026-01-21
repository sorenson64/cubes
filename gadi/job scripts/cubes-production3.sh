#!/bin/bash

#PBS -q normal
#PBS -P qk9
#PBS -l ncpus=3840
#PBS -l mem=2000GB
#PBS -l jobfs=10GB
#PBS -l walltime=05:00:00
#PBS -l storage=scratch/qk9
#PBS -l wd

module load openmpi
mpirun code/parmain
