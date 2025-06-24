#!/bin/bash
#PBS -q q-s
#PBS -l select=ncpus=48:mem=184gb
#PBS -N TEST_JOB
#PBS -o out.dat
#PBS -j oe

source $SELECT_PE INTEL
export OMP_NUM_THREADS=48
export KMP_AFFINITY=disabled
cd $PBS_O_WORKDIR
dplace 