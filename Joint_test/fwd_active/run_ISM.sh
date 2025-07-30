#!/bin/bash
# generate forward result of ACTIVE and MT for test joint inversion test
#source /opt/intel/oneapi/setvars.sh
## ACTIVE 2 sources ####################################################

#PBS -q q-m
#PBS -l select=1:ncpus=192:mem=740gb:mpiprocs=1
#PBS -N TEST_JOB
#PBS -o out.dat
#PBS -j oe

source $SELECT_PE INTEL
module load intelmkl/2023.2.0
cd $PBS_O_WORKDIR

src=../../src/solver
cd $src
make clean -f Makefile_ISM
make -f Makefile_ISM
cd -

OMP_NUM_THREADS=192
time ${src}/ebfem_bxyz.exe <<EOF > active.log 2>&1
active.ctl
EOF
