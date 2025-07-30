#!/bin/bash

#source /opt/intel/oneapi/setvars.sh

#PBS -q q-m
#PBS -l select=1:ncpus=192:mem=740gb:mpiprocs=6
#PBS -N TEST_JOB
#PBS -o out.dat
#PBS -j oe

source $SELECT_PE INTEL
export OMP_NUM_THREADS=32
export I_MPI_PIN_DOMAIN=omp
export KMP_AFFINITY=compact
cd $PBS_O_WORKDIR

#Intel MKLモジュールのロード
module load intelmkl/2023.2.0

src=../../src/src_inv_joint
cd $src
#make clean
make
cd -

echo "Start n_ebfem_joint !"
time mpirun ${src}/ebfem_inv_joint.exe <<EOF
3
active.ctl
mt.ctl
joint.ctl
EOF
