#!/bin/bash

#source /opt/intel/oneapi/setvars.sh

#PBS -q q-test
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

src=../../src/src_3DMT_mpi
cd $src
make clean -f Makefile_ISM
make -f Makefile_ISM
cd -

echo "Start 3DMT_mpi !"
./clean.sh
time mpirun ${src}/n_ebfem_3DMT.exe > result_homo/3DMT.log 2>&1 <<EOF
3DMT_homo.ctl
EOF
