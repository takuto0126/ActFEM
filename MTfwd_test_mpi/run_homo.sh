#!/bin/bash

source /opt/intel/oneapi/setvars.sh --force

src=../src/src_3DMT_mpi
cd $src
make clean
make
cd -

export OMP_NUM_THREADS=6

time mpiexec -n 4 ${src}/n_ebfem_3DMT.exe <<EOF
3DMT_homo.ctl
EOF
