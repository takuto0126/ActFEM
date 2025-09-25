#!/bin/bash

source /opt/intel/oneapi/setvars.sh --force

src=../../src/src_3DMT
cd $src
make clean -f Makefile_ISM
make       -f Makefile_ISM
cd -

export OMP_NUM_THREADS=8

time ${src}/n_ebfem_3DMT.exe <<EOF
3DMT_homo.ctl
EOF
