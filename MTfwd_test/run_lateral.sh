#!/bin/bash

source /opt/intel/oneapi/setvars.sh --force

src=../src/src_3DMT
cd $src
make clean
make
cd -

export OMP_NUM_THREADS=8

time ${src}/n_ebfem_3DMT.exe <<EOF
3DMT_lateral.ctl
EOF
