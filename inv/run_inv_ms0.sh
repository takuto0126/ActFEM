#!/bin/bash

source /opt/intel/oneapi/setvars.sh --force

SRC=../src/src_inv_ap
cd ${SRC}
make clean
make
cd -

#cat tmp.ctl
OMP_NUM_THREADS=7
mpiexec -n 2 ${SRC}/ebfem_inv_ap.exe << EOF &> ./result_ms0/ms0.log
aso_inv.ctl
inv201405_ap.ctl
EOF
