#!/bin/bash
# generate forward result of ACTIVE and MT for test joint inversion test
source /opt/intel/oneapi/setvars.sh
## ACTIVE 2 sources ####################################################
src=../../src/solver
cd $src
#make clean
make
cd -

OMP_NUM_THREADS=8
time ${src}/ebfem_bxyz.exe <<EOF #> result_fwd/active.log
active.ctl
EOF
