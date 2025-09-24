#!/bin/bash
source /opt/intel/oneapi/setvars.sh --force

## run forward
src=../../src/solver
cd $src
#make clean
make -f Makefile
cd -

OMP_NUM_THREADS=12
time ${src}/ebfem_bxyz.exe <<EOF #> result/fwd.log
active.ctl
EOF

#time mpiexec_mpt -np  dplace -s1 ./ebfem_bxyz.exe
