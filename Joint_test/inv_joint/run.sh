#!/bin/bash

source /opt/intel/oneapi/setvars.sh

src=../../src/src_inv_joint
cd $src
#make clean
make
cd -

export OMP_NUM_THREADS=8

time mpiexec -n 6 ${src}/ebfem_inv_joint.exe <<EOF
3
active_fwd.ctl
mt_fwd.ctl
joint.ctl
EOF
