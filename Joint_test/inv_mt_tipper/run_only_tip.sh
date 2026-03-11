#!/bin/bash

source /opt/intel/oneapi/setvars.sh --force

src=../../src/src_inv_joint
cd $src
#make clean
make
cd -

export OMP_NUM_THREADS=4

./clean_tip.sh
# 1: only ACTIVE, 2:only MT, 3: both ACTIVE and MT
time mpiexec -n 6 ${src}/ebfem_inv_joint.exe <<EOF | tee result_inv_only_tip/inv.log # 2025.07.31
2
mt.ctl
inv_mt_only_tip.ctl
EOF
