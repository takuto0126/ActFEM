#!/bin/bash

source /opt/intel/oneapi/setvars.sh --force

src=../../src/src_inv_joint
cd $src
make clean
make
cd -

export OMP_NUM_THREADS=4
./clean.sh

time mpiexec -n 6 ${src}/ebfem_inv_joint.exe <<EOF | tee result_inv/joint.log # 2025.07.31
1
active.ctl
inv_act.ctl
EOF
