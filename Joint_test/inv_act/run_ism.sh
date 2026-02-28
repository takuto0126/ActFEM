#!/bin/bash


###################  q-m  ########################
#PBS -q q-m
#PBS -l select=2:ncpus=192:mem=740gb:mpiprocs=3
#PBS -N inv_act
#PBS -o out.dat
#PBS -j oe

source $SELECT_PE INTEL
export OMP_NUM_THREADS=64
export I_MPI_PIN_DOMAIN=omp
export KMP_AFFINITY=compact
cd $PBS_O_WORKDIR

#Intel MKLモジュールのロード
module load intelmkl/2023.2.0

src=../../src/src_inv_joint
cd $src
make clean -f Makefile_ISM
make -f Makefile_ISM
cd -

./clean.sh
# 1: only ACTIVE, 2:only MT, 3: both ACTIVE and MT
time mpirun ${src}/ebfem_inv_joint.exe > result_inv/inv.log 2>&1 <<EOF
1
active.ctl
inv_act.ctl
EOF
