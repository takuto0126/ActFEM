#!/bin/bash

#source /opt/intel/oneapi/setvars.sh

#PBS -q q-test
#PBS -l select=ncpus=192:mem=740gb
#PBS -N TEST_JOB
#PBS -o out.dat
#PBS -j oe

source $SELECT_PE INTEL
export OMP_NUM_THREADS=192
export KMP_AFFINITY=disabled
cd $PBS_O_WORKDIR

#Intel MKLモジュールのロード
module load intelmkl/2023.2.0

src=../src/src_3DMT
cd $src
make clean -f Makefile_ISM
make -f Makefile_ISM
cd -

echo "Start 3DMT test !"
./clean.sh
time ${src}/n_ebfem_3DMT.exe > result_homo/3DMT.log 2>&1 <<EOF
3DMT_homo.ctl
EOF
