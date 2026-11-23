#!/bin/bash
#PBS -q q-test
#PBS -l select=ncpus=192:mem=740gb
#PBS -N fwd_act
#PBS -o out.dat
#PBS -j oe

source $SELECT_PE INTEL
export OMP_NUM_THREADS=192
export KMP_AFFINITY=disabled
cd $PBS_O_WORKDIR

module load intelmkl/2023.2.0

src=../../src/solver
cd $src
#make clean -f Makefile_ISM
make -f Makefile_ISM
cd -

./clean.sh
time ${src}/ebfem_bxyz.exe <<EOF > active.log 2>&1
active.ctl
EOF
