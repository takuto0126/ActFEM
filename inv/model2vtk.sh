#!/bin/bash

source /opt/intel/bin/compilervars.sh intel64
source /opt/intel/bin/debuggervars.sh intel64

src=../src/src_post
cd $src
#make clean
make
cd -

mshfile="../mesh_aso_A04/nakadake3d.msh"
connect="./result_ms0/model_connect.dat"
modelfile1="result_ms0/model12.dat"
modelfile2="result_m0/model13.dat"
#outfile="model_201405_201408.vtk"
outfile="model_201405_201408_mask.vtk"
sensmodelfile="m0sens.dat"

criterion=0.2 #2023.01.02

${src}/model2vtk.exe <<EOF
$mshfile
$connect
$modelfile1
$modelfile2
resistivity_201405Ohm.m
resistivity_201408Ohm.m
log10_res_201405-o-res_201408
1    # iunit
$outfile
1    # imask 2023.01.02#  0:no mask, 1:mask, followed by sensmodel, criterion,three values for masked part
$sensmodelfile
$criterion
10000000000.
10000000000.
10.0
EOF
