# generate model.msh from model**.dat
#!/bin/bash

source /opt/intel/bin/compilervars.sh intel64
source /opt/intel/bin/debuggervars.sh intel64

src=../src/src_post
cd $src
make clean
make
cd -

mshfile="../mesh_aso_A04/nakadake3d.msh"
connect="./result_ms0/model_connect.dat"
modelfile="result_ms0/model12.dat"
outfile="model.dat"

${src}/n_model2cond.exe << EOF
$mshfile
$connect
$modelfile
0
$outfile
EOF

cat $mshfile $outfile > model.msh

