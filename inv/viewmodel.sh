# Coded on 2020.12.04
#!/bin/bash

ite=12
folder=result_ms0
mshfile=../mesh_aso_A04/nakadake3d.msh

modelfile=./${folder}/model${ite}.dat
connect=./${folder}/model_connect.dat
outfile=m0.dat

awk '{if (NR == 1) print;else print( log($1) /log(10) )}' $modelfile > model.dat

src=../src/src_post
cd $src
make clean
make -f Makefile_gfort
cd -

${src}/n_model2cond.exe << EOF
$mshfile
$connect
model.dat
0
$outfile
EOF

cat $mshfile $outfile > model.msh

exit

# for masking
#${src}/n_model2cond.exe << EOF
#$mshfile
#$connect
#model0.dat
#1
#dmsens.dat
#0.06
#$outfile
#EOF

