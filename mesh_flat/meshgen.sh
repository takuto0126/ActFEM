# Coded on 2017.02.21
#!/bin/bash
source /opt/intel/oneapi/setvars.sh intel64

SRC=../src/src_mesh
cd $SRC
#make clean
make -f Makefile_gfort
cd -

#ctlfile=aso.ctl
ctlfile=aso_ECP.ctl
#ctlfile=aso_rotate.ctl
#ctlfile=aso_ECP_rotate.ctl

#![1]##
${SRC}/meshgen1.exe <<EOF
$ctlfile
EOF

#![2]##
gmsh -2 nakadake2d.geo -bgm nakadake2d.pos -format msh2

#![3]##
${SRC}/meshgen2.exe << EOF
$ctlfile
EOF

#![4]##
#gmsh -3 nakadake3d.geo -bgm nakadake3d.pos -format msh2 > nakadake3d.log
gmsh -3 nakadake3d.geo  -format msh2 > nakadake3d.log

#![5]## mkline
${SRC}/mkline.exe <<EOF
$ctlfile
EOF

#![6]## mkface
#${SRC}/mkface.exe <<EOF
#$ctlfile
#EOF
