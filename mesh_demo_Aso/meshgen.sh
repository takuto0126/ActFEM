# Coded on 2017.02.21
#!/bin/bash
source /opt/intel/oneapi/setvars.sh 
SRC=../src/src_mesh
cd $SRC
#make clean
make -f Makefile_gfort
#make 
cd -

#ctlfile="aso.ctl"       # UTM map projection
#ctlfile=aso_ECP.ctl    # ECP (Equidistant Cylindrical Projection)
#ctlfile=aso_rotate.ctl
ctlfile=aso_ECP_rotate.ctl

#![1]##
${SRC}/meshgen1.exe <<EOF
$ctlfile
EOF

#![2]##
gmsh -2 nakadake2d.geo -bgm nakadake2d.pos -format msh2

#![3]##
${SRC}/meshgen2.exe <<EOF
$ctlfile
EOF

#![4]##
gmsh -3 nakadake3d.geo -bgm nakadake3d.pos -format msh2 > nakadake3d.log

#![5]## mkline
${SRC}/mkline.exe <<EOF
$ctlfile
EOF

exit

#![6]## mkface
${SRC}/mkface.exe <<EOF
$ctlfile
EOF
