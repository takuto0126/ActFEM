# Coded on 2017.02.21
#!/bin/bash
source /opt/intel/bin/compilervars.sh intel64
source /opt/intel/bin/debuggervars.sh intel64
SRC=../src/src_mesh
cd $SRC
#make clean
make
cd -

ctlfile="aso_add.ctl"

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

exit

#![5]## mkline
${SRC}/mkline.exe < aso.ctl

#![6]## mkface
${SRC}/mkface.exe < aso.ctl
