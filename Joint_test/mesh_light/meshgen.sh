# Coded on 2017.02.21
#!/bin/bash
#source /opt/intel/bin/compilervars.sh intel64

SRC=../../src/src_mesh
cd $SRC
make clean -f Makefile_gfort
make -f Makefile_gfort
cd -

#![1]##
${SRC}/meshgen1.exe <<EOF
aso.ctl
EOF

#![2]##
gmsh -2 -string "Mesh.RandomSeed=1;" nakadake2d.geo -bgm nakadake2d.pos -format msh2 # random seed is fixed 2026.11.25 CC

#![3]##
${SRC}/meshgen2.exe << EOF
aso.ctl
EOF

#![4]##
gmsh -3 -string "Mesh.RandomSeed=1;" nakadake3d.geo -bgm nakadake3d.pos -format msh2 > nakadake3d.log # random seed is fixed 2026.11.25 CC

#![5]## mkline
${SRC}/mkline.exe <<EOF
aso.ctl
EOF

#exit # commented out 2026.11.25 CC (to generate faceinfo.dat by mkface.exe)

#![6]## mkface
${SRC}/mkface.exe <<EOF
aso.ctl
EOF
