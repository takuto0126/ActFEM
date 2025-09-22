#!/bin/bash
cat <<EOF > cond_homo.ctl
!------10!-------20!----
input 3d mshfile   !../mesh_joint/nakadake3d.msh
0homo,1cond,2model !0
homo resistivity   !300.0
output cond        !./cond_homo.msh
0:elevation,1:depth!0
# of cuboid        !0
EOF

SRC=../../src/src_post
cd $SRC
make clean -f Makefile_gfort
make -f Makefile_gfort
cd -

${SRC}/change_model2cond.exe < cond_homo.ctl

rm cond_homo.ctl
