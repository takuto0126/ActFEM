# Coded on 2017.02.22
#!/bin/bash
#source /opt/intel/bin/compilervars.sh intel64
#source /opt/intel/bin/debuggervars.sh intel64
SRC=../src/src_mesh
cd $SRC
make clean
make
cd -

#![1]##
${SRC}/gridmdl2msh.exe < cond_init.ctl
cat ../mesh_aso_A04/nakadake3d.msh cond_init.msh > init.msh

#[homo]##
${SRC}/gridmdl2msh.exe < cond_homo.ctl
#cat ../mesh_aso/nakadake3d.msh cond_homo.msh > homo.msh

#[kanda]##
#${SRC}/gridmdl2msh.exe < cond_kanda.ctl
#cat ../mesh_aso/nakadake3d.msh cond_kanda.msh > kanda.msh
