# Coded on 2017.02.22
#!/bin/bash
source /opt/intel/bin/compilervars.sh intel64
source /opt/intel/bin/debuggervars.sh intel64
SRC=../src/src_mesh
cd $SRC
#make clean
make
cd -

#![1]##
${SRC}/gridmdl2msh.exe < cond.ctl
cat nakadake3d.msh cond.msh > tmp.msh
