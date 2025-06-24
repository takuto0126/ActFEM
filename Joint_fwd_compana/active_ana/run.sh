#!/bin/bash

source /opt/intel/bin/compilervars.sh intel64

SRC=../../src/ana

cd ${SRC}
#make clean
make
cd -

time ${SRC}/active_ana.exe <<EOF
ana.ctl
EOF


