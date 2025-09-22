#!/bin/bash

source /opt/intel/oneapi/setvars.sh

SRC=../../src/ana

cd ${SRC}
make clean
make
cd -

time ${SRC}/active_ana.exe <<EOF
ana.ctl
EOF


