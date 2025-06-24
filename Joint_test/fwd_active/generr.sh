#!/bin/bash


cat <<EOF > err.ctl
        10       20| header
# of frequencies   |5
input folder       |./result/
output folder      |./data/
# of observatories |4
obsname 1          |A02
obsname 2          |A03
obsname 3          |A04
obsname 4          |A01
# of sources       |2
source name 1      |S1
source name 2      |S2
error/amp ratio    |0.01  
EOF

cd src
gfortran generr_active.f90 -o generr_active.exe
cd -

./src/generr_active.exe < err.ctl
