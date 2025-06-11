#!/bin/bash

cat ./topo127_134_29_36.xyz | awk -F" " '{print($1,$2,0.0)}' > topo_aso_flat.dat
