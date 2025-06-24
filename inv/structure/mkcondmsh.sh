#!/bin/bash

cat << EOF > cond.msh
\$MeshFormat
2.2 0 8
\$EndMeshFormat
EOF

cat ../mesh_aso/cond.msh >> cond.msh
