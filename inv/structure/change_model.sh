# Takuto Minami Feb 5, 2021
#!/bin/bash
source /opt/intel/bin/compilervars.sh intel64

src=../src/src_post
cd $src
make clean
make -f Makefile_gfort
cd -

${src}/change_model.exe <<EOF
!------10!-------20!----
input 3d mshfile   !../mesh_aso_A04/nakadake3d.msh
model connect file !./result_ms0/model_connect.dat
input model        !./result_ms0/model12.dat
output model       !./model12_v2.dat
# of cuboid        !2
1  xminmax [km]    !-0.500         -0.300
1  yminmax [km]    ! -0.2           0.200
1  zminmax [km]    !  0.5           0.75
1  rho    [Ohm.m]  !100.0
2  xminmax [km]    ! 0.500          0.80
2  yminmax [km]    !-0.100          0.45
2  zminmax [km]    ! 0.000          0.50
2  rho    [Ohm.m]  !  1.0
EOF
