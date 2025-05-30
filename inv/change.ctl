!------10!-------20!----
input 3d mshfile   !../mesh_aso_A04/nakadake3d.msh
model connect file !./result_ms0/model_connect.dat
input model        !./result_ms0/model12.dat
output model       !./model12_v2.dat
# of cuboid        !2
1  xminmax [km]    !-0.500         1.000
1  yminmax [km]    ! 0.45          1.000
1  zminmax [km]    !-0.500          1.200
1  rho    [Ohm.m]  !100.0
2   xminmax [km]   !-0.500          0.000
2   yminmax [km]   !-0.200          0.45
2   zminmax [km]   ! 0.000          1.120
2   rho    [Ohm.m] !100.0
