!------10!-------20!----
input 3d mshfile   !../../mesh_joint/nakadake3d.msh
0homo,1cond,2model !0
homo resistivity   !100.0
output cond        !./cond_G3.msh
0:elevation,1:depth!0
# of cuboid        !4
C22 xminmax [km]   !-0.150         0.1500
C22 yminmax [km]   !-0.1500        0.1500
C22 zminmax [km]   ! 0.5000        1.05
C22 rho    [Ohm.m] !10.0
C20 xminmax [km]   !-0.0500         0.0500
C20 yminmax [km]   !-0.05000        0.0500
C20 zminmax [km]   ! 0.500         1.05
C20 rho    [Ohm.m] !1.0
C21 xminmax [km]   !-0.200          0.0500
C21 yminmax [km]   !-0.0500         0.0500
C21 zminmax [km]   ! 0.9500         1.05
C21 rho    [Ohm.m] !1.0
C21 xminmax [km]   !-0.2000         -0.100
C21 yminmax [km]   !-0.0500         0.050
C21 zminmax [km]   ! 0.9500         1.15
C21 rho    [Ohm.m] !1.0
