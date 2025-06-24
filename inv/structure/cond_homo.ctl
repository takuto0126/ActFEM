!------10!-------20!----
input 3d mshfile   !../mesh_aso_A04/nakadake3d.msh
output cond.msh    !./cond_homo.msh
rho_n=[ratio]*rho  !1.0
# of cuboid        !1
C1 xminmax [km]    !-500.0          500.0
C1 yminmax [km]    !-500.0          500.0
C1 zminmax [km]    !-500.0          2.000
C1 rho    [Ohm.m]  !100.0
C1 ratio           !1.0
