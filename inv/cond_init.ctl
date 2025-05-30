!------10!-------20!----
input 3d mshfile   !../mesh_aso_A04/nakadake3d.msh
output cond.msh    !./cond_init.msh
rho_n=[ratio]*rho  !1.0
# of cuboid        !10
C1 xminmax [km]    !-0.500         1.000
C1 yminmax [km]    ! 0.45          1.000
C1 zminmax [km]    !-0.500          1.200
C1 rho    [Ohm.m]  !100.0
C1 ratio           !0.5
C21 xminmax [km]   !-0.500          0.000
C21 yminmax [km]   !-0.200          0.45
C21 zminmax [km]   ! 0.000          1.120
C21 rho    [Ohm.m] !100.0
C21 ratio          !0.0
C22 xminmax [km]   ! 0.000          0.400
C22 yminmax [km]   !-0.200          0.400
C22 zminmax [km]   ! 0.000          1.120
C22 rho    [Ohm.m] !3000.0
C22 ratio          !0.5
C23 xminmax [km]   ! 0.000          0.400
C23 yminmax [km]   !-0.60          -0.200
C23 zminmax [km]   ! 0.000          1.12
C23 rho    [Ohm.m] !3000.0
C23 ratio          !0.5
C24 xminmax [km]   !-0.50           0.0000
C24 yminmax [km]   !-10.00         -0.200
C24 zminmax [km]   ! 0.000          1.12
C24 rho    [Ohm.m] !3000.0
C24 ratio          !0.5
C3 xminmax [km]    !-0.100          0.300
C3 yminmax [km]    !-0.600         -0.200
C3 zminmax [km]    ! 1.1200         1.250
C3 rho    [Ohm.m]  !100.0
C3 ratio           !0.5
C4 xminmax [km]    !-0.400          0.300
C4 yminmax [km]    !-0.200          0.200
C4 zminmax [km]    ! 1.1200         1.200
C4 rho    [Ohm.m]  !100.0
C4 ratio           !0.5
C5 xminmax [km]    !-0.400          0.300
C5 yminmax [km]    ! 0.200          0.400
C5 zminmax [km]    ! 1.1200         1.250
C5 rho    [Ohm.m]  !100.0
C5 ratio           !0.5
C6 xminmax [km]    ! 0.000         10.000
C6 yminmax [km]    !-10.0           0.30
C6 zminmax [km]    ! 0.000          1.25
C6 rho    [Ohm.m]  !3000.0
C6 ratio           !0.5
C7 xminmax [km]    !-10.00         -0.5000
C7 yminmax [km]    !-0.6           0.700
C7 zminmax [km]    ! 0.000          1.25
C7 rho    [Ohm.m]  !3000.0
C7 ratio           !0.5





