## lines starting with "##" work as lines for comments ! 2020.09.28
## this file is forward control file
##-----10!-------20!
itopofile 0 or 1   !1
# of topofile      !1
topofile           !../../mesh_joint_flat/topo_aso_flat.dat
lon lat shift      !0.0         0.0
mesh file          !../../mesh_joint_flat/nakadake3d.msh
2d triangle z file !../../mesh_joint_flat/nakadake2dz.msh
local line file    !../../mesh_joint_flat/lineinfo.dat
output folder      !./result/
header2d  (a50)    !nakadake2d
header3d  (a50)    !nakadake3d
# of frequency     !5
Frequency [Hz]     !1.d0
Frequency [Hz]     !3.d0
Frequency [Hz]     !11.d0
Frequency [Hz]     !33.d0
Frequency [Hz]     !99.d0
west bound         !-1.7
east bound         !1.7
south bound        !-1.5
north bound        !1.5
lenout [km]        !50.0
upz in [km]  (>0)  !1.3
downz in [km](<0)  !-1.1
zmax   [km]        !50.0
zmin   [km]        !-50.0
sizein [km]        !0.15
sizebo [km]        !10.0
sigma_obs [km]     !0.4
A_obs     [km]     !0.01
dlen_source [km]   !0.1
sigma_src [km]     !0.3
A_src     [km]     !0.005
# of observatory   !4
lonlat(1),xyz (2)  !1
UTM ZONE           |52S
lonlatorigin       !131.084782   32.884882
1  Name            !A02
1  xyz             !131.083411   32.886706  -0.001
2  Name            !A04
2  xyz,sigma,A[km] !131.081939   32.884808  -0.001
3  Name            !A01
3  xyz,sigma,A[km] !131.083367   32.882725  -0.001
4  Name            !A03
4  xyz,sigma,A[km] !131.086847   32.881981  -0.001
ixyflg 0:no,1:surfv!0
# of sources       !2
Source Name        !S1
source start point !131.086490   32.877356   -0.001
source end   point !131.090840   32.877658   -0.001
Source Name        !S2
source start point !131.0784333  32.8908028  -0.001
source end   point !131.0814639  32.8912333  -0.001
Elcetric current[A]! 1.0
sigma_air    [S/m] !1.e-8
condflag 0:home,1: !0
nvolume            !1
cond               !0.01
##condfile           !./cond.msh
