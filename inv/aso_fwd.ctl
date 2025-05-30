!------10!-------20!
0: no, 1:topofile  !1
# of topo file     !2
topofile           !./aso.dat
lon lat shift      !0.0            0.0
topofile           !../topo/topo127_134_29_36.xyz
mesh file          !../mesh_aso_A04/nakadake3d.msh
2d triangle z file !../mesh_aso_A04/nakadake2dz.msh
local line file    !../mesh_aso_A04/lineinfo.dat
output folder      !./result_fwd/
header2d  (a50)    !nakadake2d
header3d  (a50)    !nakadake3d
# of frequency     !3
Frequency [Hz]     !1.d0
Frequency [Hz]     !11.d0
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
A_obs     [km]     !0.008
dlen_source [km]   !0.3
sigma_src [km]     !0.32
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
ixyflag 0 or 1     !0
# of sources       !1
Source Name        !S1
source start point !131.0784333  32.8908028  -0.001
source end   point !131.0814639  32.8912333  -0.001
Elcetric current[A]! 1.0
sigma_air    [S/m] !1.e-8
condflag 0:home,1: !1
condfile           !./cond_init.msh
!------10!-------20!
source start point !131.0784333  32.8908028  1.1591617
source end   point !131.0814639  32.8912333  1.1848705
!
0.0 -1.00 0.0
0.0 -0.95 0.0
0.0 -0.90 0.0
0.0 -0.85 0.0
0.0 -0.80 0.0
0.0 -0.75 0.0
0.0 -0.70 0.0
0.0 -0.65 0.0
0.0 -0.60 0.0
0.0 -0.55 0.0
0.0 -0.50 0.0
0.0 -0.45 0.0
0.0 -0.40 0.0
0.0 -0.35 0.0
0.0 -0.30 0.0
0.0 -0.25 0.0
0.0 -0.20 0.0
0.0 -0.15 0.0
0.0 -0.10 0.0
0.0 -0.05 0.0
