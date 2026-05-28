## lines starting with "##" work as lines for comments
##-----10!-------20!
## itopofile controls whether topo file is provided or not
## itopofile = 0   : no topo file and z=0 surface is assumed
##                   obs and src coordinates should be provided by cartesian position
## itopofile = 1   : topofile(s) is(are) provided, iflag_map should be provided
   itopofile       !1
## iflag_map controls the type of Map projection
## iflag_map = 1 (ECP) : lon lat topo + Equidistant Cylindrical Projection  (ECP)
##                 x=planetrad*cos(latorigin)*(lon-lonorigin), y=planetrad*(lat-latorigin)
##                 one additional line: lonlat origin (2f15.7) 
## iflag_map = 2 (UTM) : lon lat topo + Universal Transvers Mercatol (UTM) projection
##                 two additional line: UTM zone like 52S (a3) and lonlat origin (2f15..7)
iflag_map          !2
UTM zone           !52S
lonlatorigin       !131.084782     32.884882
# of topo file     !2
topofile           !../../topo/aso_map.data
lon lat shift      !0.0            0.0
topofile           !../../topo/topo127_134_29_36.xyz
lon lat shift      !0.0            0.0
mesh file          !../mesh_joint/nakadake3d.msh
2d triangle z file !../mesh_joint/nakadake2d_z.msh
local mesh file    !../mesh_joint/lineinfo.dat
angle              !0.0
water level        !0
output folder      !-
header2d  (a50)    !nakadake2d
header3d  (a50)    !nakadake3d
# of frequency     !3
Frequency [Hz]     !1.d0
Frequency [Hz]     !10.d0
Frequency [Hz]     !100.d0
west bound         !-1.8
east bound         !1.7
south bound        !-1.8
north bound        !1.7
lenout [km]        !5.0
upz in [km]        !1.4
downz in [km]      !0.3
zmax   [km]        !25.0
zmin   [km]        !-25.0
sizein [km]        !0.11
sizebo [km]        !2.0
sigma_obs [km]     !0.25
A_obs     [km]     !0.01
dlen_source [km]   !0.5
sigma_src [km]     !0.32
A_src     [km]     !0.01
# of observatory   !12
lonlat(1),xyz (2)  !1
1  Name            !A02
1  xyz             !131.0834110     32.8867060     -0.0010000
2  Name            !A04
2  xyz,sigma,A[km] !131.0819390     32.8848080     -0.0010000
3  Name            !A01
3  xyz,sigma,A[km] !131.0833670     32.8827250     -0.0010000
4  Name            !A03
4  xyz,sigma,A[km] !131.0868470     32.8819810     -0.0010000
5  Name            !DUM
5  xyz,sigma,A[km] !131.0847820     32.8848820     -0.1000000
5  Name            !DM1
5  xyz,sigma,A[km] !131.0847820     32.8848820     -0.4000000
6  Name            !DM2
6  xyz,sigma,A[km] !131.0833350     32.8848820     -0.1000000
7  Name            !DM3
7  xyz,sigma,A[km] !131.0796720     32.8849450     -0.1000000
8  Name            !MT1
8  xyz,sigma,A[km] !131.0745870     32.8929300     -0.001000000
9  Name            !MT2
9  xyz,sigma,A[km] !131.0945430     32.8929300     -0.001000000
10 Name            !MT3
10 xyz,sigma,A[km] !131.0745870     32.8763960     -0.001000000
11 Name            !MT4
11 xyz,sigma,A[km] !131.0945430     32.8763960     -0.001000000
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
condflag 0:home,1: !1
condfile           !../mesh_aso/cond.msh
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
