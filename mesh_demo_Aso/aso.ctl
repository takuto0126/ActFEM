## lines starting with "##" work as lines for comments ! 2020.09.28
## this file is mesh / forward control file
## In ActFEM, x is east, y is north, z is upward
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
topofile           !../topo/aso_map.data
lon lat shift      !0.0            0.0
topofile           !../topo/topo127_134_29_36.xyz
lon lat shift      !0.0            0.0
mesh file          !../mesh_demo_Aso/nakadake3d.msh
2d triangle z file !../mesh_demo_Aso/nakadake2d_z.msh
local mesh file    !../mesh_demo_Aso/lineinfo.dat
angle              !0.0
water level        !0
output folder      !-
header2d  (a50)    !nakadake2d
header3d  (a50)    !nakadake3d
# of frequency     !3
Frequency [Hz]     !1.d0
Frequency [Hz]     !10.d0
Frequency [Hz]     !100.d0
west bound         !-1.7
east bound         !1.7
south bound        !-1.3
north bound        !1.5
lenout [km]        !300.0
upz in [km]  (>0)  !1.3
downz in [km](<0)  !-1.1
zmax   [km]        !300.0
zmin   [km]        !-300.0
sizein [km]        !0.10
sizebo [km]        !70.0
sigma_obs [km]     !0.32
A_obs     [km]     !0.006
dlen_source [km]   !0.3
sigma_src [km]     !0.32
A_src     [km]     !0.005
# of observatory   !7
## lonlatflag = 1: src/obs coordinates are provided by (lon[deg],lat[deg],alt from the surface [km])
## lonlatflag = 2: src/obs coordinates are provided by (x[km],y[km],z[km]) 
lonlatflag         !1
1  Name            !A02
1  xyz             !131.083411   32.886706  -0.001
2  Name            !A04
2  xyz,sigma,A[km] !131.081939   32.884808  -0.001
3  Name            !A01
3  xyz,sigma,A[km] !131.083367   32.882725  -0.001
4  Name            !A03
4  xyz,sigma,A[km] !131.086847   32.881981  -0.001
5  Name            !DUM
5  xyz,sigma,A[km] !131.084782   32.884882  -0.100
6  Name            !DM2
6  xyz,sigma,A[km] !131.083335   32.884882  -0.100
7  Name            !DM3
7  xyz,sigma,A[km] !131.079672   32.884945  -0.100
ixyflg 0:no,1:surfv!0
# of sources       !1
Source Name        !S2
source start point !131.0784333  32.8908028  -0.001
source end   point !131.0814639  32.8912333  -0.001
Elcetric current[A]! 1.0
sigma_air    [S/m] !1.e-8
condflag 0:home,1: !1
condfile           !../mesh_aso/cond.msh
