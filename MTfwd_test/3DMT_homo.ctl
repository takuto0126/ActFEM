## lines starting with "##" work as lines for comments ! 2021.06.08
## this file is forward control file
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
iflag_map          !1
lonlatorigin       !131.084782     32.884882
# of topofile      !1
topofile           !../mesh_flat/topo_aso_flat.dat
lon lat shift      !0.0         0.0
mesh file          !../mesh_flat/nakadake3d.msh
2d triangle z file !../mesh_flat/nakadake2dz.msh
local line file    !../mesh_flat/lineinfo.dat
rotate angle[NdegE]!0.0
output folder      !./result_homo/
header2d  (a50)    !nakadake2d
header3d  (a50)    !nakadake3d
# of frequency     !2
Frequency [Hz]     !1.d0
Frequency [Hz]     !3.d0
west bound         !-1.7
east bound         !1.7
south bound        !-1.5
north bound        !1.5
lenout [km]        !25.0
upz in [km]  (>0)  !1.3
downz in [km](<0)  !-1.1
zmax   [km]        !50.0
zmin   [km]        !-50.0
sizein [km]        !0.10
sizebo [km]        !2.0
sigma_obs [km]     !0.4
A_obs     [km]     !0.01
dlen_source [km]   !0.1
sigma_src [km]     !0.3
A_src     [km]     !0.005
# of observatory   !4
lonlat(1),xyz (2)  !1
1  Name            !A02
1  xyz             !131.083411   32.886706  -0.001
2  Name            !A04
2  xyz,sigma,A[km] !131.081939   32.884808  -0.001
3  Name            !A01
3  xyz,sigma,A[km] !131.083367   32.882725  -0.001
4  Name            !A03
4  xyz,sigma,A[km] !131.086847   32.881981  -0.001
ixyflg 0:no,1:surfv!0
sigma_air    [S/m] !1.e-8
## case for homogeneous conductivity 
condflag 0:homo,1: !0
Volume             !1
conductivity[S/m]  !0.01
##condflag 0:homo,1: !1
##condfile           !cond_lateral_contrast.msh
