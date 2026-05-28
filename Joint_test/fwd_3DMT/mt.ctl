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
iflag_map          !2
UTM zone           !52S
lonlatorigin       !131.084782     32.884882
# of topofile      !1
topofile           !../mesh_joint/topo_aso_flat.dat
lon lat shift      !0.0         0.0
mesh file          !../mesh_joint/nakadake3d.msh
2d triangle z file !../mesh_joint/nakadake2dz.msh
local line file    !../mesh_joint/lineinfo.dat
angle              !0.0
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
# of observatory   !9
lonlat(1),xyz (2)  !1
1  Name            !MT1
1  xyz             !131.079780      32.887800      -0.001000000
2  Name            !MT2
2  xyz             !131.090989      32.887800      -0.001000000
3  Name            !MT3
3  xyz             !131.079780      32.880280      -0.001000000
4  Name            !MT4
4  xyz             !131.090989      32.880280      -0.001000000
5  Name            !MT5
5  xyz             !131.085386      32.884040      -0.001000000
6  Name            !MT6
6  xyz             !131.079780      32.884040      -0.001000000
7  Name            !MT7
7  xyz             !131.090989      32.884040      -0.001000000
8  Name            !MT8
8  xyz             !131.085386      32.887800      -0.001000000
9  Name            !MT9
9  xyz             !131.085386      32.880280      -0.001000000
ixyflg 0:no,1:surfv!0
sigma_air    [S/m] !1.e-8
condflag 0:homo,1: !1
##Volume             !1
##cond               !0.01
condfile           !../structure/cond_G3.msh

