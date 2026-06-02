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
angle              !0.0
water level        !0
header2d  (a50)    !nakadake2d
header3d  (a50)    !nakadake3d
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
ixyflg 0:no,1:surfv!0
# of sources       !2
Source Name        !S1
source start point !131.086490   32.877356   -0.001
source end   point !131.090840   32.877658   -0.001
Source Name        !S2
source start point !131.0784333  32.8908028  -0.001
source end   point !131.0814639  32.8912333  -0.001