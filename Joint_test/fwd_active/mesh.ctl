[MESH]
## lines starting with "##" work as lines for comments ! 2020.09.28
## this file is forward control file
## iflag_map = 1or 2
## = 1 (ECP) : lon lat topo + Equidistant Cylindrical Projection  (ECP)
##                 x=planetrad*cos(latorigin)*(lon-lonorigin), y=planetrad*(lat-latorigin)
##                 one additional line: lonlat origin (2f15.7) 
## = 2 (UTM) : lon lat topo + Universal Transvers Mercatol (UTM) projection
##                 two additional line: UTM zone like 52S (a3) and lonlat origin (2f15..7)
iflag_map          !2
UTM zone           !52S
lonlatorigin       !131.084782     32.884882
mesh file          !../mesh_joint/nakadake3d.msh
2d triangle z file !../mesh_joint/nakadake2dz.msh
local line file    !../mesh_joint/lineinfo.dat
angle              !0.0
water level        !0
output folder      !./result/
condflag 0:home,1: !1
##nvolume            !1
##cond               !0.01
condfile           !../structure/cond_G3.msh
