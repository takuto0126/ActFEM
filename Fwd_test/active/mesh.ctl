## lines starting with "##" work as lines for comments
## base.ctl: mesh files + conductivity (common to active / MT / tidal)
## NO geometry bounds, mesh-gen sizes, or obs/src parameters
##-----10!-------20!
[BASE]
   itopofile       !1
iflag_map          !2
UTM zone           !52S
lonlatorigin       !131.084782     32.884882
# of topo file     !1
topofile           !../mesh_flat/topo_aso_flat.dat
lon lat shift      !0.0         0.0
mesh file          !../mesh_flat/nakadake3d.msh
2d triangle z file !../mesh_flat/nakadake2dz.msh
local line file    !../mesh_flat/lineinfo.dat
angle              !0.0
water level        !0
output folder      !./result/
header2d  (a50)    !nakadake2d
header3d  (a50)    !nakadake3d
sigma_air    [S/m] !1.e-8
condflag 0:homo,1: !0
nvolume            !1
cond               !0.01
