## obs_active.ctl
## settings for forward modelling for ACTIVE
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
# of sources       !2
Source Name        !S1
source start point !131.086490   32.877356   -0.001
source end   point !131.090840   32.877658   -0.001
Source Name        !S2
source start point !131.078109   32.889686   -0.001
source end   point !131.082206   32.891490   -0.001
Elcetric current[A]! 1.0
sigma_air    [S/m] !1.e-8
condflag 0:home,1: !1
##nvolume            !1
##cond               !0.01
condfile           !../structure/cond_G3.msh
