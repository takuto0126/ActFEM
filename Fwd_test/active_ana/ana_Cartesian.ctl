## lines starting with ## wrok as comment lines
##-----10!-------20!
output folder      !./result_Cartesian/
# of frequency     !11
Frequency [Hz]     !1.d0
Frequency [Hz]     !33.d0
Frequency [Hz]     !10.d0
Frequency [Hz]     !33.d0
Frequency [Hz]     !100.d0
Frequency [Hz]     !   330.d0
Frequency [Hz]     !  1000.d0
Frequency [Hz]     !  3300.d0
Frequency [Hz]     ! 10000.d0
Frequency [Hz]     ! 33000.d0
Frequency [Hz]     !100000.d0
# of observatory   !5
lonlat(1),xyz (2)  !2
1  Name            !A01
1  xyz             !0.0       0.001         -0.001
2  Name            !A03
2  xyz,sigma,A[km] !0.0       0.003         -0.001
3  Name            !A10
3  xyz,sigma,A[km] !0.0       0.010         -0.001
4  Name            !A30
4  xyz,sigma,A[km] !0.0       0.030         -0.001
4  Name            !A100
4  xyz,sigma,A[km] !0.0       0.100         -0.001
ixyflg 0:no,1:surfv!0
# of sources       !1
Source Name        !S1
source start point !-0.1         0.0         -0.001
source end   point !0.1          0.0         -0.001
Elcetric current[A]! 1.0
ds division len [m]! 1.0
## conductivity structure
sigma_air    [S/m] !1.e-8
nalyer             !2
depth1   [km]      !1.0
cond1              !0.01
cond2              !0.01
