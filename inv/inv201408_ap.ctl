!------10!-------20! inv201408_ap.ctl
face info file     !../mesh_IUGG/faceinfo.dat
1:cond,2:model     !1
ref  cond file     !./cond_homo.msh
1:cond,2:model     !1
init cond file     !./cond_init.msh
output folder      |./result_m0/
error floor  [0-1] !0.010
Roughness type     |1
1:L,2:Cl,3:Mi,4:Gr !2
alpha init         !1.0
iflag_replace 0,1  |1
Bo 0:off,1:on,2:Tr !0
# of srces for inv !1
Bx,By,Bz,Ex,Ey(5i2)!0 0 1 0 0
index of source(S1)!1
# of obesrvatories !4
A02 amp data       !   1 ../inv2014_15_base/result_201408_ap/A02_amp.dat
A04 amp data       !   2 ../inv2014_15_base/result_201408_ap/A04_amp.dat
A01 amp data       !   3 ../inv2014_15_base/result_201408_ap/A01_amp.dat
A03 amp data       !   4 ../inv2014_15_base/result_201408_ap/A03_amp.dat
A02 pha data       !   1 ../inv2014_15_base/result_201408_ap/A02_pha.dat
A04 pha data       !   2 ../inv2014_15_base/result_201408_ap/A04_pha.dat
A01 pha data       !   3 ../inv2014_15_base/result_201408_ap/A01_pha.dat
A03 pha data       !   4 ../inv2014_15_base/result_201408_ap/A03_pha.dat
icombine:0,1,2:fix !1
10
-1.5
-1.2
-0.75
-0.45
-0.15
0.15
0.45
0.75
1.2
1.5
10
-1.5
-1.2
-0.75
-0.45
-0.15
0.15
0.45
0.75
1.2
1.5
21
0.0
0.25
0.5
0.6
0.7
0.8
0.85
0.9
0.95
1.0
1.025
1.05
1.075
1.1
1.125
1.15
1.175
1.2
1.225
1.25
1.50
ioutlevel:0,1:Jacob|1
final rms          |1.0
