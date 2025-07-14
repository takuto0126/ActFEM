## lines starting with "##" work as lines for comments ! 2020.09.28
## this file is joint inversion control file
##-----10!-------20! inv201408_ap.ctl
face info file     !../mesh_joint/faceinfo.dat
1:cond,2:model     !1
ref  cond file     !./structure/cond_homo.msh
1:cond,2:model     !1
init cond file     !./structure/cond_homo.msh
output folder      |./result_inv/
Roughness type     |1
1:L,2:Cl,3:Mi,4:Gr !2
alpha init         !1
iflag_replace 0,1  |1
## iboundflag can set upper and lower limit of the conductivity value in inversion
## iboundflag = 0 : off : no boundary for conductivity value
## iboundflag = 1 : simple upper and lower limit will be specified by cutting
## iboundflag = 2 : transformed model variable that automatically satisfies upper
##                  and lower limit are used. (see Grayver et al. 2013)
## Note that if you use Smooth constraint, 0 is recommended.
 iboundflag =      !0
## ACTIVE data #########################################################
# of srces for inv !2
Bx,By,Bz,Ex,Ey(5i2)!0 0 1 0 0
Act err floor [0-1]!0.010
index of source(S1)!1
# of obesrvatories !4
A02 amp data       !   1 ./fwd_active_check/data/A02_S1_amp.dat
A04 amp data       !   2 ./fwd_active_check/data/A04_S1_amp.dat
A01 amp data       !   3 ./fwd_active_check/data/A01_S1_amp.dat
A03 amp data       !   4 ./fwd_active_check/data/A03_S1_amp.dat
A02 pha data       !   1 ./fwd_active_check/data/A02_S1_pha.dat
A04 pha data       !   2 ./fwd_active_check/data/A04_S1_pha.dat
A01 pha data       !   3 ./fwd_active_check/data/A01_S1_pha.dat
A03 pha data       !   4 ./fwd_active_check/data/A03_S1_pha.dat
index of source(S1)!2
# of obesrvatories !4
A02 amp data       !   1 ./fwd_active_check/data/A02_S2_amp.dat
A04 amp data       !   2 ./fwd_active_check/data/A04_S2_amp.dat
A01 amp data       !   3 ./fwd_active_check/data/A01_S2_amp.dat
A03 amp data       !   4 ./fwd_active_check/data/A03_S2_amp.dat
A02 pha data       !   1 ./fwd_active_check/data/A02_S2_pha.dat
A04 pha data       !   2 ./fwd_active_check/data/A04_S2_pha.dat
A01 pha data       !   3 ./fwd_active_check/data/A01_S2_pha.dat
A03 pha data       !   4 ./fwd_active_check/data/A03_S2_pha.dat
## MT data #############################################################
## initial version support only impedance data 2021.12.13
imp: 0, 1:amp,pha  !0
MT errfloor SSQ*   !0.010
# of observatories |4
MT1 impedance      !   1 ./fwd_3DMT_check/result/MT1_MT_imp.dat
MT2 impedance      !   2 ./fwd_3DMT_check/result/MT2_MT_imp.dat
MT3 impedance      !   3 ./fwd_3DMT_check/result/MT3_MT_imp.dat
MT4 impedance      !   4 ./fwd_3DMT_check/result/MT4_MT_imp.dat
MT1 impedance err  !   1 ./fwd_3DMT_check/mt_err/MT1_MT_imp_err.dat
MT2 impedance err  !   2 ./fwd_3DMT_check/mt_err/MT2_MT_imp_err.dat
MT3 impedance err  !   3 ./fwd_3DMT_check/mt_err/MT3_MT_imp_err.dat
MT4 impedance err  !   4 ./fwd_3DMT_check/mt_err/MT4_MT_imp_err.dat
########################################################################
icombine:0,1,2:fix !0
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

