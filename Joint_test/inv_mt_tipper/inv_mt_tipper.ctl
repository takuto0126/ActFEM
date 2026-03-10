## lines starting with "##" work as lines for comments ! 2020.09.28
## this file is joint inversion control file
##-----10!-------20! inv201408_ap.ctl
face info file     !../../mesh_joint/faceinfo.dat
1:cond,2:model     !1
ref  cond file     !../structure/cond_homo.msh
1:cond,2:model     !1
init cond file     !../structure/cond_homo.msh
output folder      |./result_inv/
Roughness type     |1
1:L,2:Cl,3:Mi,4:Gr !2
alpha init         !100.
factor(10^factor)  !-1.0
iflag_replace 0,1  |0
## iboundflag can set upper and lower limit of the conductivity value in inversion
## iboundflag = 0 : off : no boundary for conductivity value
## iboundflag = 1 : simple upper and lower limit will be specified by cutting
## iboundflag = 2 : transformed model variable that automatically satisfies upper
##                  and lower limit are used. (see Grayver et al. 2013)
## Note that if you use Smooth constraint, 0 is recommended.
 iboundflag =      !0
## ACTIVE data #########################################################
## of srces for inv !2
##Bx,By,Bz,Ex,Ey(5i2)!0 0 1 0 0
##Act err floor [0-1]!0.01
##index of source(S1)!1
### of obesrvatories !4
##A02 amp data       !   1 ../fwd_active/data/A02_S1_amp.dat
##A04 amp data       !   2 ../fwd_active/data/A04_S1_amp.dat
##A01 amp data       !   3 ../fwd_active/data/A01_S1_amp.dat
##A03 amp data       !   4 ../fwd_active/data/A03_S1_amp.dat
##A02 pha data       !   1 ../fwd_active/data/A02_S1_pha.dat
##A04 pha data       !   2 ../fwd_active/data/A04_S1_pha.dat
##A01 pha data       !   3 ../fwd_active/data/A01_S1_pha.dat
##A03 pha data       !   4 ../fwd_active/data/A03_S1_pha.dat
##index of source(S2)!2
### of obesrvatories !4
##A02 amp data       !   1 ../fwd_active/data/A02_S2_amp.dat
##A04 amp data       !   2 ../fwd_active/data/A04_S2_amp.dat
##A01 amp data       !   3 ../fwd_active/data/A01_S2_amp.dat
##A03 amp data       !   4 ../fwd_active/data/A03_S2_amp.dat
##A02 pha data       !   1 ../fwd_active/data/A02_S2_pha.dat
##A04 pha data       !   2 ../fwd_active/data/A04_S2_pha.dat
##A01 pha data       !   3 ../fwd_active/data/A01_S2_pha.dat
##A03 pha data       !   4 ../fwd_active/data/A03_S2_pha.dat
## MT data #############################################################
## initial version support only impedance data 2021.12.13
## imp_ap_flag=-1 : no imp/ap data, 0: imp data, 1: amp, phase data
imp: 0, 1:amp,pha  !0
unit:1Ohm,2mV/km/nT!2
MT errfloor SSQ*   !0.1
# of observatories |5
MT1 impedance      !   1 ../fwd_3DMT/result/MT1_MT_imp.dat
MT2 impedance      !   2 ../fwd_3DMT/result/MT2_MT_imp.dat
MT3 impedance      !   3 ../fwd_3DMT/result/MT3_MT_imp.dat
MT4 impedance      !   4 ../fwd_3DMT/result/MT4_MT_imp.dat
MT5 impedance      !   5 ../fwd_3DMT/result/MT5_MT_imp.dat
MT1 impedance err  !   1 ../fwd_3DMT/mt_err/MT1_MT_imp_err.dat
MT2 impedance err  !   2 ../fwd_3DMT/mt_err/MT2_MT_imp_err.dat
MT3 impedance err  !   3 ../fwd_3DMT/mt_err/MT3_MT_imp_err.dat
MT4 impedance err  !   4 ../fwd_3DMT/mt_err/MT4_MT_imp_err.dat
MT5 impedance err  !   5 ../fwd_3DMT/mt_err/MT5_MT_imp_err.dat
########################################################################
iflag_tipper       !1
errorfloor         !0.005
MT1 tipper         !   1 ../fwd_3DMT/mt_err/MT1_TIP_err.dat
MT2 tipper         !   2 ../fwd_3DMT/mt_err/MT2_TIP_err.dat
MT3 tipper         !   3 ../fwd_3DMT/mt_err/MT3_TIP_err.dat
MT4 tipper         !   4 ../fwd_3DMT/mt_err/MT4_TIP_err.dat
MT5 tipper         !   5 ../fwd_3DMT/mt_err/MT5_TIP_err.dat
## icombine = 0: normal
## icombine = 1: integrate the outside blocks to one and assign one model parameter
## icombine = 2: integrate the outside blocks to one and fix the modelparameter with given cond
## icombine = -1: inner outer mode
icombine:0,1,2:fix !0
19
-1.5
-1.2
-0.9
-0.75
-0.6
-0.45
-0.3
-0.2
-0.1
0.0
0.1
0.2
0.3
0.45
0.6
0.75
0.9
1.2
1.5
19
-1.5
-1.2
-0.9
-0.75
-0.6
-0.45
-0.3
-0.2
-0.1
0.0
0.1
0.2
0.3
0.45
0.6
0.75
0.9
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
##
ioutlevel:0,1:Jacob|0
final rms          |1.1

