## lines starting with "##" work as lines for comments ! 2020.09.28
## this file is joint inversion control file
##-----10!-------20! inv201408_ap.ctl
face info file     !../../mesh_light/faceinfo.dat
1:cond,2:model     !1
ref  cond file     !../structure_light/cond_homo.msh
1:cond,2:model     !1
init cond file     !../structure_light/cond_homo.msh
output folder      |./result_inv/
Roughness type     |1
1:L,2:Cl,3:Mi,4:Gr !2
alpha init         !100.
factor(10^factor)  !-0.50
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
## Case for only tipper data
## imp = -1: no impedance data, 0: impedance data, 1: amp, phase data (not supported)
imp: 0, 1:amp,pha  !-1
# of observatories |9
########################################################################
iflag_tipper       !1
errorfloor         !0.05
MT1 tipper         !   1 ../fwd_3DMT_light/mt_err/MT1_TIP_err.dat
MT2 tipper         !   2 ../fwd_3DMT_light/mt_err/MT2_TIP_err.dat
MT3 tipper         !   3 ../fwd_3DMT_light/mt_err/MT3_TIP_err.dat
MT4 tipper         !   4 ../fwd_3DMT_light/mt_err/MT4_TIP_err.dat
MT5 tipper         !   5 ../fwd_3DMT_light/mt_err/MT5_TIP_err.dat
MT6 tipper         !   6 ../fwd_3DMT_light/mt_err/MT6_TIP_err.dat
MT7 tipper         !   7 ../fwd_3DMT_light/mt_err/MT7_TIP_err.dat
MT8 tipper         !   8 ../fwd_3DMT_light/mt_err/MT8_TIP_err.dat
MT9 tipper         !   9 ../fwd_3DMT_light/mt_err/MT9_TIP_err.dat
## icombine = 0: normal
## icombine = 1: integrate the outside blocks to one and assign one model parameter
## icombine = 2: integrate the outside blocks to one and fix the modelparameter with given cond
## icombine = -1: inner outer mode
icombine:0,1,2:fix !0
3
-1.0
0.0
1.0
3
-1.0
0.0
1.0
3
-3.0
0.0
0.5
##
ioutlevel:0,1:Jacob|0
final rms          |1.0

