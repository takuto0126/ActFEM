# modified on 2020.12.08
# coded on 2018.02.27
# to draw a slice
####################################################################################
## Explanation of slice.exe: (ActFEMv1.0/src/src_post/slice.f90)
## slice.exe can generate two types of slices
## one is vertical slice, the other is horizontal slice
## You can choose either of them by setting three coefficients, A, B, and D
##
## [1] For vertical slices
##     slice plane is designed by Ax + By + D = 0
##     x is eastward, y is northward, z is upward, sea level is z=0, x,y,z are in [km]
##     For example, A=1,B=0,D=-1 generate x = 1 km plane.
##     Note that A should be > 0 in slice.exe
##
## [2] For horizontal slices
##     You shold set A=0, B=0, and D not equal to 0.
##     Then slice is generated for the plane of z + D = 0 (i.e., z = - D [km] plane)
##     z is upward, sea level is z=0
##
##  In both above two cases, four boudary values should be set after setting of A,B,D
##  1. left   boundary (lb) [km]
##  2. right  boundary (rb) [km]
##  3. bottom boundary (bb) [km]
##  4. top    boundary (tb) [km]
##
##  lb [km]     rb [km]
##  v           v
##  |-----------| <- tb [km] (elevation)
##  |  cross    |
##  |  section  |
##  |-----------| <- bb [km] (elevation, 0 means sea level)
##
##  These four parameters control the width and height of rectangles cut out of the slice plane
##  The origin in the slice plane is assumed to be the point
##  where the origin of the 3-D mesh is normally prohjected to the slice plane
##  The vertical axis in the slice plane is upward, the same as z axis in 3-D mesh
##  The horizontal axis in the slice plane is rightward when you see the plane
##  from the origin of the original 3-D mesh
##
######################################################################################

#!/bin/bash
#source /opt/intel/bin/compilervars.sh intel64

ite=09
folder=result_ms0

modelfile=./${folder}/model${ite}.dat
connectfile=./${folder}/model_connect.dat

ID=ms0
igendat=1 # 0: no dat file generation, 1: dat file generation

# CODE

#[1]## generate control file
cat << EOF > tmp0.ctl
##--------------- control file for slice.exe starts here!! ----------------------
## 2020.12.10 Takuto Minami
##
## Note that lines starting with "##" work for comments ---
##
##             20->|
mshfile            !../mesh_aso_A04/nakadake3d.msh
## You should chose either of two types of resistivity structure file format,
## cond.msh type or model**.dat type
## 0 for cond.msh type, 1 for model**.dat type
##  case for cond.msh type
##0:cond,1:model     !0
##ncondfile          !1
##polygonhead        !$ID
##condfile           !./result_ms0/cond${ite}.msh
## case for model**.dat type
0:cond,1:model     !1
model connect      !$connectfile
nmodelfile         !1
## polygonhead: header of output file
polygonhead        !$ID
modelfile          !$modelfile
## ibound 0: no bound, 1:lower, 2: upper bound
## If you want not to show some resistivity values larger / smaller than some threshold
## you can set by ibound =1 for lower bound, and 2 for upper bound of resistivity
## no bound case
ibound             !0
## lower bound case
##iboud              !1
##lower bo.log10[O.m]!1.0
## upper bound case
##iboud              !2
##upper bo.log10[O.m]!2.8
## nslice: number of slices you want to generate
nslice             !2
##---- Details for the first slice ------
## First slice file will be ${ID}_01.dat
## slice y = 0 (WE vertical slice)
A                  !0.0
B                  !1.0
D                  !0.0
boundary left      !-2.0
boundary right     !2.0
boundary bottom    !-0.5
boundary top       !1.5
##--- Details for the second slice -----
## second slice file will be ${ID}_02.dat
## slice x = 0 (SN vertical slice)
A                  !1.0
B                  !0.0
D                  !0.0
boundary left      !-2.0
boundary right     !2.0
boundary bottom    !-0.5
boundary top       !1.5
##------------------ control file for slice.exe ends here !! ----------------------
EOF

#[2]## generate tmp.dat
if [ $igendat -eq 1 ];then
SRC="../src/src_post"
cd $SRC
make clean
#make   # for the use of intel fortran compiler
make -f Makefile_gfort # for the use of gfortran
cd -

${SRC}/slice.exe <<EOF1
tmp0.ctl
EOF1

fi

#[3]## Draw picture
SL1="WE"
SL2="SN"

CPT="slice.cpt"
range=-2.0/2.0/-0.5/1.5
#bb=a0.25/a0.25:"Distance\(km\)":WeSn
bb=a0.25/a0.25wesn
scl=8/5
range2=0/8/0/5

outfile=slice_ms0.ps
outpng=slice_ms0.png
outpdf=slice_ms0.pdf

bb1=a0.25/a0.25:Elevation[km]:Wesn
bb2=a1f0.25:Distance[km]:/a0.25:Elevation[km]:WeSn
xshift=2.5   ; yshift=2.5 # x=2.5 is default

gmt makecpt -Crainbow -T0.2/2.8/0.01 -I -Dred,purple > ${CPT}

echo "-10 -10" | gmt psxy -JX$scl -R$range -Sc0.1  -K  > $outfile

scl2=20/20
range3=0/20/0/20

echo $ID
yy=`echo " scale=1; 2.5 + 5.3" | bc`
if [ -e "./${ID}_01.dat" ];then
gmt psxy "${ID}_01.dat" -JX$scl -R$range -B${bb1} -C${CPT} -L -K -O -X${xshift} -Y${yy} >> $outfile
else
gmt psbasemap -JX$scl -R$range -B${bb1}  -K -O -X${xshift} -Y${yy} >> $outfile
fi
gmt pstext -JX$scl -R$range2  -G255 -K -O <<EOF >> $outfile
0.9  4.7 14 0 4 LM $ID $ite
0.15  4.6 12 0 4 LM W
7.85  4.6 12 0 4 RM E
EOF
if [ -e  "${ID}_02.dat" ];then
gmt psxy "${ID}_02.dat" -JX$scl -R$range -B${bb2} -C${CPT} -L -K -O -Y-5.5 >> $outfile
else
gmt psbasemap -JX$scl -R$range -B${bb2[$i]}  -K -O -Y-3.0 >> $outfile
fi
gmt pstext -JX$scl -R$range2  -G255 -K -O <<EOF >> $outfile
0.9  4.7 14 0 4 LM $ID
0.15  4.6 12 0 4 LM S
7.85  4.6 12 0 4 RM N
EOF

gmt psscale -Dx9.1/4.5/10/0.3 -B0.2 -G0.2/2.8 -C$CPT -E -K -O >> $outfile

gmt pstext -R0/20/0/20 -JX20/20 -O << EOF >> $outfile
9.5 11.0 13 0 4 CM log@-10@- @%30%@~r@~@%%
9.5 10.3 13 0 4 CM (@~W@~m)
EOF

#gv $outfile &
rm gmt.history tmp.ctl $CPT # 2018.03.16
convert -density 200 $outfile $outpdf
open $outpdf
