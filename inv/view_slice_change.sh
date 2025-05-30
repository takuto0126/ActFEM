# coded on Feb 5, 2021
# see view_slice.sh for details of slice.exe
#!/bin/bash
#source /opt/intel/bin/compilervars.sh intel64

connectfile=./result_ms0/model_connect.dat
model1=./result_ms0/model12.dat # before
model2=./model12_v2.dat         # after
ID1=model1
ID2=model2

igendat=1 # 0: no dat file generation, 1: dat file generation
# CODE

#[1]## generate control file
cat << EOF > tmp1.ctl
mshfile            !../mesh_aso_A04/nakadake3d.msh
0:cond,1:model     !1
model connect      !$connectfile
nmodelfile         !1
polygonhead        !$ID1
modelfile          !$model1
ibound             !0
nslice             !1
A                  !0.0
B                  !1.0
D                  !0.0
boundary left      !-2.0
boundary right     !2.0
boundary bottom    !-0.5
boundary top       !1.5
EOF
cat << EOF > tmp2.ctl
mshfile            !../mesh_aso_A04/nakadake3d.msh
0:cond,1:model     !1
model connect      !$connectfile
nmodelfile         !1
polygonhead        !$ID2
modelfile          !$model2
ibound             !0
nslice             !1
A                  !0.0
B                  !1.0
D                  !0.0
boundary left      !-2.0
boundary right     !2.0
boundary bottom    !-0.5
boundary top       !1.5
EOF


#[2]## generate tmp.dat
if [ $igendat -eq 1 ];then
SRC="../src/src_post"
cd $SRC
make clean
#make   # for the use of intel fortran compiler
make -f Makefile_gfort # for the use of gfortran
cd -

# slice 1
${SRC}/slice.exe <<EOF1
tmp1.ctl
EOF1
# slice 2
${SRC}/slice.exe <<EOF
tmp2.ctl
EOF

fi

#[3]## Draw picture
SL1="WE"
SL2="WE"

CPT="slice.cpt"
range=-2.0/2.0/-0.5/1.5
#bb=a0.25/a0.25:"Distance\(km\)":WeSn
bb=a0.25/a0.25wesn
scl=8/5
range2=0/8/0/5

outfile=slice_change.ps
outpng=slice_change.png
outpdf=slice_change.pdf

bb1=a0.25/a0.25:Elevation[km]:Wesn
bb2=a1f0.25:Distance[km]:/a0.25:Elevation[km]:WeSn
xshift=2.5   ; yshift=2.5 # x=2.5 is default

gmt makecpt -Crainbow -T0.2/2.8/0.01 -I -Dred,purple > ${CPT}

echo "-10 -10" | gmt psxy -JX$scl -R$range -Sc0.1  -K  > $outfile

scl2=20/20
range3=0/20/0/20

echo $ID1
yy=`echo " scale=1; 2.5 + 5.3" | bc`
if [ -e "./${ID1}_01.dat" ];then
gmt psxy "${ID1}_01.dat" -JX$scl -R$range -B${bb1} -C${CPT} -L -K -O -X${xshift} -Y${yy} >> $outfile
else
gmt psbasemap -JX$scl -R$range -B${bb1}  -K -O -X${xshift} -Y${yy} >> $outfile
fi
gmt pstext -JX$scl -R$range2  -G255 -K -O <<EOF >> $outfile
0.9  4.7 14 0 4 LM $ID1 $ite
0.15  4.6 12 0 4 LM W
7.85  4.6 12 0 4 RM E
EOF
if [ -e  "${ID2}_01.dat" ];then
gmt psxy "${ID2}_01.dat" -JX$scl -R$range -B${bb2} -C${CPT} -L -K -O -Y-5.5 >> $outfile
else
gmt psbasemap -JX$scl -R$range -B${bb2[$i]}  -K -O -Y-3.0 >> $outfile
fi
gmt pstext -JX$scl -R$range2  -G255 -K -O <<EOF >> $outfile
0.9  4.7 14 0 4 LM $ID2
0.15  4.6 12 0 4 LM W
7.85  4.6 12 0 4 RM E
EOF

gmt psscale -Dx9.1/4.5/10/0.3 -B0.2 -G0.2/2.8 -C$CPT -E -K -O >> $outfile

gmt pstext -R0/20/0/20 -JX20/20 -O << EOF >> $outfile
9.5 11.0 13 0 4 CM log@-10@- @%30%@~r@~@%%
9.5 10.3 13 0 4 CM (@~W@~m)
EOF

#gv $outfile &
rm gmt.history tmp.ctl $CPT # 2018.03.16
convert -density 200 $outfile $outpdf
#open $outpdf
