# coded on 2025.06.09
#!/bin/bash

folder=result_lateral
sitefile=${folder}/site_mt.dat
infile1=${folder}/A02_TIP.dat
infile2=${folder}/A04_TIP.dat
infile3=${folder}/A01_TIP.dat
infile4=${folder}/A03_TIP.dat

gmt begin ${foler} pdf

gmt basemap -JX10/10 -R-1/1/-1/1 -Bxag0.5+l"Easting[km]" -Byag0.5+l"Nrothing [km]" -BWeSn
gmt makecpt -Crainbow -T0/1/0.01
# A02
x=`head -1 $sitefile | awk '{print($2)}'`
y=`head -1 $sitefile | awk '{print($3)}'`
echo $x $y
# -Sv reflects the angle from east to north
# Here X is eastward, Y is northward, Z is upward
#Bz = Tx*Bx + Ty*By. Here, (Real(Tx),Real(Ty)) points conductor without negative sign

awk -v x="$x" -v y="$y" '{print(x,y,atan2($4,$2)*180/3.1415,5*sqrt($2*$2+$4*$4),log($1)/log(10))}' $infile1 | gmt plot -Sv10p+e+p0.5 -C -Z -G+z

# A04
x=`head -2  $sitefile |tail -1| awk '{print($2)}'` # 1Hz
y=`head -2 $sitefile |tail -1| awk '{print($3)}'`  # 3Hz
echo $x $y
awk -v x="$x" -v y="$y" '{print(x,y,atan2($4,$2)*180/3.1415,5*sqrt($2*$2+$4*$4),log($1)/log(10))}' $infile2 | gmt plot -Sv10p+e+p0.5 -C -Z -G+z
# A01
x=`head -3 $sitefile |tail -1| awk '{print($2)}'`
y=`head -3 $sitefile |tail -1| awk '{print($3)}'`
echo $x $y
awk -v x="$x" -v y="$y" '{print(x,y,atan2($4,$2)*180/3.1415,5*sqrt($2*$2+$4*$4),log($1)/log(10))}' $infile3 | gmt plot -Sv10p+e+p0.5 -C -Z -G+z
# A03
x=`head -4 $sitefile |tail -1| awk '{print($2)}'`
y=`head -4 $sitefile |tail -1| awk '{print($3)}'`

# unit vector
echo "-0.8 -0.8  0 5" | gmt plot -Sv10p+e+p0.5 -G0
echo "-0.3 -0.7 1.0 " | gmt text -F+f14,Helvetica,black
echo $x $y
awk -v x="$x" -v y="$y" '{print(x,y,atan2($4,$2)*180/3.1415,5*sqrt($2*$2+$4*$4),log($1)/log(10))}' $infile4 | gmt plot -Sv10p+e+p0.5 -C -Z -G+z

gmt colorbar -Dx11/0+w10/0.3 -B+l"log10 Frequency" -C -B1  
gmt end show