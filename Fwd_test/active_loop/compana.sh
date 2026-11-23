#!/bin/bash

caldir=result/
anadir=../active_ana/result/
outps=compana.ps
outpdf=compana.pdf

gnuplot <<EOF
set terminal postscript enhanced color
set output "$outps"
set multiplot
set size 0.5,0.5
set lmargin 9
set rmargin 5

# S1
# Amp
set origin 0,0.5
set logscale x
set logscale y
set yrange [0.01:1]
set xrange [0.9:105]
set ylabel "Bz Amp [nT/A]"
set xlabel "Frequency [Hz]"
plot "${anadir}A02_S1.dat" u 1:6 w l lw 2 title "A02 S1",\
     "${caldir}A02_S1.dat" u 1:6 w p lw 2 notitle ,\
     "${anadir}A04_S1.dat" u 1:6 w l lw 2 title "A04 S1",\
     "${caldir}A04_S1.dat" u 1:6 w p lw 2 notitle,\
     "${anadir}A01_S1.dat" u 1:6 w l lw 2 title "A01 S1",\
     "${caldir}A01_S1.dat" u 1:6 w p lw 2 notitle,\
     "${anadir}A03_S1.dat" u 1:6 w l lw 2 title "A03 S1",\
     "${caldir}A03_S1.dat" u 1:6 w p lw 2 notitle

# Phase
unset logscale y
set origin 0,0
set yrange [-60:20]
set key left bottom
set ylabel "Bz Phase [deg]"
plot "${anadir}A02_S1.dat" u 1:7 w l lw 2 title "A02 S1",\
     "${caldir}A02_S1.dat" u 1:7 w p lw 2 notitle ,\
     "${anadir}A04_S1.dat" u 1:7 w l lw 2 title "A04 S1",\
     "${caldir}A04_S1.dat" u 1:7 w p lw 2 notitle,\
     "${anadir}A01_S1.dat" u 1:7 w l lw 2 title "A01 S1",\
     "${caldir}A01_S1.dat" u 1:7 w p lw 2 notitle,\
     "${anadir}A03_S1.dat" u 1:7 w l lw 2 title "A03 S1",\
     "${caldir}A03_S1.dat" u 1:7 w p lw 2 notitle

# S2
# Amp
set origin 0.5,0.5
set logscale y
set yrange [0.01:1]
set ylabel "Bz Amp [nT/A]"
set key top right
plot "${anadir}A02_S2.dat" u 1:6 w l lw 2 title "A02 S2",\
     "${caldir}A02_S2.dat" u 1:6 w p lw 2 notitle ,\
     "${anadir}A04_S2.dat" u 1:6 w l lw 2 title "A04 S2",\
     "${caldir}A04_S2.dat" u 1:6 w p lw 2 notitle,\
     "${anadir}A01_S2.dat" u 1:6 w l lw 2 title "A01 S2",\
     "${caldir}A01_S2.dat" u 1:6 w p lw 2 notitle,\
     "${anadir}A03_S2.dat" u 1:6 w l lw 2 title "A03 S2",\
     "${caldir}A03_S2.dat" u 1:6 w p lw 2 notitle

# Phase
unset logscale y
set origin 0.5,0
set yrange [120:200]
set key left bottom
set ylabel "Bz Phase [deg]"
plot "${anadir}A02_S2.dat" u 1:7 w l lw 2 title "A02 S2",\
     "${caldir}A02_S2.dat" u 1:7 w p lw 2 notitle ,\
     "${anadir}A04_S2.dat" u 1:7 w l lw 2 title "A04 S2",\
     "${caldir}A04_S2.dat" u 1:7 w p lw 2 notitle,\
     "${anadir}A01_S2.dat" u 1:7 w l lw 2 title "A01 S2",\
     "${caldir}A01_S2.dat" u 1:7 w p lw 2 notitle,\
     "${anadir}A03_S2.dat" u 1:7 w l lw 2 title "A03 S2",\
     "${caldir}A03_S2.dat" u 1:7 w p lw 2 notitle

EOF
convert -density 300 -rotate 90 $outps $outpdf
open $outpdf
