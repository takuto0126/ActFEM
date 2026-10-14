#!/bin/bash
# Plot Bz amplitude vs. distance from the dipole source (log-log),
# comparing Primary / Secondary / Total fields at several frequencies.
#
# Data source: ./result_Cartesian/A**_S1[.dat|_pr.dat|_se.dat]
#   (produced by src/ana/n_active_ana.f90, subroutine OUTOBSFILESFWD)
#   no suffix = Total field, "_pr" = Primary field, "_se" = Secondary field
#
# Column layout of each file (see forward() in n_active_ana.f90):
#   col1 =freq[Hz]
#   col2,3 =Bx amp[nT],phase[deg]
#   col4,5 =By amp[nT],phase[deg]
#   col6,7 =Bz amp[nT],phase[deg]   <-- Bz amplitude/phase are columns 6,7
#   col8,9 =Ex amp[uV/m],phase[deg]
#   col10,11=Ey amp[uV/m],phase[deg]
#
# Station name A** encodes the distance from the dipole in meters
# (A01=1m, A03=3m, A10=10m, A30=30m, A100=100m; see xyz in ana_Cartesian.ctl).

set -e

RESDIR=./result_Cartesian
OUTNAME=bz_vs_distance
FREQS="1 10 100 1000 10000 100000"
STATIONS="A01 A03 A10 A30 A100"
BZCOL=6   # Bz amplitude column (1-indexed)

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

gmt begin ${OUTNAME} pdf
gmt set FONT_ANNOT_PRIMARY 9p,Helvetica,black FONT_LABEL 10p,Helvetica,black \
        FONT_TAG 11p,Helvetica-Bold,black MAP_FRAME_TYPE plain MAP_GRID_PEN 0.25p,gray80

gmt subplot begin 2x3 -Fs6.5c/6c -M1.5c/0.9c -R0.7/200/1e-5/1e3 -JX6.5cl/6cl \
    -T"Bz amplitude vs. distance from dipole (Source S1)"

panel=0
for f in $FREQS; do

  for kind in pr se total; do
    case $kind in
      pr)    suffix="_pri" ;;
      se)    suffix="_sec" ;;
      total) suffix=""    ;;
    esac
    dat="$TMPDIR/${kind}.d"
    : > "$dat"
    for st in $STATIONS; do
      dist=$((10#${st#A}))   # A030 -> 30 (strip leading zeros)
      file="$RESDIR/${st}_S1${suffix}.dat"
      amp=$(awk -v f="$f" -v c="$BZCOL" '$1==f{print $c; exit}' "$file")
      echo "$dist $amp" >> "$dat"
    done
    sort -n -o "$dat" "$dat"
  done

  gmt subplot set $panel
  gmt basemap -Bxa1f3g3+l"Distance from dipole [m]" -Bya1f3g3+l"Bz amplitude [nT]" -BWeSn+t"${f} Hz"

  gmt plot "$TMPDIR/pr.d"    -W1p,blue    -Sc0.22c -Gblue
  gmt plot "$TMPDIR/se.d"    -W1p,red     -St0.22c -Gred
  gmt plot "$TMPDIR/total.d" -W1p,green     -Sd0.22c -Ggreen

  if [ $panel -eq 0 ]; then
    gmt legend -DjML+w3.2c+o0.15c/0.15c -F+gwhite+p0.5p <<- EOF
	S 0.3c c 0.22c blue  1p,blue  0.6c Primary
	S 0.3c t 0.22c red   1p,red   0.6c Secondary
	S 0.3c d 0.22c green 1p,green 0.6c Total
	EOF
  fi

  panel=$((panel+1))
done

gmt subplot end
gmt end show
