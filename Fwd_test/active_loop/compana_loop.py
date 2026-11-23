#!/usr/bin/env python3
# Compare FEM Bz for polyline (multi-dipole) sources with the analytical solution.
# The analytical code (src/ana/active_ana.exe) accepts only single-dipole sources,
# so each dipole of the polyline is given as an independent source and the complex
# Bz responses are summed (superposition) for each polyline source.
#
# usage (in Fwd_test/active_loop/, after ./run.sh):
#   python3 compana_loop.py
# requires: ../../src/ana/active_ana.exe  (build with "make" in src/ana or ../active_ana/run.sh)
import cmath, math, os, subprocess, sys

anactl  = "../active_ana/ana.ctl"      # base analytical control file (obs, freq, conductivity)
srcs    = ["S1", "S2"]                 # source names (source point files are <name>.dat)
sites   = ["A02", "A04", "A01", "A03"]
femdir  = "result/"
anadir  = "result_ana/"
anaexe  = "../../src/ana/active_ana.exe"

#[1] make analytical control file with one source per dipole
base = open(anactl).read()
head = base[:base.index("# of sources")]
tail = base[base.index("Elcetric current"):]
head = head.replace("!./result/", "!./" + anadir)
dip = {}
lines = []
for s in srcs:
    L = [l for l in open(s + ".dat").read().split("\n") if l.strip()]
    n = int(L[0]); pts = [l.split() for l in L[1:n + 1]]
    dip[s] = []
    for j in range(n - 1):
        name = "%sd%02d" % (s, j + 1)
        dip[s].append(name)
        lines += ["Source Name        !" + name,
                  "source start point !" + " ".join(pts[j]),
                  "source end   point !" + " ".join(pts[j + 1])]
nd = sum(len(v) for v in dip.values())
open("ana_loop.ctl", "w").write(head + "# of sources       !%d\n" % nd + "\n".join(lines) + "\n" + tail)
os.makedirs(anadir, exist_ok=True)

#[2] run analytical code
with open(anadir + "ana.log", "w") as log:
    r = subprocess.run([anaexe], input="ana_loop.ctl\n", text=True, stdout=log, stderr=subprocess.STDOUT)
if r.returncode != 0:
    sys.exit("active_ana.exe failed, see " + anadir + "ana.log")

#[3] sum complex Bz over dipoles and compare with FEM (column 6: Bz amp [nT/A], 7: Bz phase [deg])
def rd(f):
    return [list(map(float, l.split())) for l in open(f) if l.strip()]

print("%-4s %-3s %6s %11s %11s %6s %8s %8s %6s" %
      ("site", "src", "f[Hz]", "|Bz|ana", "|Bz|FEM", "ratio", "pha_ana", "pha_FEM", "dpha"))
for s in srcs:
    for site in sites:
        tot = None
        for name in dip[s]:
            z = [r[5] * cmath.exp(1j * math.radians(r[6])) for r in rd(anadir + site + "_" + name + ".dat")]
            tot = z if tot is None else [a + b for a, b in zip(tot, z)]
        out = open(anadir + site + "_" + s + ".dat", "w")
        for r, z in zip(rd(femdir + site + "_" + s + ".dat"), tot):
            pa = math.degrees(cmath.phase(z))
            dp = (r[6] - pa + 180.) % 360. - 180.
            out.write("%15.7e %15.7e %15.7f\n" % (r[0], abs(z), pa))
            print("%-4s %-3s %6g %11.4e %11.4e %6.3f %8.2f %8.2f %6.2f" %
                  (site, s, r[0], abs(z), r[5], r[5] / abs(z), pa, r[6], dp))
        out.close()
