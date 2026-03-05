# generated on 2026.03.01
#!/bin/bash

gfortran generr_mt_tip.f90  -o generr_tip.exe

cat <<EOF > err_tip.ctl
        10     20->|header
# of obs files     |5
# of frequencies   |5
input folder       |./result/
output folder      |./mt_err/
in file1 ***.dat   |MT1_TIP.dat
errfile1           |MT1_TIP_err.dat
in file2 ***.dat   |MT2_TIP.dat
errfile2           |MT2_TIP_err.dat
in file3 ***.dat   |MT3_TIP.dat
errfile3           |MT3_TIP_err.dat
in file4 ***.dat   |MT4_TIP.dat
errfile4           |MT4_TIP_err.dat
in file4 ***.dat   |MT5_TIP.dat
errfile4           |MT5_TIP_err.dat
ratio              !0.01
EOF

./generr_tip.exe < err_tip.ctl

rm err_tip.ctl

