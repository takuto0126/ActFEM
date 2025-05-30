mshfile            !../mesh_aso_A04/nakadake3d.msh
0:c,1:m,2:dc,3:dm  !3
connectfile        !./result_m0/model_connect.dat
nmodelfile         !1
polygonhead        !m0sens_m
modelfile          !m0sens.dat
ibound 0,1:l,2:u   !2
ubound             !0.2
nslice             !5
A                  !0.2588
B                  !-0.9659
D                  !0.6
boundary left      !-1.5
boundary right     !1.5
boundary bottom    !-0.25
boundary top       !1.75
A                  !0.2588
B                  !-0.9659
D                  !0.3
boundary left      !-1.5
boundary right     !1.5
boundary bottom    !-0.25
boundary top       !1.75
A                  !0.2588
B                  !-0.9659
D                  !0.0
boundary left      !-1.5
boundary right     !1.5
boundary bottom    !-0.25
boundary top       !1.75
A                  !0.2588
B                  !-0.9659
D                  !-0.3
boundary left      !-1.5
boundary right     !1.5
boundary bottom    !-0.25
boundary top       !1.75
A                  !0.2588
B                  !-0.9659
D                  !-0.6
boundary left      !-1.5
boundary right     !1.5
boundary bottom    !-0.25
boundary top       !1.75
