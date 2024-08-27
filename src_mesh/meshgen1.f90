! Coded on 2015.09.17
! to generate .geo file for simple 2D mesh with suitable .pos file
! from topo file around Nakadake
!
! Input file
! aso_caldera.dat
! topo127_134_29_36.xyz
!
! Output file
!  header.geo
!  header.pos
!
! Note x is eastward, y is northward, z is upward
program meshgen1
use param
use horizontalresolution
use water_level ! 2023.08.31
use topo_tool   ! 2023.08.31
implicit none
type(param_forward)          :: g_param
type(param_source)           :: s_param
type(coast_data),allocatable :: g_coast(:),h_coast(:)  ! 2023.08.31
type(poly_data)              :: i_poly    ! see m_topo_tool.f90
type(poly_data), allocatable :: g_poly(:),h_poly(:) ! 2023.09.04
type(bound_data)             :: g_bound    ! ncmax, zlabel(ncmax), ibelong(ncmax)
real(8),dimension(2)         :: xy1, xy2, xy, v12, v0, dv
integer(4)                   :: ifile=1, ndat,iline
integer(4)                   :: ipoint,i,ii,j,ishift, ispline, ilineloop
character(50)                :: posfile,header2d ! 2017.09.28
real(8)                      :: elev  ! size of triangles [km], elev is water level hight [km]
real(8),dimension(4,2)       :: lonlat4
integer(4),allocatable       :: line12(:,:)
type(param_water_level)      :: g_param_water ! 2023.08.31
integer(4)                   :: node0,ncmax,inum_water_level=0,inwl
integer(4),   dimension(0:3) :: iflg,ilflg ! 2023.09.01 copied from TMTGEMv1.3
type(grid_data)              :: g_grd
real(8),dimension(0:3)       :: xcorner,ycorner ! 2019.02.20
integer(4),parameter         :: lpmax=30 ! max total number of polygons 2024.08.27
! iflag_coast_required: true -> ocean-land boundary,flase -> lake boundary in land
logical                      :: iflag_coast_required(lpmax) ! 2023.08.31

![0]## read param
 call readparam(g_param,s_param) ! 2018.06.14

!#[1]# set
 header2d = g_param%header2d

!#[1]## [water level option]# 2023.08.31
!#[1-1]# Case for water_level.ctl exist
 open(1,file="water_level.ctl",status='old',err=100)
   call readwaterlevelctl(g_param,g_param_water,iunit=1)
   write(*,*) " 'water_level.ctl' exists!"
   write(*,*) " water_level_elev(1)",g_param_water%water_level_elev(1),"[km]"
 close(1)

 !#[1-1-1]# read topo # the followings are copied from coastline.f90 in TMTGEMv1.3 2023.08.31
   CALL topo2grdxyz(g_param,g_grd,g_param_water) ! 2023.08.31
   !# make lake outline
   ncmax = g_grd%node / 10 ! 2023.09.01
   inum_water_level = g_param_water%inum_water_level
   write(*,*) "inum_water_level =", inum_water_level
   inwl = inum_water_level ! 2024.08.27
   allocate(g_poly(inwl),h_poly(inwl),g_coast(inwl),h_coast(inwl)) ! sometime, inwl=0 2024.08.27

 do i=1,inum_water_level
   !#[1-1-2]## compose polygon for land and oceans
     elev = g_param_water%water_level_elev(i)
     write(*,*) "loop for water level =",elev
     CALL allocatecoast(g_coast(i),ncmax)               ! ncmax is max # of points on coastlines
     CALL makecoast5(g_grd,g_coast(i),elev)             ! make coast nodes and output coastfile
     call allocatepoly(g_poly(i),lpmax,ncmax,lpmax)     ! nclose for g_poly is unknow here
     call makepolygon7(g_coast(i),g_grd,g_poly(i),elev) ! lpoly0= # of polygons
     call outpolygon8( g_poly(i),g_coast(i),header2d,0) ! [for check] commented out on2019.03.07

   !#[1-1-3]## 
     call identifycraterlake(iflag_coast_required,g_grd,g_coast(i),g_poly(i),g_param_water,lpmax,ncmax)
     call landpolygon10(iflag_coast_required,g_poly(i),g_coast(i),lpmax) ! ind_r is modified
     !   call outpolygon8(g_poly,g_coast,header2d,1) ! [for check] commented out on 2019.03.07

   !#[1-1-4]## smoothen the coastlines by spline and background size information
     call smoothen10_5(g_poly(i),h_poly(i),g_coast(i),h_coast(i),g_param) ! h_poly is generated
     call deallocatecp(g_coast(i),h_coast(i),g_poly(i))

 end do
  write(*,*) "loop end"

 !#[1-2]# Case for no water lake 
  goto 101
  100 continue
   write(*,*) " 'water_level.ctl' doesn't exist!"
   write(*,*) " No lake surface will be generated"
  101 continue

  !#[2]## Finalize polygons
  call integratepoly(h_poly,g_param,inum_water_level,i_poly,ncmax) ! 2024.08.27
  ! call outpolygon8(i_poly,h_coast(1),header2d,1) ! [for check] commented out on 2019.03.07
   call outpolygon12(i_poly,header2d,1)           ! commented out 2018.11.09

  call outpolygeo13(i_poly,g_param,g_param_water)

! call avoidintersections(h_poly)                ! m_intersection.f90
! call avoidintersections_interpolygons(h_poly)  ! m_intersection.f90 2019.02.27
!  call outpolygon12(h_poly,header2d,2)           ! commented out 2018.11.09

!#[9]## add corner nodes and integrate nclose polygons
! call allocatebound(g_bound,ncmax)
! call addcorner11(h_poly,i_poly,g_grd,g_bound,iflg,ilflg,xcorner,ycorner) ! i_poly is generated
! call outpolygon12(i_poly,header2d,4) ! [for check] commented out on2019.03.07

![3]## generate .pos file
 posfile=trim(header2d)//".pos"
 call calobsr(s_param,g_param) ! calculate sigma_r, A_r, nobsr
 call outbgmesh2d(posfile,g_param) ! outbgmesh2d.f90 2017.09.28
! call outbgfield2d(g_param) ! outbgfield2d.f90 2017.09.28

end program meshgen1
!###########################################################
subroutine deallocatecp(g_coast,h_coast,g_poly) ! 2024.08.27
use topo_tool   ! 2023.08.31
implicit none
type(coast_data),intent(inout) :: g_coast,h_coast  ! 2023.08.31
type(poly_data), intent(inout) :: g_poly    ! see m_topo_tool.f90
     call deallocatecoast(g_coast)
     call deallocatecoast(h_coast)
     call deallocatepoly(g_poly)
return
end
!###########################################################
subroutine outpolygeo13(i_poly,g_param,g_param_water)
use param
use topo_tool
use water_level
implicit none
type(param_forward),    intent(in) :: g_param
type(param_water_level),intent(in) :: g_param_water
type(poly_data),        intent(in) :: i_poly
integer(4)    :: i,j,ii,jj,k,is,lpoly0,inum_water_level,lpmax,ncmax
character(70) :: geofile
character(50) :: header2d
character(2)  :: num
integer(4),allocatable,dimension(:) :: lpoly
real(8),allocatable,dimension(:,:) :: xpoly,ypoly
real(8) :: si,sb

 si   = g_param%sizein
 sb   = g_param%sizebo

header2d = g_param%header2d
lpoly0 = i_poly%lpoly0
lpoly  = i_poly%lpoly
lpmax  = i_poly%lpmax
ncmax  = i_poly%ncmax
allocate(xpoly(lpmax,ncmax),ypoly(lpmax,ncmax))
xpoly(1:lpmax,1:ncmax) = i_poly%xypoly(1,1:lpmax,1:ncmax)
ypoly(1:lpmax,1:ncmax) = i_poly%xypoly(2,1:lpmax,1:ncmax)
inum_water_level = g_param_water%inum_water_level

!# topo/lake surface geo file
geofile=trim(header2d)//".geo"
 write(*,*) "Start write geo file", trim(geofile)
 open(1,file=geofile,status='replace')
 jj=0
 write(1,*) "lc=",sb,";" ! lc is the characteristic length for element size

!#[1]# generate points
 do i=1,lpoly0
   do j=1,lpoly(i)
     jj=jj+1
     write(1,*)"Point(",jj,")={",xpoly(i,j),",",ypoly(i,j),",0.0,lc} ;"
   end do
 end do
 !# corner points

!# [2] ### generate lines
 jj=0;is=0
 do i=1,lpoly0
  write(*,*) "i",i,"is",is,"lpoly(i)",lpoly(i)
  is=is+lpoly(i)
  write(*,*) "i",i,"is",is,"lpoly(i)",lpoly(i)
  do j=1,lpoly(i)-1
   jj=jj+1
   write(1,*) "Line(",jj,")={",jj,",",jj+1,"} ;"
  end do
  jj=jj+1
  write(1,*) "Line(",jj,")={",jj,",",is-lpoly(i)+1,"} ;" ! close the polygon
 end do

 !# [3] ## Line loop
  is = 0
  do i=1,lpoly0
   is=is+lpoly(i);k=is-lpoly(i)
   write(1,*) "Line Loop(",i,")={",(k+j,",",j=1,lpoly(i)-1),k+lpoly(i),"} ;"  
  end do

 !# [4] ## Plane surface
  do i=1,lpoly0 -1
   write(1,*) "Plane Surface(",i,")={",i,",",i+1,"};" ! unclosed part
  end do
   write(1,*) "Plane Surface(",lpoly0,")={",lpoly0,"};" ! lake bottom

   do i=1,inum_water_level
    jj = i + inum_water_level + 1
    write(1,*) "Plane Surface(",jj,")={",i+1,"};" ! lake surface
   end do
close(1)

write(*,*) "### OUTPOLYGEO13 END!! ###" ! 2024.08.27

return
end
!###########################################################
subroutine integratepoly(h_poly,g_param,inum_water_level,i_poly,ncmax)
use topo_tool
use param
implicit none
type(param_forward), intent(in)    :: g_param
type(poly_data),     intent(inout) :: h_poly(inum_water_level) ! sometimes inum_water_level=0 
integer(4),          intent(in)    :: inum_water_level,ncmax
type(poly_data),     intent(out)   :: i_poly
integer(4) :: i,nclose,lpmax,jj
real(8)    :: x1,x2,y1,y2,x01,x02,y01,y02, x0,y0, z0,lc

 x01 = g_param%xbound(1)
 x1  = g_param%xbound(2)
 x2  = g_param%xbound(3)
 x02 = g_param%xbound(4)
 y01 = g_param%ybound(1)
 y1  = g_param%ybound(2)
 y2  = g_param%ybound(3)
 y02 = g_param%ybound(4)
 lc  = g_param%sizebo

nclose = 0 ! # number of unclosed polygon
lpmax = inum_water_level + 1

i_poly%lpmax  = inum_water_level + 1
i_poly%ncmax  = ncmax
i_poly%lpoly0 = inum_water_level + 1 ! # of polygon including calculation boundary polygon
i_poly%nclose = 0 

!# allocate output polygon i_poly
allocate(i_poly%lpoly(lpmax))
allocate(i_poly%xypoly(2,lpmax,ncmax))

!# calculation boundary polygon
 i_poly%lpoly(1) = 4 ! land polygon
 i_poly%xypoly(1:2,1,1) = (/x02,y02/) ! top right
 i_poly%xypoly(1:2,1,2) = (/x02,y01/) ! top left
 i_poly%xypoly(1:2,1,3) = (/x01,y01/) ! bottom left
 i_poly%xypoly(1:2,1,4) = (/x01,y02/) ! bottom right

!# 2 to (inum_water_level + 1) polygons
do i=1,inum_water_level   ! skipped when inum_water_level = 0
  i_poly%lpoly(i+1) = h_poly(i)%lpoly(1)
  i_poly%xypoly(1:2,i+1,1:ncmax) = h_poly(i)%xypoly(1:2,1,1:ncmax)
end do

!# deallocate h_poly     
do i=1,inum_water_level ! skipped when inum_water_level = 0
  call deallocatepoly(h_poly(i))
end do
write(*,*) "### INTEGRATE POLY END!! ###" ! 2024.08.27

return
end
!########################################################### smoothen10_5
subroutine smoothen10_5(g_poly,h_poly,g_coast,h_coast,g_param)
use topo_tool ! 2023.09.01
use param
implicit none
type(param_forward),   intent(in)       :: g_param
type(coast_data),      intent(inout)    :: g_coast
type(poly_data),       intent(in)       :: g_poly
type(poly_data),       intent(out)      :: h_poly
type(coast_data),      intent(out)      :: h_coast
integer(4)                              :: lpoly0,nclose,nclose2, ncmax
integer(4),allocatable,dimension(:)     :: lpoly,lpoly2
integer(4),allocatable,dimension(:,:,:) :: ind2
real(8),   allocatable,dimension(:,:)   :: xpoly,ypoly,xpoly2,ypoly2
integer(4)                              :: i,j,k,lpn,lpn1,lp,lp1,lp0,lp0n
integer(4)                              :: iclose,lpoly0pre,inum
integer(4)                              :: ii,iind,ncoast
! iind is the accumulated number of nodes for ind modify
real(8)                                 :: space,tn,dist,prespace,x1,y1,alpha
real(8),   allocatable,dimension(:)     :: t,X,X2,Y,Y2
real(8),   allocatable,dimension(:)     :: tnn
real(8)                                 :: space_lake_interval=0.004
integer(4),allocatable,dimension(:,:)   :: ind
real(8),   allocatable,dimension(:,:)   :: loc ! 2019.02.25
logical                                 :: iflag_topright_land ! 2019.02.28
write(*,*) "smoothen start!"

!#[1]## set input
  ncoast = g_coast%ncoast
  ncmax  = g_coast%ncmax
  nclose = g_poly%nclose
  lpoly0 = g_poly%lpoly0
  write(*,*) "lpoly0",lpoly0,"ncmax",ncmax
  allocate( tnn(ncoast+1000),ind(ncmax,2))! 2021.11.08 +1000 is added
  allocate( lpoly(lpoly0), lpoly2(lpoly0))
  allocate( xpoly( lpoly0,ncmax), ypoly( lpoly0,ncmax))
  allocate( xpoly2(lpoly0,ncmax), ypoly2(lpoly0,ncmax))
  ind(1:ncmax,1:2) = g_coast%ind_r(1:ncmax,1:2)
  xpoly(1:lpoly0,1:ncmax)  = g_poly%xypoly(2,1:lpoly0,1:ncmax)
  ypoly(1:lpoly0,1:ncmax)  = g_poly%xypoly(1,1:lpoly0,1:ncmax)
  lpoly(1:lpoly0)  = g_poly%lpoly(1:lpoly0)
  if ( nclose .ge. 1) then       ! 2019.03.06
   allocate( ind2(nclose,2,2) )  ! 2019.03.06
   allocate( loc(nclose,2) )     ! 2019.02.25
   loc    = g_poly%loc           ! 2019.02.25
   ind2(:,:,:)=0.d0              ! 2019.03.06
  end if
  iflag_topright_land = g_poly%iflag_topright_land ! 2019.02.28

write(*,*) "2"
!#
! scale : length from the origin that permits sparce meshes
 xpoly2(:,:)=0.d0 ; ypoly2(:,:)=0.d0 ; lpoly2(:)=0
 lp0     = lpoly0
 lp0n    = 0      ! new lpoly0, which excludes small islands
 nclose2 = nclose ! # of unclosed polygons
 iind    = 0      ! iind is for accumulated lpoly
 do i=1,lp0 ! polygon loop start!
  iind = iind + lpoly(i)
  !write(*,*) "investigate",i,"-th polygon"
  !write(*,'(a,i6,a,i6)') ("ind(1)=",ind(k,1),"ind(k,2)",ind(k,2),k=iind-lpoly(i)+1,iind)
  !do i=1,5 ! polygon loop start!
  !  write(*,*) "i=",i,"loop start! lpoly(i)=",lpoly(i)
  lp0n = lp0n+1
  lp   = lpoly(i)
  lp1  = lp ! if polygon is not closed, lp1 should be lp
  if ( i .gt. nclose ) lp1=lp+1  ! if pooygon is closed, count the start node twice and refer it to both start and end
  allocate(t(lp1),Y(lp1),Y2(lp1),X(lp1),X2(lp1))
  X(1:lp)=xpoly(i,1:lp)
  Y(1:lp)=ypoly(i,1:lp)
  if (i .gt. nclose) then ! when the polygon is closed
    X(lp+1)=X(1)
    Y(lp+1)=Y(1)
  end if
  t(:)=0.d0
  t(1)=0.d0
  do j=2,lp1  ! set t as a distance
    t(j)=t(j-1)+dsqrt((X(j)-X(j-1))**2.d0+(Y(j)-Y(j-1))**2.d0)
  end do
  call SPLINE(t, X,lp1,X2) !make spline coefficients
  call SPLINE(t, Y,lp1,Y2)
  !  write(*,'(a2,f15.7,1x,a4,2f15.7)')("t=",t(j),"X,Y=",X(j),Y(j),j=1,lp1)
  lpn1=1
  tnn(1)=0.d0
  ! make tnn
  do while (tnn(lpn1) .lt. t(lp1))
    lpn1=lpn1+1
!    if (ncoast .lt. lpn1) then
!      write(*,*) "GEGEGE ncoast=",ncoast,"is less than lpn1=",lpn1 ! 2021.11.07
!      stop
!    end if
    call SPLINT(t,X,X2,lp1,tnn(lpn1-1),x1) ! interpolate
    call SPLINT(t,Y,Y2,lp1,tnn(lpn1-1),y1) ! interpolate
!    call sizebgele(x1,y1,pos,space)
    space = space_lake_interval
    tnn(lpn1)=tnn(lpn1-1)+space ! smaller node spacing on coastlines against the bgmesh lead to not adding extra nodes on coastlines
    !  if (lpn1 .eq. 2) tnn(lpn1)=space*0.7 ! in order to constrain the size of polygon in the vicinity of calculation boundaries, commented out on Oct. 25, 2015
    if (lpn1 .eq. 2 .and. i .le. nclose) tnn(lpn1)=space*0.7 ! inserted on Oct. 25, 2015
  end do
  alpha=t(lp1)/tnn(lpn1) ! shrink the length of tnn into that of t, delete comment out on Oct. 24, 2015
  tnn(1:lpn1)=alpha*tnn(1:lpn1) ! tnn is created, delete comment out on Oct. 24, 2015
  !tnn(lpn1)=t(lp1)   !    comment out on Oct. 24, 2015
  !tnn(lpn1-1)=(tnn(lpn1-2)+tnn(lpn1))/2.d0  ! comment out on Oct. 24, 2015
  lpn=lpn1                                 ! when the polygon is not closed, lpn should be lpn1
  if ( i .gt. nclose ) lpn=lpn1-1 ! when the polygon is closed, lpn should exclude the last point (identical to start node)
  lpoly2(lp0n)=lpn
  ! if (lpn .le. 3) then ! Do not add islands consisting of less than 4 points, commented out on Oct. 25, 2015
  if (lpn .le. 4) then ! Do not add islands consisting of less than 4 points, inserted on Oct. 25, 2015
    xpoly2(lp0n,1:lpn)=0.d0
    ypoly2(lp0n,1:lpn)=0.d0
    lpoly2(lp0n)=0
    lp0n=lp0n-1
    if ( i .le. nclose )   nclose2=nclose2-1 ! reduce the # of unclosed polygon
  else
    ! make xpoly2 and ypoly2
    !write(*,*) "poygon # i=",i,"alpha=",alpha
    do j=1,lpn
      call SPLINT(t,X,X2,lp1,tnn(j),xpoly2(lp0n,j)) ! interpolate
      call SPLINT(t,Y,Y2,lp1,tnn(j),ypoly2(lp0n,j))
    end do
    if ( i .le. nclose ) then !! ASSIGN ind if ( i .le. nclose), here is skipped when nclose = 0
      ! between starting and end nodes, we cannot define ind.
      ind2(lp0n,1,:)=ind(iind-lpoly(i)+1,:) ! starting node of lp0n-th polygon
      ind2(lp0n,2,:)=ind(iind,:) ! and end node of unclosed polygons, we can still assign ind()
      loc(lp0n,1:2) = loc(i,1:2) ! 2019.02.25
      write(*,*) "lp0n=",lp0n
      !   write(*,*) "ind2(lp0n,1,1-2)=",ind2(lp0n,1,1),ind2(lp0n,1,2)
      !   write(*,*) "ind2(lp0n,1,1-2)=",ind2(lp0n,2,1),ind2(lp0n,2,2)
    end if
    !      write(*,*) "lpn=",lpn,"is too small to create new polygon."
  end if
  deallocate(t,X,X2,Y,Y2)
end do  ! i loop end

!#[3]## set output h_poly and h_coast
!#[3-1]## h_poly
 call allocatepoly(h_poly,lp0n+1,g_poly%ncmax,nclose2) ! +1 is for when nclose = 0
 write(*,*) "check3 lp0n=",lp0n
 h_poly%lpoly0              = lp0n
 h_poly%xypoly(2,1:lp0n,:)  = xpoly2(1:lp0n,1:g_poly%ncmax)
 h_poly%xypoly(1,1:lp0n,:)  = ypoly2(1:lp0n,1:g_poly%ncmax)
 h_poly%lpoly(1:lp0n)       = lpoly2(1:lp0n)
 h_poly%nclose              = nclose2
 h_poly%iflag_topright_land = iflag_topright_land  ! 2019.03.07
 if ( nclose2 .ge. 1 ) then ! 2019.03.06
  h_poly%loc(1:nclose2,1:2)  = loc(1:nclose2,1:2)  ! 2019.02.25
  h_poly%ind2(1:nclose2,:,:) = ind2(1:nclose2,:,:) ! nclose2 < nclose
 end if                     ! 2019.03.06
 lpoly0pre                  = lpoly0

!#[3-2]## h_coast
 ncoast=0
 do i=1,lp0n
  ncoast = ncoast + h_poly%lpoly(i)
 end do
 call allocatecoast(h_coast,ncoast)
 h_coast%ncoast = ncoast
 write(*,'(a,i8)') " # of points on coast is", h_coast%ncoast
 ii=0
 do i=1,lp0n
  do j=1,h_poly%lpoly(i)
   ii=ii+1
   h_coast%cxy(1:2,ii)=h_poly%xypoly(1:2,i,j)
  end do
 end do

 write(*,*) "new lpoly0 is reduced to",lp0n, "from",lpoly0pre
 write(*,*) "nclose",nclose,"is reduced to nclose2",nclose2

 if ( iflag_topright_land .and. nclose2 .eq. 0 ) then ! 2019.02.28
  write(*,*) "GEGEGE! There are no oceans!"
  write(*,'(a,l3)') " iflag_topright_land",iflag_topright_land
  stop
 end if
 write(*,*) "### smoothen10_5 end!###"

 return
end subroutine smoothen10_5


!########################################################### identifycraterlakepolygon 2023.09.01
subroutine identifycraterlake(iflag_coast_required,g_grd,g_coast,g_poly,g_param_water,lpmax,ncmax)
  use topo_tool
  use param
  use water_level
  implicit none
  logical,                intent(out)    :: iflag_coast_required(lpmax) ! 2021.05.29
  integer(4),             intent(in)     :: lpmax,ncmax
  type(coast_data),       intent(in)     :: g_coast
  type(grid_data),        intent(in)     :: g_grd
  type(poly_data),        intent(in)     :: g_poly
  type(param_water_level),intent(in)     :: g_param_water
  real(8),allocatable,    dimension(:)   :: lon_grd,lat_grd
  integer(4),allocatable                 :: ind(:,:), lpoly(:)
  real(8)                                :: cxy(2),clonlat(2)
  integer(4)                             :: ic,jc,i,j,k,nsouth,neast,node,lpoly0,nclose,pid,pid_tmp1,pid_tmp2
  integer(4),allocatable,dimension(:,:)  :: iflgclosedpolygon_pid

  !#[1]# set
  iflag_coast_required = .false.    ! initial value
  node        = g_grd%node         ! node = neast * nsouth 2019.02.25
  cxy         = g_param_water%xy_lake_center(1:2)
  clonlat     = g_param_water%lonlat_lake_center(1:2)
  neast       = g_grd%neast
  nsouth      = g_grd%nsouth        ! 2019.02.25
  ind         = g_coast%ind_r
  nclose      = g_poly%nclose      ! 2019.02.25 # of polygons reaching boundaries
  lpoly       = g_poly%lpoly
  lpoly0      = g_poly%lpoly0
  write(*,*) "neast",neast
  write(*,*) "nsouth",nsouth
  write(*,*) "size lon",size(g_grd%lon)
  write(*,*) "size lat", size(g_grd%lat)
  write(*,*) "clonlat",clonlat
  lon_grd     = g_grd%lon(1:neast)  ! 2023.09.01  
  lat_grd     = g_grd%lat(1:nsouth) ! 2023.09.01

  if (.false.) then
   open(1,file="cxy.dat")
     do i=1,g_coast%ncoast
       write(1,*) g_grd%xyz(1:2,g_coast%ind_r(i,1)),i 
     end do
   close(1)
  end if

  !#[2]# look for ic,jc for 
  do i=1,neast-1
   if ( lon_grd(i) .le. clonlat(1) .and. clonlat(1) .lt. lon_grd(i+1) ) ic =i
  end do
  do i=1,nsouth-1
   !    write(*,*) i,"lat i+1,i", lat_grd(i+1),lat_grd(i)
    if ( lat_grd(i+1) .le. clonlat(2) .and. clonlat(2) .lt. lat_grd(i) ) jc =i
   end do
   write(*,*) "ic=",ic
   write(*,*) "jc=",jc

   !# [3] ## generate iflagclosedpolygon_pid
   allocate( iflgclosedpolygon_pid(node,2) ) ! 2019.02.27
   iflgclosedpolygon_pid = 0 ! initialize
   if ( nclose .eq. 0 ) then
     k=0
   else
    write(*,*) "GEGEGE nclose is not zero. THis program does not fit"
    stop
   end if
   do i=nclose+1,lpoly0
     do j=1,lpoly(i)
       k=k+1
       if (ind(k,2) .eq. 1)  iflgclosedpolygon_pid(ind(k,1),1)=i  ! i-th polygon exist east of ind(k,1) node 
       if (ind(k,2) .eq. 2)  iflgclosedpolygon_pid(ind(k,1),2)=i  ! i-th polygon exist south of ind(k,1) node 
     end do
   end do

  !#[4]## look for the closest polygon to (ic,jc)
  !#                      ic   ic +1
  !#             jc     >  .     .
   !#    lake center -  >     *　
  !#             jc + 1 >  .     .
  !# only search eastward
  do i=ic,neast
    pid_tmp1 = iflgclosedpolygon_pid( neast*(jc-1) + i,1)
    pid_tmp2 = iflgclosedpolygon_pid( neast*(jc-1) + i,2)
!    write(*,*) "i,j",i,j,"lon(i)",lon_grd(i),"pid_tmp1",pid_tmp1 
!    write(*,*) "i,j",i,j,"lon(i)",lon_grd(i),"pid_tmp2",pid_tmp2 
    if ( pid_tmp2 .ne. 0 .or. pid_tmp1 .ne. 0 ) then
      if ( pid_tmp1 .ne. 0 ) pid = pid_tmp1 
      if ( pid_tmp2 .ne. 0 ) pid = pid_tmp2 
      goto 100
    end if
  end do
  100 continue
  write(*,*) "Corresponding lake polygon is pid=",pid
  iflag_coast_required(pid) = .true.

  return
  end

  !########################################################### judgeocean9
 subroutine judgeocean9(iflag_coast_required,g_coast,g_poly,g_grd,lpmax,g_param) ! 2019.02.28
 !###        judge if polygons are in land or oceans        ###
 ! iflgsealand : 0 -> not determined, 1 -> in ocean, 2 -> in land
 ! iflgcoastpresence : 0 -> no coastline, i -> i-th polygon coastline exists
 ! iflgcoastpresence(i,1) : right side of the node, (i,2) : downside of the node
 use topo_tool ! 2023.08.31
 use param     ! 2023.08.31
 implicit none
 type(param_forward),intent(in)     :: g_param                     ! 2023.08.31
 type(poly_data),    intent(inout)  :: g_poly
 type(coast_data),   intent(inout)  :: g_coast
 type(grid_data),    intent(in)     :: g_grd
 integer(4),         intent(in)     :: lpmax
 logical,            intent(out)    :: iflag_coast_required(lpmax) ! 2021.05.29
 integer(4)                       :: irow,icolm,nclose           ! 2019.02.25
 integer(4),       allocatable    :: lpoly(:)
 integer(4),       allocatable    :: ind(:,:)
 logical,          allocatable    :: iflgcoastpresence(:,:)      ! 2019.02.25
 integer,          allocatable    :: iflgclosedpolygon_pid(:,:)  ! 2019.02.27
 logical                          :: iflag_topright_land         ! 2019.02.25
 real(8),          allocatable    :: loc(:,:)                    ! 2019.02.25
 integer(4)                       :: ln(2)                       ! 2019.02.25
 real(8)                          :: loc_current, loc_end
 real(8),          allocatable    :: xypoly(:,:,:)               ! 2019.02.28
 logical,          allocatable    :: cross_closed_once(:)        ! 2019.02.27
 logical                          :: cross_closed_once_exist     ! 2019.02.27
 real(8)                          :: len_threshold,len_coast     ! 2019.02.28
 real(8)                          :: xbound(4),ybound(4),dl      ! 2019.02.28
 real(8),          allocatable    :: len_poly(:)                 ! 2019.02.28
 integer(4) :: k,i,j,istart,ii,jj,lpoly0,ncmax,neast,node,nsouth,pid
 integer(4) :: ishift ! 2019.02.28

!#[0]## set input
 allocate(lpoly(lpmax))          ! 2018.08.28
 node       = g_grd%node         ! node = neast * nsouth 2019.02.25
 ncmax      = g_coast%ncmax      ! 2018.08.28
 allocate( ind(ncmax,2))         ! 2018.08.28
 neast      = g_grd%neast
 nsouth     = g_grd%nsouth       ! 2019.02.25
 ind        = g_coast%ind_r
 lpoly      = g_poly%lpoly
 lpoly0     = g_poly%lpoly0
 nclose     = g_poly%nclose      ! 2019.02.25 # of polygons reaching boundaries
 xbound     = g_param%xbound
 ybound     = g_param%ybound
 len_threshold = min(xbound(4)-xbound(1),ybound(4)-ybound(1)) ! 2019.02.28
 write(*,'(a,f15.7,a)') "len_threshold=",len_threshold,"[km]" ! 2022.03.01
 iflag_topright_land = g_grd%iflag_topright_land ! 2019.02.25

 xypoly     = g_poly%xypoly ! allocate and fill xypoly(2,lpmax,ncmax) 2022.03.01
 iflag_coast_required(1:lpmax) = .false. ! all are not reuquired initially 2019.02.28
 write(*,*) "## JUDGEOCEAN9 START!1 ##" ! 2019.03.06

 if ( nclose .eq. 0 ) then   !== 2019.03.06 ================================================= nclose = 0 ===========
     write(*,*) "## nclose = 0 and no polygons reaching calculation boundaries ##" ! 2019.03.06

 else                        !== 2019.03.06 ================================================= nclose >= 1 ==========
     !# [1] ### set loc(1:nclose,1:2)
     allocate(loc(nclose,2)   )     ! 2018.08.28
     loc(:,:) = 0.d0
     ii       = 0
     do i=1,nclose ! nclose is # of unclosed polygons
         ii=ii+lpoly(i)    ! ii is accumulated number of unclosed polygon points, Oct. 24, 2015
         ln(1)=ii-lpoly(i)+1; ln(2)=ii ! ln(1) is the first node id, ln(2) is the end node id, for i-th polygon
         do j=1,2  ! set loc for the start and end points of i-th polygons
             jj = ind(ln(j),1) ! jj is node id in topofile for ln(1)
             !           icolm      order in the topo file (Gebco)
             !       ------|-----     -------->
             !       |     v    |   | -------->
             ! irow --->  jj    |   | -------->
             !       |          |   | -------->
             !       ------------   v -------->
             irow  = (jj - 1)/neast + 1    ! 2019.02.20
             icolm = jj - (irow - 1)*neast ! 2019.02.20
             ! set loc, location on the boundary
             if ( icolm .eq. neast )                      loc(i,j)=(irow -1 + 0.5 )/(nsouth-1)                   ! right  0< loc<1
             if ( irow  .eq. nsouth )                     loc(i,j)=((neast -1.) - (icolm - 1.+0.5))/(neast-1)+1. ! bottom 1<loc<2
             if ( icolm .eq. 1 .and. ind(ln(j),2) .eq. 2) loc(i,j)=((nsouth-1.) - (irow  - 1.+0.5))/(nsouth-1)+2.! left   2<loc<3
             if ( irow  .eq. 1 .and. ind(ln(j),2) .eq. 1) loc(i,j)=(icolm -1+0.5)/(neast-1.)+3.                  ! top    3<loc<4
         end do
         if (loc(i,1) .eq. 0.d0 .or. loc(i,2) .eq. 0.d0) goto 999
         write(*,*) "polygon #", i, "loc(i,1)=",loc(i,1)
         write(*,*) "polygon #", i, "loc(i,2)=",loc(i,2)
      end do
      write(*,'(a)') "## set loc end for nclose polygon!! in judgeocean9 ##"

      !# [2] ## set len_poly(1:lpoly0) 2019.02.28
      allocate(len_poly(lpoly0))
      len_poly = 0.d0
      do i=1,lpoly0
          do j=1,lpoly(i) - 1
             dl = sqrt((xypoly(1,i,j+1) - xypoly(1,i,j))**2. + (xypoly(2,i,j+1) - xypoly(2,i,j))**2. )
             len_poly(i)=len_poly(i)+ dl
          end do
      end do
      write(*,'(a)') "## set len_poly end for all polygons!! in judgeocean9 ##"
   
      !# [3] ## set iflag_coast_required(1:nclose) by deleaneating ocean 2019.02.25
      istart = 1
      15 continue
      loc_current = 0.d0
      loc_end     = 4.0
      len_coast   = 0.d0 ! km
      write(*,'(a,l3)') " iflag_topright_land",iflag_topright_land
      do i=istart,nclose ! unclosed polygons
         !   write(*,*) "loc_current",loc_current,"loc_end",loc_end
         if ( i .eq. istart .and. iflag_topright_land ) then
             loc_end     = loc(i,2)
             iflag_coast_required(i)=.true.
             loc_current = loc(i,1)
             len_coast   = len_poly(1)
             goto 10
         end if
         if ( loc(i,1) .gt. loc_current .and. loc(i,1) .lt. loc_end) then
             iflag_coast_required(i)=.true.
             loc_current = loc(i,2)
             len_coast   = len_coast + len_poly(i)
         end if
         10 continue
         !  write(*,'(a,i,a,l3)') "polygon #",i," iflag_coast_required =",iflag_coast_required(i) ! 2019.02.27
      end do
      !# check coastline is longer than len_threshold 2019.02.28
      if ( len_coast .lt. len_threshold ) then
          iflag_coast_required = .false.
          istart = istart + 1
          if ( istart .gt. nclose ) then
             write(*,*) "no coastline is identified in nclose polygons with len_threshold",len_threshold
             goto 16
          end if
          goto 15
      end if
      16 continue !2019.02.28
      write(*,*) "istart/nclose",istart,"/",nclose
   
      !#[4]## reverse the node order for first boundary polygon when the top right corner is in land 2019.02.28
      if ( iflag_topright_land ) then
         xypoly(1:2,istart,1:lpoly(istart)) = xypoly(1:2,istart,lpoly(istart):1:-1)
         ishift=0
         do i=1,istart-1
              ishift = ishift + lpoly(i)
          end do
         ind(ishift+1:ishift+lpoly(istart),1:2) = ind(ishift+lpoly(istart):ishift+1:-1,1:2)
         loc(istart,1:2) = loc(istart,2:1:-1)
      end if
  end if ! nclose >= 1 end ! 2019.03.06

!# [3] ## set iflag_coast_required(nclose+1:lpoly0)

 !# [3-1] ## set iflgcoastpresence
 allocate( iflgcoastpresence(node,2) ) ! 2018.08.28
 iflgcoastpresence = .false.           ! initial value
 k = 0
 do i=1,nclose                         ! this loop is skipped when nclose = 0
   do j=1,lpoly(i)
     k=k+1
     if ( iflag_coast_required(i) ) then ! when polygons are ocean-land boundary
       if (ind(k,2) .eq. 1) iflgcoastpresence(ind(k,1),1)=.true. ! right side exist
       if (ind(k,2) .eq. 2) iflgcoastpresence(ind(k,1),2)=.true. ! down  side exist
     end if
   end do
 end do

 !# [3-2] ## closed polygons
 allocate( iflgclosedpolygon_pid(node,2) ) ! 2019.02.27
 iflgclosedpolygon_pid = 0
 do i=nclose+1,lpoly0
  do j=1,lpoly(i)
   k=k+1
    if (ind(k,2) .eq. 1)  iflgclosedpolygon_pid(ind(k,1),1)=i  ! 2019.02.27
    if (ind(k,2) .eq. 2)  iflgclosedpolygon_pid(ind(k,1),2)=i  ! 2019.02.27
  end do
 end do

!# [3-3] ## set iflag_coast_required
 j=0
 allocate( cross_closed_once(lpoly0) ) ! 2019.02.27
 do i=1,nclose      ! this loop is skipped when nclose = 0
  j = j + lpoly(i)
 end do
 do i=nclose + 1,lpoly0                 ! i polygon loop for closed polygons
  cross_closed_once = .false.           ! 2019.02.27
  k=0
  j=j+lpoly(i)
  istart = ind(j,1) ! node id in original gebco file for j-th coastline node id
  irow   = (istart - 1)/neast + 1       ! 2019.02.20
  icolm  = istart - (irow - 1)*neast    ! 2019.02.20

  !# upward check between ii and ii + 1 th row
  do ii=irow,1,-1 ! upward ; begin with jj = istart
   jj=neast*(ii-1)+icolm  ! 2019.02.20
   if ( iflgcoastpresence(jj,2) ) then
    k=k+1 ! how many times coastlines are crossed
!    write(*,*) "k",k,"cross ii -",ii
   end if
   if ( iflgclosedpolygon_pid(jj,2) .gt. 0 .and. iflgclosedpolygon_pid(jj,2) .ne. i ) then ! 2019.02.27
    pid = iflgclosedpolygon_pid(jj,2)
    if ( cross_closed_once(pid) ) then
     cross_closed_once(pid) = .false. ! second cross
    else
     cross_closed_once(pid) = .true.  ! first cross
    end if
   end if
  end do ! ii loop end

  !# rightward check between jj and jj + 1 th row
  do jj=icolm,neast-1               ! rightward 2019.02.20
   if ( iflgcoastpresence(jj,1) ) then
    k=k+1
!    write(*,*) "k",k,"cross jj |",jj
   end if
   if ( iflgclosedpolygon_pid(jj,1) .gt. 0 .and. iflgclosedpolygon_pid(jj,2) .ne. i ) then ! 2019.02.27
    pid = iflgclosedpolygon_pid(jj,1)
    if ( cross_closed_once(pid) ) then
     cross_closed_once(pid) = .false. ! second cross
    else
     cross_closed_once(pid) = .true.  ! first cross
    end if
   end if
  end do ! jj loop end

  !#  if top right is in ocean even k survive   2019.02.25
  !#  if top right is in land  odd  k survive
  if ( mod(k,2) .eq. 0 ) then ! k is even
   if ( .not. iflag_topright_land ) iflag_coast_required(i) = .true.
  else                        ! k is odd
   if (       iflag_topright_land ) iflag_coast_required(i) = .true.
  end if
  !# check whether a closed polygon is whithin another closed polygons
  cross_closed_once_exist=.false.
  do jj=nclose+1,lpoly0
   if ( cross_closed_once(jj) ) cross_closed_once_exist=.true.
  end do
  if ( cross_closed_once_exist ) iflag_coast_required(i) = .false.

  !write(*,'(a,i5,a,l3,a,i5,a,2i5)') "polygon #",i," iflag_coast_required(i)",iflag_coast_required(i)," k",k," irow,icolm",irow,icolm

 end do ! i polygon loop end

 do i=1,nclose ! skipped when nclose = 0
   write(*,'(i5,a,l3,a,2g15.7)') i,"iflag_coast_required(i)",iflag_coast_required(i)," loc(i,1:2)",loc(i,1:2)
 end do


 !# set loc to h_poly 2019.02.25
 if ( nclose .ge. 1 ) g_poly%loc = loc ! 2019.03.06
 g_poly%xypoly  = xypoly      ! 2019.02.28
 g_coast%ind_r  = ind         ! 2019.02.28

 write(*,*)"### judgeocean9 end!  ###"

return

999 continue ! moved to here on Feb 20, 2019
    write(*,*) "GEGEGE loc is not set!!"
    write(*,*) "i=",i,"loc(i,1)=",loc(i,1),"loc(i,2)=",loc(i,2)
    write(*,*) "ln(1)=",ln(1),"ln(2)=",ln(2)
    write(*,*) "ind(ln(1),1)=",ind(ln(1),1),"ind(ln(2),1)",ind(ln(2),1)
    write(*,*) "ind(ln(1),2)=",ind(ln(1),2),"ind(ln(2),2)",ind(ln(2),2)
    write(*,*) "lpoly(i)=",lpoly(i)
    stop

end subroutine judgeocean9

!########################################################### landpolygon10
subroutine landpolygon10(iflag_coast_required,g_poly,g_coast,lpmax)  ! 2019.02.25
!# eliminate polygons with iflag_coast_required = .false.
use topo_tool ! 2023.08.31
implicit none
type(coast_data),      intent(inout)    :: g_coast
type(poly_data),       intent(inout)    :: g_poly
integer(4),            intent(in)       :: lpmax
logical,               intent(in)       :: iflag_coast_required(lpmax)      ! 2019.02.25
integer(4)                              :: lpoly0, lpoly0n, nboundn, ncoastn
integer(4)                              :: ncoast,nclose,nbound,ncmax, nclosen ! 2019.02.25
integer(4)                              :: i,ii,jj,k,kk
integer(4),allocatable,dimension(:,:)   :: ind,indn
integer(4),            dimension(lpmax) :: lpolyn,lpoly
real(8),   allocatable,dimension(:,:)   :: xpoly, ypoly, xpolyn, ypolyn
real(8),   allocatable,dimension(:,:)   :: loc,locn                         ! 2019.02.25

!#[1]## set input
 ncmax  = g_coast%ncmax ! 2018.08.28
 allocate( ind(ncmax,2), indn(ncmax,2))              ! 2018.08.28
 allocate( xpoly( lpmax,ncmax), ypoly( lpmax,ncmax)) ! 2018.08.28
 allocate( xpolyn(lpmax,ncmax), ypolyn(lpmax,ncmax)) ! 2018.08.28
 ncoast = g_coast%ncoast ! total # of coastline nodes
 nbound = g_coast%nbound ! # of coastline nodes at calculation boundaries
 ind    = g_coast%ind_r
 lpoly0 = g_poly%lpoly0
 lpoly  = g_poly%lpoly
 xpoly  = g_poly%xypoly(2,:,:) ! north
 ypoly  = g_poly%xypoly(1,:,:) ! east
 nclose = nbound/2 ! # of unclosed polygons, which touch calculation boundaries 2019.02.25
 if (nclose .ge. 1 ) then ! 2019.03.06
  allocate( loc(nclose,2), locn(nclose,2) )
  loc    = g_poly%loc     ! 2019.02.25
 end if                   ! 2019.03.06

!#[2]##
 xpolyn(:,:)=0.d0;ypolyn(:,:)=0.d0
 indn(:,:)=0;lpolyn(:)=0;ncoastn=0
 ii=0;jj=0;kk=0
 nboundn = 0
 nclosen = 0   ! 2019.02.25
 write(*,*) "lpoly0=",lpoly0
 do i=1,lpoly0
  if ( iflag_coast_required(i) ) then ! 2019.02.25
   jj=jj+1
   xpolyn(jj,1:lpoly(i))      = xpoly(i,1:lpoly(i))
   ypolyn(jj,1:lpoly(i))      = ypoly(i,1:lpoly(i))
   indn(ii+1:ii+lpoly(i),1:2) = ind(kk+1:kk+lpoly(i),1:2)
   lpolyn(jj)    = lpoly(i)
   lpoly0n       = jj
   ncoastn       = ncoastn + lpoly(i)
   ii            = ii + lpoly(i)
   if ( i .le. nclose ) then      ! 2019.02.25
     nclosen      = nclosen + 1
     nboundn      = nboundn + 2
     locn(jj,1:2) = loc(i,1:2) ! 2019.02.25
   end if
  end if
  kk = kk + lpoly(i)
 end do

! set new polygon data
 g_coast%ind_r(:,:)   = indn(:,:)
 g_coast%nbound       = nboundn
 g_coast%ncoast       = ncoastn
!
 g_poly%lpoly0        = lpoly0n
 g_poly%xypoly(2,:,:) = xpolyn(:,:) ! north
 g_poly%xypoly(1,:,:) = ypolyn(:,:) ! east
 g_poly%lpoly(:)      = lpolyn(:)
 g_poly%nclose        = nclosen            ! new nclose is set 2019.02.25
 if ( nclosen .lt. nclose) then     ! when nclose is reduced 2019.02.25, here is skipped nclosen=nclose=0
  deallocate( g_poly%loc )
  allocate(   g_poly%loc(nclosen,2) )
  g_poly%loc(1:nclosen,1:2) = locn(1:nclosen,1:2)
 end if

 write(*,*) "Total number of polygon is reduced to",lpoly0
 write(*,*) "# of coastal node is reduced to",ncoast
 write(*,*) "# of coastal nodes on the calculation boundaries is",nbound
 write(*,*) "# of boundary polygon (nclose) is",nclosen
 write(*,*) "### landpolygon10 end!###"

return

end subroutine landpolygon10
!############################################ topo2grdxyz
!Copied from extrude.f90 and modified on Sep. 17, 2015
!# modify to interpolate in lat and lon, then convert to UTM coordinate on Dec. 17, 2015
!# modified for multiple topofraphy files on 2017.09.26
!# this program assumes eariler files are narrow compared to later ones 2017.09.26
subroutine topo2grdxyz(g_param,g_grd,g_param_water)
!use mesh_type
use param
use topo_tool   ! 2023.08.31
use water_level ! 2023.08.31
implicit none
type(param_forward),    intent(in)    :: g_param
type(grid_data) ,       intent(out)   :: g_grd ! only for the first topofile (only for the focus area)
type(param_water_level),intent(in)    :: g_param_water ! 2023.08.31
real(8)       :: lon_c, lat_c, x_c, y_c
real(8)       :: rlat,rlon,x_utm,y_utm
real(8),      allocatable,dimension(:,:) :: lon1,lat1, z1
real(8),      allocatable,dimension(:)   :: x,y
character(50),allocatable,dimension(:)   :: files           ! 2017.09.25
integer(4),   allocatable,dimension(:)   :: nt,nsouth,neast ! 2017.09.25
real(8),      allocatable,dimension(:,:) :: lonlatshift     ! 2017.09.26
real(8)       :: cx0_utm,cy0_utm
real(8)       :: lonorigin, latorigin, elev
character(3)  :: UTM
character(2)  :: num ! 20200728
character(120):: line
character(50) :: header2d,header3d,pregeo,topofile
integer(4)    :: ifile, ntmax
logical       :: iflag_topright_land ! 2023.08.31 iflag_topright_land
integer(4)    :: ntopo, i, j, ii, k, neast1,nsouth1,nfile
real(8)       :: zz1,zz2,zz3,zz4             ! 2017.09.26
real(8)       :: lonw,lone,latn,lats         ! 2017.09.26

!#[0]## set
nfile       = g_param%nfile           ! 2017.09.25
allocate(files(nfile))                ! 2017.09.25
allocate(nt(nfile),neast(nfile),nsouth(nfile)) ! 2017.09.25
allocate(lonlatshift(2,nfile) )       ! 2017.09.26
files       = g_param%topofile        ! 2017.09.25
lonlatshift = g_param%lonlatshift     ! 2017.09.26
UTM         = g_param%UTM   ! zone
lonorigin   = g_param%lonorigin
latorigin   = g_param%latorigin

!# [1] ### count lines of gebcofile and allocate
ntmax = 0     ! 2017.09.25
do ifile=1,nfile
  ntopo=0
  write(*,*) "Topo File is ",files(ifile)
  open(1,file=files(ifile))
  do while (ntopo .ge. 0)
    read(1,*,end=99)
    ntopo=ntopo+1 ! ntopo is # of lines in gebcofile
  end do
  99 nt(ifile)=ntopo
  close(1)
  write(*,*) "  nt(ifile)=",nt(ifile)
  ntmax=max(ntmax, nt(ifile)) ! 2017.09.25
end do
allocate( lon1(ntmax,nfile),lat1(ntmax,nfile)) ! 2017.09.25
allocate( z1(ntmax,nfile) ) ! x1(i,j) is x coord of j-th node of i-th file

!# [2] ### read coordinates in gebcofile
do ifile=1,nfile
  open(1,file=files(ifile))
  write(*,'(a,2f15.7,a,i5)') " lonlatshift=",lonlatshift(1:2,ifile)," ifile",ifile ! 2021.09.29
  do i=1,nt(ifile)
    read(1,*) lon1(i,ifile),lat1(i,ifile), z1(i,ifile)   ! 2017.09.25
    lon1(i,ifile) = lon1(i,ifile) + lonlatshift(1,ifile) ! 2017.09.25
    lat1(i,ifile) = lat1(i,ifile) + lonlatshift(2,ifile) ! 2017.09.25
  end do
  z1(:,ifile)=z1(:,ifile)/1000.d0 ! [m] -> [km]
  close(1)
end do

!# [3] ### change order and measure # of nodes in horizontal and vertical directions
!# reorder the values to ((west -> east), north -> south )
do i=1,nfile ! file loop
!  write(*,*) "lon1(1:2,ifile)=",lon1(1:2,i),"ifile",i
  call changeorder(lon1(1:nt(i),i),lat1(1:nt(i),i),z1(1:nt(i),i), nt(i), nsouth(i), neast(i))
  write(num,'(i2.2)') i                                       ! 20200728
!  open(1,file="reordered_topo"//num//".dat")                  ! 20200728
!  write(1,'(3f15.7)') (lon1(j,i),lat1(j,i),z1(j,i),j=1,nt(i)) ! 20200728
!  close(1)                                                    ! 20200728
end do

!# [4] ### calculate depth "z" at all the horizontal nodes
rlon=lonorigin !copy because lonorigin and latorigin are parameters
rlat=latorigin
CALL UTMGMT(rlon,rlat,x_utm,y_utm,UTM,0) ! 0:LONLAT2UTM , 1:UTM2LONLAT
cx0_utm=x_utm
cy0_utm=y_utm
write(*,'(a,f15.7,a,f15.7)') " (cx0_utm,cy0_utm)=",cx0_utm,"  ",cy0_utm ! 2021.09.30

allocate(x(nt(1)),y(nt(1)))
CALL UTMGMT_N(nt(1),lon1(1:nt(1),1),lat1(1:nt(1),1),x,y,UTM,0) ! see m_param.f90
x = x - x_utm
y = y - y_utm

!# [4] ## set g_grd only for the first topo file
allocate(g_grd%xyz(3,nt(1)))
g_grd%xyz(1,:) = x*1.d-3 ! [km]
g_grd%xyz(2,:) = y*1.d-3 ! [km]
g_grd%xyz(3,:) = z1(1:nt(1),1) ! [km]
g_grd%node   = nt(1)
g_grd%neast  = neast(1)
g_grd%nsouth = nsouth(1)
allocate(g_grd%lon(neast(1)),g_grd%lat(nsouth(1))) ! 2023.09.01
g_grd%lon = lon1(1:neast(1),1) ! 2023.09.01
do i=1,nsouth(1) ! 2023.09.01
  g_grd%lat(i) = lat1((i-1)*neast(1) + 1,1) ! 2023.09.01
end do        ! 2023.09.01

!#[5]## iflag_topright_land ### 2023.08.23
elev = g_param_water%water_level_elev(1) 
 iflag_topright_land=.false. ! 2023.08.31 inversed from TMTGEMv1.3
 if ( g_grd%xyz(3,g_grd%neast) - elev .ge. 0.d0 ) iflag_topright_land=.true. ! 2023.08.31 above water level of interest
 g_grd%iflag_topright_land = iflag_topright_land              ! 2023.08.31

write(*,*) "### TOPO2GRDXYZ END!! ###"

 if ( .true. ) then ! 2019.03.08
  open(1,file="topo.xyz")
  write(1,'(3g15.7)') (g_grd%xyz(1,i),g_grd%xyz(2,i),g_grd%xyz(3,i),i=1,g_grd%node)
  close(1)
 end if

return
999 continue
write(*,*) "GEGEGE when node # is k=",k,"x(k),y(k)=",x(k),y(k)
write(*,*) "neast1=",neast1,"nsouth1=",nsouth1
stop
end subroutine topo2grdxyz

!########################################################### makecoast5 copied from TMTGEMv1.0 2023.08.31
!# way for search of coastline points 2019.02.19
!# 1. right  boundary
!# 2. bottom boundary
!# 3. left   boundary
!# 4. top    boundary
!# 5. inside boundaries
subroutine makecoast5(g_grd,g_coast,elev) ! elev is added 2023.08.31
use topo_tool ! 2023.08.31
implicit none     ! make points on coastlines
type(grid_data),       intent(in)     :: g_grd
real(8),               intent(in)     :: elev   ! [km] above sea level
type(coast_data),      intent(inout)  :: g_coast
integer(4)                            :: node,neast,nsouth,ncmax
integer(4)                            :: i,ii,j,jj,i1,i2,i3,nbound,ncoast
real(8),   allocatable,dimension(:)   :: x,y,h ! x (northward), y(eastward), z(upward)
real(8),   allocatable,dimension(:)   :: cx,cy
integer(4),allocatable,dimension(:,:) :: ind0
real(8),               dimension(2)   :: cxy


!#[1]## set input
neast     = g_grd%neast
nsouth    = g_grd%nsouth
node      = g_grd%node
ncmax     = g_coast%ncmax

allocate( x(node),y(node),h(node))
allocate( cx(ncmax),cy(ncmax)    )
allocate( ind0(ncmax,2)         ) ! 2018.08.28
x(1:node) = g_grd%xyz(2,1:node)    ! north    [km]
y(1:node) = g_grd%xyz(1,1:node)    ! east     [km]
h(1:node) = g_grd%xyz(3,1:node)    ! altitude [m]

jj=0
!!## NOTE THAT the top right corner of the calculation space is in the OCEAN ##!!!!!
!# right boundary
j=neast
do i=1,nsouth-1 ! downward
  ii=neast*(i-1)+j
  i1=ii;i3=ii+neast
  if ((h(i1)-elev)*(h(i3)-elev) .lt. 0) then ! when coastline exists between i1 and i3 ! elev is added 2023.08.31
    jj=jj+1
    call interpolate(h(i1),h(i3),x(i1),x(i3),y(i1),y(i3),cxy(1:2),elev)
!    write(2,*) cxy(1), cxy(2),ii," 2"
    cx(jj)=cxy(1)
    cy(jj)=cxy(2)
    ind0(jj,1:2)=(/ii,2/) ! gebco node ii is upside of jjthe coast node 2019.02.19
    end if
end do
!# bottom boundary
i=nsouth ! leftward
do j=neast-1,1,-1
 ii=neast*(i-1)+j
 i1=ii;i2=ii+1
! write(*,*) "i1",i1,"i2",i2,"h(i1)",h(i1),"h(i2)",h(i2),"nsouth",nsouth,"neast",neast
 if ((h(i1)-elev)*(h(i2)-elev) .lt. 0) then ! 2023.08.31
  jj=jj+1
  call interpolate(h(i1),h(i2),x(i1),x(i2),y(i1),y(i2),cxy(1:2),elev)
  !write(2,*) cxy(1), cxy(2),ii," 1"
  cx(jj)=cxy(1)
  cy(jj)=cxy(2)
  ind0(jj,1:2)=(/ii,1/) ! gebco node ii is left of jjthe coast node 2019.02.19
 end if
end do
!# left boundary
j=1
do i=nsouth-1,1,-1 ! upward
 ii=neast*(i-1)+j
 i1=ii;i3=ii+neast
 if ((h(i1)-elev)*(h(i3)-elev) .lt. 0) then ! 2023.08.31
  jj=jj+1
  call interpolate(h(i1),h(i3),x(i1),x(i3),y(i1),y(i3),cxy(1:2),elev)
  !write(2,*) cxy(1), cxy(2),ii," 2"
  cx(jj)=cxy(1)
  cy(jj)=cxy(2)
  ind0(jj,1:2)=(/ii,2/)
 end if
end do
!# top boundary
i=1 ! rightward
do j=1,neast-1
 ii=neast*(i-1)+j
 i1=ii;i2=ii+1
 if ((h(i1)-elev)*(h(i2)-elev) .lt. 0) then ! 2023.08.31
  jj=jj+1
  call interpolate(h(i1),h(i2),x(i1),x(i2),y(i1),y(i2),cxy(1:2),elev)
  !write(2,*) cxy(1), cxy(2),ii," 1"
  cx(jj)=cxy(1)
  cy(jj)=cxy(2)
  ind0(jj,1:2)=(/ii,1/)
 end if
end do
nbound=jj
!# inside the boundaries
!# "1"
do i=2,nsouth-1
 do j=1,neast-1
  ii=neast*(i-1)+j
  i1=ii;i2=ii+1
! 1
  if ((h(i1)-elev)*(h(i2)-elev) .lt. 0) then ! 2023.08.31
   jj=jj+1
   call interpolate(h(i1),h(i2),x(i1),x(i2),y(i1),y(i2),cxy(1:2),elev)
   !write(2,*) cxy(1), cxy(2),ii," 1"
   cx(jj)=cxy(1)
   cy(jj)=cxy(2)
   ind0(jj,1:2)=(/ii,1/)
  end if
 end do
end do
!# "1"
do i=1,nsouth-1
 do j=2,neast-1
  ii=neast*(i-1)+j
  i1=ii;i3=ii+neast
  ! 2
  if ((h(i1)-elev)*(h(i3)-elev) .lt. 0) then ! 2023.08.31
   jj=jj+1
   call interpolate(h(i1),h(i3),x(i1),x(i3),y(i1),y(i3),cxy(1:2),elev)
   !write(2,*) cxy(1), cxy(2),ii," 2"
   cx(jj)=cxy(1)
   cy(jj)=cxy(2)
   ind0(jj,1:2)=(/ii,2/)
  end if
 end do
end do
ncoast=jj
write(*,*) "ncoast (# of nodes on coastline is)",ncoast
if ( g_coast%ncmax .lt. ncoast +4   ) then ! +4 is required for addcorner.f90
 write(*,*) "GEGEGE! g_coast%ncmax=",g_coast%ncmax,"is smaller than ncoast+4",ncoast+4
 stop
end if
write(*,*) "### makecoast5 end! ###"
!close(2)

!#[3]## set output
 g_coast%ind=ind0
 g_coast%cxy(1,1:ncmax) = cy(1:ncmax) ! east
 g_coast%cxy(2,1:ncmax) = cx(1:ncmax) ! north
 g_coast%ncoast = ncoast
 g_coast%nbound = nbound ! # of coastline nodes at calculation boundaries

!#[4]## output coast node
if ( .true. ) then
 open(2,file="coast.dat",status='replace')
  do i=1,g_coast%ncoast
!  write(*,*) "i",i,g_coast%cxy(1:2,i),g_coast%ind(i,1:2)
  write(2,*) g_coast%cxy(1:2,i),g_coast%ind(i,1:2)
  end do
  close(2)
 end if

return
end subroutine makecoast5

!########################################################### makepolygon7 copied from coastline.f90 in TMTGEMv1.3
subroutine makepolygon7(g_coast,g_grd,g_poly,elev) ! lpoly0= # of polygons
!###            make polygons and output "polygon.dat"           ###
!label(num_coastnode,1) : polygon number, default = 0, meaning not belonging to any polygon
!label(num_coastnode,2) : node number in n-th polygon nodes, which begins with 1 in each polygon
!nuwd : up wind direction, 1 -> from top or left, 2 -> from bottom or right
use topo_tool ! 2023.08.31
implicit none
type(coast_data),      intent(inout)   :: g_coast ! label is calculated here
type(grid_data),       intent(in)      :: g_grd
type(poly_data),       intent(inout)   :: g_poly
real(8),               intent(in)      :: elev ! 2023.08.31
integer(4)                             :: lpoly0
integer(4),allocatable,dimension(:)    :: lpoly
real(8),   allocatable,dimension(:,:)  :: xpoly,ypoly
integer(4),allocatable,dimension(:,:)  :: ind0, label, ind
real(8),   allocatable,dimension(:)    :: cx,cy
integer(4)                             :: nbound,node,neast,i,l,lpmax, ncmax, ncoast
real(8),   allocatable,dimension(:)    :: h
logical                                :: iflag_topright_land ! 2019.02.25

!#[0]## set input
 ! g_coast
 ncoast = g_coast%ncoast
 ncmax  = g_coast%ncmax
 node   = g_grd%node
 lpmax  = g_poly%lpmax
 write(*,*) "ncoast",ncoast
 write(*,*) "ncmax",ncmax
 write(*,*) "node",node
 write(*,*) "lpmax",lpmax
 allocate( cx(ncmax),cy(ncmax))
 allocate( ind0(ncmax,2), label(ncmax,2), ind(ncmax,2))
 allocate( h(node) )
 allocate( xpoly(lpmax,ncmax) )
 allocate( ypoly(lpmax,ncmax) )
 allocate( lpoly(lpmax)       )
 nbound = g_coast%nbound
 write(*,*) "nbound",nbound
 cx(:)  = g_coast%cxy(2,:) ! north
 cy(:)  = g_coast%cxy(1,:) ! east
 ind0   = g_coast%ind
 neast  = g_grd%neast
 h(:)   = g_grd%xyz(3,:) ! altitude [m]
 label  = 0
 lpoly0 = 0
 iflag_topright_land = g_grd%iflag_topright_land ! 2019.02.25

!#[1]## unclosed polygon, which touches boundaries
 do i=1,nbound
  if (label(i,1) .eq. 0) then ! when i-th coastline node doesn't belong to any polygons 2019.02.19
   call loopelement(i,lpoly0,node,neast,ncmax,ncoast,h,label,ind0,cx,cy,lpoly,xpoly,ypoly,ind,lpmax,elev)
  end if
 end do
 write(*,*) "nbound loop end"

!#[2]## including inside the boundaries
 do i=nbound+1,ncoast
!  write(*,*) "label(i,1)",label(i,1),"i/ncoast",i,"/",ncoast
  if (label(i,1) .eq. 0) then
   call loopelement(i,lpoly0,node,neast,ncmax,ncoast,h,label,ind0,cx,cy,lpoly,xpoly,ypoly,ind,lpmax,elev)
  end if
 end do

!#[4]## set label to g_coast
!  g_coast%label = label ! commented out 2019.02.25
  g_coast%ind_r = ind
  g_poly%lpoly0 = lpoly0
  g_poly%lpoly  = lpoly
  g_poly%xypoly(1,:,:) = ypoly ! east
  g_poly%xypoly(2,:,:) = xpoly ! north
  g_poly%nclose = nbound/2      ! 2019.02.25
  g_poly%iflag_topright_land = iflag_topright_land ! 2019.02.28

 write(*,*) "Total number of polygon is",lpoly0
 write(*,*) "### makepolygon7 end!###"

 return
end subroutine makepolygon7

!############################################### outpolygon8
 subroutine outpolygon8(g_poly,g_coast,head,inum)
 use topo_tool ! 2023.08.31
 implicit none
 type(poly_data),       intent(in)     :: g_poly
 type(coast_data),      intent(in)     :: g_coast
 character(50),         intent(in)     :: head
 integer(4),            intent(in)     :: inum
 integer(4)                            :: ncoast, nbound, ncmax
 integer(4)                            :: nclose,iclose,ii,jj,ij,i,j
 integer(4),allocatable,dimension(:,:) :: ind
 character(1)                          :: num
 character(70)                         :: polygonfile
 character(70)                         :: polygon_numfile ! 2019.02.27
 real(8)                               :: x_sum,y_sum     ! 2019.02.27

 write(num,'(i1.1)') inum
 polygonfile=head(1:len_trim(head))//num//".dat"
 polygon_numfile=head(1:len_trim(head))//num//"_num.dat" ! 2019.02.27
 open(1,file=polygonfile)
 open(2,file=polygon_numfile) ! 2019.02.27

!#[1]## set input
  ncmax  = g_coast%ncmax
  ncoast = g_coast%ncoast
  nbound = g_coast%nbound
  allocate( ind(ncmax,2)) ! 2018.08.28
  ind    = g_coast%ind_r
  nclose = g_poly%nclose  ! 2019.02.25

!# polygon.dat read
 write(*,*) "nclose=",nclose
 ii=0
 do i=1,g_poly%lpoly0
   iclose=1
   if ( i .le. nclose ) iclose=0
   !# iclose=0 -> polygon is not closed because reaching calculation boundaries
   write(1,*) i,g_poly%lpoly(i),iclose
!  if ( iclose .eq. 0 ) nclose=nclose+1
   x_sum = 0.d0 ! 2019.02.27
   y_sum = 0.d0 ! 2019.02.27
   do j=1,g_poly%lpoly(i)
     ii=ii+1      !# ii : new ncoast node number
!      write(1,'(i8,2g15.7,2i10)') j, g_poly%xypoly(2,i,j),g_poly%xypoly(1,i,j),ind(ii,1),ind(ii,2)!2019.02.19
      write(1,'(i8,2g15.7)') j, g_poly%xypoly(2,i,j),g_poly%xypoly(1,i,j) !2019.02.19      x_sum = x_sum + g_poly%xypoly(1,i,j)/g_poly%lpoly(i) ! 2019.02.27
      y_sum = y_sum + g_poly%xypoly(2,i,j)/g_poly%lpoly(i) ! 2019.02.27
   end do
   write(2,*) x_sum,y_sum, i ! 2019.02.27
 end do
 close(1)
 close(2) ! 2019.02.27
 !# Now the ncoast coastline nodes are ordered
 !    to be {1 2 3},{4,5,6,7,8,9}, where {} means one polygon
 if (ncoast .ne. ii) then
  write(*,*) "GEGEGE!!! ncoast is not equal to ii!!ncoast=",ncoast,"ii=",ii
  stop
 end if

!# polygon.dat read end!
 if ( nbound/2 .ne. nclose) then
  write(*,*) "GEGEGE nbound/2 is not equal to nclose!!! nclose=",nclose,"nbound=",nbound
  stop
 end if
 write(*,*) "# of coastal nodes on the calculation boundaries is",nbound
 write(*,*) "# of boundary polygon (nclose) is",nclose
 write(*,'(a,a)') " # out put ",trim(polygonfile) ! 2019.03.05
 write(*,*) "### outpolygon8 end!###"

return

end subroutine outpolygon8
!######################################################### interpolate
subroutine interpolate(h1,h2,x1,x2,y1,y2,cxy,elev) ! 2023.08.31
implicit none
real(8)              :: h1,h2,x1,x2,y1,y2,a,elev ! 2023.08.31 elev is added
real(8),dimension(2) :: cxy

a=dabs(h1-elev)/(dabs(h1-elev)+dabs(h2-elev)) ! 2023.08.31
cxy(1)=x1*(1.d0-a)+x2*a
cxy(2)=y1*(1.d0-a)+y2*a
if (cxy(2) > 1.d+4) then
write(*,*) "a,1-a",a,1.-a
write(*,*) "h1,h2=",h1,h2
write(*,*) "x1,x2=",x1,x2
write(*,*) "y1,y2=",y1,y2
write(*,*) "cxy(1,2)=", cxy(1),cxy(2)
stop
end if

return
end

!######################################################### loopelement
! iii : start node # of polygon
! label(num_coastnode,1) : polygon number, default = 0, meaning not belonging to any polygon
! label(num_coastnode,2) : node number in n-th polygon nodes, which begins with 1
! iclose : 0 -> polygon is not closed, i.e. touching boundaries, 1 -> polygon is closed
subroutine loopelement(iii,l,node,neast,ncmax,ncoast,h,label,ind0,&
&                      cx,cy,lpoly,xpoly,ypoly,ind,lpmax,elev)
implicit none
integer(4),                       intent(in)    :: ncmax,node,ncoast,neast      !      2019.02.19
integer(4)                   ,    intent(in)    :: iii ! starting node id in topo grid 2019.02.19
integer(4)                   ,    intent(inout) :: l   ! current polygon #             2019.02.19
integer(4)                                      :: nuwd,nuwd2,npoly,iclose,lpmax
integer(4)                                      :: i,ii,j,k
real(8),                          intent(in)    :: elev ! 2023.08.31
integer(4),dimension(ncmax,2),    intent(in)    :: ind0
integer(4),dimension(ncmax,2),    intent(inout) :: label,ind
integer(4),dimension(lpmax),      intent(inout) :: lpoly
real(8),   dimension(lpmax,ncmax),intent(inout) :: xpoly,ypoly
integer(4),dimension(ncoast)                    :: polygon
real(8),   dimension(node)                      :: h
real(8),   dimension(ncmax)                     :: cx,cy

 i=iii  ! i is the starting node
 l=l+1  ! increase polygon number
 j=1
! i-th node on coast nods
 label(i,1)=l ! set polygon number
 label(i,2)=j ! set node number in l-th polygon
 polygon(1)=i
!write(*,*) "polygon #=",l,"start!!"
! set first nuwd
 nuwd=0
 if (ind0(i,2) .eq. 1) then ! gebco node is left side of i-th coastline node
    nuwd=1 ! from top
    if ( ind0(i,1) .gt. node - neast ) nuwd = 3 ! from bottom
 end if
 if (ind0(i,2) .eq. 2) then
    nuwd=2 ! from left
    if ( mod(ind0(i,1),neast) .eq. 0 ) nuwd = 4 ! from right
 end if
!   write(*,*) "ind(i,1)=",ind(i,1),"ind(i,2)=",ind(i,2),"nuwd=",nuwd
!##  polygon loop start
100 continue
 call findnextnode(i,nuwd,ncoast,ncmax,node,neast,h,ind0,ii,nuwd2,elev) ! elev is added 2023.08.31
!write(*,*) "i,ind(i,1),ind(i,2),nuwd=",i,ind(i,1),ind(i,2)
!write(*,*) "ii,ind(ii,1),ind(ii,2),nuwd2=",ii,ind(ii,1),ind(ii,2),nuwd2
!### if the circle is closed ###
if ( label(ii,1) .eq. l ) then
    !   write(*,*) "polygon is closed"
    iclose=1
    goto 200
end if
 label(ii,1)=l
 j=j+1
 label(ii,2)=j
 polygon(j)=ii ! ii is the next node of i in j-th polygon
!### if reach calculation boundaries ###
 if ( (ind0(ii,1) .le. neast .and. ind0(ii,2) .eq. 1 ) .or. &  !top boundary
    &     (ind0(ii,1) .gt. node-neast .and. ind0(ii,2) .eq. 1 ) .or. & !bottom boundary
    &     (mod(ind0(ii,1),neast) .eq. 0 .and. ind0(ii,2) .eq. 2 ) .or. & !right boundary
    &     (mod(ind0(ii,1),neast) .eq. 1 .and. ind0(ii,2) .eq. 2 )) then ! left boundary
!	write(*,*) "reach calculation boundary"
    iclose=0
   goto 200
 end if
! setting for next step
 i=ii
 nuwd=nuwd2
 goto 100
!## polygon loop end
200 continue
npoly=j ! npoly is the # of nodes which constitute l-th polygon
!write(ifile,'(3i7)') l,npoly,iclose
!write(ifile,'(i7,2g15.7,2i7)') (j,cx(polygon(j)),cy(polygon(j)),&
!&                       ind0(polygon(j),1),   ind0(polygon(j),2),j=1,npoly)
!write(*,*) "polygon #",l,"end! # of corner is",npoly
k=0
do j=1,l-1 ! l is the polygon #
   k=k+lpoly(j) ! k is summation of nodes of polygons from 1-st to l-th polygons
end do
lpoly(l)=npoly ! # of nodes in l-th polygon is npoly
do j=1,npoly ! npoly is the # of nodes which constitue l-th polygon
   xpoly(l,j)=cx(polygon(j)) ! polygon(j) is the node # of j-th node of l-th polygon
   ypoly(l,j)=cy(polygon(j))
   k=k+1
   ind(k,1:2)=ind0(polygon(j),1:2)
   ! convert ind0 to ind, which is ordered in polygon node order
end do
return
end
!
!############################################## findnextnode
! findnextnode return the next node on coastline
!# nuwd stands for up wind direction in looking for on-coastline nodes
! nuwd = 1 : from the top or left side
! nuwd = 2 : from the bottom or right side
! then we can make tabele as follows
!  nuwd    1         2              3             4
!         from top  from left      from bottom   from right
!# i is the i-th node in ncoast nodes on coastlines
!# h are altitude data from gebco original file [m]
!# ind(ncoast,2) points to left/top original gebco grid nodes
! ind(i,1) : grid node # of original gebco file on left or up side
! ind(i,2) : 1 -> ind(1) indicates left node, 2 -> ind(2) is top side node
subroutine findnextnode(i,nuwd,ncoast,ncmax,node,neast,h,ind0,ii,nuwd2,elev) ! elev added on 2023.08.31
implicit none
integer(4)                    :: i,nuwd,neast,ncoast,node,ii,j,ij,jj,nuwd2,ncmax
real(8),   intent(in)         :: elev ! 2023.08.31
real(8),   dimension(node)    :: h
integer(4),dimension(ncmax,2) :: ind0
integer(4),dimension(3,2)     :: ind1
integer(4),dimension(3)       :: id3,iflag
integer(4),dimension(3,4)     :: nuwd34
data nuwd34(1,1:4) /4,3,4,3/
data nuwd34(2,1:4) /2,1,2,1/
data nuwd34(3,1:4) /1,2,3,4/
! look for on-coastline nodes on the other three edges
! determine three (ind(1), ind(2)) pairs
! set ind1
ind1(:,:)=0
if ( nuwd .eq. 1) then
ind1(1,1)=ind0(i,1)       ; ind1(1,2)=2
ind1(2,1)=ind0(i,1)+1     ; ind1(2,2)=2
ind1(3,1)=ind0(i,1)+neast ; ind1(3,2)=1
end if
if ( nuwd .eq. 2) then
ind1(1,1)=ind0(i,1)       ; ind1(1,2)=1
ind1(2,1)=ind0(i,1)+neast ; ind1(2,2)=1
ind1(3,1)=ind0(i,1)+1     ; ind1(3,2)=2
end if
if ( nuwd .eq. 3) then
ind1(1,1)=ind0(i,1)-neast   ; ind1(1,2)=2
ind1(2,1)=ind0(i,1)-neast+1 ; ind1(2,2)=2
ind1(3,1)=ind0(i,1)-neast   ; ind1(3,2)=1
end if
if ( nuwd .eq. 4) then
ind1(1,1)=ind0(i,1)-1       ; ind1(1,2)=1
ind1(2,1)=ind0(i,1)-1+neast ; ind1(2,2)=1
ind1(3,1)=ind0(i,1)-1       ; ind1(3,2)=2
end if
if ( ind1(1,1) .eq. 0 ) then
write(*,*) "GEGEGE nuwd was not set from 1 to 4!!"
stop
end if
!##       ind1 set end      ####
!
!## look for next node on coastline ##
ij=0
id3(1:3)=0
iflag(1:3)=0
do j=1,ncoast
! write(*,*) "j,ind(j,1),ind(j,2)",j,ind(j,1),ind(j,2)
 do jj=1,3
  if (ind0(j,1) .eq. ind1(jj,1) .and. ind0(j,2) .eq. ind1(jj,2)) then
  ij=ij+1
  id3(jj)=j
  iflag(jj)=1
  end if
 end do
end do
if (ij .eq. 0 ) then
write(*,*) "GEGEGE there are no next coastline nodes!! ij=",ij
write(*,*) h(ind0(i,1)-1),h(ind0(i,1)),h(ind0(i,1)+1)
write(*,*) h(ind0(i,1)-1+neast),h(ind0(i,1)+neast),h(ind0(i,1)+1+neast)
stop
end if
!## when ij = 1, 2, and 3 ##
if ( ij .eq. 1 ) then
 do j=1,3
  if ( iflag(j) .eq. 1) then
  ii=id3(j)
  nuwd2=nuwd34(j,nuwd)
  end if
 end do
end if
!## ij = 2
if ( ij .eq. 2 ) then
 write(*,*) "GEGEGE there are two next coast line nodes!! node # on coasatline is",i
 write(*,*) "nuwd=",nuwd
 if (ind0(i,1)-1-neast .ge.1) then
  write(*,*) h(ind0(i,1)-1-neast),h(ind0(i,1)-neast),h(ind0(i,1)+1-neast)
 end if
 write(*,*) h(ind0(i,1)-1),h(ind0(i,1)),h(ind0(i,1)+1)
 if (ind0(i,1)+1+neast .le. node) then
  write(*,*) h(ind0(i,1)-1+neast),h(ind0(i,1)+neast),h(ind0(i,1)+1+neast)
 end if
stop
end if
if ( ij .eq. 3 ) then
! if ( h(ind0(i,1)) > 0 ) then ! commented out 2019.03.31 "keep nallow water roads"
 if ( h(ind0(i,1)) -elev < 0 ) then ! 2019.03.31 avoid nallow water roads
 ii=id3(1)
 nuwd2=nuwd34(1,nuwd)
 else
 ii=id3(2)
 nuwd2=nuwd34(2,nuwd)
 end if
end if
return
end
!############################################################## SPLINT
SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
IMPLICIT REAL*8 (A-H,O-Z)
DIMENSION XA(N), YA(N), Y2A(N)
KLO = 1
KHI = N
1 IF ( KHI - KLO .GT. 1 ) THEN
K = ( KHI + KLO )/2
IF ( XA(K) .GT. X ) THEN
KHI = K
ELSE
KLO = K
END IF
GO TO 1
END IF
H = XA(KHI) - XA(KLO)
IF ( H .EQ. 0.D0 ) PRINT *,'BAD XA INPUT'
A = ( XA(KHI) - X )/H
B = ( X - XA(KLO) )/H
Y = A*YA(KLO) + B*YA(KHI) + ( ( A**3 - A )*Y2A(KLO) + ( B**3 - B ) &
&    *Y2A(KHI) )*(H**2)/6.D0
RETURN
END
!############################################################## SPLINE
SUBROUTINE SPLINE(X,Y,N,Y2)
IMPLICIT REAL*8 (A-H,O-Z)
!DIMENSION X(N), Y(N), Y2(N), U(2048)
DIMENSION X(N), Y(N), Y2(N), U(N)
! Interpolation by Natural Spline: Y2(1)=U(1)=QN=UN=0.d0
Y2(1) = 0.D0
U(1) = 0.D0
DO 10 I = 2, N - 1
SIG = ( X(I) - X(I-1) )/( X(I+1) - X(I-1) )
P = SIG*Y2(I-1) + 2.D0
Y2(I) = ( SIG - 1.D0 )/P
U(I) = ( 6.D0*( ( Y(I+1) - Y(I) )/( X(I+1) - X(I) ) - ( Y(I) - &
& Y(I-1) )/( X(I) - X(I-1) ) )/( X(I+1) - X(I-1) ) - SIG*U(I-1) )/P
10 CONTINUE
!
QN = 0.D0
UN = 0.D0
!
Y2(N) = ( UN - QN*U(N-1) )/( QN*Y2(N-1) + 1.D0 )
!
DO 20 K = N - 1, 1, -1
Y2(K) = Y2(K)*Y2(K+1) + U(K)
20 CONTINUE
!
RETURN
END




!########################################################### outpolygon12
 subroutine outpolygon12(g_poly,head,inum)
 use topo_tool
 implicit none
 type(poly_data),intent(in)  :: g_poly
 character(50),  intent(in)  :: head
 integer(4),     intent(in)  :: inum
 integer(4)                  :: i,k,iclose
 character(1)                :: num
 character(70)               :: polygonfile
!# iclose : 0 -> not closed polygon, 1 -> closed polygon

!#[1]## setting
 write(num,'(i1)') inum
 polygonfile=head(1:len_trim(head))//num//".dat"

!#[2]## output
 open(1,file=polygonfile)
 do i=1,g_poly%lpoly0
  iclose=1 ! closed
  if ( i .le. g_poly%nclose) iclose=0 ! unclosed
  write(1,'(3i7)') i,g_poly%lpoly(i),iclose
  write(1,'(i7,2g15.7)')(k,g_poly%xypoly(2,i,k),g_poly%xypoly(1,i,k),k=1,g_poly%lpoly(i))
 end do
 close(1)

 write(*,*) "### outpolygon12 end! ###"
 return

 end subroutine outpolygon12
