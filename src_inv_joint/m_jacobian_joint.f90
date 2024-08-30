!# Modified on 2017.09.03 for multisource inversion
!# Jacobian matrix for ACTIVE responses
!# coded on May 13, 2017
!#     | J2 |
!# J = | J3 | [Ndat,Nmodel]
!#     | vv |
!#     | JN | Ji is for i-th frequency
!#
!# Ji [Nobs,Nmodel]; j = obs index, k = model index
!# Ji[j,k]  = {d/dmk(|bz_i,j|)|bz_1,j| - |bz_i,j|d/dmk(|bz_1,j|)}/|bz_1,j|^2
!#
!# d/dmk(|bz_i,j|) = |bz_i,j| Re{(1/PT*Al)*PT*d/dmk(Al) }
!#
!# Au = P
!# P [nedge, nobs]
!# u [nedge, nobs]
module jacobian_joint
use matrix
implicit none

!# 2017.06.07
 type amp_phase_dm
  integer(4) :: nobs          ! for each freqency
  integer(4) :: nmodelactive  ! 2018.06.25
  integer(4) :: nsr_inv  ! # of sources,     2017.08.31
!  real(8),  allocatable,dimension(:) :: ampbz
  type(real_crs_matrix),allocatable,dimension(:,:) :: dampdm ! damp/dm (5,nsr) 2018.10.05
  type(real_crs_matrix),allocatable,dimension(:,:) :: dphadm ! dpha/dm (5,nsr) 2018.10.05
 end type

 !# 2022.01.05
 type mt_dm
  integer(4) :: nobs_mt          ! for each freqency
  integer(4) :: nmodelactive     ! 2022.01.05
  type(complex_crs_matrix) :: dzxxdm ! [nobs,nmodelactive]  
  type(complex_crs_matrix) :: dzxydm ! 
  type(complex_crs_matrix) :: dzyxdm ! 
  type(complex_crs_matrix) :: dzyydm ! 
 end type

 !# 2023.12.22 Jacoain for tipper
 type tip_dm
  integer(4) :: nobs_mt             ! for each freqency
  integer(4) :: nmodelactive        ! 2023.12.22
  type(complex_crs_matrix) :: dtxdm ! [nobs,nmodelactive]  
  type(complex_crs_matrix) :: dtydm ! 
 end type

contains

!###########################
!# nmode -> nmodelactive on 2018.06.25
!# generated on 2017.08.31 for multiple sources
subroutine allocateapdm(nobs,nfreq,nmodelactive,nsr_inv,g_apdm) ! 2018.06.25
implicit none
integer(4),        intent(in)    :: nobs,nfreq
integer(4),        intent(in)    :: nmodelactive    ! 2017.09.04
integer(4),        intent(in)    :: nsr_inv         ! 2017.07.14
type(amp_phase_dm),intent(inout) :: g_apdm(nfreq)
integer(4) :: i

do i=1,nfreq
 g_apdm(i)%nobs         = nobs
 g_apdm(i)%nsr_inv      = nsr_inv        ! 2017.07.14
 g_apdm(i)%nmodelactive = nmodelactive   ! 2018.06.25
 allocate( g_apdm(i)%dampdm(5,nsr_inv) ) ! 2018.10.05
 allocate( g_apdm(i)%dphadm(5,nsr_inv) ) ! 2018.10.05
end do

return
end
!###########################
!# nmode -> nmodelactive on 2018.06.25
!# generated on 2017.08.31 for multiple sources
subroutine allocatemtdm(nobs_mt,nfreq_mt,nmodelactive,g_mtdm) ! 2018.06.25
  implicit none
  integer(4),        intent(in)    :: nobs_mt,nfreq_mt
  integer(4),        intent(in)    :: nmodelactive    ! 2017.09.04
  type(mt_dm),    intent(inout)    :: g_mtdm(nfreq_mt)
  integer(4) :: i
  
  do i=1,nfreq_mt
   g_mtdm(i)%nobs_mt      = nobs_mt
   g_mtdm(i)%nmodelactive = nmodelactive   ! 2018.06.25
  end do
  
  return
  end
!########################### Transxtomodel
!# coded on 2018.01.22
!# iflag = 1 : model to X
!# iflag = 2 : X to model
subroutine TRANSMODELX(h_model,g_param_joint,iflag)
use param_jointinv       ! see m_param_jointinv.f90
use caltime
implicit none
integer(4)             , intent(in)     :: iflag
type(param_joint),       intent(in)     :: g_param_joint! see m_param_jointinv.f90
type(model),             intent(inout)  :: h_model
integer(4) :: nmodel,i
real(8),allocatable,dimension(:) :: logrho_model
real(8) :: p_param
real(8) :: logrho_upper,logrho_lower
type(watch) :: t_watch

!#[1]## set
 call watchstart(t_watch)
 nmodel = h_model%nmodel
 allocate(logrho_model(nmodel))
 logrho_model = h_model%logrho_model
 p_param      = g_param_joint%p_param
 logrho_upper = g_param_joint%logrho_upper
 logrho_lower = g_param_joint%logrho_lower

!#[2]## transform
 if ( iflag .eq. 1 ) then ! model -> X
  do i=1,nmodel
  logrho_model(i) = 1./p_param*&
  & log( (logrho_model(i)-logrho_lower)/(logrho_upper-logrho_model(i)))
  end do
 end if

 if ( iflag .eq. 2 ) then ! X to model
  do i=1,nmodel
  logrho_model(i) = (logrho_lower + logrho_upper*exp(p_param*logrho_model(i)))&
                                      & /(1.d0 + exp(p_param*logrho_model(i)))
  end do
 end if

!#[3]## set output
 h_model%logrho_model = logrho_model

 call watchstop(t_watch) ! 2018.01.23
 write(*,'(a,g15.7,a)') "### TRANSMODELX       END!! Time=",t_watch%time,"[min]"
return
end

!########################### Transferjacob 2018.01.22
!# 2018.01.22
subroutine transjacob(g_param_joint,JJ,h_model)
use caltime           ! see m_caltime.f90 2017.09.06
use param_jointinv       ! see m_param_jointinv.f90
use modelpart            !
implicit none
type(real_crs_matrix),   intent(inout)  :: JJ          ! 2017.06.07
type(param_joint), intent(in)     :: g_param_joint! see m_param_jointinv.f90
type(model),             intent(in)     :: h_model
real(8),    allocatable,dimension(:)    :: logrho_model,dmdx
integer(4)                              :: i,ii,ntot,nrow,ncolm,nmodel ! 2018.06.25
integer(4)                              :: nmodelactive ! 2018.06.25
integer(4), allocatable,dimension(:)    :: iactive   ! 2018.06.25
integer(4), allocatable,dimension(:)    :: ptrnmodel ! 2018.06.25
real(8)     :: p_param, mi
real(8)     :: logrho_upper
real(8)     :: logrho_lower
type(watch) :: t_watch

!#[0]## set
 call watchstart(t_watch) ! 2018.01.23
 nmodel       = h_model%nmodel
 nmodelactive = h_model%nmodelactive ! 2018.06.25
 allocate( iactive(nmodel))          ! 2018.06.25
 iactive      = h_model%iactive      ! 2018.06.25
 allocate( ptrnmodel(nmodelactive))  ! 2018.06.25
 allocate(logrho_model(nmodel),dmdx(nmodel))
 logrho_model = h_model%logrho_model
 p_param      = g_param_joint%p_param
 logrho_upper = g_param_joint%logrho_upper
 logrho_lower = g_param_joint%logrho_lower
 ntot         = JJ%ntot
 nrow         = JJ%nrow
 ncolm        = JJ%ncolm
 if ( ncolm .ne. nmodelactive ) then ! 2018.06.25
 write(*,*) "GEGEGE! nrow",nrow,"nmodelactive",nmodelactive,"nmodel",nmodel ! 2018.06.26
 stop
 end if

!#[1]## set ptrnmodel 2018.06.25
 ii = 0
 do i=1,nmodel
  if ( iactive(i) .eq. 1 ) then
   ii=ii+1
   ptrnmodel(ii) = i
  end if
 end do

!#[2]## gen dm/dx
 do i=1,nmodel
  mi   = logrho_model(i)
  dmdx(i) = p_param*(logrho_upper - mi)*(mi-logrho_lower)/(logrho_upper - logrho_lower)
 end do

!#[3]## J-> J'
 do i=1,ntot
  JJ%val(i)=JJ%val(i)*dmdx(ptrnmodel(JJ%item(i))) ! 2018.06.25
 end do

 call watchstop(t_watch) ! 2018.01.23
 write(*,'(a,g15.7,a)') "### TRANSJACOB END!! ### Time=",t_watch%time,"[min]"

return
end

!########################### gen jacobian
!# nmodelactive is introduced on 2018.06.22
!# Modified on 2017.09.03
!# Coded on 2017.05.16
subroutine genjacobian1(nobs,nline,nsr_inv,ut,bs,PT,g_model,h_mesh,l_line,omega,&
                      & g_apdm,g_param_joint,ip,np) !2020.09.18
use matrix            ! see m_matrix.f90
use modelpart         ! see m_modelpart.f90
use line_type         ! see m_line_type.f90
use constants         ! see m_constants.f90
use outerinnerproduct ! see m_outerinnerproduct.f90
use fem_util          ! for volume, intv, (see m_fem_utiil.f90 )
use fem_edge_util     ! see fem_edge_util.f90
use mesh_type         ! see m_mesh_type.f90
use caltime           ! see m_caltime.f90 2017.09.06
use param_jointinv       ! see m_param_jointinv.f90  2018.10.05
implicit none
integer(4),              intent(in)     :: ip,np            ! 2020.09.18
type(param_joint),       intent(in)     :: g_param_joint    ! 2018.10.05
integer(4),              intent(in)     :: nobs,nline
integer(4),              intent(in)     :: nsr_inv          ! # of sources 2017.09.03
type(complex_crs_matrix),intent(in)     :: ut(5)            ! [nobs,nline]  2018.10.05
type(complex_crs_matrix)                :: utt              ! 2018.10.05
complex(8),              intent(in)     :: bs(nline,nsr_inv)!        2017.09.03
type(real_crs_matrix),   intent(in)     :: PT(5)            ! [nobs,nline] 2018.10.05
type(model),             intent(in)     :: g_model
type(line_info),         intent(in)     :: l_line
type(mesh),              intent(in)     :: h_mesh
real(8),                 intent(in)     :: omega
type(amp_phase_dm),      intent(inout)  :: g_apdm           ! 2017.06.07
type(real_crs_matrix)                   :: dampdm(nsr_inv)  ! 2017.09.03 (->g_apdm)
type(real_crs_matrix)                   :: dphadm(nsr_inv)  ! 2017.09.03 (->g_apdm)
complex(8),     dimension(nobs,nsr_inv) :: be               ! 2018.10.05
integer(4)                              :: nmodel,nphys2
integer(4)                              :: imodel,i,j,k,l,m,n,ii,jj,kk!2018.06.21
integer(4)                              :: iele,idirection(6),n6line(6)
integer(4),allocatable,dimension(:)     :: stack,item
real(8),   allocatable,dimension(:)     :: logrho_model
complex(8),allocatable,dimension(:,:)   :: AAL,dbeobs       ! 2017.09.03
real(8),   allocatable,dimension(:,:,:) :: dampdmfull       ! 2017.09.03
real(8),   allocatable,dimension(:,:,:) :: dphadmfull       ! 2017.09.03
real(8)                                 :: threshold =1.d-10
integer(4)                              :: nmodelactive     ! 2018.06.21
integer(4),allocatable,dimension(:)     :: iactive          ! 2018.06.21
integer(4),            dimension(5)     :: iflag_comp       ! 2018.10.05
integer(4)                              :: icomp            ! 2018.10.05
!# for dA/dm
complex(8)                              :: S1(6,6)
real(8)                                 :: elm_xyz(3,4),yy,gn(3,4),v
real(8)                                 :: RM,sigma         ! model = log10(rho)
complex(8)                              :: dBBdm, iunit=(0.d0,1.d0)
real(8),   parameter                    :: L0=1.d+3         ! [m]  scale length
complex(8),allocatable,dimension(:,:)   :: utfull           ! 2017.09.03
complex(8)                              :: z                ! 2017.09.03
type(watch)                             :: t_watch          ! see m_caltime.f90 2017.09.06

 call watchstart(t_watch) ! see m_caltime.f90 2017.12.22

!#[0]## set
 nmodel       = g_model%nmodel
 nmodelactive = g_model%nmodelactive ! 2018.06.21
 allocate( iactive(nmodel) )         ! 2018.06.21
 iactive      = g_model%iactive      ! 2018.06.21
 nphys2       = g_model%nphys2
 allocate(stack(0:nmodel),item(nphys2))
 allocate(logrho_model(nmodel))
 stack        = g_model%model2ele%stack
 item         = g_model%model2ele%item  ! element id for whole element space
 logrho_model = g_model%logrho_model
 iflag_comp   = g_param_joint%iflag_comp ! 2018.10.05
 allocate(utfull(nline,nobs))        ! 2018.10.05
 allocate( AAL(6,nsr_inv) )
 allocate( dbeobs(nobs,nsr_inv)            )       ! 2018.10.05
 allocate( dampdmfull(nobs,nmodelactive,nsr_inv) ) ! 2018.06.22
 allocate( dphadmfull(nobs,nmodelactive,nsr_inv) ) ! 2018.06.22

!#[1]## component loop start          2018.10.05
 do icomp = 1,5                     ! 2018.10.05
  if (iflag_comp(icomp) .eq. 0)cycle! 2018.10.05

!#[1]## cal bx,by,bz,ex,ey dependent on iflag_comp
  utt    = ut(icomp)                 ! 2018.10.05
  utfull = 0.d0                      ! 2017.09.03
  do   i = 1,utt%nrow                ! 2017.09.03
   do j=utt%stack(i-1)+1,utt%stack(i)! 2017.09.03
    utfull(utt%item(j),i) = utt%val(j) ! 2017.09.03
   end do                            ! 2017.09.03
  end do                             ! 2017.09.03
  if ( icomp .ge. 4 ) utfull = -iunit*omega*utfull ! only for Ex,Ey, 2018.10.05

!#[1]## cal either of bx,by,bz,ex,ey dependent on icomp
  do i=1,nsr_inv                     ! 2017.09.03
   call mul_matcrs_cv(PT(icomp),bs(:,i),nline,be(:,i)) ! 2018.10.05
  end do
  if ( icomp .ge. 4 ) be     = -iunit*omega*be        ! only for Ex,Ey, 2018.10.05

!#[2]## cal ut*dA/dm*Al

 kk=0 ! 2018.06.21
 do imodel=1,nmodel !=========================== model loop start
   if ( iactive(imodel) .ne. 1 ) cycle ! 2018.06.22
   kk = kk + 1                         ! 2018.06.22
   AAL    = 0.d0                       ! 2017.09.03
   dbeobs = 0.d0                       ! 2018.10.05 either of bx,by,bz,ex,ey
   do jj=stack(imodel-1)+1,stack(imodel)! element loop 
     iele=item(jj)

     !#[2-1]## ! check the direction of edge, compared to the defined lines
     idirection(1:6)=1
     n6line(1:6) = l_line%n6line(iele,1:6)
     do j=1,6
       if ( n6line(j) .lt. 0 ) idirection(j)=-1
       n6line(j) = n6line(j)*idirection(j)
       end do
     do j=1,4
       elm_xyz(1:3,j)=h_mesh%xyz(1:3,h_mesh%n4(iele,j)) ! [km]
       end do
     call gradnodebasisfun(elm_xyz,gn,v) ! see fem_util.f90

     !#[2-2]## Second term from i * omega * mu * sigma * int{ sigma w cdot w }dv {Bsl}
     !# [4-1] ## assemble coefficient for i * omega*
     if ( h_mesh%n4flag(iele,1) .lt. 2 ) goto 99 ! in the case of not in land
     RM    = g_model%logrho_model(imodel) ! log10(rho)
     sigma = 10**(-RM)                     ! sigma = 1/10**(log10(rho)) = 10**(-M)
     !   write(*,*) "imodel=",imodel,"sigma=",sigma,"RM=",RM

     !#[2-3]##
     !# BB     = iunit*omega*dmu*10^(-M)  (where M=log10(rho))
     !# dBB/dm = iunit*omega*dmu*{10^-M)}*(-log10)
     dBBdm =iunit*omega*dmu*sigma*L0**2.d0*(-log(10.d0))

     !#[2-4] ## assemble scheme No.2 ( analytical assembly)
     S1(:,:)=(0.d0,0.d0)
     do i=1,6
      do j=1,6
       k=kl(i,1);l=kl(i,2) ; m=kl(j,1) ; n=kl(j,2) ! gn*gn [km^-2], intv [km^3], yy[km]
       yy =     intv(k,m,v)*inner(gn(:,l), gn(:,n))   & ! first term
       &	-  intv(k,n,v)*inner(gn(:,l), gn(:,m))   & ! second term
       &      -  intv(l,m,v)*inner(gn(:,k), gn(:,n))   & ! third term
       &      +  intv(l,n,v)*inner(gn(:,k), gn(:,m))     ! forth term
       S1(i,j)= yy*idirection(i)*idirection(j)*dBBdm  ! S1 [km*rad/s*S/m]
      end do
     end do
     AAL = 0.d0  ! 2017.09.03
     do i=1,6
       ii = n6line(i)
       do j=1,6
         !     write(*,101) "iele",iele,"i,j,",i,j,"bs(n6line)=",bs(n6line(j)),"S1(i,j)",S1(i,j)
         do k=1,nsr_inv    ! 2017.09.03
           AAL(i,k)=AAL(i,k) + S1(i,j)*bs(n6line(j),k) ! (dA/dm)*Al 2017.09.03
          end do            ! 2017.09.03
       end do
     end do
     !# uT * (-AAL)
     do i=1,utt%nrow  ! 2018.10.05
       do k=1,nsr_inv  ! 2017.09.03
	       do j=1,6        ! 2017.09.03
	         ii = n6line(j) ! 2017.09.03
	         dbeobs(i,k) = dbeobs(i,k) + utfull(ii,i)*(-AAL(j,k)) ! 2017.07.21
	       end do          ! 2017.09.03
	     end do          ! 2017.09.03
     end do           ! 2017.09.03
  end do ! element loop (iele) for i-thmodel

  call watchstart(t_watch) ! 2017.09.06
  !#[2-4]## cal ut*dA/dm*Al
  !# d/dm(|bz(iobs)|) = |bz(iobs)| Re{1/bz(iobs)*uT*(-dA/dm)*Al}

  do k=1,nsr_inv  ! 2017.09.03
   do i=1,nobs ! kk : model index
    dampdmfull(i,kk,k)=1./log(10.d0)*dreal(dbeobs(i,k)/be(i,k)) ! 2018.10.05
    dphadmfull(i,kk,k)=180./pi      *dimag(dbeobs(i,k)/be(i,k)) ! 2018.10.05
    !   write(*,100) "imodel=",imodel,"iobs=",i,"dbeobs=",dbeobs(i), &
   !    & "bz=",bz(i),"dbzdmfull", dbzdmfull(i,imodel)
   end do
  end do       ! nsr loop, 2017.09.03

 end do  ! model loop   (imodel)

!#[3]## full to real_crs_matrix
 do k=1,nsr_inv   ! 2017.09.03
  call conv_full2crs(dampdmfull(:,:,k),nobs,nmodelactive,dampdm(k),threshold)!2017.06.22
  call conv_full2crs(dphadmfull(:,:,k),nobs,nmodelactive,dphadm(k),threshold)!2018.06.22
 end do           ! 2017.09.03

! write(*,*) "dbzdm:" ! realcrs [nobs,nmodel]
! call realcrsout(dbzdm)

!#[4]## set output  2017.06.07
 g_apdm%nobs     = nobs
 do k=1,nsr_inv ! 2017.09.03
  g_apdm%dampdm(icomp,k) = dampdm(k) ! 2018.10.05
  g_apdm%dphadm(icomp,k) = dphadm(k) ! 2018.10.05
 end do

 end do ! comp loop end 2018.10.05

 call watchstop(t_watch)  ! see m_caltime.f90 2017.12.12

!write(*,'(a,i2,a,i2,a,f9.4,a)') " ### GENJACOBIAN1   END !! ###  ip =",ip," /",np," Time =",t_watch%time," [min]" !2020.09.18
write(*,'(a,i2,a,i2)') " ### GENJACOBIAN1     END !! ###  ip =",ip !" /",np," Time =",t_watch%time," [min]" !2020.09.18
return

!# error #
99 continue
   write(*,*) "GEGEGE air region is included in nmodel. imodel=",imodel
   write(*,*) "iele=",iele,"h_mesh%n4flag(iele,1)=",h_mesh%n4flag(iele,1)
   do j=1,4
    write(*,*) "node=",h_mesh%n4(iele,j),"elm_xyz=",elm_xyz(1:3,j)
   end do
stop

!# format #
100 format(a,i3,a,i3,a,2g15.7,a,2g15.7,a,g15.7)
101 format(a,i3,a,2i3,a,2g15.7,a,2g15.7)

end

!##################################################### gen jacobian
subroutine genjacobian1_mt(nobs_mt,nline,ut_mt,bs_mt,PT_mt,g_model,h_mesh,&
 &                         l_line,omega,g_mtdm,g_tipdm,g_param_joint,ip,np) !2023.12.26
 !# g_tipdm is added on 2023.12.26
 !# Explanation
   !#  Ax = f
   !#  dA/dm*x + Adx/dm = 0
   !#  dx/dm =  - A^-1 [dA/dm*x]
   !#  p^T dx/dm = - u^T [dA/dm*x] (Au = p)
 !# modules
   use matrix            ! see m_matrix.f90
   use modelpart         ! see m_modelpart.f90
   use line_type         ! see m_line_type.f90
   use constants         ! see m_constants.f90
   use outerinnerproduct ! see m_outerinnerproduct.f90
   use fem_util          ! for volume, intv, (see m_fem_utiil.f90 )
   use fem_edge_util     ! see fem_edge_util.f90
   use mesh_type         ! see m_mesh_type.f90
   use caltime           ! see m_caltime.f90 2017.09.06
   use param_jointinv       ! see m_param_jointinv.f90  2018.10.05
   implicit none
 !# declaration
   integer(4),              intent(in)     :: ip,np          ! 2020.09.18
   type(param_joint),       intent(in)     :: g_param_joint  ! 2018.10.05
   integer(4),              intent(in)     :: nobs_mt,nline
   type(complex_crs_matrix),intent(in)     :: ut_mt(4)       ! [nobs,nline]  
   complex(8),              intent(in)     :: bs_mt(nline,2) !  2017.09.03
   type(real_crs_matrix),   intent(in)     :: PT_mt(4)       ! [nobs,nline] 2018.10.05
   type(model),             intent(in)     :: g_model
   type(line_info),         intent(in)     :: l_line
   type(mesh),              intent(in)     :: h_mesh
   real(8),                 intent(in)     :: omega
   type(mt_dm),             intent(inout)  :: g_mtdm          ! 2022.01.05
   type(tip_dm),            intent(inout)  :: g_tipdm         ! 2023.12.25
 !# declaration for internal variables 
   type(complex_crs_matrix)                :: utt   ! 2018.10.05
   type(complex_crs_matrix)                :: dzxx  ! 2017.09.03 (->g_mtdm)
   type(complex_crs_matrix)                :: dzxy  ! 2017.09.03 (->g_mtdm)
   type(complex_crs_matrix)                :: dzyx  ! 2017.09.03 (->g_mtdm)
   type(complex_crs_matrix)                :: dzyy  ! 2017.09.03 (->g_mtdm)
   type(complex_crs_matrix)                :: dtx   ! 2024.08.30
   type(complex_crs_matrix)                :: dty   ! 2024.08.30
   complex(8),    dimension(nobs_mt,2)       :: be               ! 2022.01.05
   integer(4)                              :: nmodel,nphys2
   integer(4)                              :: imodel,i,j,k,l,m,n,ii,jj,kk!2018.06.21
   integer(4)                              :: iele,idirection(6),n6line(6)
   integer(4),allocatable,dimension(:)     :: stack,item
   real(8),   allocatable,dimension(:)     :: logrho_model
   complex(8),allocatable,dimension(:,:)   :: AAL,dbeobs       ! 2017.09.03
   complex(8),allocatable,dimension(:,:)   :: dzxxdm       ! 2017.09.03
   complex(8),allocatable,dimension(:,:)   :: dzxydm       ! 2017.09.03
   complex(8),allocatable,dimension(:,:)   :: dzyxdm       ! 2017.09.03
   complex(8),allocatable,dimension(:,:)   :: dzyydm       ! 2017.09.03
   complex(8),allocatable,dimension(:,:)   :: dtxdm       ! 2024.08.29
   complex(8),allocatable,dimension(:,:)   :: dtydm       ! 2024.08.29
   real(8)                                 :: threshold =1.d-10
   integer(4)                              :: nmodelactive     ! 2018.06.21
   integer(4),allocatable,dimension(:)     :: iactive          ! 2018.06.21
   integer(4),            dimension(5)     :: iflag_comp       ! 2018.10.05
   integer(4)                              :: icomp            ! 2018.10.05
 !# declaration for dA/dm
   complex(8)                              :: S1(6,6)
   real(8)                                 :: elm_xyz(3,4),yy,gn(3,4),v
   real(8)                                 :: RM,sigma         ! model = log10(rho)
   complex(8)                              :: dBBdm, iunit=(0.d0,1.d0)
   real(8),   parameter                    :: L0=1.d+3         ! [m]  scale length
   complex(8),allocatable,dimension(:,:)   :: utfull           ! 2017.09.03
   complex(8)                              :: z                ! 2017.09.03
   type(watch)                             :: t_watch          ! see m_caltime.f90
   complex(8) :: D,Dinv,Dinv2,hx1,hx2,hy1,hy2,ex1,ex2,ey1,ey2 ! 2022.01.05
   complex(8) :: dhx1,dhx2,dhy1,dhy2,dex1,dex2,dey1,dey2
   complex(8),dimension(nobs_mt,2 )        :: exo,eyo,hxo,hyo
   complex(8),allocatable,dimension(:,:,:) :: dexo,deyo,dhxo,dhyo ! 2022.12.06
   integer(4),parameter :: nsr_mt=2
 !# for addition of Tip
   logical                                 :: TIP=.false.      ! 2023.12.26
   integer(4)                              :: ncomp

 !#[-1]## allocation and watch start
   call watchstart(t_watch) ! see m_caltime.f90 2017.12.22
   allocate(dexo(nobs_mt,2,g_model%nmodelactive)) ! 2022.12.06
   allocate(deyo(nobs_mt,2,g_model%nmodelactive)) ! 2022.12.06
   allocate(dhxo(nobs_mt,2,g_model%nmodelactive)) ! 2022.12.06
   allocate(dhyo(nobs_mt,2,g_model%nmodelactive)) ! 2022.12.06
   if ( g_param_joint%iflag_tipper .eq. 1 ) TIP=.true. ! 2023.12.26

 !#[0]## set integers and allocate
   nmodel       = g_model%nmodel
   nmodelactive = g_model%nmodelactive ! 2018.06.21
   allocate( iactive(nmodel) )         ! 2018.06.21
   iactive      = g_model%iactive      ! 2018.06.21
   nphys2       = g_model%nphys2
   allocate(stack(0:nmodel),item(nphys2))
   allocate(logrho_model(nmodel))
   stack        = g_model%model2ele%stack
   item         = g_model%model2ele%item  ! element id for whole element space
   logrho_model = g_model%logrho_model
   allocate(utfull(nline,nobs_mt))        ! 2018.10.05
   allocate( AAL(6,nsr_mt) )
   allocate( dbeobs(nobs_mt,2)            )       ! 2018.10.05
   allocate( dzxxdm(nobs_mt,nmodelactive) ) ! 2018.06.22
   allocate( dzxydm(nobs_mt,nmodelactive) ) ! 2018.06.22
   allocate( dzyxdm(nobs_mt,nmodelactive) ) ! 2018.06.22
   allocate( dzyydm(nobs_mt,nmodelactive) ) ! 2018.06.22

 !#[1]## cal bx,by,ex,ey
 ncomp = 4 ! bx,by,ex,ey       2023.12.26
 if (TIP) ncomp = 5 ! 5 is bz  2023.12.26
 do icomp = 1,4 
   utt    = ut_mt(icomp)                 ! 2018.10.05
   utfull = 0.d0                         ! 2017.09.03
   do   i = 1,utt%nrow                   ! 2017.09.03
     do j=utt%stack(i-1)+1,utt%stack(i)  ! 2017.09.03
       utfull(utt%item(j),i) = utt%val(j)! 2017.09.03
     end do                              ! 2017.09.03
   end do                                ! 2017.09.03
   if ( icomp .ge. 3 ) utfull = -iunit*omega*utfull ! only for Ex, Ey

   !#[1-1]## cal either of bx,by,ex,ey dependent on icomp
     be = 0.d0
     do i=1,nsr_mt                     ! nsr_mt = 2 2022.01.12
       call mul_matcrs_cv(PT_mt(icomp),bs_mt(:,i),nline,be(:,i))!be[nobs_mt,2] 2018.10.05
     end do
     if ( icomp .ge. 3 ) be    = -iunit*omega*be        ! only for electric field
  
   !#[1-2]## cal ut*dA/dm*Al
     kk=0
     do imodel=1,nmodel
       if ( iactive(imodel) .ne. 1 ) cycle ! 2018.06.22
       kk = kk + 1                         ! 2018.06.22
       AAL    = 0.d0                       ! 2017.09.03
       dbeobs = 0.d0                       ! 2018.10.05 either of bx,by,bz,ex,ey
       do jj=stack(imodel-1)+1,stack(imodel) ! element loop for i-th model
         iele=item(jj)

         !#[1-2-1]## ! check the direction of edge, compared to the defined lines
           idirection(1:6)=1
           n6line(1:6) = l_line%n6line(iele,1:6)
           do j=1,6
             if ( n6line(j) .lt. 0 ) idirection(j)=-1
             n6line(j) = n6line(j)*idirection(j)
           end do
           do j=1,4
             elm_xyz(1:3,j)=h_mesh%xyz(1:3,h_mesh%n4(iele,j)) ! [km]
           end do
           call gradnodebasisfun(elm_xyz,gn,v) ! see fem_util.f90

         !#[1-2-2]## Second term from i * omega * mu * sigma * int{ sigma w cdot w }dv {Bsl}
           !# [4-1] ## assemble coefficient for i * omega*
           if ( h_mesh%n4flag(iele,1) .lt. 2 ) goto 99 ! in the case of not in land
           RM    = g_model%logrho_model(imodel) ! log10(rho)
           sigma = 10**(-RM) ! sigma = 1/10**(log10(rho)) = 10**(-M)
           !   write(*,*) "imodel=",imodel,"sigma=",sigma,"RM=",RM
  
         !#[1-2-3]## dBB/dm
           !# BB     = iunit*omega*dmu*10^(-M)  (where M=log10(rho))
           !# dBB/dm = iunit*omega*dmu*{10^-M)}*(-log10)
           dBBdm =iunit*omega*dmu*sigma*L0**2.d0*(-log(10.d0))

         !#[2-4] ## assemble scheme No.2 ( analytical assembly)
           S1(:,:)=(0.d0,0.d0)
           do i=1,6
             do j=1,6
               k=kl(i,1);l=kl(i,2) ; m=kl(j,1) ; n=kl(j,2) ! gn*gn [km^-2], intv [km^3], yy[km]
               yy =      intv(k,m,v)*inner(gn(:,l), gn(:,n))   & ! first term
               &	    -  intv(k,n,v)*inner(gn(:,l), gn(:,m))   & ! second term
               &      -  intv(l,m,v)*inner(gn(:,k), gn(:,n))   & ! third term
               &      +  intv(l,n,v)*inner(gn(:,k), gn(:,m))     ! forth term
               S1(i,j)= yy*idirection(i)*idirection(j)*dBBdm  ! S1 [km*rad/s*S/m]
             end do
           end do
           AAL = 0.d0  ! 2017.09.03
           do i=1,6
             ii = n6line(i)
             do j=1,6
               do k=1,nsr_mt    ! 2017.09.03
                 AAL(i,k)=AAL(i,k) + S1(i,j)*bs_mt(n6line(j),k) ! (dA/dm)*Al 2017.09.03
               end do            ! 2017.09.03
             end do
           end do
           !# uT * (-AAL)
           do i=1,utt%nrow  ! 2018.10.05 obs loop
             do k=1,nsr_mt  ! 2022.01.05  src loop
               do j=1,6        ! 2017.09.03  inside tetrahedron
                 ii = n6line(j) ! 2017.09.03
                 dbeobs(i,k) = dbeobs(i,k) + utfull(ii,i)*(-AAL(j,k)) ! i: obs, k: MT pol. for model kk
               end do          ! 2017.09.03
             end do          ! 2017.09.03
           end do           ! 2017.09.03           
       end do ! element loop (iele) for i-thmodel elemnt loop

       if (     icomp .eq. 1) then ! bx, dbx
         hxo(:,1:2)     = be(:,1:2) 
         dhxo(:,1:2,kk) = dbeobs(:,1:2)
       else if (icomp .eq. 2) then ! by, dby
         hyo(:,1:2)     = be(:,1:2)
         dhyo(:,1:2,kk) = dbeobs(:,1:2)
       else if (icomp .eq. 3) then ! ex, dex
         exo(:,1:2)     = be(:,1:2)
         dexo(:,1:2,kk) = dbeobs(:,1:2)
       else if (icomp .eq. 4) then ! ey, dey
         eyo(:,1:2)     = be(:,1:2)
         deyo(:,1:2,kk) = dbeobs(:,1:2)
       end if
     end do ! model loop end
 end do ! icomp end

 call watchstart(t_watch) ! 2017.09.06

 !#[2-4]## cal ut*dA/dm*Al
 !# d/dm(|bz(iobs)|) = |bz(iobs)| Re{1/bz(iobs)*uT*(-dA/dm)*Al}
 if (.false.) then
  write(*,*) "hx 1, 2"
  do i=1,nobs_mt
   write(*,'(2g15.7)') hxo(i,1),hxo(i,2) 
  end do
  write(*,*) "hy 1, 2"
  do i=1,nobs_mt
   write(*,'(2g15.7)') hyo(i,1),hyo(i,2) 
  end do
 end if

 kk=0
 do imodel = 1,nmodel
   if ( iactive(imodel) .ne. 1 ) cycle ! 2018.06.22
   kk = kk + 1                         ! 2018.06.22
   do i=1,nobs_mt
     hx1 = hxo(i,1) ; hy1 = hyo(i,1)
     hx2 = hxo(i,2) ; hy2 = hyo(i,2)
     ex1 = exo(i,1) ; ey1 = eyo(i,1)
     ex2 = exo(i,2) ; ey2 = eyo(i,2)
     dhx1 = dhxo(i,1,kk) ; dhy1 = dhyo(i,1,kk)
     dhx2 = dhxo(i,2,kk) ; dhy2 = dhyo(i,2,kk)
     dex1 = dexo(i,1,kk) ; dey1 = deyo(i,1,kk)
     dex2 = dexo(i,2,kk) ; dey2 = deyo(i,2,kk)
     if ( i .eq. 1 .and. imodel .eq. 1 .and. ip .eq. 4 .and. .false. ) then
       !write(*,*) "hx1",hx1 ! 0
       write(*,*) "hx2",hx2 !
       write(*,*) "hy1",hy1 !
       !write(*,*) "hy1",hy1 ! 0
       write(*,*) "ex1",ex1
       !write(*,*) "ex2",ex2 ! 0
       !write(*,*) "ey1",ey1 ! 0
       write(*,*) "ey2",ey2
       write(*,*) "Zxy",ex1/hy1
       write(*,*) "Zyx",ey2/hx2
       write(*,*) ""
       !write(*,*) "dhx1",dhx1 ! 0
       write(*,*) "dhx2",dhx2
       write(*,*) "dhy1",dhy1
       !write(*,*) "dhy2",dhy2 ! 0
       write(*,*) "dex1",dex1
       !write(*,*) "dex2",dex2 ! 0
       !write(*,*) "dey1",dey1 ! 0
       write(*,*) "dey2",dey2
     end if
     !# generate dzxx, dzxy, dzyx, dzyy

     D = hx1*hy2-hx2*hy1
     !   write(*,*) "D",D,"hx1",hx1,"hy2",hy2
     Dinv=1./D
     Dinv2=Dinv**2.

     dzxxdm(i,kk)     =    Dinv*(hy2*dex1 + ex1*dhy2 - hy1*dex2 - ex2*dhy1) &
     & -(ex1*hy2-ex2*hy1)*Dinv2*(hy2*dhx1 + hx1*dhy2 - hy1*dhx2 - hx2*dhy1) ! ok

     dzxydm(i,kk)     =    Dinv*(hx1*dex2 + ex2*dhx1 - hx2*dex1 - ex1*dhx2) &
     & -(ex2*hx1-ex1*hx2)*Dinv2*(hy2*dhx1 + hx1*dhy2 - hy1*dhx2 - hx2*dhy1) ! ok

     dzyxdm(i,kk)     =     Dinv*(hy2*dey1 + ey1*dhy2 - hy1*dey2 - ey2*dhy1)&
     & -(ey1*hy2-ey2*hy1)*Dinv2*(hy2*dhx1 + hx1*dhy2 - hy1*dhx2 - hx2*dhy1) ! ok

     dzyydm(i,kk)     =     Dinv*(hx1*dey2 + ey2*dhx1 - hx2*dey1 - ey1*dhx2) &
     & -(ey2*hx1-ey1*hx2)*Dinv2*(hy2*dhx1 + hx1*dhy2 - hy1*dhx2 - hx2*dhy1) ! ok
   
    if (TIP) then

    end if

   end do
 end do ! imodel loop

 !#[3]## full to complex_crs_matrix
   call conv_full2crs_complex(dzxxdm,nobs_mt,nmodelactive,dzxx,threshold)!2022.01.05
   call conv_full2crs_complex(dzxydm,nobs_mt,nmodelactive,dzxy,threshold)!2022.01.05
   call conv_full2crs_complex(dzyxdm,nobs_mt,nmodelactive,dzyx,threshold)!2022.01.05
   call conv_full2crs_complex(dzyydm,nobs_mt,nmodelactive,dzyy,threshold)!2022.01.05
   if (TIP) then !2023.12.26
     call conv_full2crs_complex(dtxdm,nobs_mt,nmodelactive,dtx,threshold)!2023.12.26
     call conv_full2crs_complex(dtydm,nobs_mt,nmodelactive,dty,threshold)!2023.12.26
   end if        !2023.12.26

 !#[4]## set output
   g_mtdm%nobs_mt = nobs_mt
   g_mtdm%dzxxdm = dzxx
   g_mtdm%dzxydm = dzxy
   g_mtdm%dzyxdm = dzyx
   g_mtdm%dzyydm = dzyy
   if (TIP) then              ! 2023.12.26
    g_tipdm%nobs_mt = nobs_mt ! 2023.12.26
   end if                     ! 2023.12.26

 call watchstop(t_watch)

 !write(*,'(a,i2,a,i2,a,f9.4,a)') " ### GENJACOBIAN1   END !! ###  ip =",ip," /",np," Time =",t_watch%time," [min]" !2020.09.18
 write(*,'(a,i2,a,i2)') " ### GENJACOBIAN1_MT  END !! ###  ip =",ip !" /",np," Time =",t_watch%time," [min]" !2020.09.18

 return

 !# error #
 99 continue
 write(*,*) "GEGEGE air region is included in nmodel. imodel=",imodel
 write(*,*) "iele=",iele,"h_mesh%n4flag(iele,1)=",h_mesh%n4flag(iele,1)
 do j=1,4
  write(*,*) "node=",h_mesh%n4(iele,j),"elm_xyz=",elm_xyz(1:3,j)
 end do
 stop

 !# format #
   100 format(a,i3,a,i3,a,2g15.7,a,2g15.7,a,g15.7)
   101 format(a,i3,a,2i3,a,2g15.7,a,2g15.7)

 end
!############################################################# genjacobian2
!# modified on 2018.10.05 for multiple components
!# nmodelactive is introduced on 2018.06.25
!# modified on 2017.09.04 for multiple source inversion
!# Coded on 2017.05.17
subroutine genjacobian2(g_param_joint,nfreq,g_apdm,JJ) ! 2018.06.25
use param_jointinv ! 2017.09.04
use matrix
use caltime     ! 2017.12.22
implicit none
type(param_joint),   intent(in)    :: g_param_joint   ! 2017.09.04
integer(4),                intent(in)    :: nfreq           ! 2017.09.04
type(amp_phase_dm),        intent(in)    :: g_apdm(nfreq)   ! 2017.06.07
type(real_crs_matrix),     intent(out)   :: JJ              ! Jacobian matrix [ndat,nmodel]
type(real_crs_matrix),allocatable,dimension(:,:,:,:)   :: dapdm       !(2,5,nsr,nfreq) [nobs,nmodel]
logical,              allocatable,dimension(:,:,:,:,:) :: data_avail  ! 2018.10.05
integer(4)                               :: ii,i,nmodel,ntotr,nsr_inv, nobs    ! 2017.09.04
integer(4)                               :: j,k,l,nsft,lsft,ncount,ndat        ! 2017.09.04
integer(4)                               :: nmodelactive,icomp    ! 2018.10.05
integer(4),                dimension(5)  :: iflag_comp            ! 2018.10.05
type(watch) :: t_watch  ! see m_caltime.f90 2017.12.22

 call watchstart(t_watch) ! see m_caltime.f90 2017.12.22

!#[0]## set
 ndat           = g_param_joint%ndat          ! 2017.09.04
 nobs           = g_param_joint%nobs          ! 2017.09.04
 nsr_inv        = g_param_joint%nsr_inv       ! 2017.09.04
 nmodelactive   = g_apdm(1)%nmodelactive      ! 2018.06.25
 iflag_comp     = g_param_joint%iflag_comp    ! 2018.10.05
 allocate(dapdm(2,5,nsr_inv,nfreq))           ! 2018.10.05
 allocate(data_avail(2,5,nfreq,nobs,nsr_inv)) !2018.10.05
 data_avail = g_param_joint%data_avail      ! 2017.09.04
 do k=1,nsr_inv ! 2017.09.04
  do i=1,nfreq  ! 2017.06.07
   do icomp = 1,5 ! 2018.10.05
    if ( iflag_comp(icomp) .eq. 0 ) cycle          ! 2018.10.05
    dapdm(1,icomp,k,i) = g_apdm(i)%dampdm(icomp,k) ! 2017.10.05
    dapdm(2,icomp,k,i) = g_apdm(i)%dphadm(icomp,k) ! 2018.10.05
   end do         ! 2018.10.05
  end do
 end do

!#[1]## count ntotr; total number of elements in JJ
 ntotr = 0
 do i = 1,nfreq
  do k= 1, nsr_inv  ! 2017.09.04
   do j = 1,nobs    ! 2017.09.04
    do icomp = 1,5  ! 2018.10.05
     do l = 1,2     ! 2017.09.04
      if (data_avail(l,icomp,i,j,k)) then    ! 2018.10.05
       ntotr = ntotr + dapdm(l,icomp,k,i)%stack(j) - dapdm(l,icomp,k,i)%stack(j-1) ! 2018.10.05
      end if
     end do
    end do          ! 2018.10.05
   end do
  end do
 end do

!#[2]## assemble J
 JJ%nrow  = ndat         ! 2017.09.04
 JJ%ntot  = ntotr        ! 2017.09.04
 JJ%ncolm = nmodelactive ! 2018.06.25
 allocate(JJ%stack(0:JJ%nrow),JJ%item(ntotr),JJ%val(ntotr))
 JJ%stack(:)=0
 ii = 0            ! 2017.09.04
 do i =1,nfreq
  do k=1,nsr_inv   ! 2017.09.04
   do j = 1,nobs   ! 2017.09.04
    do icomp = 1,5 ! 2018.10.5
     do l = 1,2    ! 2017.09.04 amp phase
      if ( data_avail(l,icomp,i,j,k) ) then  ! 2018.10.05
       ii = ii + 1 ! current row index
	 ncount = dapdm(l,icomp,k,i)%stack(j) - dapdm(l,icomp,k,i)%stack(j-1)
	 JJ%stack(ii) = JJ%stack(ii-1) + ncount ! 2017.09.04
	 nsft = JJ%stack(ii-1)                  ! 2017.09.04
	 lsft = dapdm(l,icomp,k,i)%stack(j-1)         ! 2017.09.04
	 JJ%item(nsft+1:nsft+ncount) = dapdm(l,icomp,k,i)%item(lsft+1:lsft+ncount)
	 JJ%val( nsft+1:nsft+ncount) = dapdm(l,icomp,k,i)%val( lsft+1:lsft+ncount)
     end if
    end do ! l     loop 2017.09.04
    end do ! icomp loop 2018.10.05
   end do  ! j     loop 2017.09.04
  end do   ! k     loop 2017.09.04
 end do    ! i     loop 2017.09.04
 if ( ii .ne. ndat) then
  write(*,*) "GEGEGE ii",ii,"should be ndat",ndat,"!!!"
  stop
 end if

!#[3]## set output
 call watchstop(t_watch)  ! see m_caltime.f90 2017.12.12
 write(*,'(a,f8.4,a)') " ### GENJACOBIAN2 END!! ### Time =",t_watch%time," [min]"!2020.09.18
return
end
!############################################################# genjacobian2
!# modified on 2018.10.05 for multiple components
!# nmodelactive is introduced on 2018.06.25
!# modified on 2017.09.04 for multiple source inversion
!# Coded on 2017.05.17
subroutine genjacobian2_MT(g_param_joint,nfreq_mt,g_mtdm,JJ_mt) ! 2018.06.25
  use param_jointinv ! 2017.09.04
  use matrix
  use caltime     ! 2017.12.22
  implicit none
  type(param_joint),       intent(in)    :: g_param_joint   ! 2017.09.04
  integer(4),              intent(in)    :: nfreq_mt        ! 2017.09.04
  type(mt_dm),             intent(in)    :: g_mtdm(nfreq_mt)! 2017.06.07
  type(real_crs_matrix),   intent(out)   :: JJ_mt           ! Jacobian matrix [ndat,nmodel]
  !
  type(complex_crs_matrix),dimension(4,nfreq_mt) :: dz       !(2,5,nsr,nfreq) [nobs,nmodel]
  logical,              allocatable,dimension(:,:,:,:) :: data_avail_mt  ! 2018.10.05
  integer(4)                               :: ii,i,nmodel,ntotr,nsr_inv, nobs_mt   ! 2017.09.04
  integer(4)                               :: j,k,l,nsft,lsft,ncount,ndat_mt     ! 2017.09.04
  integer(4)                               :: nmodelactive,icomp    ! 2018.10.05
  integer(4) :: ireal,iimag,nsft1,nsft2,lsft1,lsft2
  type(watch) :: t_watch  ! see m_caltime.f90 2017.12.22
  
   call watchstart(t_watch) ! see m_caltime.f90 2017.12.22
  
  !#[0]## set
   ndat_mt        = g_param_joint%ndat_mt     ! 2017.09.04
   nobs_mt        = g_param_joint%nobs_mt     ! 2017.09.04
   nmodelactive   = g_mtdm(1)%nmodelactive    ! 2018.06.25
   write(*,*) "genjacobian2_mt ndat_mt",ndat_mt !           2022.01.05
   write(*,*) "genjacobian2_mt nobs_mt",nobs_mt !           2022.01.05
   write(*,*) "genjacobian2_mt nmodelactive",nmodelactive ! 2022.01.05

   allocate(data_avail_mt(2,4,nobs_mt,nfreq_mt)) !2018.10.05
   data_avail_mt  = g_param_joint%data_avail_mt   ! 2017.09.04
   do i=1,nfreq_mt  ! 2017.06.07
      dz(1,i) = g_mtdm(i)%dzxxdm ! 2017.10.05
      dz(2,i) = g_mtdm(i)%dzxydm ! 2018.10.05
      dz(3,i) = g_mtdm(i)%dzyxdm ! 2017.10.05
      dz(4,i) = g_mtdm(i)%dzyydm ! 2018.10.05
   end do
  
   open(1,file="dzxxdm1Hz.dat")
    call complexcrsout(dz(1,1),1) !m_matrix.f90
   close(1)
   open(1,file="dzxydm1Hz.dat")
    call complexcrsout(dz(2,1),1) !m_matrix.f90
   close(1)
   open(1,file="dzyxdm1Hz.dat")
    call complexcrsout(dz(3,1),1) !m_matrix.f90
   close(1)

  !#[1]## count ntotr; total number of elements in JJ
   ntotr = 0
   do i = 1,nfreq_mt
     do j = 1,nobs_mt    ! 2017.09.04
       do icomp=1,4    ! zxx,zxy,zyx,zyy
           if ( data_avail_mt(1,icomp,j,i) ) then    ! 2018.10.05
             ntotr = ntotr + 2*(dz(icomp,i)%stack(j) - dz(icomp,i)%stack(j-1)) ! 2 for real and imag
           end if
        end do
      end do
    end do          ! 2018.10.05
  write(*,*) "genjacobian2_mt ntotr",ntotr

  !#[2]## assemble J
   JJ_mt%nrow  = ndat_mt      ! 2017.09.04
   JJ_mt%ntot  = ntotr        ! 2017.09.04
   JJ_mt%ncolm = nmodelactive ! 2018.06.25
   allocate(JJ_mt%stack(0:JJ_mt%nrow),JJ_mt%item(ntotr),JJ_mt%val(ntotr))
   JJ_mt%stack(:)=0
   ii = 0            ! 2017.09.04
   do i =1,nfreq_mt
     do j = 1,nobs_mt   ! 2017.09.04
      do icomp = 1,4 ! 2018.10.5
        if ( data_avail_mt(1,icomp,j,i) ) then  ! 2018.10.05
         ireal=ii+1
         iimag=ii+2
         ncount = dz(icomp,i)%stack(j) - dz(icomp,i)%stack(j-1)
         JJ_mt%stack(ireal) = JJ_mt%stack(ireal-1) + ncount ! 2017.09.04
         JJ_mt%stack(iimag) = JJ_mt%stack(iimag-1) + ncount ! 2017.09.04
         nsft1 = JJ_mt%stack(ireal-1)      ! 2022.01.05
         nsft2 = JJ_mt%stack(ireal)        ! 2022.01.05
         lsft  = dz(icomp,i)%stack(j-1)    ! 2022.01.05
         !
         JJ_mt%item(nsft1+1:nsft1+ncount) = dz(icomp,i)%item(lsft+1:lsft+ncount) ! ireal
         JJ_mt%item(nsft2+1:nsft2+ncount) = dz(icomp,i)%item(lsft+1:lsft+ncount)
         !
         JJ_mt%val( nsft1+1:nsft1+ncount) = real(dz(icomp,i)%val( lsft+1:lsft+ncount) )
         JJ_mt%val( nsft2+1:nsft2+ncount) = imag(dz(icomp,i)%val( lsft+1:lsft+ncount) )

         ii = ii + 2 ! current row index
       end if
      end do ! icomp loop 2018.10.05
     end do  ! j     loop 2017.09.04
   end do    ! i     loop 2017.09.04
   if ( ii .ne. ndat_mt) then
    write(*,*) "GEGEGE ii",ii,"should be ndat_mt",ndat_mt,"!!!"
    stop
   end if
  
  open(1,file="JJ_mt_in.dat")
  call realcrsout(JJ_mt,1)
  close(1)

  !#[3]## set output
   call watchstop(t_watch)  ! see m_caltime.f90 2017.12.12
   write(*,'(a,f8.4,a)') " ### GENJACOBIAN2_MT END!! ### Time =",t_watch%time," [min]"!2020.09.18
  return
  end
  
!######################################################### OUTJACOB
!# modified for multiple components on 2018.10.05
!# coded on 2018.06.25
subroutine OUTJACOB(g_param_joint,g_data_ap,JJ,ite,g_model,g_param,sparam) ! 2022.11.01
use param_jointinv
use matrix
use modelpart
use outresp  ! for free_unit ./common/m_outresp.f90
use param
implicit none
type(real_crs_matrix),    intent(in)        :: JJ
type(param_joint),  intent(in)              :: g_param_joint
type(data_vec_ap),        intent(in)        :: g_data_ap     ! 2022.11.01
type(param_source),       intent(in)        :: sparam
integer(4),               intent(in)        :: ite
type(model),              intent(in)        :: g_model
type(param_forward),      intent(in)        :: g_param
character(50)                               :: head,outfile,site,sour
character(2)                                :: num,nf
integer(4)                                  :: i,j,k,l,m,n,ii ! 2018.06.25
logical,   allocatable,dimension(:,:,:,:,:) :: data_avail        ! 2018.10.05
integer(4)                                  :: nfreq,nobs,nsr_inv
integer(4)                                  :: nmodel,nmodelactive
integer(4),allocatable,dimension(:)         :: iactive      ! 2018.06.25
integer(4),allocatable,dimension(:)         :: ptrnmodel    ! 2018.06.25
integer(4)                                  :: idev,nsi,nso ! 2018.10.05
integer(4),allocatable,dimension(:,:,:,:,:) :: icount       ! 2018.10.05
real(8),   allocatable,dimension(:,:,:)     :: JP           ! 2018.06.25
real(8),   allocatable,dimension(:)         :: freq
integer(4),allocatable,dimension(:)         :: srcindex
integer(4)                                  :: icomp        ! 2018.10.05
integer(4),            dimension(5)         :: iflag_comp   ! 2018.10.05
real(8),   allocatable,dimension(:)         :: modelvolume ! 2022.11.01
real(8),   allocatable,dimension(:)         :: error       ! 2022.11.01
character(5)   :: fre     ! 2022.11.01
character(100) :: asdfile ! 2022.11.01
real(8) ::   er1,er2,asd  ! 2022.11.01

!#[0]## set
 nfreq        = g_param_joint%nfreq           ! 2018.06.25
 nobs         = g_param_joint%nobs            ! 2018.06.25
 nsr_inv      = g_param_joint%nsr_inv         ! 2018.06.25
 head         = g_param_joint%outputfolder    ! 2018.06.25
 nmodel       = g_model%nmodel                ! 2018.06.25
 nmodelactive = g_model%nmodelactive          ! 2018.06.25
 iactive      = g_model%iactive               ! iactive(nmodel)  allocate and fill
 write(num,'(i2.2)') ite                      ! 2017.07.14
 srcindex     = g_param_joint%srcindex        ! srcindex(nsr_inv) allocate and fill
 freq         = g_param%freq                  ! freq(nfreq) allocate and fill
 data_avail   = g_param_joint%data_avail      ! 2018.06.25 (2,5,nfreq,nobs,nsr_inv) allocate and fill
 modelvolume  = g_model%modelvolume           ! (nmodel) 2022.11.01
 error        = g_data_ap%error               ! (ndat)   2022.11.01
 if ( nmodelactive .ne. JJ%ncolm ) then
  write(*,*) "GEGEGE nmodelactive",nmodelactive,"JJ%ncolm",JJ%ncolm
 end if

!#[1]## set ptractive
 allocate( ptrnmodel(nmodelactive))           ! 2018.06.25
 ii=0              ! 2018.06.25
 do i=1,nmodel     ! 2018.06.25
  if ( iactive(i) .eq. 1 ) then ! 2018.06.25
   ii = ii + 1     ! 2018.06.25
   ptrnmodel(ii)=i ! 2018.06.25
  end if           ! 2018.06.25
 end do            ! 2018.06.25
 iflag_comp = g_param_joint%iflag_comp ! 2018.10.05

!#[2]## output
 allocate(JP(2,nfreq,nmodel))   ! 2018.06.25

!#[3]## generate icount
 ii = 0 ! 2018.06.25
 allocate(icount(2,5,nfreq,nobs,nsr_inv)) ! 2018.10.05
 icount = g_data_ap%idata    ! 2022.11.01
   !do i=1,nfreq
   !do k=1,nsr_inv
   !do l=1,nobs
   ! do icomp = 1,5                      ! 2018.10.05
   ! do m=1,2                            ! 2018.06.25
   ! if (data_avail(m,icomp,i,l,k)) then ! 2018.10.05
   !  ii = ii + 1
   !  icount(m,icomp,i,l,k) = ii         ! 2018.10.05
   ! end if
   !end do                               ! 2018.06.25
   !end do                               ! 2018.10.05
   !end do
   !end do
   !end do


!#[4]## generate JP and output
do icomp = 1,5                       ! 2018.10.05
  if (iflag_comp(icomp) .eq. 0) cycle ! 2018.10.05

  do k=1,nsr_inv
   do l=1,nobs
     site    = g_param%obsname(l)
     sour    = sparam%sourcename(srcindex(k))
     outfile = trim(head)//"Jacob"//trim(site)//"_"//trim(sour)//"_"//comp(icomp)//num(1:2)//".dat"!2018.10.05
     idev    = free_unit()   ! see m_outresp.f90

     !# gen JP
     JP = 0.d0 ! 2018.06.25
     do m=1,2  ! 1:amp, 2:phase
       do i=1,nfreq
         if ( data_avail(m,icomp,i,l,k) ) then        ! 2018.10.05
           ii = icount(m,icomp,i,l,k)                  ! 2018.10.05
           do j=JJ%stack(ii-1)+1,JJ%stack(ii)
             JP(m,i,ptrnmodel(JJ%item(j))) = JJ%val(j)  ! 2018.06.25
           end do
         end if
       end do
     end do    ! 2018.06.25


     !# output JP
     if ( .false. ) then ! 2022.11.01
       open(idev,file=outfile)
         write(idev,'(8g15.7)') freq(1:nfreq)
         write(idev,*) nmodel
         do m=1,2      ! 2018.06.25 1:amp, 2:phase
           do n=1,nmodel
             write(idev,'(8g15.7)') JP(m,1:nfreq,n)
           end do
         end do        ! 2018.06.25
       close(idev)
       !
     else               ! 2022.11.01
       !# output ASD      2022.11.01
       do i=1,nfreq
         write(fre,'(i5.5)') int(freq(i))       
         asdfile = trim(head)//"ASD"//trim(site)//"_"//trim(sour)//"_"//comp(icomp)//trim(fre)//"Hz_"//num(1:2)//".dat"
       
         open(idev,file=asdfile)
           write(idev,*) nmodel
           er1 = error(icount(1,icomp,i,l,k)) ! amp error
           er2 = error(icount(2,icomp,i,l,k)) ! phase error
           do n=1,nmodel
             asd = (abs(JP(1,i,n))/er1 + abs(JP(2,i,n))/er2)/modelvolume(n) ![km]^-3
             write(idev,*) asd
           end do
           !
         close(idev)
       end do
      end if

    end do ! obs loop
  end do  ! src loop
end do   ! icomp loop   2018.10.05

 write(*,*) "### OUTJACOB END !! ###" ! 2022.11.01
return
end
!########################################################## getnewmodel_joint
!# nmodelactive is introduced on 2018.06.25
!# coded on 2017.05.17
!#  [C]beta = X
!#  where
!#  [C] = [J*Cm*JT + alpha*Cd]
!#   X  = d_obs - d(M) - J(M - M_ref)
!# ---
!#  M_ite+1 = M_ref + Cm*JT*beta
!# 
!#
subroutine getnewmodel_joint(JJ_ac,JJ_mt,g_model_ref,h_model,g_data,h_data,&
  &                          g_data_mt,h_data_mt,CM,CD_ac,CD_mt,alpha,g_param_joint)
  use modelpart
  use matrix
  use param_jointinv ! 2017.06.09
  use caltime ! 2017.12.22
  implicit none
  type(param_joint),    intent(in)  :: g_param_joint
  real(8),              intent(in)    :: alpha
  type(real_crs_matrix),intent(in)    :: JJ_ac       ! [ndat,nmodel]          2022.01.05
  type(real_crs_matrix),intent(in)    :: JJ_mt       ! [ndat_mt,nmodel]       2022.01.05
  type(model),          intent(in)    :: g_model_ref ! reference model
  type(model),          intent(inout) :: h_model     ! old -> new
  type(data_vec_ap),    intent(in)    :: g_data      ! observed
  type(data_vec_ap),    intent(in)    :: h_data      ! calculated
  type(data_vec_mt),    intent(in)    :: g_data_mt   ! observed               2022.01.05
  type(data_vec_mt),    intent(in)    :: h_data_mt   ! calculated             2022.01.05
  type(real_crs_matrix),intent(in)    :: CM          ! model covariance matrix
  type(real_crs_matrix),intent(in)    :: CD_ac       ! data covariance matrix 2022.01.05
  type(real_crs_matrix),intent(in)    :: CD_mt       ! data covariance matrix 2022.01.05
  type(real_crs_matrix)               :: C,JCMJT,ACD,crsout      ! 2018.01.23
  type(real_crs_matrix)               :: CMJT, JT    ! 2018.01.23
  type(real_crs_matrix)               :: JJ,CD       !                        2022.01.05
  real(8),allocatable,dimension(:)    :: X, dobs, dcal,JM,beta   ! [ndat]
  real(8),allocatable,dimension(:)    :: dobs_mt, dcal_mt        ! [ndat_mt]  2022.01.05
  real(8),allocatable,dimension(:)    :: dobs_ac, dcal_ac        ! [ndat_ac]  2022.01.05
  real(8),allocatable,dimension(:)    :: model1,model_ref,dmodel ! [nmodel]
  integer(4)                          :: ndat,nmodel,i,ii ! 2018.06.25
  integer(4)                          :: ndat_mt,ndat_ac !                    2022.01.05
  integer(4)                          :: nmodelactive ! 2018.06.25
  integer(4),allocatable,dimension(:) :: iactive      ! 2018.06.25
  integer(4)                          :: iboundflag,ijoint,ioutlevel ! 2022.10.14
  integer(4)                          :: icheck = 1
  type(model)                         :: ki,kiref
  type(model)                         :: m0,m_ref
  logical                             :: MT, ACT,TIPPER ! 2023.12.23
  type(watch) :: t_watch ! 2017.12.22
  type(watch) :: t_watch1 ! 2018.01.23
  
   call watchstart(t_watch) ! 2017.12.22
   write(*,'(a,f9.4)') " alpha =",alpha," is adopted in GETNEWMODEL"  ! 2020.09.29
  
  !#[0]## set
    call watchstart(t_watch1)
    ijoint       = g_param_joint%ijoint     ! 2022.10.14 1:ACTIVE,2:MT,3:Joint
    call setnec(ijoint,g_param_joint,ACT,MT,TIPPER) ! 2023.12.23 see m_param_joint.f90
    ioutlevel    = g_param_joint%ioutlevel  ! 2022.10.14
    ndat_ac      = g_data%ndat              !            2022.01.05
    ndat_mt      = g_data_mt%ndat_mt        !            2022.01.05
    nmodel       = g_model_ref%nmodel
    nmodelactive = g_model_ref%nmodelactive ! 2018.06.25
    allocate( iactive(nmodel) )             ! 2018.06.25
    iactive      = g_model_ref%iactive      ! 2018.06.25
    allocate(dobs_ac(ndat_ac),dcal_ac(ndat_ac)) ! 2022.01.05
    allocate(dobs_mt(ndat_mt),dcal_mt(ndat_mt)) !        2022.01.05
    dobs_ac      = g_data%dvec
    dcal_ac      = h_data%dvec
    dobs_mt      = g_data_mt%dvec_mt        !            2022.01.05
    dcal_mt      = h_data_mt%dvec_mt        !            2022.01.05
    allocate( model1(   nmodelactive) )     ! 2018.06.25
    allocate( model_ref(nmodelactive) )     ! 2018.06.25
    allocate( dmodel(   nmodelactive) )     ! 2018.06.25
    !#
    iboundflag = g_param_joint%iboundflag   ! 0:off, 1:on, 2:transfer

    call watchstop(t_watch1)
    if (icheck .eq. 1) write(*,10) " ### GETNEWMODEL [0]   END!! Time =",t_watch1%time," [min]"!2020.09.18
  
  !#[1]## combine JJ_ac and JJ_mt -> JJ     !            2022.01.05
        ! combine CD_ac and CD_mt -> CD     !            2022.01.05
    if ( ACT .and. ioutlevel == 1) then ! 2022.10.14
      open(1,file="JJ_ac.dat")
         call realcrsout(JJ_ac,1)
      close(1)
    end if ! 2022.10.14
    if (MT .and. ioutlevel == 1) then ! 2022.10.14
      open(1,file="JJ_mt.dat")
       call realcrsout(JJ_mt,1)
      close(1)
    end if ! 2022.10.14

    if (ijoint == 1 ) JJ = JJ_ac ! 2022.10.14
    if (ijoint == 2 ) JJ = JJ_mt ! 2022.10.14
    if (ijoint == 3 ) call  vertical_combine_real_crs(JJ_ac,JJ_mt,JJ) ! Joint       2022.10.14

    write(*,*) "combine JJ end"
    if (ioutlevel == 1 .and. .false.) then ! 2022.10.14
     open(1,file="JJ.dat")
      call realcrsout(JJ,1)
     close(1)
    end if

    ! As for combining CD, the following works for all ijoint = 1,2,3 2022.10.14
    CD%ntot  = CD_ac%ntot  + CD_mt%ntot 
    CD%ncolm = CD_ac%ncolm + CD_mt%ncolm 
    CD%nrow  = CD_ac%nrow  + CD_mt%nrow 
    ndat = ndat_ac + ndat_mt

    allocate(CD%stack(0:ndat),CD%item(ndat),CD%val(ndat))
    CD%stack(0)=0
    do i=1,ndat
      CD%stack(i)= i
      CD%item(i) = i
      if ( i .le. ndat_ac)      CD%val(i) = CD_ac%val(i)
      if ( i .gt. ndat_ac)      CD%val(i) = CD_mt%val(i-ndat_ac)
    end do
      write(*,*) "combine CD end"

  !#[2]## combine data vector
    ndat = ndat_ac + ndat_mt
    allocate( dobs(ndat),dcal(ndat))         !           2022.01.05
    allocate( X(ndat),JM(ndat) )             !           2022.01.05
    if ( ijoint == 1) then
     dobs = dobs_ac
     dcal = dcal_ac
    else if (ijoint == 2) then
     dobs = dobs_mt
     dcal = dcal_mt
    else if (ijoint == 3) then
    dobs(        1:ndat_ac) = dobs_ac        !           2022.01.05
    dobs(ndat_ac+1:ndat   ) = dobs_mt        !           2022.01.05
    dcal(        1:ndat_ac) = dcal_ac        !           2022.01.05
    dcal(ndat_ac+1:ndat   ) = dcal_mt        !           2022.01.05
    end if
    write(*,*) "combine data vector end"
    
  !#[2-1]## when transformation is used 2018.01.22
   call watchstart(t_watch1)
   if ( iboundflag .eq. 2 ) then   ! 2018.06.25
    ki    = h_model                ! 2018.06.25
    kiref = g_model_ref            ! 2018.06.25
    call TRANSMODELX(ki,   g_param_joint,1) ! model -> X
    call TRANSMODELX(kiref,g_param_joint,1) ! model -> X
    m0    = ki                     ! 2018.06.25
    m_ref = kiref                  ! 2018.06.25
   else
    m0    = h_model                ! 2018.06.25
    m_ref = g_model_ref            ! 2018.06.25
   end if
   call watchstop(t_watch1)
  if (icheck .eq. 1)  write(*,10) " ### GETNEWMODEL [1]   END!! Time =",t_watch1%time," [min]"!2020.09.18
  
  !#[1-2]## set model1 and model_ref
    ii=0
    do i=1,nmodel                           ! 2018.06.25
     if (iactive(i) .eq. 1) then            ! 2018.06.25
      ii=ii+1                               ! 2018.06.25
      model1(ii)    = m0%logrho_model(i)    ! 2018.06.25
      model_ref(ii) = m_ref%logrho_model(i) ! 2018.06.25
     end if                                 ! 2018.06.25
    end do                                  ! 2018.06.25
    !#
    if ( ii .ne. nmodelactive) then
     write(*,*) "GEGEGE ii",ii,"nmodelactive",nmodelactive
     stop
    end if
  
  !#[2]## assemble C
  !#[2-1]# gen JT
    call watchstart(t_watch1)
    call trans_crs2crs(JJ,JT) ! 2018.01.23 m_matrix.f90
    call watchstop(t_watch1)
  if (icheck .eq. 1) write(*,10) " ### GETNEWMODEL [2-1] END!! Time =",t_watch1%time," [min]"!2020.09.18
  
  !#[2-2]# cal CMJT = CM*JT
    call watchstart(t_watch1)
    call mulreal_crs_crs_crs(CM,JT,CMJT)
    call watchstop(t_watch1)
  if (icheck .eq. 1) write(*,10) " ### GETNEWMODEL [2-2] END!! Time =",t_watch1%time," [min]"!2020.09.18
  
  !#[2-3]# cal JCMJT = JJ * CMJT
    call watchstart(t_watch1)
    call mulreal_crs_crs_crs(JJ,CMJT,JCMJT)
    call watchstop(t_watch1)
  if (icheck .eq. 1) write(*,10) " ### GETNEWMODEL [2-3] END!! Time =",t_watch1%time," [min]"!2020.09.18

  open(1,file="cd.dat")
  call realcrsout(CD,1)
  close(1)
  open(1,file="jcmjt.dat")
  call realcrsout(JCMJT,1)
  close(1)

  !#[2-4]# cal C = JCMJT + alpha*CM
    call watchstart(t_watch1)
    ACD = CD
    ACD%val = Cd%val * alpha
    call add_crs_crs_crs2(JCMJT,ACD,C) ! C = JCMJT + alpha*CM 2022.01.05
    call watchstop(t_watch1)
  if (icheck .eq. 1) write(*,10) " ### GETNEWMODEL [2-4] END!! Time =",t_watch1%time," [min]"!2020.09.18

  open(1,file="c.dat")
  call realcrsout(c,1)
  close(1)
  
  !#[3]## gen X
  !#   X  = d_obs - d(M) + J(M - M_ref)
  !#[3-1]##
    call watchstart(t_watch1)
    dmodel = model1 - model_ref
    call mul_matcrs_rv(JJ,dmodel,nmodel,JM)
    X = dobs - dcal + JM           ! 2017.07.20
    call watchstop(t_watch1)
  if (icheck .eq. 1) write(*,10) " ### GETNEWMODEL [3-1] END!! Time =",t_watch1%time," [min]"!2020.09.18
  
  !#[4]## cal beta
  !# [C]beta = X
  !# M_ite+1 = M_ref + Cm*JT*beta
  !#[4-1]## solve beta
   call watchstart(t_watch1)
   allocate(beta(ndat))                     ! [ndat,1]
   call solvebeta(C,X,ndat,beta)            ! [ndat,ndat]*[ndat,1] = [ndat,1]
   call watchstop(t_watch1)
  if (icheck .eq. 1) write(*,10) " ### GETNEWMODEL [4-1] END!! Time =",t_watch1%time," [min]"!2020.09.18
  
  !#[4-2]## generate dmodel (crs format)
   call watchstart(t_watch1)
   call mul_matcrs_rv(CMJT,beta,ndat,dmodel)! [nmodel,ndat]*[ndat,1]=[nmodel,1]
   call watchstop(t_watch1)
  if (icheck .eq. 1)  write(*,10) " ### GETNEWMODEL [4-2] END!! Time =",t_watch1%time," [min]"!2020.09.18
  
  ii = 0 ! 2018.06.25
   do i=1,nmodel
    if ( iactive(i) .eq. 1) then ! 2018.06.25
     ii = ii + 1
     h_model%logrho_model(i) = dmodel(ii) + model_ref(ii)  ! [nmodel,1] 2018.06.25
     if (.false.) write(*,*) i,h_model%logrho_model(i),dmodel(ii),model_ref(ii)
    end if
   end do
  
  !#[5-1]## cut the model
    if ( iboundflag .eq. 1 ) then  ! 2018.01.18
     call BOUNDMODEL(h_model,g_param_joint)      ! 2018.01.18
    end if
  
  !#[5-2]## transform X to model
    if ( iboundflag .eq. 2 ) then                !  2018.01.22
     call watchstart(t_watch1)
     !# iflag=1: model -> x, iflag = 2: X to model  2018.01.22
     call TRANSMODELX(h_model,g_param_joint,2)   !  2018.01.22 X -> Model
     call watchstop(t_watch1)
  if (icheck .eq. 1) write(*,10) " ### GETNEWMODEL [5-2] END!! Time =",t_watch1%time," [min]"!2020.09.18
    end if                                       !  2018.01.22
  
   call watchstop(t_watch) ! 2017.12.22
   write(*,10) " ### GETNEWMODEL_JOINT END!! ###   Time =",t_watch%time," [min]"!2020.09.18
  
  return
  
  10 format(a,f8.4,a) ! 2020.09.18
  end
!#################################################### BOUNDMODEL
!# coded on 2018.01.18
subroutine BOUNDMODEL(h_model,g_param_joint)
use modelpart
use param_jointinv
implicit none
type(param_joint),intent(in)    :: g_param_joint
type(model),            intent(inout) :: h_model
integer(4) :: i
real(8)    :: logrho,logrho_upper,logrho_lower
integer(4) :: icount_low,icount_up

!#[1]## set
logrho_upper = g_param_joint%logrho_upper
logrho_lower = g_param_joint%logrho_lower
icount_low = 0
icount_up  = 0

!#[2]## cut off
do i=1,h_model%nmodel
 logrho=h_model%logrho_model(i)
 if ( logrho .gt. logrho_upper ) then
  logrho = logrho_upper
  icount_up = icount_up + 1
 else if ( logrho .lt. logrho_lower ) then
  logrho = logrho_lower
  icount_low = icount_low + 1
 else
  goto 100
 end if
  h_model%logrho_model(i) = logrho
 100 continue
end do

write(*,*) "icount_low=",icount_low," icount_up=",icount_up
write(*,*) "### COUNDMODEL END!! ###"

return
end

end module
