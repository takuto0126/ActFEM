! Coded by T. MINAMI based on
! static_linear/set_iccg_var.f90 on 2015.08.15
! This is module file for n_ebfem.f90
! aiccg_u, aiccg_l, diag, x_vec, b_vec are complex
!# Modified on Dec. 22, 2015
!# only for node FEM
module iccg_var_takuto
implicit none
type global_matrix
!# matrix for induction equation
 integer(4),dimension(:), pointer :: IAU, IAL ! coumn id of entries
 integer(4),dimension(:), pointer :: INU, INL ! accumulated # of row entries
 complex(8),dimension(:), pointer :: AU, AL
 complex(8),dimension(:), pointer :: D
 integer(4) :: doftot,intdoftot
 integer(4) :: iau_tot  ! size of au, iau
 integer(4) :: ial_tot  ! size of al, ial
end type

contains
!######################################## INITIALIZE 2021.09.15
!moved from forward_3DMT.f90 on 2021.09.15
subroutine initialize(A,b_vec,bs,doftot,nsr) ! 2017.07.11
!use iccg_var_takuto ! global_matrix
implicit none
type(global_matrix)             :: A
integer(4),         intent(in)  :: doftot, nsr
complex(8),         intent(out) :: b_vec(doftot,nsr),bs(doftot,nsr) ! 2017.07.11

A%D(:)     = (0.d0,0.d0)
A%AU(:)    = (0.d0,0.d0)
A%AL(:)    = (0.d0,0.d0)
bs(:,:)    = (0.d0,0.d0)
b_vec(:,:) = (0.d0,0.d0)

return
end

!############################################# sup_iccg 2021.09.16
! moved from forward_3DTM.f90 on 2021.09.16
! Coppied from static_linear/sup_iccg.f90 on June 10, 2015
! here elm_k, duag, aiccg_u are complex
subroutine  sup_iccg ( &
&      elm_k, table_dof_elm, edof,                                       &
&      diag , istack_u   , item_u , aiccg_u, doftot, item_u_tot     )
!    superpose symmetric element stiffness matrix to gloval stiffness
!    for iccg solver (lower element, upper whole stiffness)
implicit  none
integer(4),intent(in) :: edof, doftot, item_u_tot
complex   (8),   intent(in)   ::  elm_k(edof,edof) !  element stiffness matrix
integer(4), intent(in)   :: table_dof_elm(edof)
!    digree of freedom number table at a element table_dof_elm(i)
!      where i    : dof order of 'elm_k'
complex(8), intent(inout) ::  diag(doftot)                  ! diagonal value   (i-th dof)
integer(4), intent(in) ::  istack_u(0:item_u_tot) !  last term count at each freedom  (i-th dof)
integer(4), intent(in) ::  item_u(item_u_tot)       ! freedom number at each term      (i-th term)
complex(8), intent(inout) ::  aiccg_u(item_u_tot)      ! upper triangular k value         (i-th term)
integer(4)  ::  i_gl, j_gl, i_elm, j_elm
integer(4)  ::  num , i_pos
logical        ::  found
! * superpose element stiffness matrix
!      store diagonal and upper term for packed whole stiffness matrix
!       ( element stiffness matrix stored only lower triangle )
!
!    k(i_gl,j_gl) = k(i_gl,j_gl) + elm_k(i_elm,j_elm)
!
!   . make diagonal and upper triangle
!
      do i_elm=1,edof
         i_gl = table_dof_elm(i_elm)
         do j_elm=1,edof
            j_gl = table_dof_elm(j_elm)
            if      ( i_gl  == j_gl  ) then            ! diagonal
              diag(i_gl) = diag(i_gl) + elm_k(i_elm,j_elm)
            else if ( i_gl  <  j_gl  ) then            ! upper triangle
              found = .false.
!                                                 .. search term pos.
              do num=istack_u(i_gl-1)+1,istack_u(i_gl)
                if ( j_gl == item_u(num) ) then
                  i_pos =  num
                  found = .true.
                  exit
                endif
              enddo
              if ( .not.found ) then
                write(*,*) ' ***** error on sup_iccg : dof not found  ', i_gl, j_gl
!                call geofem_abort(34,'internal logic error')
                stop
              endif
!                                                 .. superpose the term
              if ( i_elm >  j_elm )  then              ! elm. lower
                aiccg_u(i_pos) = aiccg_u(i_pos) + elm_k(i_elm,j_elm)
              else                                     ! elm. upper(sym)
                aiccg_u(i_pos) = aiccg_u(i_pos) + elm_k(j_elm,i_elm)
              endif
            endif
         enddo
      enddo
      return
      end subroutine  sup_iccg

!############################################### COPY_UL_ICCG12 2021.09.15
! copied from forward_3DMT.f90 on 2021.09.16
! Copied from static_linear/copy_ul_iccg.f90 on June 10, 2015
! Converted aiccg_u, aiccg_l, diag to complex on Aug. 21, 2015
!      subroutine  copy_ul_iccg12( doftot,item_u_tot,item_l_tot,istack_u,item_u,aiccg_u,istack_l,item_l,aiccg_l,ip)
     subroutine  copy_ul_iccg12( A,ip)
!    copy upper triangle term value to lower for iccg solver
      implicit  none
    type(global_matrix),intent(inout) :: A
    integer(4),intent(in) :: ip ! 2021.10.04
      integer(4 ) ::  i_gl, j_gl
      integer(4)  ::  num , num2, i_pos
      logical     ::  found
    logical,allocatable,dimension(:) :: assigned
    integer(4),allocatable,dimension(:) :: order
    complex(8) :: test0=0.d0, test1=0.d0
    A%AL(:)=(0.d0,0.d0)
    allocate (assigned(A%ial_tot),order(A%ial_tot))
    assigned(:)=.false.
!  * copy upper triangle to lower triangle  (i_gl,j_gl) -> (j_gl,i_gl)
!
      do i_gl=1,A%doftot
         do num=A%INU(i_gl-1)+1,A%INU(i_gl)
              j_gl=A%IAU(num)
              found = .false.
!                                                 .. search term pos.
!              do num2=istack_l(j_gl-1)+1,istack_l(j_gl)
               do num2=A%INL(j_gl-1)+1,A%INL(j_gl)
!                if ( i_gl == item_l(num2) ) then
                if ( i_gl == A%IAL(num2) ) then
                  i_pos =  num2
                  found = .true.
                  exit
                endif
              enddo
              if ( .not.found ) then
                write(*,*) ' ***** error on copy_ul_iccg: dof not found', i_gl, j_gl
!                call geofem_abort(34,'internal logic error')
                stop
              endif
!                                                 .. set the term
!              aiccg_l(i_pos) = aiccg_u(num)
               A%AL(i_pos) = A%AU(num)
           assigned(i_pos)=.true.
           order(num)=i_pos
!           test0=test0+A%AU(num)*dconjg(A%AU(num))
!           test1=test1+A%AL(i_pos)*dconjg(A%AL(i_pos))
               if (cdabs(A%AL(i_pos)) .ne. cdabs(A%AU(num)) ) then
                write(*,*) "GEGEGE"
            stop
           end if
        enddo
      enddo
      do i_gl=1,A%ial_tot
     if ( .not. assigned(i_gl) ) then
      write(*,*) "GEGEGE not assigned i_gl=",i_gl
      stop
     end if
!     test0=test0+A%AU(i_gl)
!     test1=test1+A%AL(i_gl)
      end do
!    write(*,*) "test0=",test0
!    write(*,*) "test1=",test1
!      if (ip .eq. 0 ) write(*,*) "### COPY_UL_ICCG12 END!! ### ip=",ip ! commented out 2022.01.02
      return
      end

!####################################################################################
! coded on 2021.09.16 by Takuto Minami
! simplify set_iccg_var7_dofn assuming dofn=1, and deleting table_dof
!
! In case of edge-base, interpret as
! nodtot    -> edgetot
! intnodtot -> intedgetot
! enod      -> # of edges per element
! dofn      -> order of edge element, 1 for 1st order, where 1 dof is assigned per edge
!
!    subroutine  set_iccg_var7_dofn(dofn,enod,nodtot,intnodtot,doftot,elmtot,elmid,table_dof,A,ip)
    subroutine  set_size_of_A(           enod,nodtot,intnodtot,doftot,elmtot,elmid,          A,ip)
    implicit  none
    type(global_matrix) :: A
!    integer(4),intent(in) :: ip, dofn, enod, table_dof(nodtot,dofn),doftot
    integer(4),intent(in) :: ip, enod, doftot
    integer(4),intent(in) :: nodtot, intnodtot
    integer(4),intent(in) :: elmtot, elmid(elmtot,enod) !  total element count on this PE
    integer(4), parameter :: max_retry    = 10                 !   max. retry count
    integer(4), parameter :: dofn = 1               ! 2021.09.16
    integer(4)            :: table_dof(nodtot,dofn) ! 2021.09.16
    integer(4), dimension(max_retry)         ::  max_inum_table  !  max_inum retry table
    integer(4), dimension(:,:), allocatable  ::  index_u_tmp ! freedom number at each term (upper, temporaty)
    integer(4), dimension(:,:), allocatable  ::  index_l_tmp ! freedom number at each term (lower, temporaty)
    integer(4), dimension(:),   allocatable  ::  inum_u_tmp ! term count at each freedom  (i-th dof)
    integer(4), dimension(:),   allocatable  ::  inum_l_tmp  ! term count at each freedom  (i-th dof)
    integer(4)     ::  i_dof, num, max_inum,  retry_count   ! retry count
    logical        ::  inum_overflow                                         !  table overflow status key
    data max_inum_table / 98,120,160,200,300,400,500,600,800,1000/
!     data max_inum_table / 2000, 3000, 4000, 8000, 5*10000/
    do i_dof=1,doftot         ! 2021.09.16
     table_dof(i_dof,1)=i_dof ! 2021.09.16
    end do                    ! 2021.09.16
    A%doftot=doftot           ! for example, doftot is nline 2021.10.14
    A%intdoftot=intnodtot*dofn
    allocate (  inum_u_tmp ( doftot),  inum_l_tmp ( doftot ) ) ! how many terms in i-th dof row
!      write(*,*) "check0"

! 1. allocate temporary array 'index_u_tmp', 'index_l_tmp','inum_u_tmp' , 'inum_l_tmp'
    do retry_count=1,max_retry
        inum_overflow = .false.
        max_inum = max_inum_table(retry_count)  !  max. term count at a d.o.f.
        allocate ( index_u_tmp ( doftot, max_inum), index_l_tmp ( doftot, max_inum) ) ! index of dof colmun

        call C_set_dof2(dofn,enod,elmtot,elmid,table_dof,nodtot,doftot,max_inum, inum_overflow, &
      &                        inum_u_tmp, index_u_tmp, inum_l_tmp, index_l_tmp)

        if ( .not. inum_overflow ) exit ! exit the loop
      write(*,*) 'set_iccg_val : retry # ', retry_count, max_inum
        deallocate( index_u_tmp, index_l_tmp )   !   .. next try
      enddo
      if (inum_overflow)  goto 99
!      write(*,*) "check1"

! 3. allocate and set 'istack_u','istack_l' (last term count at each line)
      allocate ( A%INU ( 0:doftot ),  A%INL ( 0:doftot ) )
      A%INU(0) = 0  ;  A%INL(0) = 0
      do i_dof=1,doftot
        A%INU(i_dof) = A%INU(i_dof-1) + inum_u_tmp(i_dof)
        A%INL(i_dof) = A%INL(i_dof-1) + inum_l_tmp(i_dof)
      enddo
!      write(*,*) "check3"

! 4. set array length  'item_u_tot', 'item_l_tot'
      A%iau_tot = A%INU(doftot)
      A%ial_tot = A%INL(doftot)
!    write(*,*) "check4"

! 5. allocate and set 'item_u','item_l'
!       (freedom number at each term)
      allocate ( A%IAU ( A%iau_tot ),  A%IAL  ( A%ial_tot ) )
    write(*,'(3(a,i10,2x))') "iau_tot=",A%iau_tot,"ial_tot=",A%ial_tot,"ip=",ip
      A%IAU = 0
      A%IAL = 0
!    set item_u and item_l
      do i_dof=1,doftot
        do num=1,inum_u_tmp(i_dof)
          A%IAU(A%INU(i_dof-1)+num) = index_u_tmp(i_dof,num)
        enddo
        do num=1,inum_l_tmp(i_dof)
          A%IAL(A%INL(i_dof-1)+num) = index_l_tmp(i_dof,num)
        enddo
      enddo
!    write(*,*) "check5"

! 6. deallocate temporary array 'index_u_tmp', 'index_l_tmp',  'inum_u_tmp', 'inum_l_tmp'
      deallocate ( index_u_tmp )
      deallocate ( index_l_tmp )
      deallocate ( inum_u_tmp )
      deallocate ( inum_l_tmp )
! 7. allocate other iccg variables
      allocate (     A%D (   doftot          ) )
      allocate (     A%AU (   A%iau_tot         ) )
      allocate (     A%AL (   A%ial_tot         ) )

     write(*,'(a,i10)') "### SET_SIZE_OF_A END!! ### ip=",ip
      return
99       write(*,*) 'set_iccg_val : max_inum exceed in retry ',max_retry, max_inum
        stop
      end subroutine set_size_of_A ! 2021.09.16

!####################################################################################
! In case of edge-base, interpret as
! nodtot    -> edgetot
! intnodtot -> intedgetot
! enod      -> # of edges per element
! dofn      -> order of edge element, 1 for 1st order, where 1 dof is assigned per edge
!
    subroutine  set_iccg_var7_dofn(dofn,enod,nodtot,intnodtot,doftot,elmtot,elmid,table_dof,A,ip)
    implicit  none
	type(global_matrix) :: A
	integer(4),intent(in) :: ip, dofn, enod, table_dof(nodtot,dofn),doftot
	integer(4),intent(in) :: nodtot, intnodtot
	integer(4),intent(in) :: elmtot, elmid(elmtot,enod) !  total element count on this PE
    integer(4), parameter   ::  max_retry    = 10                 !   max. retry count
    integer(4), dimension(max_retry) ::  max_inum_table  !  max_inum retry table
    integer(4), dimension(:,:), allocatable  ::  index_u_tmp ! freedom number at each term (upper, temporaty)
    integer(4), dimension(:,:), allocatable  ::  index_l_tmp ! freedom number at each term (lower, temporaty)
    integer(4), dimension(:),   allocatable  ::  inum_u_tmp ! term count at each freedom  (i-th dof)
    integer(4), dimension(:),   allocatable  ::  inum_l_tmp  ! term count at each freedom  (i-th dof)
    integer(4)  ::  i_dof, num, max_inum,  retry_count   ! retry count
    logical        ::  inum_overflow                                         !  table overflow status key
    data max_inum_table / 98,120,160,200,300,400,500,600,800,1000/
!     data max_inum_table / 2000, 3000, 4000, 8000, 5*10000/
	A%doftot=doftot
	A%intdoftot=intnodtot*dofn
	allocate (  inum_u_tmp ( doftot),  inum_l_tmp ( doftot ) ) ! how many terms in i-th dof row
!      write(*,*) "check0"

! 1. allocate temporary array 'index_u_tmp', 'index_l_tmp','inum_u_tmp' , 'inum_l_tmp'
	do retry_count=1,max_retry
        inum_overflow = .false.
        max_inum = max_inum_table(retry_count)  !  max. term count at a d.o.f.
        allocate ( index_u_tmp ( doftot, max_inum), index_l_tmp ( doftot, max_inum) ) ! index of dof colmun

        call C_set_dof2(dofn,enod,elmtot,elmid,table_dof,nodtot,doftot,max_inum, inum_overflow, &
	  &                        inum_u_tmp, index_u_tmp, inum_l_tmp, index_l_tmp)

        if ( .not. inum_overflow ) exit ! exit the loop
	  write(*,*) 'set_iccg_val : retry # ', retry_count, max_inum
        deallocate( index_u_tmp, index_l_tmp )   !   .. next try
      enddo
      if (inum_overflow)  goto 99
!      write(*,*) "check1"

! 3. allocate and set 'istack_u','istack_l' (last term count at each line)
      allocate ( A%INU ( 0:doftot ),  A%INL ( 0:doftot ) )
      A%INU(0) = 0  ;  A%INL(0) = 0
      do i_dof=1,doftot
        A%INU(i_dof) = A%INU(i_dof-1) + inum_u_tmp(i_dof)
        A%INL(i_dof) = A%INL(i_dof-1) + inum_l_tmp(i_dof)
      enddo
!      write(*,*) "check3"

! 4. set array length  'item_u_tot', 'item_l_tot'
      A%iau_tot = A%INU(doftot)
      A%ial_tot = A%INL(doftot)
!	write(*,*) "check4"

! 5. allocate and set 'item_u','item_l'
!       (freedom number at each term)
      allocate ( A%IAU ( A%iau_tot ),  A%IAL  ( A%ial_tot ) )
!	write(*,'(3(a,i10,2x))') "iau_tot=",A%iau_tot,"ial_tot=",A%ial_tot,"ip=",ip
      A%IAU = 0
      A%IAL = 0
!    set item_u and item_l
      do i_dof=1,doftot
        do num=1,inum_u_tmp(i_dof)
          A%IAU(A%INU(i_dof-1)+num) = index_u_tmp(i_dof,num)
        enddo
        do num=1,inum_l_tmp(i_dof)
          A%IAL(A%INL(i_dof-1)+num) = index_l_tmp(i_dof,num)
        enddo
      enddo
!	write(*,*) "check5"

! 6. deallocate temporary array 'index_u_tmp', 'index_l_tmp',  'inum_u_tmp', 'inum_l_tmp'
      deallocate ( index_u_tmp )
      deallocate ( index_l_tmp )
      deallocate ( inum_u_tmp )
      deallocate ( inum_l_tmp )
! 7. allocate other iccg variables
      allocate (     A%D (   doftot          ) )
      allocate (     A%AU (   A%iau_tot         ) )
      allocate (     A%AL (   A%ial_tot         ) )

!	write(*,'(a,i10)') "### SET_ICCG_VAR7_node END!! ### ip=",ip
      return
99 	  write(*,*) 'set_iccg_val : max_inum exceed in retry ',max_retry, max_inum
        stop
	  end subroutine set_iccg_var7_dofn
	  !##############################################
      subroutine  C_set_dof2( dofn,enod,elmtot,elmid,table_dof,nodtot,doftot,max_inum, inum_overflow, &
	&  inum_u_tmp, index_u_tmp, inum_l_tmp, index_l_tmp)
      implicit none
	integer(4),intent(in)  :: dofn,enod,elmtot, elmid(elmtot,enod)
	integer(4),intent(in)  :: nodtot,doftot, max_inum, table_dof(nodtot,dofn)
	integer(4),intent(out) :: inum_u_tmp(doftot),  inum_l_tmp(doftot) ! # of items in the row in concern
	integer(4),intent(out) ::  index_u_tmp(doftot, max_inum) ! band matrix of column info
	integer(4),intent(out) ::  index_l_tmp(doftot,  max_inum) ! band matrix of column info
	logical,intent(inout)  :: inum_overflow
!
      integer(4) :: elm_dof(enod*dofn)
	integer(4) :: kf_e, nod, nod_e, num_dof
	integer(4) :: elm, kf, lf, i_dof, j_dof, num, i_pos
	logical :: exist
       inum_u_tmp = 0
      index_u_tmp = 0

       inum_l_tmp = 0
      index_l_tmp = 0

      do elm=1,elmtot
!                             .. set element parameter table 'lem_char'
!        call  set_echa ( elmtyp(elm) )
!                             .. set element freedom table 'elm_dof'
        kf_e     = 0
        elm_dof  = 0
        do nod_e=1,enod
          nod = abs(elmid(elm,nod_e)) ! abs for n6line
          do kf=1,dofn
            num_dof = table_dof(nod,kf)
            if ( num_dof /= 0 ) then
              kf_e = kf_e + 1
              elm_dof(kf_e) = num_dof
            endif
          enddo
        enddo
!                             .. check element freedom table
        do kf=1,kf_e
          i_dof  = elm_dof(kf)
          do lf=1,kf_e
            j_dof  = elm_dof(lf)
!                                  .. upper triangle (i_dof=j_dof : diag)
            if (i_dof<j_dof) then
!                         k(i_dof,j_dof).. search term already exist (upper)
              exist = .false.
              do num=1,inum_u_tmp(i_dof)
                if      ( j_dof == index_u_tmp(i_dof,num) ) then
                  exist = .true.
                  exit
                else if ( j_dof <  index_u_tmp(i_dof,num) ) then
                  exit
                endif
              enddo
!                                       .. add term when non-existent
              if ( .not.exist ) then
!                                                 .. add 'inum_u_tmp'
                inum_u_tmp(i_dof) = inum_u_tmp(i_dof) + 1

                if ( inum_u_tmp(i_dof) > max_inum     ) then
                  inum_overflow = .true.
                  return
                endif
!                                                 .. set 'index_u_tmp'
!                                         . search insert position
                i_pos =  inum_u_tmp(i_dof)
                do num=1,inum_u_tmp(i_dof) - 1
                  if    ( j_dof <  index_u_tmp(i_dof,num) ) then
                    i_pos = num
                    exit
                  endif
                enddo
!                                         . shift index
                do num = inum_u_tmp(i_dof),i_pos+1,-1
                  index_u_tmp(i_dof,num  ) = index_u_tmp(i_dof,num-1)
                enddo
                  index_u_tmp(i_dof,i_pos) = j_dof

              endif
!                                  .. lower triangle
            else if (i_dof>j_dof) then
!                         k(i_dof,j_dof).. search term already exist (lower)
              exist = .false.
              do num=1,inum_l_tmp(i_dof)
                if      ( j_dof == index_l_tmp(i_dof,num) ) then
                  exist = .true.
                  exit
                else if ( j_dof <  index_l_tmp(i_dof,num) ) then
                  exit
                endif
              enddo
!                                       .. add term when non-existent
              if ( .not.exist ) then
!                                                 .. add 'inum_l_tmp'
                inum_l_tmp(i_dof) = inum_l_tmp(i_dof) + 1

                if ( inum_l_tmp(i_dof) > max_inum     ) then
                  inum_overflow = .true.
                  return
                endif
!                                                 .. set 'index_l_tmp'
!                                         . search insert position
                i_pos =  inum_l_tmp(i_dof)
                do num=1,inum_l_tmp(i_dof) - 1
                  if    ( j_dof <  index_l_tmp(i_dof,num) ) then
                    i_pos = num
                    exit
                  endif
                enddo
!                                         . shift index
                do num = inum_l_tmp(i_dof),i_pos+1,-1
                  index_l_tmp(i_dof,num  ) = index_l_tmp(i_dof,num-1)
                enddo
                  index_l_tmp(i_dof,i_pos) = j_dof

              endif

            endif

          enddo
        enddo

      enddo

      return
	end subroutine  C_set_dof2

end module iccg_var_takuto




