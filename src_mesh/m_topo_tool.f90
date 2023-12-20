module topo_tool
implicit none

type grid_data
 integer(4) :: node
 integer(4) :: neast,nsouth
 real(8),dimension(:,:),allocatable :: lonlatalt ! lon [deg], lat [deg], altitude [m]
 ! eastward[km], northward[km], upward[km]
  real(8),dimension(:,:),allocatable :: xyz
  logical :: iflag_topright_land     ! 2019.02.25
  real(8),allocatable,dimension(:) :: lon,lat ! 2023.09.01
end type grid_data

type coast_data
 integer(4) :: ncmax  ! max # of points on coastline
 integer(4) :: ncoast ! # of points on coastline
 integer(4) :: nbound ! # of boundary points
! cxy(1,i)  :  eastward, cxy(2,i) : northward
! ind(i,1)  :  grid node id of original gebco file on left or up side
! ind(i,2)  :  1 -> ind(1) indicates left node, 2 -> ind(2) is top side node
 real(8),   allocatable,dimension(:,:) :: cxy
 integer(4),allocatable,dimension(:,:) :: ind ! ind(ncmax,2)
 integer(4),allocatable,dimension(:,:) :: ind_r ! ind_r(ncmax,2) reordered ind
 integer(4),allocatable,dimension(:,:) :: label ! label(ncmax,2)
end type coast_data
! ind and ind_r are different in the order of node_num.
! former is in the same manner as topofile, but latter is in polygon manner

type poly_data
integer(4) :: lpmax  ! max # of polygons
integer(4) :: ncmax  ! max # of points for each polygon
integer(4) :: lpoly0 ! # of polygons
integer(4) :: nclose ! # of un-closed polygons, =nbound/2
integer(4),allocatable,dimension(:,:,:) :: ind2 ! ind2(nclose,2,2)
! ind2(i,1,:) : ind for start node of polygon
! ind2(i,2,:) : ind for end node of polygon
! ind2 is defined only for unclosed polygons because the start and end node is always
! assumed on the calculation boundary
integer(4),allocatable,dimension(:)     :: lpoly  ! # of points in each polygon
real(8),   allocatable,dimension(:,:,:) :: xypoly
! xypoly(1:2,1:lpoly0,1:ncmax) x,y component of
logical   :: iflag_topright_land       ! 2019.02.21
real(8),allocatable,dimension(:,:)      :: loc ! 2019.02.21 see addcorner.f90
end type

type bound_data
 integer(4) :: ncmax
 integer(4),allocatable,dimension(:) :: zlabel
 integer(4),allocatable,dimension(:) :: ibelong
end type bound_data

contains
!######################################## ALLOCATEBOUND
subroutine allocatebound(g_bound, ncmax)
type(bound_data),intent(out) :: g_bound
integer(4),      intent(in)  :: ncmax

 g_bound%ncmax=ncmax
 allocate(g_bound%zlabel(ncmax))
 allocate(g_bound%ibelong(ncmax))

return
end subroutine allocatebound

!######################################## ALLOCATEPOLY
subroutine allocatepoly(g_poly,lpmax,ncmax,nclose)
implicit none
type(poly_data),intent(out) :: g_poly
integer(4),     intent(in)  :: ncmax, lpmax, nclose

g_poly%lpmax=lpmax
g_poly%ncmax=ncmax
g_poly%lpoly0=0

allocate(g_poly%lpoly(lpmax))
allocate(g_poly%xypoly(2,lpmax,ncmax))
allocate(g_poly%ind2(nclose,2,2))
allocate(g_poly%loc(nclose,2)) ! 2019.02.25

write(*,*) "### ALLOCATEPOLY END!! ###"
return
end
!######################################## ALLOCATECOAST
subroutine allocatecoast(g_coast,ncmax)
implicit none
type(coast_data),intent(out) :: g_coast
integer(4),      intent(in)  :: ncmax

 g_coast%ncmax = ncmax
 allocate( g_coast%cxy(2,ncmax)  )
 allocate( g_coast%ind(ncmax,2)  )
 allocate( g_coast%ind_r(ncmax,2))
 allocate( g_coast%label(ncmax,2))

write(*,*) "### ALLOCATECOAST END!! 3##"
return
end subroutine allocatecoast

!######################################## DEALLOCATECOAST
subroutine deallocatecoast(g_coast)
implicit none
type(coast_data),intent(inout) :: g_coast

if (allocated(g_coast%cxy))   deallocate( g_coast%cxy)
if (allocated(g_coast%ind))   deallocate( g_coast%ind)
if (allocated(g_coast%ind_r)) deallocate( g_coast%ind_r)
if (allocated(g_coast%label)) deallocate( g_coast%label)
 g_coast%ncmax=0  ! max # of points on coastline
 g_coast%ncoast=0 ! # of points on coastline
 g_coast%nbound=0 ! # of boundary points

write(*,*) "### DEALLOCATE COAST END!! ###"

return
end subroutine deallocatecoast

!######################################## DEALLOCATEPOLY
subroutine deallocatepoly(g_poly)
implicit none
type(poly_data),intent(inout) :: g_poly

 if (allocated(g_poly%lpoly))  deallocate(g_poly%lpoly)
 if (allocated(g_poly%xypoly)) deallocate(g_poly%xypoly)
 if (allocated(g_poly%ind2))   deallocate(g_poly%ind2)
 if (allocated(g_poly%loc) )   deallocate(g_poly%loc)
 g_poly%lpmax=0  ! max # of polygons
 g_poly%ncmax=0  ! max # of points for each polygon
 g_poly%lpoly0=0 ! # of polygons
 g_poly%nclose=0 ! # of un-closed polygons, =nbound/2
 g_poly%iflag_topright_land=.false.

write(*,*) "### DEALLOCATE POLY END!! ###"

return
end subroutine deallocatepoly

!#######################################################
! modified on 20200728
! Coded on Sep. 17, 2015
subroutine changeorder(x1,y1,z1,ntopo,n_south, n_east)
implicit none
integer(4),  intent(in)    :: ntopo
real(8),     intent(inout) :: x1(ntopo),y1(ntopo),z1(ntopo)
integer(4),  intent(out)   :: n_south, n_east
real(8),allocatable,dimension(:) :: x2,y2,z2
integer(4) :: i, j, ishift, jshift
allocate(x2(ntopo),y2(ntopo),z2(ntopo))
!write(*,*) "changeorder start!" ! commented out on 2021.09.29

if ( x1(1) .eq. x1(2) ) then ! [if #1 start] if the order is based on colmun
  if ( y1(1) .lt. y1(2) ) then ! [if #2 start] ! upward (south to north)
   i=1
   do while ( y1(i) .lt. y1(i+1) )
     i=i+1
   end do
   n_south=i ; n_east=ntopo/n_south ! ok
   do i=1,n_south
      ishift=(i-1)*n_east
      do j=1, n_east
	  jshift=(j-1)*n_south
        x2(ishift+j)=x1(jshift + n_south - (i-1) )
	  y2(ishift+j)=y1(jshift + n_south - (i-1) )
	  z2(ishift+j)=z1(jshift + n_south - (i-1) )
	end do
    end do
!
  else if ( y1(2) .lt. y1(1) ) then ! downward (north to south)
   i=1
   do while ( y1(i) .gt. y1(i+1) )
     i=i+1
   end do
   n_south=i ; n_east=ntopo/n_south
   do i=1,n_south
      ishift=(i-1)*n_east
      do j=1, n_east
	  jshift=(j-1)*n_south
      x2(ishift+j)=x1(jshift + i )
	  y2(ishift+j)=y1(jshift + i )
	  z2(ishift+j)=z1(jshift + i )
	end do
   end do
!
  else
    write(*,*) "GEGEGE! x1(1) = x1(2) and y1(1)=y2(2)"
    stop
  end if ! [if #2 end] upward or downward end when column based

  x1=x2
  y1=y2
  z1=z2
  i=1
  do while ( x1(i) .lt. x1(i+1) )
   i=i+1
  end do
  n_east=i ; n_south=ntopo/n_east

elseif ( y1(1) .eq. y1(2) ) then ![if #1 ]  if the order is based on row 20200728
!# check upward or downward
  i=1
  do while ( y1(i) .eq. y1(i+1) )
    i=i+1
  end do
  n_east=i ; n_south = ntopo/n_east

  if ( y1(n_east) .lt. y1(n_east+1) ) then ![if #3] upward when row based
  !# change to downward
  do i=1,n_south
   x2(n_east *(n_south - i) + 1 : n_east *(n_south - i) + n_east) = x1((i-1)*n_east + 1 : i*n_east )
   y2(n_east *(n_south - i) + 1 : n_east *(n_south - i) + n_east) = y1((i-1)*n_east + 1 : i*n_east )
   z2(n_east *(n_south - i) + 1 : n_east *(n_south - i) + n_east) = z1((i-1)*n_east + 1 : i*n_east )
  end do
  x1=x2
  y1=y2
  z1=z2

  else !![if #3] downward when row based
   ! #nothing to do
  end if !![if #3 end] if upward and downward when the order is row based

end if   ! [if #1 end] if the order is based on column or row

return

end subroutine changeorder

end module topo_tool