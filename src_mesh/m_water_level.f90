!2023.08.31
module water_level
implicit none

type param_water_level
 integer(4) :: inum_water_level
 real(8),allocatable,dimension(:) :: water_level_elev
 real(8)    :: lonlat_lake_center(2)
 real(8)    :: xy_lake_center(2)
 !#
!# water level info 2023.09.05
 ! Case for inum_water_lake = 3
 !iflag_nodez_2d : 0:normal topo, i: i-th lake surface node
 ! surface 1 (flag = 0)
 !------------- ___5___ ------------ 
 !             \ __6__ / < surface 2  <- flag = 0
 !              \ _7_ /  < surface 3  <- flag = 0
 !               \ _ /   < surface 4  <- flag = 0
 ! surface 7 is not included in nakadake2d.msh 
 ! surface 4 is copied to generate surface 7 first
 ! surface 5 (flag = 1)
 ! surface 6 (flag = 2)
 ! surface 7 (flag = 3) <- water_level_elev will be assigned               
 integer(4) :: node2d
 integer(4),allocatable,dimension(:) :: iflag_nodez_2d ! (node2d)

end type

contains
!#################################################################
subroutine readwaterlevelctl(g_param,g_param_water,iunit)
use param
implicit none
type(param_forward),    intent(in)  :: g_param
type(param_water_level),intent(out) :: g_param_water
integer(4),intent(in) :: iunit
integer(4) :: inum_water_level,i
real(8)    :: lonorigin,latorigin,xy(2),lonlat(2)

read(iunit,*) ! header
read(iunit,'(20x,i10)') inum_water_level
allocate(g_param_water%water_level_elev(inum_water_level))
 read(iunit,'(20x,g15.7)') g_param_water%lonlat_lake_center(1) ! 2023.09.01
 read(iunit,'(20x,g15.7)') g_param_water%lonlat_lake_center(2) ! 2023.09.01
do i=1,inum_water_level
 read(iunit,'(20x,g15.7)') g_param_water%water_level_elev(i)
end do

!# get xy_lake_center 2023.09.01
lonlat(1:2) = g_param_water%lonlat_lake_center(1:2) !! allocate and fill
lonorigin = g_param%lonorigin
latorigin = g_param%latorigin

 call UTMXY(lonlat,lonorigin,latorigin,xy,g_param%UTM)
 g_param_water%xy_lake_center(1:2)=xy(1:2)
 g_param_water%inum_water_level = inum_water_level ! 2023.09.04

open(1,file="xy_center_lake.dat")
 write(1,*) xy(1:2)
close(1)

return
end subroutine readwaterlevelctl

end module water_level