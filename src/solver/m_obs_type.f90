
! coded on May 20, 2016
!
module obs_type
use matrix
implicit none

type obs_info
 integer(4)                          :: nobs
 character(4)                        :: name
 real(8),allocatable,dimension(:,:)  :: xyz_obs
 type (real_crs_matrix)              :: coeff(2,3)
 integer(4),allocatable,dimension(:) :: devnum     ! 2017.10.12
end type


end module obs_type
