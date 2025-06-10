! coded on 2016.09.03
module constants
implicit none

! earthrad, pi, d2r, r2d
real(8),parameter :: earthrad=6371.2d0, pi=4.d0*datan(1.d0),d2r=pi/180., r2d=180./pi
real(8),parameter :: planetrad=6371.2d0
real(8),parameter :: epsilon=8.85418782*1.d-12
! mu
real(8),parameter :: dmu=4.d0*pi*1.d-7

end module constants
