program test_watch
  use caltime
  implicit none

  type(watch) :: mywatch
  integer :: isleep_sec

  isleep_sec = 5  ! 5秒停止

  call watchstart(mywatch)
  print *, "Start timing..."

  call sleep(isleep_sec)  ! 5秒間待機

  call watchstop(mywatch)
  print *, "Elapsed time (min) = ", mywatch%time*60

  call watchstart(mywatch)
  print *, "Start timing...2"

  call sleep(isleep_sec)  ! 5秒間待機

  call watchstop(mywatch)
  print *, "Elapsed time (min) = ", mywatch%time*60

end program test_watch