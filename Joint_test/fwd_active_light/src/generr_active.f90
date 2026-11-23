program generr_active
implicit none
character(3),allocatable,dimension(:) :: obsname 
character(3),allocatable,dimension(:) :: srcname
character(50) :: outputfolder,inputfolder
integer(4)    :: i,j,k,l,m,nobs,nsr,nfreq
character(100) :: infile,outamp,outpha
real(8) :: bxyz(2,3),exy(2,2),f,eamp,epha,ratio,r2d,pi

pi=4.*atan(1.d0)
r2d=180./pi

!#[1]## read param
read(*,*) ! header
read(*,'(20x,i10)') nfreq
read(*,'(20x,a)') inputfolder
read(*,'(20x,a)') outputfolder
read(*,'(20x,i10)') nobs
print *, "nobs =",nobs

allocate(obsname(nobs))

do i=1,nobs
 read(*,'(20x,a)') obsname(i)
end do

read(*,'(20x,i10)') nsr
print *, "nsr =",nsr
allocate(srcname(nsr))

do i=1,nsr
read(*,'(20x,a)') srcname(i)
end do

read(*,'(20x,g15.7)') ratio
print *, "ratio",ratio
!#[2]## amp data
do i=1,nsr
 do j=1,nobs
    infile=trim(inputfolder)//trim(obsname(j))//"_"//trim(srcname(i))//".dat"
    outamp=trim(outputfolder)//trim(obsname(j))//"_"//trim(srcname(i))//"_amp.dat"
    outpha=trim(outputfolder)//trim(obsname(j))//"_"//trim(srcname(i))//"_pha.dat"
    open(1,file=infile,status='old')
    open(2,file=outamp)
    open(3,file=outpha)
    
    do k=1,nfreq
        read(1,*)  f,((bxyz(l,m),l=1,2),m=1,3),((exy(l,m),l=1,2),m=1,2)
        eamp = bxyz(1,3)*ratio
        epha = atan2(ratio,1.d0)*r2d
        write(2,*) f,"1",bxyz(1,3),eamp
        write(3,*) f,"1",bxyz(2,3),epha
    end do
    
    close(1)
    close(2)
    close(3)
 end do
end do


!#[3]## phase data

end program generr_active
