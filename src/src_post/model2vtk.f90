!# coded on 2018.11.08
program model2vtk
    use modelpart
    use param
    use mesh_type
    implicit none
    
    type(model)            :: g_model(3)   ! see m_modelpart.f90 2022.12.19
    type(param_cond)       :: g_cond(3)    ! see m_param.f90 : goal conductivity 2022.12.19
    type(mesh)             :: g_mesh
    character(50)          :: mshfile,connectfile,modelfile1,modelfile2,outfile!2022.12.19
    integer(4)             :: nphys1,nphys2, nmodel,i,j
    integer(4)             :: nlin,node,ntri,ishift,ntet
    integer(4),allocatable,dimension(:) :: id,ele2model,index
    real(8),   allocatable,dimension(:,:) :: value ! 2022.12.19
    real(8),   allocatable,dimension(:)   :: resistivity ! 2022.12.19
    !# for mask
    real(8)               :: threshold
    character(50)         :: sensfile
    integer(4)            :: iunit ! 0 : Ohm.m, 1: log10(Ohm.m) 2021.01.08
    integer(4)            :: imask ! 0:nomask, 1:mask
    real(8)               :: val_mask(3) !1 for model1, 2for model2, 3for diff
    character(50),dimension(3) :: var_label ! 2022.12.19

    
    !#[0]## read parameters
      call readfilename(mshfile,connectfile,modelfile1,modelfile2,var_label,&
    &                       iunit,imask,threshold,sensfile,val_mask,outfile)!2022.12.19
    
    !#[1]## read mesh
     write(*,*) "mshfile",mshfile
      CALL READMESH_TOTAL(g_mesh,mshfile)
    
    !#[2]## read model
    !#[2-1]## read connectfile
      open(1,file=connectfile)
      read(1,*) nphys2,nmodel
      allocate(id(nphys2),ele2model(nphys2))
      do i=1,nphys2
       read(1,'(2i10)') id(i),ele2model(i)
      end do
      close(1)
      write(*,*) "### READ CONNECTIONFILE END!! ###"
    
     !#[2-2]## read model1
     allocate(value(3,nmodel))
     open(1,file=modelfile1,status='old',err=90) ! 2017.12.21
     open(2,file=modelfile2,status='old',err=90) ! 2017.12.21
      read(1,*,err=80) nmodel
      read(2,*,err=80) nmodel
      do i=1,nmodel
       read(1,*) value(1,i) ! 2022.12.19
       read(2,*) value(2,i) ! 2022.12.19
      end do
     close(1)
     close(2) ! 2022.12.19
     write(*,*) "### READ MODELFILE END!! ###"
    
     !#
     if (imask .eq. 1 ) then
      write(*,*) "sensfile =",sensfile
      open(1,file=sensfile)
       read(1,*) nmodel
       do i=1,nmodel
        read(1,*) value(3,i) ! sensitivity values
       end do
      close(1)
     end if
     !#
     goto 100
      80 continue
      close(1)
      write(*,*) "Cannot read modelfile"
      stop
      goto 100
      90 continue
      write(*,*) "File is not exist",trim(modelfile1),trim(modelfile2)
      stop
     100 continue ! 2018.01.17
     nphys1 = g_mesh%ntet - nphys2
     write(*,*) "nphys1",nphys1
     write(*,*) "nphys2",nphys2 ! 2020.12.02
     allocate(index(nphys2))
     do i=1,nphys2
      index(i) = nphys1+i
     end do
    
     !#[2-3]## set model
      do i=1,3 ! 2023.01.02 1:model1, 2:model2, 3:model change
       g_model(i)%nmodel       = nmodel
       g_model(i)%nphys1       = g_mesh%ntet -nphys2
       g_model(i)%nphys2       = nphys2
       g_model(i)%ele2model    = ele2model
       g_model(i)%index        = index
       if ( i == 3 ) cycle ! 2023.01.02
       if ( iunit .eq. 0 ) then             ! 2021.01.08
         g_model(i)%logrho_model = log10(value(i,:)) ! 2021.01.08
       elseif (iunit .eq. 1) then           ! 2021.01.08
         g_model(i)%logrho_model = value(i,:)        ! 2021.01.08
       else                                 ! 2021.01.08
         write(*,*) "GEGEGE iunit",iunit     ! 2021.01.08
         write(*,*) "input should be 0: [Ohm.m] or 1:(log10(Ohm.m))" ! 2021.01.08
         stop  ! 2021.01.08
       end if ! 2021.01.08

      end do ! 2022.12.19

    !#[2-4]## generate change model 2023.01.02
     allocate(g_model(3)%logrho_model(nmodel))
     do j=1,nmodel
!      write(*,*) "j",j
      g_model(3)%logrho_model(j) = log10(value(2,j)) - log10(value(1,j))
     end do

    !#[2-5]## apply mask
      if ( imask .eq. 1) then
        do j=1,nmodel

         if ( value(3,j) .lt. threshold ) then
           do i=1,3
             g_model(i)%logrho_model(j) = val_mask(i) ! replaced by mask value
           end do
         end if

        end do
       end if


    !#[3]## generate cond for 2 models
      do i=1,3

       !#[3-1]## initial g_cond
       g_cond(i)%condflag  = 0 ! 0:homogeneous
       g_cond(i)%sigma_air = 1.d-8 ! cond for air
       g_cond(i)%nvolume   = 1  ! # of physical volume in land
       allocate(g_cond(i)%sigma_land(1))
       g_cond(i)%sigma_land(1)=0.01
    
       !#[4]## generate cond
       if (iunit .eq. 0 ) then ! unit [Ohm.m] 2021.01.08
         call model2cond(g_model(i),g_cond(i),0) ! 0 means g_cond%rho = 10**g_model%logrho_model
       elseif (iunit .eq. 1) then ! unit log10([Ohm.m]) 2021.01.08
         call model2cond(g_model(i),g_cond(i),1) ! 1 means g_cond%rho = g_model%logrho_model
       end if

      end do ! 2022.12.19
    
    !#[5]## output
     node = g_mesh%node
     nlin = g_mesh%nlin
     ntri = g_mesh%ntri
     ntet = g_mesh%ntet
     open(1,file=outfile)
     write(1,'(a)') "# vtk DataFile Version 3.0"
     write(1,'(a)') "vtk output"
     write(1,'(a)') "ASCII"
     write(1,'(a)') "DATASET UNSTRUCTURED_GRID"
     write(1,'(a,i10,a)') "POINTS",node," double"
     do i=1,node
       write(1,'(3g15.7)') g_mesh%xyz(1:3,i) 
     end do
!     write(1,'(a,2i8)') "CELLS",ntet,ntet*5
     write(1,'(a,2i8)') "CELLS",nphys2,nphys2*5
     do i=nphys1+1,ntet
      write(1,*) 4,g_mesh%n4(i,1:4)-1 ! vtk use index starting with 0
     end do
     write(1,'(a,i10)') "CELL_TYPES",nphys2
     write(1,'(10000000i3)') (10,i=1,nphys2)
     write(1,'(a,i10)') "CELL_DATA",nphys2
     allocate(resistivity(nphys2))
!     resistivity(1:nphys1) = 1.d+8 ! air

     !#[1]## resistivity 1
     write(1,'(a)') "SCALARS "//trim(var_label(1))//" double"
     write(1,'(a)') "LOOKUP_TABLE default"
     resistivity(1:nphys2)=g_cond(1)%rho(1:nphys2) ! earth
     write(1,101) (resistivity(i),i=1,nphys2)

     !#[2]## resistivity 2
     write(1,'(a)') "SCALARS "//trim(var_label(2))//" double"
     write(1,'(a)') "LOOKUP_TABLE default"
     resistivity(1:nphys2)=g_cond(2)%rho(1:nphys2) ! earth
     write(1,101) (resistivity(i),i=1,nphys2)

    !#[3]## resistivity change
     write(1,'(a)') "SCALARS "//trim(var_label(3))//" double"
     write(1,'(a)') "LOOKUP_TABLE default"
!     resistivity(1:nphys1) = 0.d0 ! air
     resistivity(1:nphys2)=g_cond(3)%rho(1:nphys2) ! resistivity change 2023.01.02
     write(1,101) (resistivity(i),i=1,nphys2)

     99 continue
    close(1)
    
    101 format(10000000g15.7)
    end program
    
    !#######################################################################
    subroutine readfilename(mshfile,connectfile,modelfile1,modelfile2,var_label,&
    &                       iunit,imask,threshold,sensfile,val_mask,outfile)!2021.01.08 iunit is added
    implicit none
    character(50),intent(out) :: mshfile,connectfile,modelfile1,modelfile2,outfile,sensfile
    character(50),dimension(3) :: var_label
    real(8),  intent(out)  :: val_mask(3) ! 2023.01.02
    integer(4)             :: iunit,i ! 2021.01.08
    integer(4),intent(out) :: imask   ! 2021.01.02
    real(8)                :: threshold
    
    write(*,*) "input mshfile"
    read(5,'(a)') mshfile
    write(*,*) "mesh file",trim(mshfile)
    write(*,*) "connect connectfile"
    read(*,'(a)') connectfile
    write(*,'(a,a)') "connectfile: ",trim(connectfile)
    write(*,*) "input modelfile"
    read(*,'(a)') modelfile1
    write(*,'(a,a)') "modelfile1: ",trim(modelfile1)
    read(*,'(a)') modelfile2 ! 2022.12.19
    write(*,'(a,a)') "modelfile2: ",trim(modelfile2) ! 2022.12.19
    read(*,'(a)') var_label(1) !2022.12.19
    read(*,'(a)') var_label(2) !2022.12.19
    read(*,'(a)') var_label(3) !2022.12.19
    write(*,'(a,2x,a,2x,a)')  (trim(var_label(i)),i=1,3) ! 2022.12.19
    write(*,*) "parameter input : 0 for Ohm.m , 1 for log10(Ohm.m)" ! 2021.01.08
    read(*,'(i5)') iunit ! 0 for Ohm.m , 1 for log10(Ohm.m)
    write(*,*) "input output vtk file" ! 2023.01.02
    read(5,'(a)') outfile

    ! 2023.01.02
    write(*,*) "0 for no mask, 1 for mask model"
    read(*,'(i5)') imask ! 2023.01.02
    write(*,*) "imask =",imask
    if ( imask .eq. 1 ) then
     write(*,*) "input sensitivity model"
     read(*,'(a)') sensfile
     write(*,*) " input threshold"
     read(5,*)  threshold    ! sensitivity threshold
     read(*,*) val_mask(1) ! for model 1
     read(*,*) val_mask(2) ! for model 2
     read(*,*) val_mask(3) ! for difference
    end if

    return
    end
    
    