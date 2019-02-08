subroutine DVRCore()

  use constants
  use util_mod, only : printheader, printversion
  use ReadInput, only : ReadInputMain
  use InputData, only : dvr_diag, dvr_integrals, orbital_ints
  use DVRData, only:  eigen_vecs, one_e_rad_int, two_e_rad_int
  use DVRDiag, only : SetDVR, DVRDiagonalization,ReadEigValVec
  use RHFData, only: tRHF
  use DVRIntRad, only : GetRadialElements
  use DVRIntAng, only : GetAngularElements
  use OrbInts, only  : GetOrbInts
  use PrimOrbInts, only  : GetPrimOrbInts
  use RadCheck
  use DVRRHF, only: DoRHF
  use FieldIntegrals

  implicit none

  character(64) :: file_in, file_out
  real(dp)  :: start, finish

  call cpu_time(start)

  file_in = 'dvr.inp'
  file_out = 'dvr.out'

  ! open the output file
  open (iout, File=file_out, status='UNKNOWN', form= "FORMATTED")

  ! Write down the version of the code and related information at the beginning of the code
  !call printheader(iout)
  call printversion(iout)

  ! Reading the input from an input file
  call ReadInputMain(file_in)

  ! Initializing the FE-DVR basis
  call SetDVR()

  ! Doing the Diagonalization, i.e. solving the radial Schroedinger equation
  ! If the option 'dvr_diag' is not true, then it reads the eigenvalues and 
  ! eigenvectors from the dat files
  if (dvr_diag) then
    call DVRDiagonalization()
  else 
    call ReadEigValVec()
  end if

  ! Now the integrals are calculated in the primitive DVR basis. 
  ! The radial and angular part of this integrals are calculate separately. 
  ! Radial parts corresponding to both one and two electron integral are calculated.
  !The angular only correspond to  the two electron integrals 
  if (dvr_integrals) then
    call GetRadialElements()
    call GetAngularElements()
  endif

  if (tRHF) then
    if (para%split_grid) then
      call DoRHFSplit()
    else
      call DoRHF(para%ng, eigen_vecs, one_e_rad_int, two_e_rad_int)
    end if
  end if

  ! Combine the 1-e and 2-e, radial and angular parts of integrals, calculated in the
  ! primitive basis to get the integrals in the eigenbasis as obtained from solving 
  ! the radial Schroedinger equation.
  if (orbital_ints) then
    if (prim_integrals) then
      call GetPrimOrbInts()
    else
      call GetOrbInts()
    end if
  end if

  ! If any external field is applied, corresponding integrals will be requiered and 
  ! these integrals are calculated here
  if (with_field) then
    call CalcFieldInts()
  end if 

! call radial_check()

! call check_with_analytic()

  call DeallocMatrices()

  write(iout, *) 'The calculation is done. '
  
  call cpu_time(finish)
  write(iout,*) 'Time = ', finish-start, 'seconds.'
  close(iout)

end subroutine DVRCore

subroutine DoRHFSplit()

  use constants
  use util_mod, only: stop_all
  use DVRData, only : para, eigen_vecs, one_e_rad_int, two_e_rad_int
  use DVRRHF, only: DoRHF

  real(dp), allocatable :: EigVec(:,:,:), OneInts(:,:,:), TwoInts(:,:,:)
  integer :: len_1, len_2

    len_1 = para%m1*para%nl
    len_2 = para%m2*para%nl - 1

    ! Case I
    if (para%diagtype == 'only_inner') then

      allocate(EigVec(len_1, len_1, para%l+1))
      allocate(OneInts(len_1, len_1, 2*para%l+1))
      allocate(TwoInts(len_1, len_1, 2*para%l+1))

      do i = 1, len_1
        do j = 1, len_1
          EigVec(i,j,:) = eigen_vecs(i,j,:)
          OneInts(i,j,:) = one_e_rad_int(i,j,:)
          TwoInts(i,j,:) = two_e_rad_int(i,j,:)
        end do
      end do

      call DoRHF(len_1, EigVec, OneInts, TwoInts)
    
      eigen_vecs = 0.0d0
      do i = 1, len_1
        do j = 1, len_1
          eigen_vecs(i,j,:) = EigVec(i,j,:)
        end do
      end do

      deallocate(EigVec, OneInts, TwoInts)

    ! Case II
    elseif (para%diagtype == 'only_outer') then

      allocate(EigVec(len_2, len_2, para%l+1))
      allocate(OneInts(len_2, len_2, 2*para%l+1))
      allocate(TwoInts(len_2, len_2, 2*para%l+1))

      do i = 1, len_2
        do j = 1, len_2
          EigVec(i,j,:) = eigen_vecs(i+len_1,j+len_1,:)
          OneInts(i,j,:) = one_e_rad_int(i+len_1,j+len_1,:)
          TwoInts(i,j,:) = two_e_rad_int(i+len_1,j+len_1,:)
        end do
      end do

      call DoRHF(len_2, EigVec, OneInts, TwoInts)
    
      eigen_vecs = 0.0d0
      do i = 1, len_2
        do j = 1, len_2
          eigen_vecs(i+len_1,j+len_1,:) = EigVec(i,j,:)
        end do
      end do

      deallocate(EigVec, OneInts, TwoInts)

    ! Case III
    elseif ((para%diagtype == 'both').or.(para%diagtype == 'dont')) then

      call DoRHF(para%ng, eigen_vecs, one_e_rad_int, two_e_rad_int)

    else
       
      call stop_all('DoRHFSplit','Not a possible option to run RHF')

    end if

end subroutine DoRHFSplit

subroutine DeallocMatrices()
   
  use constants, only : iout
  use DVRData, only : eigen_vals, eigen_vecs, one_e_rad_int, &
    &                 two_e_rad_int, integrals_ang, sph_harm
  use OrbData, only : SpatialOrbInd, OneEInts, TwoEInts, TwoERadOrbInts
  
  if (allocated(eigen_vals)) deallocate(eigen_vals)
  if (allocated(eigen_vecs)) deallocate(eigen_vecs)
  if (allocated(one_e_rad_int)) deallocate(one_e_rad_int)
  if (allocated(two_e_rad_int)) deallocate(two_e_rad_int)
  if (allocated(integrals_ang)) deallocate(integrals_ang)

  if (allocated(SpatialOrbInd)) deallocate(SpatialOrbInd)
  if (allocated(OneEInts)) deallocate(OneEInts)
  if (allocated(TwoEInts)) deallocate(TwoEInts)
  if (allocated(TwoERadOrbInts)) deallocate(TwoERadOrbInts)
  
  if (allocated(sph_harm%n_m)) deallocate(sph_harm%n_m)
  if (allocated(sph_harm%m_init)) deallocate(sph_harm%m_init)

  write(iout, *) 'Deallocated all the main matrices.'

end subroutine DeallocMatrices
